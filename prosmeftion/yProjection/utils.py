import logging
from fractions import Fraction
from typing import List

import sage.all
import sage.matrix as mx
from sage.rings.rational_field import QQ
from .index import Index
from .indices import Indices_Operator
from .operators import Tensors
from .read_write import get_terms, get_type
from .summand import Summand
# from .get_FORM.form_read_in import Term_Model, Term_s
from .term import Term, TermType

from autoeft.base.tensors import SUNTableau  # TODO: Check autoeft.invariants -> autoeft.base.tensors
# , field_projection_operator, symmetrize_tensors
# FIXME: field_projection_operator, symmetrize_tensors not defined in new AutoEFT
# from autoeft.model import SUNGroup
# from autoeft.sun_projection import tensor_projection
from .. import FORM_PATH, FORM_GENERAL_PATH, get_antisymEps, op_config, bosons, fermions, tensors, \
    run_form

logger_autoeft = logging.getLogger("autoeft.projection")
logger = logger_autoeft.getChild(__name__)

def equalize_indices(ref_term: Summand, eq_term: Summand) -> (List[Index], List[Index]):
    """
    Renames indices of the fields for a given Summand eq_term in the same way as in the reference Summand ref_term
    and replaces them accordingly in the tensors of the eq_term.
    Parameters
    ----------
    ref_term
        Reference Summand.
    eq_term
        Summand for which the indices are equalized.
    Returns
    -------
    """
    # Note the order of indices on a field is very specific and specified in the op_config.yml file.
    # Nevertheless, the order of derivative indices is not specific and also derivative indices might stand
    # before or behind the field indices.
    #
    for ref_field, eq_field in zip(ref_term.fields, eq_term.fields):
        assert ref_field.name == eq_field.name
        skip_field_indices = []  # Those are the indices of the non derivative indices on a field which are already replaced -> Important if there are 2 identical non derivative indices in a field.
        for i, ref_index in enumerate(ref_field.indices):
            index_tensor_found = False
            if der := ref_index.derIndex:
                for j, eq_index in enumerate(eq_field.indices):
                    # Assure that the correct derivative index is changed, i.e. the index of the correct derivative of the correct type.
                    if eq_index.derIndex == der and eq_index.typ == ref_index.typ:
                        for tensor in eq_term.tensors:
                            for k, tensor_index in enumerate(tensor.indices):
                                # find the tensor in the eq_term with which the field is contracted
                                if tensor_index.indname == eq_index.indname:
                                    index_tensor_found = True
                                    tensor.indices[k] = ref_index.dual_index
                                    break
                            if index_tensor_found: break
                    # Replace Field index afterwards to ensure that contraction stay the same.
                    if index_tensor_found:
                        eq_field.indices[j] = ref_index
                        break
            else:
                for j, eq_index in enumerate(eq_field.indices):
                    if j in skip_field_indices:
                        # This index is already replaced.
                        continue
                    # Assure that it is not a derivative index
                    if not eq_index.derIndex and eq_index.typ == ref_index.typ:
                        for tensor in eq_term.tensors:
                            for k, tensor_index in enumerate(tensor.indices):
                                # find the tensor in the eq_term with which the field is contracted
                                if tensor_index.indname == eq_index.indname:
                                    index_tensor_found = True
                                    tensor.indices[k] = ref_index.dual_index
                                    break
                            if index_tensor_found: break
                    # Replace Field index afterwards to ensure that contraction stay the same.
                    if index_tensor_found:
                        eq_field.indices[j] = ref_index
                        skip_field_indices.append(j)
                        break
                # for tensor in eq_term.tensors:
                #     for j, tensor_index in enumerate(tensor.indices):
                #         # find the tensor in the eq_term with which the field is contracted
                #         if tensor_index.indname == eq_field.indices[i].indname:
                #             index_tensor_found = True
                #             # print(f"{index} -> {index.dual_index}")
                #             tensor.indices[j] = ref_index.dual_index
                #             break
                #     if index_tensor_found: break
                # # Replace Field index afterwards to ensure that contraction stay the same.
                # eq_field.indices[i] = ref_index
    #  Contractions within tensors like for example between 2 Yukawa matrices are still allowed, but occur also
    #  only between Yukawa matrices.
    ref_tensor_indices = Indices_Operator([])  # All indices occuring in the tensors of ref_term
    for ref_tensor in ref_term.tensors:
        ref_tensor_indices += ref_tensor.indices
    eq_tensor_indices = Indices_Operator([])  # All indices occuring in the tensors of eq_term
    for eq_tensor in eq_term.tensors:
        eq_tensor_indices += eq_tensor.indices

    ref_tensor_indices = list(ref_tensor_indices.indices)
    eq_tensor_indices = list(eq_tensor_indices.indices)
    # Remove double indices in eq_tensor_indices and ref_tensor_indices
    for index in ref_tensor_indices.copy():
        while index.is_in(ref_tensor_indices) > 1:
            ref_tensor_indices.remove(index)
    for index in eq_tensor_indices.copy():
        while index.is_in(eq_tensor_indices) > 1:
            eq_tensor_indices.remove(index)
    # Remove all indices from eq_tensor_indices which are already in ref_tensor_indices, so that only non renamed
    # indices in eq_tensor_indices remain:
    eq_tensor_indices_copy = eq_tensor_indices.copy()
    for index in eq_tensor_indices_copy:
        if index.is_in(ref_tensor_indices):
            eq_tensor_indices.remove(index)
    for index in ref_tensor_indices.copy():
        if index.is_in(eq_tensor_indices_copy):
            ref_tensor_indices.remove(index)
    eq_tensor_indices = list(eq_tensor_indices)
    ref_tensor_indices = list(ref_tensor_indices)
    # Worst case:
    # In  eq_term: [ye+](flav2765,flav3459)*ye(flav3459,flav2991)*[ye+](flav2991,flav3648)
    # In ref_term: [ye+](flav2765,flav2917)*ye(flav2917,flav3140)*[ye+](flav3140,flav3648)
    # -> Sum over contractions within tensors in FORM, so that FORM will generate dummy indices
    # and the naming isn't important anymore.
    # In this example: sum flav3459, flav2991, flav2917, flav3140;
    # Read in output and replace dummy indices by newly generated flavor indices, which will be generated by
    # a generator in python, respecting already set indices.
    return ref_tensor_indices, eq_tensor_indices


def equalize_field_indices(single_terms):
    """
    Equalize indices in expression of same fieldstructure, by comparison of the field indices.
    Parameters
    ----------
    single_terms

    Returns
    -------
        List of in tensors fully contracted and therefore not equalized indices.
    """

    for type in single_terms.values():
        for term_mass_dim in type.values():
            term_with_specific_field_structure = {}
            for term in term_mass_dim:
                try:
                    term_with_specific_field_structure[term.fieldstructure].append(term)
                except KeyError:
                    term_with_specific_field_structure[term.fieldstructure] = [term]
            for name_of_term, terms_specific in term_with_specific_field_structure.items():
                if len(terms_specific) > 1:
                    for i, term in enumerate(terms_specific[1:]):
                        if not i:
                            # i == 0
                            ref_tensor_indices, eq_tensor_indices = equalize_indices(terms_specific[0], term)
                        else:
                            _, eq_tensor_indices_tmp = equalize_indices(terms_specific[0], term)
                            eq_tensor_indices += eq_tensor_indices_tmp
                    sum_indices = ref_tensor_indices + eq_tensor_indices

    return single_terms


def renew_indices(ref_term: Summand, fp_exclude_indices: List[Index]) -> Summand:
    """
    Assign to all FIELDS of the given term new indices. New in this case means that they do not occur in the old expression
    AND they do not occur in the specified list of indices fp_exclude_indices.

    Note: Contractions among tensors, are not replaced.

    Parameters
    ----------
    ref_term
       Summand which gets new indices.
    exclude_indices
        List of indices which ar not chosen for the new indices.
    Returns
    -------
    """
    # First, rewrite the indices of all FIELDS in the ref_term by new indices, which never occur in any of the terms in terms_specific
    gen_index = {}
    gen_index["gauge"] = ref_term.possible_indices.generate_index("gauge",  exclude_indices=fp_exclude_indices)  # get a new, i.e. unused SU2 index with 'next(gen_index_gauge)'
    gen_index["colf"]  = ref_term.possible_indices.generate_index("colf",   exclude_indices=fp_exclude_indices)
    gen_index["flav"]  = ref_term.possible_indices.generate_index("flav",   exclude_indices=fp_exclude_indices)
    gen_index["sl"]    = ref_term.possible_indices.generate_index("Lsl",    exclude_indices=fp_exclude_indices)  # next(gen_index_sl) will generate new Lsl and new Usl index
    gen_index["sldot"] = ref_term.possible_indices.generate_index("Lsldot", exclude_indices=fp_exclude_indices)

    for ref_field in ref_term.fields:
        for i, field_index in enumerate(ref_field.indices):
            index_tensor_found = False
            for tensor in ref_term.tensors:
                for j, tensor_index in enumerate(tensor.indices):
                    if tensor_index.indname == field_index.indname:
                        index_tensor_found = True
                        #  Create new indices
                        if tensor_index.typ == "Lsl":
                            new_tensor_index, new_field_index = next(gen_index["sl"])
                        elif tensor_index.typ == "Usl":
                            new_field_index, new_tensor_index = next(gen_index["sl"])
                        elif tensor_index.typ == "Lsldot":
                            new_tensor_index, new_field_index = next(gen_index["sldot"])
                        elif tensor_index.typ == "Usldot":
                            new_field_index, new_tensor_index = next(gen_index["sldot"])
                        else:
                            new_tensor_index = next(gen_index[tensor_index.typ])
                            new_field_index = new_tensor_index.dual_index

                        new_tensor_index.derIndex = tensor_index.derIndex
                        new_field_index.derIndex  = field_index.derIndex

                        tensor.indices[j] = new_tensor_index
                        break
                if index_tensor_found: break
            # Replace Field index afterwards to ensure that contraction stay the same.
            if index_tensor_found:
                ref_field.indices[i] = new_field_index

    return ref_term


def remove_doubles(single_terms, model):
    """
    Remove terms which occur in the exact same way, i.e. only with a different coefficient, more than ones.
    Note: The equalize_field_indices function cannot be used, since the not renamed indices in 'sum_indices' had
    to be known at the time, when they are combined in FORM. It could be possible to save those indices for all types,
    but this isn't done for now.

    Parameters
    ----------
    single_terms

    Returns
    -------
    """
    merged_terms = []
    for type_name, type in single_terms.items():
        for term_mass_dim in type.values():
            term_with_specific_field_structure = {}
            for term in term_mass_dim:
                try:
                    term_with_specific_field_structure[term.fieldstructure].append(term)
                except KeyError:
                    term_with_specific_field_structure[term.fieldstructure] = [term]
            for name_of_term, terms_specific in term_with_specific_field_structure.items():
                if len(terms_specific) > 1:
                    name_form = "".join([f"{name}{nD}" for name, nD in name_of_term])
                    logger.info(f"Combine {', '.join([term.name for term in terms_specific])} of type {name_form}.")
                    ref_term = terms_specific[0]

                    excl_indices = terms_specific[1].possible_indices
                    if len(terms_specific) > 2:
                        for term in terms_specific[2:]:
                            excl_indices += term.possible_indices

                    # assign new indices to the ref_term
                    ref_term = renew_indices(ref_term, excl_indices)

                    for i, term in enumerate(terms_specific[1:]):
                        if not i:
                            # i == 0
                            ref_tensor_indices, eq_tensor_indices = equalize_indices(ref_term, term)
                        else:
                            _, eq_tensor_indices_tmp = equalize_indices(ref_term, term)
                            eq_tensor_indices += eq_tensor_indices_tmp
                    sum_indices = ref_tensor_indices + eq_tensor_indices
                    form = "Off statistics;\n"
                    form += "#include declarations_general.h # coefficient\n"
                    form += "#include declarations_general.h # indices\n"
                    form += "#include declarations_general.h # tensors\n"
                    form += "#include declarations_general.h # operators\n"
                    form += "\n"

                    exprs = [f"{term:c}"for term in terms_specific]
                    for i, expr in enumerate(exprs):
                        form += f"Local expr{i:d} = {expr};\n"#
                    form += "\n"
                    list_expr_names = [f"expr{i:d}" for i in range(len(exprs))]
                    form += f"Local expr = {' + '.join(list_expr_names)};\n"
                    form += ".sort\n"
                    form += f"Drop {', '.join(list_expr_names)};\n"
                    form += "\n"
                    if sum_indices:
                        form += f"sum {', '.join(map(str, sum_indices))};\n"
                        form += "\n"
                    form += "* Bring indices of epsilons in order:\n"
                    epss = [list(tensor["mathematica"].values())[0] for tensor_name, tensor in
                           op_config["tensors"].items() if
                           "eps" in tensor_name]
                    for eps in epss:
                        antisymeps = get_antisymEps(eps)
                        form += f"Multiply replace_({eps},{antisymeps});\n"
                        form += ".sort\n"
                        form += f"Multiply replace_({antisymeps},{eps});\n"
                        form += ".sort\n"
                    form += "\n"
                    form += f"Bracket {', '.join(tensors + bosons + fermions)};\n\n"

                    name_form = "".join([f"{name}{nD}" for name, nD in name_of_term])
                    TERM_PATH = FORM_PATH / name_form

                    form += f'#write <{TERM_PATH / "combined_term.h"}> "%E", expr\n'

                    form += "Print +ss;\n"
                    form += ".end\n"

                    TERM_PATH.mkdir(parents=True, exist_ok=True)  # Create directories if they don't exist.
                    with open(TERM_PATH / f"{name_form}.frm", "w") as file:
                        file.write(form)

                    run_form(fp_cwd=TERM_PATH, filename=f"{name_form}.frm", fp_p=FORM_GENERAL_PATH)

                    with open(TERM_PATH / "combined_term.h", "r") as file:
                        expression = file.read()
                    terms = get_terms(expression, model, name=name_form)

                    if terms.indices["dummy"]:
                        for k, summand in enumerate(terms):
                            if len(terms) > 1:
                                for j in range(len(terms)):
                                    if j != k:
                                        excl_indices = terms_specific[j].possible_indices
                                        ref_ind = j
                                        break
                                if len(terms) > 2:
                                    for term in [term for l, term in enumerate(terms) if l not in (ref_ind, k)]:
                                        excl_indices += term.possible_indices

                            dummy_indices = summand.indices["dummy"]
                            contracted = {index.expr: 0 for index in dummy_indices}
                            dummy_to_new_index = {}

                            # Note: In the following, only dummy indices in the tensors are substituted.
                            for tensor in summand.tensors:
                                for i, index in enumerate(tensor.indices):
                                    if index.is_in(dummy_indices):
                                        if contracted[index.expr] == 0:
                                            gen_index = summand.possible_indices.generate_index(index.will_be_typ, exclude_indices=excl_indices)
                                            dummy_to_new_index[index.expr] = next(gen_index)

                                            contracted[index.expr] += 1
                                            tensor.indices[i] = dummy_to_new_index[index.expr]
                                        elif contracted[index.expr] == 1:
                                            contracted[index.expr] += 1
                                            tensor.indices[i] = dummy_to_new_index[index.expr].dual_index
                                if all([True if contract == 2 else False for contract in contracted.values()]):
                                    break
                            assert all([True if contract == 2 else False for contract in contracted.values()]), "Something went wrong in the dummy index replacement."

                    merged_terms.append(terms)
                    # print("--------------------------")
                else:
                    name_form = "".join([f"{name}{nD}" for name, nD in name_of_term])
                    logger.info(f"{terms_specific[0].name} of type {name_form} is not combined")
                    # Rename Summand properly by their structure
                    terms_specific[0].name = name_form
                    merged_terms.append(Term(terms_specific, name_form))
                    # print("--------------------------")
    single_terms = get_type(merged_terms, model)
    return single_terms




