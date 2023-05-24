import logging
import sys
from fractions import Fraction
from pathlib import Path
from typing import List, Dict

import sage.all
import sage.matrix as mx
from sage.rings.rational_field import QQ
from .yProjection.basisTensors import MonBasisTensor, SymBasisTensor, MonBasisTensors, SymBasisTensors
from .yProjection.index import Index
from .yProjection.indices import Indices_Operator
from .yProjection.operator import Tensor
from .yProjection.operators import Tensors
from .yProjection.read_write import get_terms
from .yProjection.summand import Summand
# from .get_FORM.form_read_in import Term_Model, Term_s
from .yProjection.term import Term, TermType

from autoeft.invariants import SUNTableau, field_projection_operator, symmetrize_tensors
from autoeft.model import SUNGroup
from autoeft.sun_projection import tensor_projection
from . import AUTOEFT_PATH, FORM_PATH, FORM_GENERAL_PATH, model, get_antisymEps, op_config, bosons, fermions, tensors, \
    run_form, get_basis
logger_autoeft = logging.getLogger("autoeft.projection")
logger = logger_autoeft.getChild(__name__)

def get_type(terms):
    """
    Summands are sorted by their "type", i.e. field content and derivative.
    Parameters
    ----------
    terms
        List of Term objects, where each Term object contains one or more Summands.
    Returns
    -------
        Sorted terms.
    """
    logger.info("Sort fields by type, indicated by a tuple filled with integers. They specify the "
                "number of fields in the single summand in the following order: "
                f"{' '.join([field.name for field in model.fields.values()])}")
    single_terms = {}  # Ordered terms (by "type") with just a single term in it.
    for term in terms:
        if repr(term) == "0":
            continue
        for summand in term:
            typ = tuple(summand.fieldcounter.values())
            try:
                type(single_terms[typ])
                try:
                    single_terms[typ][summand.nD].append(summand)
                except KeyError:
                    single_terms[typ][summand.nD] = TermType(summand, summand.fieldcounter_stripped)
            except KeyError:
                single_terms[typ] = {summand.nD: TermType(summand, summand.fieldcounter_stripped)}

    return single_terms


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

def remove_doubles(single_terms):
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
                    terms = get_terms(expression, name=name_form)

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
    single_terms = get_type(merged_terms)
    return single_terms

def sun_internal_projector(op_type, group):
    """Construct a projector to impose the 'internal symmetries'."""
    field_id = 0 #!
    projection_operator = []
    for field, multiplicity in op_type.fields:
        for i in range(multiplicity):
            field_id += 1 #!
            field_tableau = SUNTableau.fill_tableau(
                field_id, field.sun_reprs[group].partition
            )
            projection_operator.append(d:=field_projection_operator(field_tableau))
            # print(field_id, d)
    return projection_operator

def yt_to_tensors(yt: SUNTableau, field_content: Dict, sun_group: SUNGroup, factor: Fraction):
    """
    Create Note Schreibe als gaugeF1I2,...
    Parameters
    ----------
    yt
    sun_group

    Returns
    -------
    """
    tensors = []
    if sun_group.N == 2:
        ind_prefix = "gauge"
    elif sun_group.N == 3:
        ind_prefix = "colf"
    else:
        logger.error("No index prefix defined for this group.")
        sys.exit("STOP")
    for column in zip(*yt.tableau):
        expr = f"[su{sun_group.N:d}eps]("
        expr += ",".join(f"{ind_prefix}F{index.field_num}I{index.index_num}" for index in column)
        expr += ")"
        tensor = Tensor(expr)
        for index in tensor.indices:
            index.projection = index.expr
        tensors.append(tensor)

    return MonBasisTensor(sun_group.name, field_content, Tensors(tensors), factor)

def basis_tensors(op_type, field_content, tensors, group):
    """Return monomial and symmetrized tensors in FORM-readable output."""
    internal_projector = sun_internal_projector(op_type, group)

    monomial_basis = []
    tensor_basis = []
    for i, tensor in enumerate(tensors):
        # Monomial basis tensors
        tensorI = yt_to_tensors(tensor, field_content, group, Fraction(1))
        monomial_basis.append(tensorI)
        # print(f"MonBasisTensor: {tensorI:a}")

        # Symmetrized basis tensors
        sym_tensors = symmetrize_tensors([(tensor, Fraction(1))], internal_projector)
        sym_basis_tensor = []
        for sym_tensor in sym_tensors:
            sym_basis_tensor.append(yt_to_tensors(sym_tensor[0], field_content, group, sym_tensor[1]))
        sym_basis_tensor = SymBasisTensor(sym_basis_tensor, i)
        # print(f"SymBasisTensor: {sym_basis_tensor:a}")
        tensor_basis.append(sym_basis_tensor)
    return MonBasisTensors(monomial_basis), SymBasisTensors(tensor_basis)

def get_basis_tensors(basis, field_content, derivatives, mass_dim):
    sun_projection_tensors = {}
    # get all the invariants associated with the operator
    op_class, op_subclass, op_type, operator = basis[mass_dim].get_operator(field_content, derivatives)
    for sun_group in model.sun_groups.values():
        # # get all the invariants associated with the operator
        # op_class, op_subclass, op_type, operator = basis[mass_dim].get_operator(field_content, derivatives)
        # get the SU(N) tensors of sun_group
        sun_basis = operator.sun_tensors[sun_group]
        # convert the tensors to FORM-readable output
        sun_monom, sun_tensor = basis_tensors(op_type, field_content, sun_basis, sun_group)
        # sun_monom_test, sun_tensor_test = form_basis(op_type, sun_basis, sun_group)
        if sun_monom:
            sun_projection_tensors[sun_group.name] = {"sun_monom": sun_monom, "sun_tensor": sun_tensor}
        else:
            sun_projection_tensors[sun_group.name] = False
    return sun_projection_tensors

def sun_projection(single_terms, basis, max_dim: int):
    logger.info(f"Project all terms for a specific type onto the {'-, '.join(list(model.sun_groups.keys()))}-group "
                f"basis when the term of specific type is part of the basis.")
    for type in single_terms.values():
        for term_mass_dim in type.values():
            field_content = term_mass_dim.field_content
            derivatives = term_mass_dim.nD
            mass_dim = term_mass_dim.d
            try:
                # get the SUN basis Tensors from autoeft for the projection
                term_mass_dim.sun_projection_tensors = get_basis_tensors(basis, field_content, derivatives, mass_dim)
                if not any(term_mass_dim.sun_projection_tensors.values()):
                    # This case occurs, when the term is part of the basis and thus no KeyError exception is thrown, but still doesn't has a projection.
                    # The only reason this could be is that the term doesn't contain any of the indices of all sun-groups, i.e. is a singlet of all of them.
                    assert not term_mass_dim.indices["gauge"]
                    assert not term_mass_dim.indices["colf"]
                    term_mass_dim.sun_projection_tensors = False
                    continue
                # There exists a term in the basis which matches the type and the projection can be done:
                # sun_tensors = list of all sun tensors contained in all terms of the same type.
                for sun_group, tensors in term_mass_dim.sun_projection_tensors.items():
                    if tensors:
                        N = model.sun_groups[sun_group].N
                        sun_monom  = tensors["sun_monom"]
                        sun_tensor = tensors["sun_tensor"]
                    else:
                        # Since the basis doesn't contain basis tensors for this group, also the tensor itself shouldn't
                        # and thus the next possible SUN group is tried.
                        continue
                    # collect all sun tensors of the group and the tensor in one ordered list for the projection.
                    sun_field_tensors = []
                    for term in term_mass_dim:
                        sun_field_tensors.append(term.gaugeTensorsSUN[sun_group])

                    # project all tensors simultaneously
                    sun_monomial_basis = [f"{tensor:p}" for tensor in sun_monom]
                    sun_tensor_basis = [f"{tensor:p}" for tensor in sun_tensor]
                    P_map, G_map = tensor_projection(N, sun_monomial_basis, sun_tensor_basis, [f"{tensor:p}" for tensor in sun_field_tensors])
                    dim = len(sun_monom)
                    n_projection_op = len(sun_field_tensors)
                    G = mx.constructor.matrix(QQ, dim, dim, G_map)
                    P = mx.constructor.matrix(QQ, n_projection_op, dim, P_map)
                    # P*G^(-1) gives the projection matrix
                    projection_matrix = P * G.inverse()
                    term_mass_dim.sun_projection_matrix[sun_group] = projection_matrix
                    logger.debug(f"Projection matrix of {sun_group}-group for field content {term_mass_dim.field_content}:\n{projection_matrix}.")
            except KeyError:
                # Term isn't part of the basis and thus either still contains redundancies or is (FIXME) part of the h.c. part, which isn't given at the moment.
                # In both cases their exist no basis tensor.
                term_mass_dim.sun_projection_tensors = False
    return single_terms

def replace_sun_tensors_by_projected_ones(single_terms):
    """
    Replace sun epsilon tensors by projected basistensors of autoeft.
    Parameters
    ----------
    single_terms

    Returns
    -------

    """
    for type in single_terms.values():
        for term_mass_dim in type.values():
            for sun_group, projection_matrix in term_mass_dim.sun_projection_matrix.items():
                if projection_matrix:
                    # projection_matrix exists
                    sun_basis_tensors = term_mass_dim.sun_projection_tensors[sun_group]["sun_tensor"]
                    basis_dim = len(sun_basis_tensors)
                    assert len(term_mass_dim.terms) == projection_matrix.nrows(), "Projection matrix has the wrong shape."
                    assert basis_dim == projection_matrix.ncols(), "Projection matrix has the wrong shape."
                    for i, term in enumerate(term_mass_dim.terms):
                        try:
                            term.projected_tensors[sun_group] = [Fraction(str(projection_matrix[i][j])) for j in range(basis_dim)]
                        except AttributeError:
                            term.projected_tensors = {sun_group: [Fraction(str(projection_matrix[i][j])) for j in range(basis_dim)]}


    # Combine terms in FORM:
    merged_terms = []
    tensors_for_field_content = {}
    for type in single_terms.values():
        for term_mass_dim in type.values():
            term_with_specific_field_structure = {}
            if term_mass_dim.sun_projection_tensors:
                sun_proj_tensors = {group_name: tensors["sun_tensor"] if tensors else False for group_name, tensors in term_mass_dim.sun_projection_tensors.items()}
                tensors_for_field_content[term_mass_dim.name] = sun_proj_tensors
            else:
                sun_proj_tensors = False
            # Sort terms by explicit field structure
            for term in term_mass_dim:
                term.replace_SUN_indices_by_projection_indices()
                try:
                    term_with_specific_field_structure[term.fieldstructure].append(term)
                except KeyError:
                    term_with_specific_field_structure[term.fieldstructure] = [term]

            for name_of_term, terms_specific in term_with_specific_field_structure.items():
                ref_term = terms_specific[0]
                eq_tensor_indices = []
                ref_tensor_indices = []
                for i, term in enumerate(terms_specific[1:]):
                    if not i:
                        # i == 0
                        ref_tensor_indices_tmp, eq_tensor_indices_tmp = equalize_indices(ref_term, term)
                        ref_tensor_indices.append(ref_tensor_indices_tmp)
                        eq_tensor_indices.append(eq_tensor_indices_tmp)
                    else:
                        ref_tensor_indices_tmp, eq_tensor_indices_tmp = equalize_indices(ref_term, term)
                        ref_tensor_indices.append(ref_tensor_indices_tmp)  # save also those reference tensor indices, because the reference term could have more yukawa matrices and therefore more flavor indices then the term which indices are equalized.
                        eq_tensor_indices.append(eq_tensor_indices_tmp)

                sum_indices = ref_tensor_indices + [index for short_list in eq_tensor_indices for index in short_list]
                n_ref = [len(indices) for indices in ref_tensor_indices]
                n = [len(indices) for indices in eq_tensor_indices]
                n_yukawa_matrices_ref = len(terms_specific[0].tensors["yukawa"])
                n_yukawa_matrices = []
                for i in range(1, len(terms_specific[1:]) + 1):
                    n_yukawa_matrices.append(len(terms_specific[i].tensors["yukawa"]))

                check_flavor_indices = [abs(n_indices - n_ref_indices) == abs(n_yuk_matr - n_yukawa_matrices_ref) for n_indices, n_ref_indices, n_yuk_matr in zip(n, n_ref, n_yukawa_matrices)]

                assert all(check_flavor_indices)  # Note: all([]) is True
                if sun_proj_tensors:
                    exprs = []
                    for term in terms_specific:
                        sun_tensors = {}
                        for sun_group in [group for group in sun_proj_tensors.keys() if sun_proj_tensors[group]]:  # term_mass_dim.sun_projection_matrix.keys():
                            del term.tensors[sun_group]
                            sun_tensors[sun_group] = [f"({str(proj_coeff)})*{sun_proj_tensors[sun_group][i]:abb}" for i, proj_coeff in enumerate(term.projected_tensors[sun_group])]  # if proj_coeff
                        all_tensors = [f"({'+'.join(sui_tensors)})" for sui_tensors in sun_tensors.values()]
                        exprs.append(f"{'*'.join(all_tensors)}*{term:c}")
                    name_form = "".join([f"{name}{nD}" for name, nD in name_of_term])
                    logger.info(f"Substitute SUN projection tensors in {', '.join([term.name for term in terms_specific])} of type {name_form}.")
                    form = "Off statistics;\n"
                    form += "#include declarations_general.h # coefficient\n"
                    form += "#include declarations_general.h # indices\n"
                    form += "#include declarations_general.h # tensors\n"
                    form += "#include declarations_general.h # operators\n"
                    form += "\n"
                    for i, expr in enumerate(exprs):
                        form += f"Local expr{i:d} = {expr};\n"  #
                    form += "\n"
                    list_expr_names = [f"expr{i:d}" for i in range(len(exprs))]
                    form += f"Local expr = {' + '.join(list_expr_names)};\n"
                    form += ".sort\n"
                    form += f"Drop {', '.join(list_expr_names)};\n"
                    form += "\n"
                    form += "* Bring indices of sl2C-epsilons in order:\n"
                    epss = [list(tensor["mathematica"].values())[0] for tensor_name, tensor in
                            op_config["tensors"].items() if
                            "sl2Ceps" in tensor_name]
                    for eps in epss:
                        antisymeps = get_antisymEps(eps)
                        form += f"Multiply replace_({eps},{antisymeps});\n"
                        form += ".sort\n"
                        form += f"Multiply replace_({antisymeps},{eps});\n"
                        form += ".sort\n"
                    form += "\n"
                    form += f"Bracket {', '.join(tensors + bosons + fermions)};\n\n"

                    TERM_PATH = FORM_PATH / name_form

                    form += f'#write <{TERM_PATH / "substituted_T-SUN_tensors.h"}> "%E", expr\n'

                    form += "Print +ss;\n"
                    form += ".end\n"

                    TERM_PATH.mkdir(parents=True, exist_ok=True)  # Create directories if they don't exist.
                    with open(TERM_PATH / f"{name_form}.frm", "w") as file:
                        file.write(form)

                    run_form(fp_cwd=TERM_PATH, filename=f"{name_form}.frm", fp_p=FORM_GENERAL_PATH)

                    with open(TERM_PATH / "substituted_T-SUN_tensors.h", "r") as file:
                        expression = file.read()

                    terms = get_terms(expression, name=name_form)

                    merged_terms.append(terms)
                else:
                    # IF there isn't a projection matrix for this group and operator there - continue.
                    name_form = "".join([f"{name}{nD}" for name, nD in name_of_term])
                    logger.info(f"{terms_specific[0].name} of type {name_form} doesn't has a projection matrix.")
                    # Rename Summand properly by their structure
                    terms_specific[0].name = name_form
                    merged_terms.append(Term(terms_specific, name_form))
    single_terms = get_type(merged_terms)
    for type in single_terms.values():
        for term_mass_dim in type.values():
            del term_mass_dim.sun_projection_matrix
            for type_name, sun_tensors in tensors_for_field_content.items():
                if term_mass_dim.name == type_name:
                    term_mass_dim.sun_projection_tensors = sun_tensors
    return single_terms