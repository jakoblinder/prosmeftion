import logging
import sys

from itertools import permutations
from copy import copy
from fractions import Fraction
from math import factorial
from typing import Iterator, List, Tuple, Dict, Union
from prosmeftion import model, op_config
from .term import Term
from .summand import Summand
from .coefficient import Factor
from .operators import Fields
from prosmeftion.sun_projection import get_type


logger_autoeft = logging.getLogger("autoeft.projection")
logger = logger_autoeft.getChild(__name__)

def permutation_sign(l: list):
    """
    Determines the permutations of sign of the given list of digits:
    For each entry, count how many smaller numbers are to the right of it. This indicates how many adjacent
    transpositions are needed to put the permutation in order.
    This is counted by 'cnt' and thus the permutation sign is given by (-1)**(cnt % 2)
    Parameters
    ----------
    l

    Returns
    -------

    """
    if len(l) <= 1:
        return 1
    assert all([isinstance(i, int) for i in l]), f"Entries have to be integers not {type(l[0])}"
    # check that all entries are different
    n_entry = [l.count(i) for i in l]
    if any([i > 1 for i in n_entry]):
        logger.warning(f"Some numbers appear at least twice in {l}.")
        return
    n = len(l)
    cnt = 0
    for i in range(n):
        for j in range(i + 1, n):
            if (l[i] > l[j]):
                cnt += 1
    return (-1) ** (cnt % 2)

def symmetrize_list(input:Tuple[List], sym_entries:List, antisym:bool=False) -> List[Tuple]:
    """

    Parameters
    ----------
    input
        (coeff, List of only different numbers).
    sym_entries
        Entries in this list which will be symmetrized.
    antisym
        If True, entries will be antisymmetrized.
    Returns
    -------
        [(coeff1, permutation1), (coeff2, permutation2), ...]
    """
    input_coeff = Fraction(input[0])
    input_list  = input[1].copy()
    assert len(input_list) == len(set(input_list)), "Some numbers occur more than ones in input_list."
    assert len(sym_entries) == len(set(sym_entries)), "Some numbers occur more than ones in sym_entries."
    sym_entries = tuple(sorted(sym_entries))  # The entries have to be sorted for an unambiguous permutations sign
    # Find positions of the numbers which should be (anti-)symmetrized:
    pos_sym_entries = []
    for i in sym_entries:
        assert i in input_list, f"{i:d} is not in input_list."
        pos_sym_entries.append(input_list.index(i))

    pos_sym_entries = tuple(pos_sym_entries)
    # Permutations of the positions:
    perms_pos = list(permutations(pos_sym_entries))

    permuted_list = []  # permutations
    for perm in perms_pos:
        new_list = input_list.copy()
        for index_sym, index_perm in enumerate(perm):
            new_list[index_perm] = sym_entries[index_sym]
        new_coeff = input_coeff * Fraction(1,factorial(len(sym_entries)))
        if antisym:
            new_coeff *= Fraction(permutation_sign(perm))
        permuted_list.append((new_coeff, new_list))

    return permuted_list

def symmetrize(input_list:List, sym_structures:List[Tuple[List]]) -> List[Tuple]:
    """

    Parameters
    ----------
    input_list
        List which will be symmetrized.
    sym_structures
        list of the following tuples:
            (sym_entries, 'S') -> Entries specified in sym_entries will be symmetrized.
            (sym_entries, 'A') -> Entries specified in sym_entries will be antisymmetrized.
    Returns
    -------

    """
    # Check that sym_structures do not intersect:
    list_all_sym_entries = [entry for sym_structure in sym_structures for entry in sym_structure[0]]
    assert len(list_all_sym_entries) == len(set(list_all_sym_entries)), "Some symmetrisations have the same entries in common."

    permuted_list = {0:[]}
    input = (Fraction(1), input_list)
    def iterative_symmetrisation(coeff_list, sym_strucs:List[Tuple[List]]):
        if sym_strucs[0][1] == "S":
            permuted = symmetrize_list(coeff_list, sym_strucs[0][0], False)
        else:
            permuted = symmetrize_list(coeff_list, sym_strucs[0][0], True)

        if len(sym_strucs) == 1:
            permuted_list[0] += permuted
        else:
            for new_list in permuted:
                iterative_symmetrisation(new_list, sym_strucs[1:])

    iterative_symmetrisation(input, sym_structures)

    return permuted_list[0]

def shift_fields(term:Summand, field_indices:Tuple[List]):
    """
    A term with e.g. 5 field has field position indices [1,2,3,4,5]. The field_indices declare to which position each
    field will be shifted. For example field_indices = ('-1/2', [2,1,3,4,5])
    would change the position of the first two field and multiply the coefficient of the term by (-1/2). This would
    correspond to the antisymmetrisation of two fermions.

    Parameters
    ----------
    term
    field_indices

    Returns
    -------
        Term with changed order of fields and adjusted coefficient.
    """
    n_fields = sum(term.fieldcounter_stripped.values())
    ref_order = list(range(1, n_fields + 1))
    new_order = field_indices[1]

    new_fields = [None for i in range(n_fields)]
    new_coeff = term.coeff.copy()
    new_coeff *= Factor(field_indices[0])

    for i, field in enumerate(term.fields):
        field_cp = field.copy()
        assert field_cp.field_pos == ref_order[i]
        field_cp.field_pos = new_order[i]
        field_cp.gaugeIndicesforProjection()
        new_fields[new_order[i] - 1] = field_cp

    return Summand(tensors=term.tensors.copy(), fields=Fields(new_fields), coeff=new_coeff, fp_name=term.name)


def rfr(single_terms):
    terms = []
    for term_type in single_terms.values():
        for term_mass_dim in term_type.values():
            if any([n_field > 1 for n_field in term_mass_dim.field_content.values()]):
                # At least one field occurs at least twice:
                logger.info(f"Terms with field content {term_mass_dim.field_content} are (anti-)symmetrized.")
                n_fields = sum(term_mass_dim.field_content.values())  # number of fields in the term
                for term in term_mass_dim.terms:
                    symmetrized_terms = []
                    field_indices = list(range(1, n_fields + 1))
                    symmetrizations = []
                    n = 1
                    for field_name, n_field in term.fieldcounter_stripped.items():
                        if n_field > 1:
                            if model.fields[field_name].ac:
                                sym_antisym = "A"
                            else:
                                sym_antisym = "S"
                            symmetrizations.append((list(range(n, n + n_field)), sym_antisym))
                        n += n_field
                    sym_field_indices = symmetrize(field_indices, symmetrizations)
                    for field_indices in sym_field_indices:
                        symmetrized_terms.append(shift_fields(term, field_indices))


                    for term_sym in symmetrized_terms:
                        name_form = "".join([f"{name}{nD}" for name, nD in term_sym.fieldstructure])
                        terms.append(Term([term_sym], name_form))

            else:
                # no field occurs at least twice
                logger.info(f"For terms with field content {term_mass_dim.field_content} no field occurs at least twice.")
                for term in term_mass_dim.terms:
                    name_form = "".join([f"{name}{nD}" for name, nD in term.fieldstructure])
                    terms.append(Term([term], name_form))

    single_terms = get_type(terms)

    return single_terms

def rfr_sun(single_terms):
    terms = []
    for term_type in single_terms.values():
        for term_mass_dim in term_type.values():
            if any([n_field > 1 for n_field in term_mass_dim.field_content.values()]):
                # At least one field occurs at least twice:
                logger.info(f"Terms with field content {term_mass_dim.field_content} are (anti-)symmetrized.")
                n_fields = sum(term_mass_dim.field_content.values())  # number of fields in the term

                for term in term_mass_dim.terms:
                    symmetrized_terms = []
                    field_indices = list(range(1, n_fields + 1))
                    symmetrizations = []
                    n = 1
                    for field_name, n_field in term.fieldcounter_stripped.items():
                        if field_name not in (op_config["bosonfields"]["H"]["autoeft"]["H"], op_config["bosonfields"]["H"]["autoeft"]["[H+]"]):
                            # Substitution only possible for Higgs fields
                            n += n_field
                            continue
                        elif not all([term.fieldstructure[i][1] == 0 for i in range(n-1, n-1 + n_field)]):
                            #  No derivatives on the corresponding fields.
                            n += n_field
                            continue
                        if n_field > 1:
                            if model.fields[field_name].ac:
                                sym_antisym = "A"
                            else:
                                sym_antisym = "S"
                            symmetrizations.append((list(range(n, n + n_field)), sym_antisym))
                        n += n_field
                    if not symmetrizations:
                        # Fields are not symmetrized due to the occurrence of derivatives.
                        logger.info(
                            f"Indices of terms with field content {term_mass_dim.field_content} are not symmetrized.")
                        for term in term_mass_dim.terms:
                            name_form = "".join([f"{name}{nD}" for name, nD in term.fieldstructure])
                            terms.append(Term([term], name_form))
                        continue
                    sym_field_indices = symmetrize(field_indices, symmetrizations)
                    for field_indices in sym_field_indices:
                        symmetrized_terms.append(shift_fields(term, field_indices))


                    for term_sym in symmetrized_terms:
                        name_form = "".join([f"{name}{nD}" for name, nD in term_sym.fieldstructure])
                        terms.append(Term([term_sym], name_form))

            else:
                # no field occurs at least twice
                logger.info(f"For terms with field content {term_mass_dim.field_content} no field occurs at least twice.")
                for term in term_mass_dim.terms:
                    name_form = "".join([f"{name}{nD}" for name, nD in term.fieldstructure])
                    terms.append(Term([term], name_form))

    single_terms = get_type(terms)

    return single_terms

# def rfr_sun(single_terms):
#     """
#     Symmetrize fields in the same way as rfr(), but keep only the gauge group tensors from the expressions so that
#     indices in gauge group tensors are symmetrized but building block order stays the same.
#
#     Otherwise, the derivative order from the ibp algorithm would be destroyed.
#
#     Parameters
#     ----------
#     single_terms
#
#     Returns
#     -------
#
#     """
#     terms = []
#     for term_type in single_terms.values():
#         for term_mass_dim in term_type.values():
#             if any([n_field > 1 for n_field in term_mass_dim.field_content.values()]):
#                 # At least one field occurs at least twice:
#                 logger.info(f"SUN tensor indices of the terms with field content {term_mass_dim.field_content} are (anti-)symmetrized.")
#                 n_fields = sum(term_mass_dim.field_content.values())  # number of fields in the term
#                 for term in term_mass_dim.terms:
#                     symmetrized_terms = []
#                     field_indices = list(range(1, n_fields + 1))
#                     symmetrizations = []
#                     n = 1
#                     sym_indices = []
#                     for field_name, n_field in term.fieldcounter_stripped.items():
#                         if field_name not in (op_config["bosonfields"]["H"]["autoeft"]["H"], op_config["bosonfields"]["H"]["autoeft"]["[H+]"]):
#                             # Substitution only possible for Higgs fields
#                             n += n_field
#                             continue
#                         elif not all([term.fieldstructure[i][1] == 0 for i in range(n-1, n-1 + n_field)]):
#                             #  No derivatives on the corresponding fields.
#                             n += n_field
#                             continue
#                         if n_field > 1:
#                             if model.fields[field_name].ac:
#                                 sym_antisym = "A"
#                             else:
#                                 sym_antisym = "S"
#                             sym_indices = list(range(n, n + n_field))
#
#                             symmetrizations.append((sym_indices, sym_antisym))
#                         n += n_field
#                     if not symmetrizations:
#                         # Fields are not symmetrized due to the occurrence of derivatives.
#                         logger.info(f"Indices of terms with field content {term_mass_dim.field_content} are not symmetrized.")
#                         for term in term_mass_dim.terms:
#                             name_form = "".join([f"{name}{nD}" for name, nD in term.fieldstructure])
#                             terms.append(Term([term], name_form))
#                         continue
#
#                     sym_field_indices = symmetrize(field_indices, symmetrizations)
#                     for field_indices in sym_field_indices:
#                         symmetrized_terms.append(shift_fields(term, field_indices))
#
#                     term.replace_SUN_indices_by_projection_indices()
#                     for term_sym in symmetrized_terms:
#                         name_form = "".join([f"{name}{nD}" for name, nD in term_sym.fieldstructure])
#                         term_sym.replace_SUN_indices_by_projection_indices()
#                         # Keep only the tensors in this symmetrized form:
#                         term_tensors_sym_only = term_sym.copy()
#                         term_tensors_sym_only.fields = term.fields
#                         # for sun_group in [group for group in term_sym.gaugeTensorsSUN.keys() if term_sym.gaugeTensorsSUN[group]]:
#                         #     term_tensors_sym_only.fields[sun_group] = term.fields[sun_group]
#                         terms.append(Term([term_tensors_sym_only], name_form))
#
#             else:
#                 # no field occurs at least twice
#                 logger.info(f"For terms with field content {term_mass_dim.field_content} no field occurs at least twice.")
#                 for term in term_mass_dim.terms:
#                     name_form = "".join([f"{name}{nD}" for name, nD in term.fieldstructure])
#                     terms.append(Term([term], name_form))
#
#     single_terms = get_type(terms)
#
#     return single_terms

if __name__ == "__main__":
    # a = symmetrize_list((2, [1,2,3,4,5]), [3,1,5], True)
    b = symmetrize([1, 2, 3, 4, 5, 6, 7, 8], [([3, 1, 5], "S"), ([2, 4], "A"), ([7, 8], "S")])
    print(b)