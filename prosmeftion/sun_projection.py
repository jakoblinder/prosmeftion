import logging
import sys
import re
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
from .yProjection.read_write import get_terms, get_type
from .yProjection.utils import equalize_indices
from .yProjection.summand import Summand
# from .get_FORM.form_read_in import Term_Model, Term_s
from .yProjection.term import Term, TermType

from autoeft.base.tensors import SUNTableau, SUNTensor  # TODO: Check autoeft.invariants -> autoeft.base.tensors
from autoeft.utils import Vector
# , field_projection_operator, symmetrize_tensors
# FIXME: field_projection_operator, symmetrize_tensors not defined in new AutoEFT
# from autoeft.model import SUNGroup
# from autoeft.sun_projection import tensor_projection
from . import AUTOEFT_PATH, FORM_PATH, FORM_GENERAL_PATH, get_antisymEps, op_config, bosons, fermions, tensors, run_form, get_basis
logger_autoeft = logging.getLogger("autoeft.projection")
logger = logger_autoeft.getChild(__name__)


def basis_tensors(sun_inv, sun_reps, sun_group):  # op_type, field_content, tensors, group):
    """Return monomial and symmetrized tensors in FORM-readable output."""
    def convert_in_basisTensors(sym_tab, sun_group, i) -> SymBasisTensors:
        mon_basis_tensors = []
        for tableau, coefficient in sym_tab.items():
            mon_basis_tensors.append(MonBasisTensor(sun_group, tableau, coefficient))

        return SymBasisTensor(mon_basis_tensors, i)

    tensor_basis = []
    for i, tab in enumerate(sun_inv):
        sym_tab = tab.symmetrize_sun_indices(sun_reps)
        sym_basis_tensor = convert_in_basisTensors(sym_tab, sun_group, i)
        # print(f"SymBasisTensor: {sym_basis_tensor:a}")
        # print(f"SymBasisTensor: {sym_basis_tensor:abb}")
        tensor_basis.append(sym_basis_tensor)

    return SymBasisTensors(tensor_basis)


def get_basis_tensors(operator):
    sun_projection_tensors = {}

    for sun, sun_group in operator.model.symmetries.sun_groups.items():
        # get the SU(N) tensor invariants of the sun_group
        sun_inv = [ybasis_element for ybasis_element, sign in operator.invariants[sun]]
        sun_reps = operator.op_type.sun_content(sun_group)

        # convert the tensors to FORM-readable output
        sym_basis_tensors = basis_tensors(sun_inv, sun_reps, sun_group)
        if not sun_inv[0]:
            sun_projection_tensors[sun_group.name] = False
        else:
            sun_projection_tensors[sun_group.name] = {"sun_monom": sun_inv, "sun_tensor": sym_basis_tensors}

    return sun_projection_tensors

def sun_projection(single_terms, basis, max_dim: int):
    model = basis[4].model
    logger.info(f"Project all terms for a specific type onto the {'-, '.join(list(model.symmetries.sun_groups.keys()))}-group "
                f"basis when the term of specific type is part of the basis.")
    for type in single_terms.values():
        for term_mass_dim in type.values():
            field_content = term_mass_dim.field_content
            derivatives = term_mass_dim.nD
            if derivatives:
                field_content["D"] = term_mass_dim.nD
            mass_dim = term_mass_dim.d
            if mass_dim < 4:
                for summand in term_mass_dim:
                    summand.sun_projection = False
                continue
            try:
                # get the SUN basis Tensors from autoeft for the projection
                operator = basis[mass_dim][field_content]
                term_mass_dim.sun_projection_tensors = get_basis_tensors(operator)
            except KeyError:
                # For some reason, no basis tensor exists for this operator type.
                if mass_dim == 4:
                    # FIXME: KeyError '1H_1H+_2D' -> Dim 4 kinetic terms should be implemented in AutoEFT.
                    for summand in term_mass_dim:
                        summand.sun_projection = False
                    logger.warning(
                        f"The operator type {field_content} corresponds to a kinetic term which is not part of the AutoEFT-Basis. The SUN part is, therefore, not projected.")
                    continue
                else:
                    logger.error("Operator should at this place be part of the basis.")
                    sys.exit("STOP")


            for sun, sun_group in model.symmetries.sun_groups.items():
                # Get invariants of the operator
                sun_reps = operator.op_type.sun_content(sun_group)

                for summand in term_mass_dim:
                    eps = f"{summand.gaugeTensorsSUN[sun]:p}"
                    if not eps:
                        continue

                    # express the eps tensors as tableaus
                    tableau = SUNTableau.parse_eps(eps)
                    # build a tensor object out of them
                    tensor = SUNTensor(Vector({tableau: 1}))
                    # FIXME: Expected: tensor = SUNTensor(tableau)

                    # project the (basis) tensor onto the basis
                    projection = tensor.reduction(term_mass_dim.sun_projection_tensors[sun]["sun_monom"], sun_group, sun_reps)
                    try:
                        summand.sun_projection[sun] = projection
                    except (NameError, AttributeError):
                        summand.sun_projection = {sun: projection}

                    # logger.debug(
                    #     f"Projection matrix of {sun}-group for field content {field_content}:\n{projection}.")
    return single_terms

def replace_sun_tensors_by_projected_ones(single_terms, model):
    """
    Replace sun epsilon tensors by projected basistensors of autoeft.
    Parameters
    ----------
    single_terms

    Returns
    -------

    """
    group_name_index = {"SU3_C": "colf", "SU2_W": "gauge"}

    # for type in single_terms.values():
    #     for term_mass_dim in type.values():
    #         for summand in term_mass_dim:
    #             if not summand.sun_projection:
    #                 continue
    #             for sun, proj_tensor in summand.sun_projection.items():


            # for sun_group, projection_mix in term_mass_dim.sun_projection_matrix.items():
            #     if projection_matrix:
            #         # projection_matrix exists
            #         sun_basis_tensors = term_mass_dim.sun_projection_tensors[sun_group]["sun_tensor"]
            #         basis_dim = len(sun_basis_tensors)
            #         assert len(term_mass_dim.terms) == projection_matrix.nrows(), "Projection matrix has the wrong shape."
            #         assert basis_dim == projection_matrix.ncols(), "Projection matrix has the wrong shape."
            #         for i, term in enumerate(term_mass_dim.terms):
            #             try:
            #                 term.projected_tensors[sun_group] = [Fraction(str(projection_matrix[i][j])) for j in range(basis_dim)]
            #             except AttributeError:
            #                 term.projected_tensors = {sun_group: [Fraction(str(projection_matrix[i][j])) for j in range(basis_dim)]}


    # Combine terms in FORM:
    merged_terms = []
    tensors_for_field_content = {}
    for type in single_terms.values():
        for term_mass_dim in type.values():
            tensors_for_field_content[term_mass_dim.name] = term_mass_dim.sun_projection_tensors
            # Sort terms by explicit field structure
            term_with_specific_field_structure = {}
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
                    ref_tensor_indices_tmp, eq_tensor_indices_tmp = equalize_indices(ref_term, term)
                    ref_tensor_indices.append(ref_tensor_indices_tmp)
                    eq_tensor_indices.append(eq_tensor_indices_tmp)

                sum_indices = [index for short_list in ref_tensor_indices for index in short_list] + [index for short_list in eq_tensor_indices for index in short_list]
                # FIXME: Do a summation over the contracted flavour indices?
                n_ref = [len(indices) for indices in ref_tensor_indices]
                n = [len(indices) for indices in eq_tensor_indices]
                n_yukawa_matrices_ref = len(terms_specific[0].tensors["yukawa"])
                n_yukawa_matrices = []
                for i in range(1, len(terms_specific[1:]) + 1):
                    n_yukawa_matrices.append(len(terms_specific[i].tensors["yukawa"]))

                check_flavor_indices = [abs(n_indices - n_ref_indices) == abs(n_yuk_matr - n_yukawa_matrices_ref) for n_indices, n_ref_indices, n_yuk_matr in zip(n, n_ref, n_yukawa_matrices)]

                assert all(check_flavor_indices)  # Note: all([]) is True
                if terms_specific[0].sun_projection:
                    exprs = []
                    for summand in terms_specific:
                        sun_tensors = {}
                        for sun, projection in summand.sun_projection.items():  # term_mass_dim.sun_projection_matrix.keys():
                            # find symmetrised SUN basis tensor:
                            sun_basis = term_mass_dim.sun_projection_tensors[sun]
                            if not sun_basis:
                                continue

                            sun_monom_tensors, sun_sym_tensors = sun_basis.values()

                            tensor_coeff = []
                            for projected_tensor, proj_coeff in projection.items():
                                for sun_monom_tensor, sun_sym_tensor in zip(sun_monom_tensors, sun_sym_tensors):
                                    if sun_monom_tensor == projected_tensor:
                                        tensor_coeff.append([f"{sun_sym_tensor:abb}", str(proj_coeff)])
                                        break

                            assert len(tensor_coeff) == len(projection)


                            del summand.tensors[sun]
                            sun_tensors[sun] = [f"({proj_coeff:s})*{proj_tensor:s}" for proj_tensor, proj_coeff in tensor_coeff]  # if proj_coeff
                        all_tensors = [f"({'+'.join(sui_tensors)})" for sui_tensors in sun_tensors.values()]
                        exprs.append(f"{'*'.join(all_tensors)}*{summand:c}")
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

                    terms = get_terms(expression, model, name=name_form)

                    merged_terms.append(terms)
                else:
                    # IF there isn't a projection matrix for this group and operator there - continue.
                    name_form = "".join([f"{name}{nD}" for name, nD in name_of_term])
                    logger.info(f"{terms_specific[0].name} of type {name_form} doesn't have a projection matrix.")
                    # Rename Summand properly by their structure
                    terms_specific[0].name = name_form
                    merged_terms.append(Term(terms_specific, name_form))
    single_terms = get_type(merged_terms, model)
    for type in single_terms.values():
        for term_mass_dim in type.values():
            term_mass_dim.sun_projection_tensors = tensors_for_field_content[term_mass_dim.name]
    return single_terms