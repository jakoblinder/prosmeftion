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
from .yProjection.read_write import get_terms, get_type
from .yProjection.utils import equalize_indices
from .yProjection.summand import Summand
# from .get_FORM.form_read_in import Term_Model, Term_s
from .yProjection.term import Term, TermType

from autoeft.base.tensors import SUNTableau  # TODO: Check autoeft.invariants -> autoeft.base.tensors
# , field_projection_operator, symmetrize_tensors
# FIXME: field_projection_operator, symmetrize_tensors not defined in new AutoEFT
# from autoeft.model import SUNGroup
# from autoeft.sun_projection import tensor_projection
from . import AUTOEFT_PATH, FORM_PATH, FORM_GENERAL_PATH, model, get_antisymEps, op_config, bosons, fermions, tensors, \
    run_form, get_basis
logger_autoeft = logging.getLogger("autoeft.projection")
logger = logger_autoeft.getChild(__name__)


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
    Create Note Schreibe als gaugeF1I2, ...
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

def get_basis_tensors(basis, field_content, derivatives, mass_dim, model):
    sun_projection_tensors = {}
    # get all the invariants associated with the operator
    op_class, op_subclass, op_type, operator = basis[mass_dim].get_operator(field_content, derivatives)
    for sun_group in model.symmetries.sun_groups.values():
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

def sun_projection(single_terms, basis, max_dim: int, model):
    logger.info(f"Project all terms for a specific type onto the {'-, '.join(list(model.symmetries.sun_groups.keys()))}-group "
                f"basis when the term of specific type is part of the basis.")
    for type in single_terms.values():
        for term_mass_dim in type.values():
            field_content = term_mass_dim.field_content
            derivatives = term_mass_dim.nD
            mass_dim = term_mass_dim.d
            try:
                # get the SUN basis Tensors from autoeft for the projection
                term_mass_dim.sun_projection_tensors = get_basis_tensors(basis, field_content, derivatives, mass_dim, model)
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
                        N = model.symmetries.sun_groups[sun_group].N
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

def replace_sun_tensors_by_projected_ones(single_terms, model):
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

                    terms = get_terms(expression, model, name=name_form)

                    merged_terms.append(terms)
                else:
                    # IF there isn't a projection matrix for this group and operator there - continue.
                    name_form = "".join([f"{name}{nD}" for name, nD in name_of_term])
                    logger.info(f"{terms_specific[0].name} of type {name_form} doesn't has a projection matrix.")
                    # Rename Summand properly by their structure
                    terms_specific[0].name = name_form
                    merged_terms.append(Term(terms_specific, name_form))
    single_terms = get_type(merged_terms, model)
    for type in single_terms.values():
        for term_mass_dim in type.values():
            del term_mass_dim.sun_projection_matrix
            for type_name, sun_tensors in tensors_for_field_content.items():
                if term_mass_dim.name == type_name:
                    term_mass_dim.sun_projection_tensors = sun_tensors
    return single_terms