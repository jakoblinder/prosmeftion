import logging
import re
from fractions import Fraction
from pathlib import Path
from typing import List, Dict, Union, Tuple
from copy import copy

import sage.all
import sage.matrix as mx
from sage.rings.rational_field import QQ
from tioc.get_FORM_refactored.term import Term, TermType
from tioc.get_FORM_refactored.summand import Summand
from tioc.get_FORM_refactored.coefficient import Factor
from tioc.get_FORM_refactored.operator import Tensor, Field
from tioc.get_FORM_refactored.operators import Tensors
from tioc.get_FORM_refactored.indices import Indices_Operator, Possible_Indices
from tioc.get_FORM_refactored.index import Index, LP_Index
from tioc.get_FORM_refactored.read_write import get_terms
from tioc.sun_projection import get_type
from tioc.general import create_procedure
from tioc.get_FORM_refactored.tableau import get_op_class, Young_Tableau
from tioc.get_FORM_refactored.lorentz import LR_Tableaux

from autoeft.io import load_basis
from autoeft.sun_projection import tensor_projection
from . import AUTOEFT_PATH, FORM_PATH, FORM_GENERAL_PATH, PROJECTION_PATH,model, get_antisymEps, op_config, bosons, fermions, tensors, bosons_non_conj, fermions_non_conj, run_form, get_SUN_name
from .general import declaration_SL2C_sets

logger_autoeft = logging.getLogger("autoeft.projection")
logger = logger_autoeft.getChild(__name__)

@create_procedure(FORM_GENERAL_PATH)
def form_sortDerivativesIBP():
    """
    Sort derivatives of the TermType with H: 1, H+: 1 and 4 derivatives so that either.
    to contracted derivatives stand on each Higgs field or 3 derivatives (in arbitrary order) stand on H
    and 1 derivative on H+.

    Returns
    -------

    """

    def id_ibp_HHdagger4(lHS: Tuple[Tuple], rHS: Tuple[Tuple], sign: str = "-"):
        """
        Constructs id statements for a special integration by parts:
        The derivatives for the id statements are specified by the number on the SL2C indices of them.
        The epsilons are number increasingly by 1,2,3,4 and the derivatives have to be adjusted accordingly.
        -> The derivatives with indices 1 & 2 and 3 & 4 are contracted.
        For example: ((1,3),(2,4)), ((2,1,3),(4,)) create on the left-hand side
        D(Lsl1?LUsl[i1],Usldot1?ULsldot[j1],D(Lsl3?LUsl[i3],Usldot3?ULsldot[j3],H(?a)))
        *D(Lsl2?LUsl[i2],Usldot2?ULsldot[j2],D(Lsl4?LUsl[i4],Usldot4?ULsldot[j4],[H+](?b)))
        and on the right-hand side
        D(Lsl2,Usldot2,D(Lsl1,Usldot1,D(Lsl3,Usldot3,H(?a))))
        *D(Lsl4,Usldot4,[H+](?b)).

        Parameters
        ----------
        lHS
        rHS

        Returns
        -------

        """
        sl2Ceps = op_config["tensors"]["[sl2Ceps]"]["mathematica"]["sl2Ceps"]  # '[sl2Ceps]'
        cov = op_config["fermionfields"]["D"]["mathematica"]["cov"]  # 'D'
        # Higgs field
        higgs = op_config["bosonfields"]["H"]["mathematica"]["H"]  # 'H'
        higgs_dagger = op_config["bosonfields"]["H"]["mathematica"]["conj[H]"]  # 'H'

        epss = []

        permutations = [(1, 2, 3, 4), (2, 1, 3, 4), (1, 2, 4, 3), (2, 1, 4, 3)]
        for p in permutations:
            for q in permutations:
                epss.append((
                    f"{sl2Ceps}(Usl{p[0]}?ULsl[i{p[0]}],Usl{p[1]}?ULsl[i{p[1]}])*{sl2Ceps}(Usl{p[2]}?ULsl[i{p[2]}],Usl{p[3]}?ULsl[i{p[3]}])"
                    f"*{sl2Ceps}(Lsldot{q[0]}?LUsldot[j{q[0]}],Lsldot{q[1]}?LUsldot[j{q[1]}])*{sl2Ceps}(Lsldot{q[2]}?LUsldot[j{q[2]}],Lsldot{q[3]}?LUsldot[j{q[3]}])",
                    f"{sl2Ceps}(Usl{p[0]},Usl{p[1]})*{sl2Ceps}(Usl{p[2]},Usl{p[3]})*{sl2Ceps}(Lsldot{q[0]},Lsldot{q[1]})*{sl2Ceps}(Lsldot{q[2]},Lsldot{q[3]})"))

        idst = ""

        def add_derivative_LHS(field: str, number: int):
            return f"{cov}(Lsl{number}?LUsl[i{number}],Usldot{number}?ULsldot[j{number}],{field})"

        def add_derivative_RHS(field: str, number: int):
            return f"{cov}(Lsl{number},Usldot{number},{field})"

        def derivatives_LHS(field: str, tuple_1_field: Tuple):
            der = f"{field:s}"
            for i in tuple_1_field[::-1]:
                der = add_derivative_LHS(der, i)
            return der

        def derivatives_RHS(field: str, tuple_1_field: Tuple):
            der = f"{field:s}"
            for i in tuple_1_field[::-1]:
                der = add_derivative_RHS(der, i)
            return der

        for eps in epss:
            idst += "id once "
            idst += f"{eps[0]}"
            idst += f"*{derivatives_LHS(f'{higgs}(?a)', lHS[0])}*{derivatives_LHS(f'{higgs_dagger}(?b)', lHS[1])}"
            idst += " = "
            idst += f"{sign}{eps[1]}"
            idst += f"*{derivatives_RHS(f'{higgs}(?a)', rHS[0])}*{derivatives_RHS(f'{higgs_dagger}(?b)', rHS[1])};\n"

        return idst

    form = ""
    form += "repeat;\n"

    form += "* D1D2H*D1D2H+ = -D1D1D2H*D2H+:\n"
    form += id_ibp_HHdagger4(((1, 3), (2, 4)), ((2, 1, 3), (4,)), sign="-")
    form += "\n"

    form += "* D1D2H*D2D1H+ = -D2D1D2H*D1H+:\n"
    form += id_ibp_HHdagger4(((1, 3), (4, 2)), ((4, 1, 3), (2,)), sign="-")
    form += "\n"

    form += "* D1H*D1D2D2H+ = -D1D1H*D2D2H+:\n"
    form += id_ibp_HHdagger4(((1,), (2, 3, 4)), ((2, 1), (3, 4)), sign="-")
    form += "\n"

    form += "* D1H*D2D1D2H+ = + D1D2D1H*D2H+:\n"
    form += id_ibp_HHdagger4(((1,), (3, 2, 4)), ((2, 3, 1), (4,)), sign="+")
    form += "\n"

    form += "* D1H*D2D2D1H+ = + D2D2D1H*D1H+:\n"
    form += id_ibp_HHdagger4(((1,), (3, 4, 2)), ((4, 3, 1), (2,)), sign="+")
    form += "\n"

    form += "* D1D2D2H*D1H+ = - D2D2H*D1D1H+:\n"
    form += id_ibp_HHdagger4(((1, 3, 4), (2,)), ((3, 4), (1, 2)), sign="-")
    form += "\n"

    form += "* H*D1D1D2D2H+ = + D1D1H*D2D2H+:\n"
    form += id_ibp_HHdagger4(((), (1, 2, 3, 4)), ((2, 1), (3, 4)), sign="+")
    form += "\n"

    form += "* H*D1D2D1D2H+ = - D1D2D1H*D2H+:\n"
    form += id_ibp_HHdagger4(((), (1, 3, 2, 4)), ((2, 3, 1), (4,)), sign="-")
    form += "\n"

    form += "* H*D1D2D2D1H+ = - D2D2D1H*D1H+:\n"
    form += id_ibp_HHdagger4(((), (1, 3, 4, 2)), ((4, 3, 1), (2,)), sign="-")
    form += "\n"

    form += "* D1D1D2D2H*H+ = + D2D2H*D1D1H+:\n"
    form += id_ibp_HHdagger4(((1, 2, 3, 4), ()), ((3, 4), (2, 1)), sign="+")
    form += "\n"

    form += "* D1D2D1D2H*H+ = - D2D1D2H*D1H+:\n"
    form += id_ibp_HHdagger4(((1, 3, 2, 4), ()), ((3, 2, 4), (1,)), sign="-")
    form += "\n"

    form += "* D1D2D2D1H*H+ = - D2D2D1H*D1H+:\n"
    form += id_ibp_HHdagger4(((1, 3, 4, 2), ()), ((3, 4, 2), (1,)), sign="-")
    form += "\n"

    form += "endrepeat;\n"
    return form

def rearrange_derivatives(single_terms):
    """
    In order to avoid that derivatives act on multiple fields and the only place where this occurs for the BS-input is
    the type nD: 4 & {"H": 1, "H+": 1}, the derivatives on these types are rearranged such that they look like
    (D^2 H) * (D^2 H+).

    Parameters
    ----------
    single_terms

    Returns
    -------

    """
    # def contracted_with_field(term: Summand, ref_index: Index) -> int:
    #     if ref_index.typ == "Usldot":
    #         tensors = term.tensors["sldot"]
    #     elif ref_index.typ == "Lsl":
    #         tensors = term.tensors["sl"]
    #     # else:
    #     #     tensors = term.tensors
    #
    #     # Find contracted tensor and save the other index on that tensor:
    #     found_index = False
    #     for tensor in tensors:
    #         for i, index in enumerate(tensor.indices):
    #             if index == ref_index.dual_index:
    #                 found_index = index
    #                 if i == 0:
    #                     other_index = tensor.indices[1]
    #                 else:
    #                     other_index = tensor.indices[0]
    #                 break
    #         if found_index: break
    #
    #     # Find field which is contracted with other_index:
    #     for field in term.fields:
    #         if other_index.dual_index.is_in(field.indices):
    #             pos = field.field_pos
    #             break
    #
    #     return pos
    #
    # def ibp_2_fields(term: Summand, before: int, afterwards: int, derIndices: List[Index]):
    #     """
    #
    #     Parameters
    #     ----------
    #     term
    #     before: int
    #         Field position on which the derivative acts before (1 or 2)
    #     afterwards
    #         Field position on which the derivative should act afterwards (1 or 2).
    #     derIndices
    #         These are the indices of the derivative which are then removed from the one field and attached to the other field.
    #     Returns
    #     -------
    #
    #     """
    #     assert len(term.fields) == 2, "This function is only allowed for a term with exactly 2 fields."
    #     assert before in (1,2) and afterwards in (1,2)
    #     assert before != afterwards
    #     assert len(derIndices) == 2, "A derivative has to have 2 indices."
    #
    #     for derIndex in derIndices:
    #         for i, field_index in enumerate(term.fields[before-1].indices):
    #             if derIndex == field_index:
    #                 del term.fields[before-1].indices[i]
    #
    #     # Recalculate the number of derivatives 'derIndex', which stand on every index.
    #     term.fields[before-1].reset_derIndex()
    #
    #     # Before appending the indices to the other field, adjust the 'derIndex' attribute of those indices:
    #     for der_index in derIndices:
    #         der_index.derIndex = term.fields[afterwards-1].nD + 1
    #
    #     # Append indices
    #     for derIndex in derIndices:
    #         term.fields[afterwards - 1].indices.append(derIndex)
    #
    #     # Multiply coefficient of the term by (-1), since we have done essentially a partial integration.
    #     term.coeff *= Factor("-1")
    #
    #     return term
    #
    # def balance_derivatives(term: Summand):
    #     """
    #     (Recursive function)
    #     Go through each derivative and check if it is contracted with a derivative on the other field.
    #     If yes -> ibp, i.e. deleted indices of the derivative on the one field, but store them (!), in order to
    #     insert the indices on the other field and multiply the coefficient by (-1).
    #     Note: After deleting the indices from the one field, the derIndex attribute of the remaining derivative
    #     indices on this field have to be set newly. When appending the indices to the other field derIndex will
    #     become just nD + 1.
    #     Parameters
    #     ----------
    #     term
    #
    #     Returns
    #     -------
    #
    #     """
    #     derFields = {}
    #     for field in term.fields:
    #         derIndices = {der: [] for der in range(1, field.nD + 1)}
    #         for index in field.indices["sl"]:
    #             if index.derIndex:
    #                 derIndices[index.derIndex].append(index)
    #         for index in field.indices["sldot"]:
    #             if index.derIndex:
    #                 derIndices[index.derIndex].append(index)
    #
    #         derIndices_contracted_crosswise = {
    #             der: all([field.field_pos != contracted_with_field(term, derIndices[der][i]) for i in range(2)]) for der
    #             in range(1, field.nD + 1)}
    #         derFields[field.field_pos] = (derIndices, derIndices_contracted_crosswise)
    #     assert len(derFields) == 2
    #     # Termination condition:
    #     if all([all(not crossed for crossed in ind_cross[1].values()) for ind_cross in derFields.values()]) and all([field.nD <= 2 for field in term.fields]):
    #         return term
    #     else:
    #         if len(derFields[1][1]) == 2 and len(derFields[2][1]) == 2:
    #             # 2 derivatives on both fields, both with crossed indices
    #             # Move the first derivative of the first field
    #             term = ibp_2_fields(term, before=1, afterwards=2, derIndices=derFields[1][0][1])
    #             return balance_derivatives(term)
    #         elif len(derFields[1][1]) == 4 and len(derFields[2][1]) == 0:
    #             # 4 derivatives on the first field.
    #             # Move the first derivative of the first field
    #             term = ibp_2_fields(term, before=1, afterwards=2, derIndices=derFields[1][0][1])
    #             return balance_derivatives(term)
    #         elif len(derFields[1][1]) == 0 and len(derFields[2][1]) == 4:
    #             # 4 derivatives on the first field.
    #             # Move the first derivative of the first field
    #             term = ibp_2_fields(term, before=2, afterwards=1, derIndices=derFields[2][0][1])
    #             return balance_derivatives(term)
    #         else:
    #             # 3 derivatives on the first field and 1 on the second or 3 derivatives on the second field and 1 on the first.
    #             # Move the one derivative from field which has 3 derivatives, where the indices are contracted cross wise.
    #             nD1 = len(derFields[1][1])  #  Number of derivatives acting on the first field
    #             nD2 = len(derFields[2][1])  #  Number of derivatives acting on the second field
    #             assert nD1 != nD2
    #             if nD1 > nD2:
    #                 der = list(derFields[1][1].values()).index(True) + 1
    #                 term = ibp_2_fields(term, before=1, afterwards=2, derIndices=derFields[1][0][der])
    #             elif nD1 < nD2:
    #                 der = list(derFields[2][1].values()).index(True) + 1
    #                 term = ibp_2_fields(term, before=2, afterwards=1, derIndices=derFields[2][0][der])
    #             return balance_derivatives(term)

    try:
        nD4_h2h_dagger2 = single_terms[(0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0)][4]
        logger.info(f"Commutate derivatives in nD=4 and H:1, H+:1.")
    except KeyError:
        logger.info("There is no term with nD=4 and H:1, H+:1 and thus no derivatives need to be rearranged.")
        return single_terms

    # for term in nD4_h2h_dagger2.terms:
    #     term = balance_derivatives(term)
    #     term.name = "".join([f"{name}{nD}" for name, nD in term.fieldstructure])

    # Write every operator as D1D1D2H*D1Hdagger or D1D2D1H*D2Hdagger
    expression = ""
    indices = Possible_Indices([])
    for term in nD4_h2h_dagger2:
        expression += f"+{term:c}\n"
        indices += term.possible_indices

    folder = f"4DHHdagger"
    TERM_PATH = FORM_PATH / f"{folder}"
    TERM_PATH.mkdir(parents=True, exist_ok=True)  # Create directories if they don't exist.

    form_sortDerivativesIBP()

    # Write SL2C and set FORM-file:
    form_SL2C = declaration_SL2C_sets(indices)
    with open(TERM_PATH / "declaration_SL2C.h", "w") as file:
        file.write(form_SL2C)

    form = "Off statistics;\n"
    form += "#include declarations_general.h # coefficient\n"
    form += "#include declarations_general.h # indices\n"
    form += "#include declarations_general.h # tensors\n"
    form += "#include declarations_general.h # operators\n"
    form += "#include declaration_SL2C.h\n"
    form += "\n"

    form += f"Local expr = {expression};\n"
    form += ".sort\n"
    form += "\n"

    form += "#call sortDerivativesIBP"

    # Symplify expressions
    form += "\n"
    form += "* Symplify expression\n"
    form += "#call simplifySL2CEps\n"
    form += "#call simplifyEpsSU2\n"
    form += "#call simplifyEpsSU3\n"

    form += ".sort\n"
    form += "#call simplifySymmetricFieldstrengthtensor\n"

    form += ".sort\n"
    form += "\n"
    form += "* Write derivatives implicit with indices inside of fields.\n"
    form += "#call derivativeasIndex\n"
    form += "* Order fields by their helicity and then alpabetically.\n"
    form += "#call sortfields\n"
    form += "* Write derivatives again outside of fields.\n"
    form += "#call indexasDerivative\n"

    form += "\n"
    form += ".sort\n"
    form += f"Bracket {', '.join(tensors + bosons + fermions)}, D2, EOM;\n\n"
    form += f'#write <{TERM_PATH / "higgs_sorted.h"}> "%E", expr\n'

    form += "Print +ss;\n"
    form += ".end\n"

    with open(TERM_PATH / f"higgs_sort_derivatives.frm", "w") as file:
        file.write(form)

    run_form(fp_cwd=TERM_PATH, filename=f"higgs_sort_derivatives.frm", fp_p=FORM_GENERAL_PATH)

    terms = get_terms(TERM_PATH / "higgs_sorted.h", as_one=False)

    terms = [summand for term in terms for summand in term.terms]
    h3Hdagger1 = TermType(terms,terms[0].fieldcounter_stripped)

    assert all([(term.fieldstructure[0][1], term.fieldstructure[1][1]) in ((3,1), (2,2)) for term in h3Hdagger1])

    sl2Ceps = op_config["tensors"]["[sl2Ceps]"]["mathematica"]["sl2Ceps"]  # '[sl2Ceps]'
    su2eps = op_config["tensors"]["[su2eps]"]["mathematica"]["su2eps"]
    cov = op_config["fermionfields"]["D"]["mathematica"]["cov"]  # 'D'
    # Higgs field
    higgs = op_config["bosonfields"]["H"]["mathematica"]["H"]  # 'H'
    higgs_dagger = op_config["bosonfields"]["H"]["mathematica"]["conj[H]"]  # 'H'
    # field strength tensors
    bl = op_config["bosonfields"]["BL"]["mathematica"]["BL"]
    br = op_config["bosonfields"]["BL"]["mathematica"]["conj[BL]"]
    wl = op_config["bosonfields"]["WL"]["mathematica"]["WL"]
    wr = op_config["bosonfields"]["WL"]["mathematica"]["conj[WL]"]
    # Coefficients
    g1 = op_config["coefficients"]["g1"]["mathematica"]["g1"]
    g2 = op_config["coefficients"]["g2"]["mathematica"]["g2"]
    # imaginary unit
    mi = op_config["coefficients"]["I"]["mathematica"]["I"]

    term_list = []  # flat list of Summands, which is later sorted by their types
    for summand in h3Hdagger1:
        @create_procedure(TERM_PATH)
        def form_commuteDerivatives():
            """
            Commute derivatives.
            """
            gen_index_gauge = summand.possible_indices.generate_index("gauge")  # get a new, i.e. unused SU2 index with 'next(gen_index_gauge)'
            # gen_index_sl = summand.possible_indices.generate_index("Lsl")  # next(gen_index_sl) will generate new Lsl and new Usl index
            # gen_index_sldot = summand.possible_indices.generate_index("Lsldot")

            gauge3 = next(gen_index_gauge)
            gauge4 = next(gen_index_gauge)

            # lsl1, usl1 = next(gen_index_sl)
            # lsl2, usl2 = next(gen_index_sl)
            # lsldot1, usldot1 = next(gen_index_sldot)
            # lsldot2, usldot2 = next(gen_index_sldot)

            epss = []  # contains epsilons of the LHS and the corresponding sign of the RHS corresponding to permutations of indices.

            permutations = [((1, 2, 3, 4), +1), ((2, 1, 3, 4), -1), ((1, 2, 4, 3), -1), ((2, 1, 4, 3), +1)]
            for p in permutations:
                for q in permutations:
                    epss.append((
                        f"{sl2Ceps}(Usl{p[0][0]}?ULsl[i{p[0][0]}],Usl{p[0][1]}?ULsl[i{p[0][1]}])*{sl2Ceps}(Usl{p[0][2]}?ULsl[i{p[0][2]}],Usl{p[0][3]}?ULsl[i{p[0][3]}])"
                        f"*{sl2Ceps}(Lsldot{q[0][0]}?LUsldot[j{q[0][0]}],Lsldot{q[0][1]}?LUsldot[j{q[0][1]}])*{sl2Ceps}(Lsldot{q[0][2]}?LUsldot[j{q[0][2]}],Lsldot{q[0][3]}?LUsldot[j{q[0][3]}])",
                        p[1] * q[1]))

            form = ""
            form += "repeat;\n"

            form += "* D1D2D1H*D2H+ = ...:\n"
            for eps in epss:
                form += "id once "
                form += f"{eps[0]}"
                form += f"*{cov}(Lsl1?LUsl[i1],Usldot1?ULsldot[j1],{cov}(Lsl3?LUsl[i3],Usldot3?ULsldot[j3],{cov}(Lsl2?LUsl[i2],Usldot2?ULsldot[j2],{higgs}(gauge1?gauge))))"
                form += f"*{cov}(Lsl4?LUsl[i4],Usldot4?ULsldot[j4],{higgs_dagger}(gauge2?gauge))"
                form += f" = ({eps[1]}) * ("
                form += f"-{sl2Ceps}(Usl1,Usl2)*{sl2Ceps}(Usl3,Usl4)*{sl2Ceps}(Lsldot1,Lsldot2)*{sl2Ceps}(Lsldot3,Lsldot4)"
                form += f"*{cov}(Lsl1,Usldot1,{cov}(Lsl2,Usldot2,{higgs}(gauge1)))*{cov}(Lsl3,Usldot3,{cov}(Lsl4,Usldot4,{higgs_dagger}(gauge2)))"

                form += f" +({mi}*{g1}/2)*{sl2Ceps}(Usl1,Usl2)*{sl2Ceps}(Usl3,Usl4)*{sl2Ceps}(Lsldot4,Lsldot2)"
                form += f"*{bl}(Lsl1,Lsl3)*{cov}(Lsl2,Usldot2,{higgs}(gauge1))*{cov}(Lsl4,Usldot4,{higgs_dagger}(gauge2))"
                form += f" +({mi}*{g2}/2)*{sl2Ceps}(Usl1,Usl2)*{sl2Ceps}(Usl3,Usl4)*{sl2Ceps}(Lsldot4,Lsldot2)"
                form += f"*{wl}(gauge1,{gauge3},Lsl1,Lsl3)*{su2eps}({gauge3},{gauge4})*{cov}(Lsl2,Usldot2,{higgs}({gauge4}))*{cov}(Lsl4,Usldot4,{higgs_dagger}(gauge2))"

                form += f" -({mi}*{g1}/2)*{sl2Ceps}(Usl4,Usl2)*{sl2Ceps}(Lsldot1,Lsldot2)*{sl2Ceps}(Lsldot3,Lsldot4)"
                form += f"*{cov}(Lsl2,Usldot2,{higgs}(gauge1))*{cov}(Lsl4,Usldot4,{higgs_dagger}(gauge2))*{br}(Usldot1,Usldot3)"
                form += f" -({mi}*{g2}/2)*{sl2Ceps}(Usl4,Usl2)*{sl2Ceps}(Lsldot1,Lsldot2)*{sl2Ceps}(Lsldot3,Lsldot4)"
                form += f"*{cov}(Lsl2,Usldot2,{higgs}({gauge4}))*{cov}(Lsl4,Usldot4,{higgs_dagger}(gauge2))*{wr}(gauge1,{gauge3},Usldot1,Usldot3)*{su2eps}({gauge3},{gauge4})"

                form += ");\n"

            form += "* D1D1D2H*D2H+ = ...:\n"
            for eps in epss:
                form += "id once "
                form += f"{eps[0]}"
                form += f"*{cov}(Lsl1?LUsl[i1],Usldot1?ULsldot[j1],{cov}(Lsl2?LUsl[i2],Usldot2?ULsldot[j2],{cov}(Lsl3?LUsl[i3],Usldot3?ULsldot[j3],{higgs}(gauge1?gauge))))"
                form += f"*{cov}(Lsl4?LUsl[i4],Usldot4?ULsldot[j4],{higgs_dagger}(gauge2?gauge))"
                form += f" = ({eps[1]}) * ("
                form += f"-{sl2Ceps}(Usl1,Usl2)*{sl2Ceps}(Usl3,Usl4)*{sl2Ceps}(Lsldot1,Lsldot2)*{sl2Ceps}(Lsldot3,Lsldot4)"
                form += f"*{cov}(Lsl1,Usldot1,{cov}(Lsl2,Usldot2,{higgs}(gauge1)))*{cov}(Lsl3,Usldot3,{cov}(Lsl4,Usldot4,{higgs_dagger}(gauge2)))"

                form += f" +({mi}*{g1}/2)*{sl2Ceps}(Usl1,Usl2)*{sl2Ceps}(Usl3,Usl4)*{sl2Ceps}(Lsldot4,Lsldot2)"
                form += f"*{bl}(Lsl1,Lsl3)*{cov}(Lsl2,Usldot2,{higgs}(gauge1))*{cov}(Lsl4,Usldot4,{higgs_dagger}(gauge2))"
                form += f" +({mi}*{g2}/2)*{sl2Ceps}(Usl1,Usl2)*{sl2Ceps}(Usl3,Usl4)*{sl2Ceps}(Lsldot4,Lsldot2)"
                form += f"*{wl}(gauge1,{gauge3},Lsl1,Lsl3)*{su2eps}({gauge3},{gauge4})*{cov}(Lsl2,Usldot2,{higgs}({gauge4}))*{cov}(Lsl4,Usldot4,{higgs_dagger}(gauge2))"

                form += f" -({mi}*{g1}/2)*{sl2Ceps}(Usl4,Usl2)*{sl2Ceps}(Lsldot1,Lsldot2)*{sl2Ceps}(Lsldot3,Lsldot4)"
                form += f"*{cov}(Lsl2,Usldot2,{higgs}(gauge1))*{cov}(Lsl4,Usldot4,{higgs_dagger}(gauge2))*{br}(Usldot1,Usldot3)"
                form += f" -({mi}*{g2}/2)*{sl2Ceps}(Usl4,Usl2)*{sl2Ceps}(Lsldot1,Lsldot2)*{sl2Ceps}(Lsldot3,Lsldot4)"
                form += f"*{cov}(Lsl2,Usldot2,{higgs}({gauge4}))*{cov}(Lsl4,Usldot4,{higgs_dagger}(gauge2))*{wr}(gauge1,{gauge3},Usldot1,Usldot3)*{su2eps}({gauge3},{gauge4})"

                form += f" +({mi}*{g1}/2)*{sl2Ceps}(Usl1,Usl2)*{sl2Ceps}(Usl3,Usl4)*{sl2Ceps}(Lsldot1,Lsldot4)"
                form += f"*{cov}(Lsl1,Usldot1,{bl}(Lsl2,Lsl3))*{higgs}(gauge1)*{cov}(Lsl4,Usldot4,{higgs_dagger}(gauge2))"
                form += f" +({mi}*{g2}/2)*{sl2Ceps}(Usl1,Usl2)*{sl2Ceps}(Usl3,Usl4)*{sl2Ceps}(Lsldot1,Lsldot4)"
                form += f"*{cov}(Lsl1,Usldot1,{wl}(gauge1,{gauge3},Lsl2,Lsl3))*{su2eps}({gauge3},{gauge4})*{higgs}({gauge4})*{cov}(Lsl4,Usldot4,{higgs_dagger}(gauge2))"

                form += f" +({mi}*{g1}/2)*{sl2Ceps}(Usl1,Usl2)*{sl2Ceps}(Usl3,Usl4)*{sl2Ceps}(Lsldot1,Lsldot4)"
                form += f"*{bl}(Lsl2,Lsl3)*{cov}(Lsl1,Usldot1,{higgs}(gauge1))*{cov}(Lsl4,Usldot4,{higgs_dagger}(gauge2))"
                form += f" +({mi}*{g2}/2)*{sl2Ceps}(Usl1,Usl2)*{sl2Ceps}(Usl3,Usl4)*{sl2Ceps}(Lsldot1,Lsldot4)"
                form += f"*{wl}(gauge1,{gauge3},Lsl2,Lsl3)*{su2eps}({gauge3},{gauge4})*{cov}(Lsl1,Usldot1,{higgs}({gauge4}))*{cov}(Lsl4,Usldot4,{higgs_dagger}(gauge2))"

                form += f" -({mi}*{g1}/2)*{sl2Ceps}(Usl1,Usl4)*{sl2Ceps}(Lsldot1,Lsldot2)*{sl2Ceps}(Lsldot3,Lsldot4)"
                form += f"*{higgs}(gauge1)*{cov}(Lsl4,Usldot4,{higgs_dagger}(gauge2))*{cov}(Lsl1,Usldot1,{br}(Usldot2,Usldot3))"
                form += f" -({mi}*{g2}/2)*{sl2Ceps}(Usl1,Usl4)*{sl2Ceps}(Lsldot1,Lsldot2)*{sl2Ceps}(Lsldot3,Lsldot4)"
                form += f"*{higgs}({gauge4})*{cov}(Lsl4,Usldot4,{higgs_dagger}(gauge2))*{cov}(Lsl1,Usldot1,{wr}(gauge1,{gauge3},Usldot2,Usldot3))*{su2eps}({gauge3},{gauge4})"

                form += f" -({mi}*{g1}/2)*{sl2Ceps}(Usl1,Usl4)*{sl2Ceps}(Lsldot1,Lsldot2)*{sl2Ceps}(Lsldot3,Lsldot4)"
                form += f"*{cov}(Lsl1,Usldot1,{higgs}(gauge1))*{cov}(Lsl4,Usldot4,{higgs_dagger}(gauge2))*{br}(Usldot2,Usldot3)"
                form += f" -({mi}*{g2}/2)*{sl2Ceps}(Usl1,Usl4)*{sl2Ceps}(Lsldot1,Lsldot2)*{sl2Ceps}(Lsldot3,Lsldot4)"
                form += f"*{cov}(Lsl1,Usldot1,{higgs}({gauge4}))*{cov}(Lsl4,Usldot4,{higgs_dagger}(gauge2))*{wr}(gauge1,{gauge3},Usldot2,Usldot3)*{su2eps}({gauge3},{gauge4})"

                form += ");\n"
            form += "endrepeat;\n"
            return form

        form_commuteDerivatives()

        summand.name = "".join([f"{name}{nD}" for name, nD in summand.fieldstructure])

        # Write SL2C and set FORM-file:
        form_SL2C = declaration_SL2C_sets(summand.possible_indices)
        with open(TERM_PATH / "declaration_SL2C.h", "w") as file:
            file.write(form_SL2C)

        form = "Off statistics;\n"
        form += "#include declarations_general.h # coefficient\n"
        form += "#include declarations_general.h # indices\n"
        form += "#include declarations_general.h # tensors\n"
        form += "#include declarations_general.h # operators\n"
        form += "#include declaration_SL2C.h\n"
        form += "\n"

        form += f"Local expr = {summand:c};\n"
        form += ".sort\n"
        form += "\n"

        form += "#call commuteDerivatives\n"
        form += ".sort\n"

        # Symplify expressions
        form += "\n"
        form += "* Symplify expression\n"
        form += "#call simplifySL2CEps\n"
        form += "#call simplifyEpsSU2\n"
        form += "#call simplifyEpsSU3\n"

        form += ".sort\n"
        form += "#call simplifySymmetricFieldstrengthtensor\n"

        form += ".sort\n"
        form += "\n"
        form += "* Write derivatives implicit with indices inside of fields.\n"
        form += "#call derivativeasIndex\n"
        form += "* Order fields by their helicity and then alpabetically.\n"
        form += "#call sortfields\n"
        form += "* Write derivatives again outside of fields.\n"
        form += "#call indexasDerivative\n"

        form += "\n"
        form += ".sort\n"
        form += f"Bracket {', '.join(tensors + bosons + fermions)}, D2, EOM;\n\n"

        form += f'#write <{TERM_PATH / "commute_derivatives.h"}> "%E", expr\n'

        form += "Print +ss;\n"
        form += ".end\n"

        with open(TERM_PATH / f"higgs_commute_derivatives.frm", "w") as file:
            file.write(form)

        run_form(fp_cwd=TERM_PATH, filename=f"higgs_commute_derivatives.frm", fp_p=FORM_GENERAL_PATH)

        terms = get_terms(TERM_PATH / "commute_derivatives.h", as_one=False)
        term_list.append(terms)

    term_list = [term for short_list in term_list for term in short_list]

    for name, type in single_terms.items():
        for i, term_mass_dim in type.items():
            if name == (0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0) and i == 4:
                # nD4_h2h_dagger2 term
                continue
            else:
                for summand in term_mass_dim:
                    term_list.append(Term([summand], name=summand.name))

    single_terms = get_type(term_list)

    return single_terms

def replace_eoms(single_terms):
    """
    Replace equation of motions first by a placeholder and then for the dimension 6 operators by the equation of motion
    of the SM-Lagrangian.
    Parameters
    ----------
    single_terms

    Returns
    -------

    """
    term_list = []  # flat list of Summands, which is later sorted by their types
    for type in single_terms.values():
        for term_mass_dim in type.values():
            if term_mass_dim.d != 6:
                term_list.append([Term(term_mass_dim.terms, name=term_mass_dim.name)])
                continue
            for summand in term_mass_dim:
                summand.name = "".join([f"{name}{nD}" for name, nD in summand.fieldstructure])
                check_derivatives_higgs = [field.nD for field in summand.fields["H"]]

                q = summand.fields[f"{op_config['fermionfields']['Q']['mathematica']['Q']}"]
                l = summand.fields[f"{op_config['fermionfields']['L']['mathematica']['L']}"]
                u = summand.fields[f"{op_config['fermionfields']['[u_C]']['mathematica']['uC']}"]
                d = summand.fields[f"{op_config['fermionfields']['[d_C]']['mathematica']['dC']}"]
                e = summand.fields[f"{op_config['fermionfields']['[e_C]']['mathematica']['eC']}"]
                bl = summand.fields[f"{op_config['bosonfields']['BL']['mathematica']['BL']}"]
                wl = summand.fields[f"{op_config['bosonfields']['WL']['mathematica']['WL']}"]

                fields_max_1_derivative = q + l + u + d + e + bl + wl
                check_derivatives_max_1_derivative = [field.nD for field in fields_max_1_derivative]

                # FIXME: For the moment, don't consider Higgs fields with more than 2 derivatives.
                if any([nD > 2 for nD in check_derivatives_higgs]):
                    term_list.append([Term([summand], name=summand.name)])
                    continue

                # FIXME: For the moment, don't consider fermions and field strength tensors with more than 1 derivatives.
                if any([nD > 1 for nD in check_derivatives_max_1_derivative]):
                    term_list.append([Term([summand], name=summand.name)])
                    continue

                # If all fermions and field strength tensors don't have any derivatives and the higgs fields have only one,
                # then there cannot occur any EOMs of the SM
                if all([nD < 2 for nD in check_derivatives_higgs]) and all([nD < 1 for nD in check_derivatives_max_1_derivative]):
                    term_list.append([Term([summand], name=summand.name)])
                    continue

                name_form = summand.name  # "".join([f"{name}{nD}" for name, nD in name_of_term])
                TERM_PATH = FORM_PATH / name_form
                TERM_PATH.mkdir(parents=True, exist_ok=True)  # Create directories if they don't exist.

                logger.info(f"Find EOMs in term {name_form}.")
                # The function which subsitute the EOMs should already be run here, so that all necessary indices are then declared afterwards.
                summand.form_higgsEOM(TERM_PATH)
                summand.form_fieldstrengthtensor(TERM_PATH)
                summand.form_fermionEOM(TERM_PATH)
                summand.form_fieldstrengthtensorEOM(TERM_PATH)

                # Write SL2C and set FORM-file:
                form_SL2C = declaration_SL2C_sets(summand.possible_indices)
                with open(TERM_PATH / "declaration_SL2C.h", "w") as file:
                    file.write(form_SL2C)

                form = "Off statistics;\n"
                form += "#include declarations_general.h # coefficient\n"
                form += "#include declarations_general.h # indices\n"
                form += "#include declarations_general.h # tensors\n"
                form += "#include declarations_general.h # operators\n"
                form += "#include declaration_SL2C.h\n"
                form += "\n"

                form += f"Local expr = {summand:c};\n"
                form += ".sort\n"
                form += "\n"

                # EOM identifications:
                form += "#call antisymDerivative\n"
                form += "#call spinorEOMidentification\n"
                form += "#call fieldstrengthtensorEOMidentification\n"
                form += "label 2;\n"
                # EOM substitutions:
                form += "#call higgsEOM\n"
                form += "label 3;\n"
                form += "#call fieldstrengthtensor\n"
                form += "#call fermionEOM\n"
                form += "#call fieldstrengthtensorEOM\n"

                form += "label 4;\n"

                # Symplify expressions
                form += "\n"
                form += "* Symplify expression\n"
                form += "#call simplifySL2CEps\n"
                form += "#call simplifyEpsSU2\n"
                form += "#call simplifyEpsSU3\n"

                form += ".sort\n"
                form += "#call simplifySymmetricFieldstrengthtensor\n"

                form += ".sort\n"
                form += "\n"
                form += "* Write derivatives implicit with indices inside of fields.\n"
                form += "#call derivativeasIndex\n"
                form += "* Order fields by their helicity and then alpabetically.\n"
                form += "#call sortfields\n"
                form += "* Write derivatives again outside of fields.\n"
                form += "#call indexasDerivative\n"

                form += "\n"
                form += ".sort\n"
                form += f"Bracket {', '.join(tensors + bosons + fermions)}, D2, EOM;\n\n"
                form += f'#write <{TERM_PATH / "term_with_less_eom.h"}> "%E", expr\n'

                form += "Print +ss;\n"
                form += ".end\n"

                with open(TERM_PATH / f"{name_form}eom.frm", "w") as file:
                    file.write(form)

                run_form(fp_cwd=TERM_PATH, filename=f"{name_form}eom.frm", fp_p=FORM_GENERAL_PATH)

                terms = get_terms(TERM_PATH / "term_with_less_eom.h", as_one=False)
                term_list.append(terms)
                # op = op_config
                print("---------------")

    term_list = [term for short_list in term_list for term in short_list]
    single_terms = get_type(term_list)

    return single_terms

def ibp_and_schouten_ids(single_terms, max_dim: int):
    """
    Apply the integration by parts and Schouten identities to the lorentz structure.
    Parameters
    ----------
    single_terms
    max_dim

    Returns
    -------

    """
    term_check = []  # This checks whether all terms are already Young tableaux or not.
    term_list  = []  # flat list of Summands, which is later sorted by their types
    for type in single_terms.values():
        for term_mass_dim in type.values():
            for summand in term_mass_dim:
                skip_schouten_ids = False
                skip_ibp_ids = False
                if not summand.nD: skip_ibp_ids = True  # No ibp necessary if no derivative is there. -> Schouten id are still necessary!

                lr = summand.lr
                if lr.l_tab.ncols() < 2 and lr.r_tab.ncols() < 2: skip_schouten_ids = True  # No Schouten ids can be applied when their are less then 2 epsilon tensors
                if skip_schouten_ids and skip_ibp_ids:
                    # Neither the Schouten nor the ibp relations need to be applied.
                    term_list.append(summand)
                    continue

                # print(f"{lr.r_tab:fp}")
                # print(f"{lr.l_tab:fp}")

                print(f"{summand:c}")

                tabs = summand.lr.remove_ibp(summand.op_class.N)

                tabs_without_Schouten = []
                if not skip_schouten_ids:
                    # Schouten identities may occur
                    for lr_with_Schouten in tabs:
                        if lr_with_Schouten[1]:
                            # EOMs might still occur
                            tabs_without_Schouten.append(lr_with_Schouten)
                        else:
                            tabs_without_Schouten += [(lr_without_Schouten, False) for lr_without_Schouten in lr_with_Schouten[0].remove_schouten()]
                else:
                    tabs_without_Schouten = tabs

                for new_lr in tabs_without_Schouten:
                    print(f"{new_lr[0]:lp}")
                    term = summand.get_term_from_lr_tabs(new_lr[0])
                    print(f"{term:c}")
                    term_list.append(term)

                for new_lr in tabs_without_Schouten:
                    if summand.d != 4 and summand.d == max_dim:  # FIXME: The part 'summand.d == max_dim' can be removed when EOMs are removed via fieldrefefinitions.
                        if new_lr[1]:
                            # EOMs might still occur
                            term_check.append(False)
                            print("False: EOM")
                        elif lTs := new_lr[0].from_lr_tableaux(summand.op_class):
                            lorentz_tab, sign = lTs[0], lTs[1]
                            term_check.append(lorentz_tab.is_ssyt)
                            if not lorentz_tab.is_ssyt:
                                print("False: No SSYT")
                        else:
                            print("False: else")
                            term_check.append(False)
                    # assert lorentz_tab.is_ssyt, f"At this point this should be a SSYT, but it is:\n{lorentz_tab:nice}"
                print("=======================================")
                # # Change the derivative structure by only changing the tableau indices (no sign change):
                # spec_derivative = object()
                # lr_copy = lr.copy()
                # del lr
                #
                # lr_copy.l_tab[0,0].lp = LP_Index(1, spec_derivative)
                # lr_copy.l_tab[0, 0].derIndex = spec_derivative  # 1
                # lr_copy.r_tab[0, 0].lp = LP_Index(1, spec_derivative)
                # lr_copy.r_tab[0, 0].derIndex = spec_derivative  # 1
                #
                # term = summand.get_term_from_lr_tabs(lr_copy)
                # # Infer from the l_tab and r_tab the derivative structure for a given field structure:
                # print("------------------------------------------")
                # print(f"{term:c}")
                # print(f"{term.lr:lp}")
                # # Integration by parts of derivative on first field
                # tabs1 = term.lr.ibp(LP_Index(1,1), summand.op_class.N)
                # # Integration by parts of derivative on second field
                # tabs2 = term.lr.ibp(LP_Index(2,1), summand.op_class.N)
                #
                # #print both terms
                # print("Integration by parts of derivative on first field:")
                # for lr in tabs1:
                #     print(f"{lr:fp}")
                #     term1 = summand.get_term_from_lr_tabs(lr)
                #     print(f"{term1:c}")
                #
                # print("Integration by parts of derivative on second field:")
                # for lr in tabs2:
                #     print(f"{lr:fp}")
                #     term1 = summand.get_term_from_lr_tabs(lr)
                #     print(f"{term1:c}")
                # print("Old term:")
                # print(f"{term:c}")
                # print("=> Old term stays the same!")

    term_list_Term = [Term([term], term.name) for term in term_list]
    single_terms = get_type(term_list_Term)

    return single_terms, all(term_check)