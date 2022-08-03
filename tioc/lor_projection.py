import logging
import re
from fractions import Fraction
from pathlib import Path
from typing import List, Dict, Union
from copy import copy

import sage.all
import sage.matrix as mx
from sage.rings.rational_field import QQ
from tioc.get_FORM_refactored.term import Term, TermType
from tioc.get_FORM_refactored.summand import Summand
from tioc.get_FORM_refactored.coefficient import Factor
from tioc.get_FORM_refactored.operators import Tensors
from tioc.get_FORM_refactored.indices import Indices_Operator
from tioc.get_FORM_refactored.index import Index, LP_Index
from tioc.get_FORM_refactored.read_write import get_terms
from tioc.sun_projection import get_type
from tioc.tableau import get_op_class, Young_Tableau

from autoeft.io import load_basis
from autoeft.sun_projection import tensor_projection
from . import AUTOEFT_PATH, FORM_PATH, FORM_GENERAL_PATH, model, get_antisymEps, op_config, bosons, fermions, tensors, bosons_non_conj, fermions_non_conj, run_form, get_SUN_name
from .general import declaration_SL2C_sets

logger_autoeft = logging.getLogger("autoeft.projection")
logger = logger_autoeft.getChild(__name__)

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
    def contracted_with_field(term: Summand, ref_index: Index) -> int:
        if ref_index.typ == "Usldot":
            tensors = term.tensors["sldot"]
        elif ref_index.typ == "Lsl":
            tensors = term.tensors["sl"]
        # else:
        #     tensors = term.tensors

        # Find contracted tensor and save the other index on that tensor:
        found_index = False
        for tensor in tensors:
            for i, index in enumerate(tensor.indices):
                if index == ref_index.dual_index:
                    found_index = index
                    if i == 0:
                        other_index = tensor.indices[1]
                    else:
                        other_index = tensor.indices[0]
                    break
            if found_index: break

        # Find field which is contracted with other_index:
        for field in term.fields:
            if other_index.dual_index.is_in(field.indices):
                pos = field.field_pos
                break

        return pos

    def ibp_2_fields(term: Summand, before: int, afterwards: int, derIndices: List[Index]):
        """

        Parameters
        ----------
        term
        before: int
            Field position on which the derivative acts before (1 or 2)
        afterwards
            Field position on which the derivative should act afterwards (1 or 2).
        derIndices
            These are the indices of the derivative which are then removed from the one field and attached to the other field.
        Returns
        -------

        """
        assert len(term.fields) == 2, "This function is only allowed for a term with exactly 2 fields."
        assert before in (1,2) and afterwards in (1,2)
        assert before != afterwards
        assert len(derIndices) == 2, "A derivative has to have 2 indices."

        for derIndex in derIndices:
            for i, field_index in enumerate(term.fields[before-1].indices):
                if derIndex == field_index:
                    del term.fields[before-1].indices[i]

        # Rearrange the number of derivatives, which stand on every index:
        nD = 1
        found = {i: 0 for i in range(1, term.fields[before-1].nD + 2)}
        for i in range(1, term.fields[before-1].nD + 2):
            for derIndexLsl in term.fields[before-1].indices['Lsl']:
                if derIndexLsl.derIndex == i:
                    derIndexLsl.derIndex = nD
                    found[i] += 1
            for derIndexUsldot in term.fields[before-1].indices['Usldot']:
                if derIndexUsldot.derIndex == i:
                    derIndexUsldot.derIndex = nD
                    found[i] += 1
            if found[i] == 2: nD += 1

        assert all([i in (0,2) for i in found.values()]), "Something went wrong in the replacement of derIndex of the derivative indices."

        # Before appending the indices to the other field, adjust the derIndex of those indices:
        for der_index in derIndices:
            der_index.derIndex = term.fields[afterwards-1].nD + 1

        # Append indices
        for derIndex in derIndices:
            term.fields[afterwards - 1].indices.append(derIndex)

        # Multiply coefficient of the term by (-1), since we have done essentially a partial integration.
        term.coeff *= Factor("-1")

        return term

    def balance_derivatives(term: Summand):
        """
        (Recursive function)
        Go through each derivative and check if it is contracted with a derivative on the other field.
        If yes -> ibp, i.e. deleted indices of the derivative on the one field, but store them (!), in order to
        insert the indices on the other field and multiply the coefficient by (-1).
        Note: After deleting the indices from the one field, the derIndex attribute of the remaining derivative
        indices on this field have to be set newly. When appending the indices to the other field derIndex will
        become just nD + 1.
        Parameters
        ----------
        term

        Returns
        -------

        """
        derFields = {}
        for field in term.fields:
            derIndices = {der: [] for der in range(1, field.nD + 1)}
            for index in field.indices["sl"]:
                if index.derIndex:
                    derIndices[index.derIndex].append(index)
            for index in field.indices["sldot"]:
                if index.derIndex:
                    derIndices[index.derIndex].append(index)

            derIndices_contracted_crosswise = {
                der: all([field.field_pos != contracted_with_field(term, derIndices[der][i]) for i in range(2)]) for der
                in range(1, field.nD + 1)}
            derFields[field.field_pos] = (derIndices, derIndices_contracted_crosswise)
        assert len(derFields) == 2
        # Termination condition:
        if all([all(not crossed for crossed in ind_cross[1].values()) for ind_cross in derFields.values()]) and all([field.nD <= 2 for field in term.fields]):
            return term
        else:
            if len(derFields[1][1]) == 2 and len(derFields[2][1]) == 2:
                # 2 derivatives on both fields, both with crossed indices
                # Move the first derivative of the first field
                term = ibp_2_fields(term, before=1, afterwards=2, derIndices=derFields[1][0][1])
                return balance_derivatives(term)
            elif len(derFields[1][1]) == 4 and len(derFields[2][1]) == 0:
                # 4 derivatives on the first field.
                # Move the first derivative of the first field
                term = ibp_2_fields(term, before=1, afterwards=2, derIndices=derFields[1][0][1])
                return balance_derivatives(term)
            elif len(derFields[1][1]) == 0 and len(derFields[2][1]) == 4:
                # 4 derivatives on the first field.
                # Move the first derivative of the first field
                term = ibp_2_fields(term, before=2, afterwards=1, derIndices=derFields[2][0][1])
                return balance_derivatives(term)
            else:
                # 3 derivatives on the first field and 1 on the second or 3 derivatives on the second field and 1 on the first.
                # Move the one derivative from field which has 3 derivatives, where the indices are contracted cross wise.
                nD1 = len(derFields[1][1])  #  Number of derivatives acting on the first field
                nD2 = len(derFields[2][1])  #  Number of derivatives acting on the second field
                assert nD1 != nD2
                if nD1 > nD2:
                    der = list(derFields[1][1].values()).index(True) + 1
                    term = ibp_2_fields(term, before=1, afterwards=2, derIndices=derFields[1][0][der])
                elif nD1 < nD2:
                    der = list(derFields[2][1].values()).index(True) + 1
                    term = ibp_2_fields(term, before=2, afterwards=1, derIndices=derFields[2][0][der])
                return balance_derivatives(term)

    nD4_h2h_dagger2 = single_terms[(0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0)][4]
    for term in nD4_h2h_dagger2.terms:
        term = balance_derivatives(term)
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

                # FIXME: For the moment, don't consider fermins and field strength tensors with more than 1 derivatives.
                if any([nD > 1 for nD in check_derivatives_max_1_derivative]):
                    term_list.append([Term([summand], name=summand.name)])
                    continue

                # If all fermions and field strength tensors don't have any derivatives and the higgs fields have only one,
                # then there cannot occur any EOMs of the SM
                if all([nD < 2 for nD in check_derivatives_higgs]) and all([nD < 1 for nD in check_derivatives_max_1_derivative]):
                    term_list.append([Term([summand], name=summand.name)])
                    continue

                a, b = bosons_non_conj, fermions_non_conj
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
                op = op_config
                print("---------------")

    term_list = [term for short_list in term_list for term in short_list]
    single_terms = get_type(term_list)

    return single_terms

def get_lr_tabs(tensors: Tensors):
    # get SL2C-eps tensors for undotted and dotted indices:
    eps_lh = tensors["sl"]
    eps_rh = tensors["sldot"]
    if eps_lh:
        eps0 = eps_lh[0]  # first epsilon tensor
        l_tab = Young_Tableau([[eps0.indices[0]], [eps0.indices[1]]])
        for eps in eps_lh[1:]:
            # Append a column for the next epsilon tensor
            l_tab.append_col(Young_Tableau([[eps.indices[0]], [eps.indices[1]]]))

    else:
        l_tab = Young_Tableau([])
    if eps_rh:
        eps0 = eps_rh[0]  # first epsilon tensor
        r_tab = Young_Tableau([[eps0.indices[0]], [eps0.indices[1]]])
        for eps in eps_rh[1:]:
            # Append a column for the next epsilon tensor
            r_tab.append_col(Young_Tableau([[eps.indices[0]], [eps.indices[1]]]))
    else:
        r_tab = Young_Tableau([])

    print(f"{eps_lh}\n->\n{l_tab:nice}\n{l_tab:lp}")
    print(f"{eps_rh}\n->\n{r_tab:nice}\n{r_tab:lp}")
    return l_tab, r_tab

def ibp_and_schouten_ids(single_terms):
    """
    Apply the integration by parts and Schouten identities to the lorentz structure.
    Parameters
    ----------
    single_terms

    Returns
    -------

    """
    def set_lp_in_tensors(tensors: Tensors, ref_index: Index, lp: LP_Index, derIndex: Union[bool, int]):
        """
        Set in the index in tensors which is conjugated with the index 'ref_index', the 'lp' attribute to lp
        and the derIndex attribute to 'derIndex'.
        Parameters
        ----------
        tensors
        ref_index
        lp
        derIndex

        Returns
        -------
            True if index successfully replace - otherwise False.
        """
        found = False
        for tensor in tensors:
            for sl_index in tensor.indices["sl"]:
                if sl_index == ref_index.dual_index:
                    found = True
                    sl_index.lp = lp
                    sl_index.derIndex = derIndex
                    break
            for sldot_index in tensor.indices["sldot"]:
                if sldot_index == ref_index.dual_index:
                    found = True
                    sldot_index.lp = lp
                    sldot_index.derIndex = derIndex
                    break
            if found: break

        if found:
            return True
        else:
            return False

    # term_list = []  # flat list of Summands, which is later sorted by their types
    for type in single_terms.values():
        for term_mass_dim in type.values():
            for summand in term_mass_dim:
                # op_class = get_op_class(summand.fieldcounter_stripped, summand.nD, summand.d)
                for field in summand.fields:
                    lp = LP_Index(field.field_pos)
                    for sl_index in field.indices["sl"]:
                        sl_index.lp = lp
                        assert set_lp_in_tensors(summand.tensors, sl_index, lp, sl_index.derIndex)
                    for sldot_index in field.indices["sldot"]:
                        sldot_index.lp = lp
                        assert set_lp_in_tensors(summand.tensors, sldot_index, lp, sldot_index.derIndex)
                # TODO: Infer from field and derivative structure the l_tab and r_tab for all summands.
                # TODO: Vice versa: Infer from the l_tab and r_tab the derivative structure for a given field structure.
                l_tab, r_tab = get_lr_tabs(summand.tensors)
                print("=========================================================================")
                l_tab = Young_Tableau([[3], [1]])
                r_tab = Young_Tableau([[3], [1]])

    return single_terms