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
from tioc.get_FORM_refactored.operator import Tensor, Field
from tioc.get_FORM_refactored.operators import Tensors
from tioc.get_FORM_refactored.indices import Indices_Operator
from tioc.get_FORM_refactored.index import Index, LP_Index
from tioc.get_FORM_refactored.read_write import get_terms
from tioc.sun_projection import get_type
from tioc.get_FORM_refactored.tableau import get_op_class, Young_Tableau
from tioc.get_FORM_refactored.lorentz import LR_Tableaux

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

        # Recalculate the number of derivatives 'derIndex', which stand on every index.
        term.fields[before-1].reset_derIndex()

        # Before appending the indices to the other field, adjust the 'derIndex' attribute of those indices:
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
        term.name = "".join([f"{name}{nD}" for name, nD in term.fieldstructure])

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
    print("=========================================================================")
    return LR_Tableaux(l_tab, r_tab, 1)

def get_term_from_lr_tabs(summand: Summand, lr: LR_Tableaux, ignore_op_class: bool=False):
    """
    Infer from the l_tab and r_tab the derivative structure for a given field structure.
    For this purpose, the field content and the other indices are necessary, which is why the term with a predominantly
    incorrect derivative structure must also be specified.
    -> In a first step, all SL2C-indices and thus also all derivatives of the fields are removed and also the
       SL2C-epsilon tensors are removed.
    -> Since each row in the l_tab and r_tab specifies uniquely an epsilon tensor and those specify uniquely the indices
       of themselves and also those of the field.
    Parameters
    ----------
    term
    lr = LR_tableaux(l_tab, r_tab, factor)
    ignore_op_class
        If True, the operator class of the Summand will not be checked, i.e. it is possible to have more than derivatives then allowed by the operator class.

    Returns
    -------

    """
    if not ignore_op_class:
        assert lr.l_tab.ncols() == summand.op_class.nl, "The number of columns in l_tab doesn't match the the one necessary for the operator class."
        assert lr.r_tab.ncols() == summand.op_class.nr, "The number of columns in r_tab doesn't match the the one necessary for the operator class."

    def get_eps_from_tab(tab: Young_Tableau):
        sl2Ceps = op_config["tensors"]["[sl2Ceps]"]["mathematica"]["sl2Ceps"]  # '[sl2Ceps]'
        epss = []
        for i in range(tab.ncols()):
            col = tab[:, i]
            tensor = Tensor(f"{sl2Ceps}({col[0, 0]},{col[1, 0]})")
            for j, index in enumerate(tensor.indices):
                index.lp = col[j, 0].lp
                index.derIndex = col[j, 0].derIndex

            epss.append(tensor)

        dual_indices = []
        for eps in epss:
            for index in eps.indices:
                dual_indices.append(index.dual_index)

        return Tensors(epss), Indices_Operator(dual_indices)

    eps_undotted, lsl_indices = get_eps_from_tab(lr.l_tab)
    eps_dotted, usldot_indices = get_eps_from_tab(lr.r_tab)
    # sort indices by fields:
    field_indices = {i: [] for i in range(1, summand.op_class.N + 1)}

    for index in lsl_indices:
        field_indices[index.lp.fp].append(index)
    for index in usldot_indices:
        field_indices[index.lp.fp].append(index)

    # Delete old SL2C-tensors
    del summand.tensors["sl2C"]
    # Append new SL2C-tensors
    for eps in eps_undotted:
        summand.tensors.append(eps)
    for eps in eps_dotted:
        summand.tensors.append(eps)

    def delete_sl2c_index(indices: Indices_Operator):
        """
        Recursively delete all SL2C-indices in 'indices'.
        Parameters
        ----------
        indices

        Returns
        -------

        """
        sl2c = [index.typ in ("Lsl", "Usl", "Lsldot", "Usldot") for index in indices]
        if any(sl2c):
            for i, index in enumerate(indices):
                if index.typ in ("Lsl", "Usl", "Lsldot", "Usldot"):
                    del field.indices[i]
                    break
            return delete_sl2c_index(indices)
        else:
            return indices

    # Delete old SL2C-indices
    for field in summand.fields:
        delete_sl2c_index(field.indices)
    # Append new SL2C-indices
    for field in summand.fields:
        for new_index in field_indices[field.field_pos]:
            field.indices.append(new_index)
        # Recalculate the number of derivatives 'derIndex', which stand on every index.
        field.reset_derIndex()

    # Adjust the coefficient
    summand.coeff *= Factor(lr.factor)

    return summand

def ibp(lr: LR_Tableaux, derivative: LP_Index,  N: int):
    """

    Parameters
    ----------
    lr = LR_tableaux(l_tab, r_tab, factor)
    derivative: LP_Index
        This index specifies the derivative, which should be integrated by parts:
    N
        Number of fields in the operator.

    Returns
    -------

    """
    assert lr.l_tab.ncols() == 1
    assert lr.r_tab.ncols() == 1

    found = 0

    l_tab_index = lr.l_tab.index_lp(derivative)
    if l_tab_index:
        l_row_index_der, colum = l_tab_index
        assert colum == 0
        del l_tab_index
        found += 1

    r_tab_index = lr.r_tab.index_lp(derivative)
    if r_tab_index:
        r_row_index_der, colum = r_tab_index
        assert colum == 0
        del r_tab_index
        found += 1
    # Note: l_row_index_der and r_row_index_der denote the row in LH and RH-tableau where the derivative Index sits.

    if found != 2:
        logger.warning("No index was found for the given LP_Index. Therefore, no ibp was made.")
        return lr
    else:
        # Indices of the derivative Indices where found and the integration by parts can be done.
        tabs = []
        spec_derivative = object()  # specify uniquely the derivative
        for i in range(1, N + 1):  # iterate over all possible field positions
            if i == derivative.fp:  # except the one which should be integrated
                continue
            else:
                # Change the field position and specify the derivative 'position' uniquely
                lr.l_tab[l_row_index_der, 0].lp = LP_Index(field_pos=i, derIndex=spec_derivative)
                lr.l_tab[l_row_index_der, 0].derIndex = spec_derivative

                lr.r_tab[r_row_index_der, 0].lp = LP_Index(field_pos=i, derIndex=spec_derivative)
                lr.r_tab[r_row_index_der, 0].derIndex = spec_derivative
                #Sign changes due to ibp
                lr *= (-1)

                tabs.append(lr)
        return tabs


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
                    # derivative index has to be set in order to be able to infer later on the derivative structure
                    # only from the tensors.
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
                skip_schouten_ids = False
                skip_ibp_ids = False
                if not summand.nD: skip_ibp_ids = True  # No ibp necessary if no derivative is there. -> Schouten id are still necessary!
                for field in summand.fields:
                    fp = field.field_pos
                    for sl_index in field.indices["sl"]:
                        lp = LP_Index(fp, sl_index.derIndex)
                        sl_index.lp = lp
                        assert set_lp_in_tensors(summand.tensors, sl_index, lp, sl_index.derIndex)
                    for sldot_index in field.indices["sldot"]:
                        lp = LP_Index(fp, sldot_index.derIndex)
                        sldot_index.lp = lp
                        assert set_lp_in_tensors(summand.tensors, sldot_index, lp, sldot_index.derIndex)

                # Infer from field and derivative structure the l_tab and r_tab, specifying the epsilon tensors.
                lr = get_lr_tabs(summand.tensors)

                if lr.l_tab.ncols() < 2 and lr.r_tab.ncols() < 2: skip_schouten_ids = True  # No Schouten ids can be applied when their are less then 2 epsilon tensors
                if skip_schouten_ids and skip_ibp_ids: continue  # Neither the Schouten nor the ibp relations need to be applied.

                print(f"{summand:c}")
                # Change the derivative structure by only changing the tableau indices:
                spec_derivative = object()
                lr.l_tab[0,0].lp = LP_Index(1, spec_derivative)
                lr.l_tab[0, 0].derIndex = spec_derivative  # 1
                lr.r_tab[0, 0].lp = LP_Index(1, spec_derivative)
                lr.r_tab[0, 0].derIndex = spec_derivative  # 1

                summand = get_term_from_lr_tabs(summand, lr)
                # Infer from the l_tab and r_tab the derivative structure for a given field structure:
                print("------------------------------------------")
                print(f"{summand:c}")
                print(f"{lr:fp}")
                tabs1 = ibp(lr, LP_Index(1,1), summand.op_class.N)
                tabs2 = ibp(lr, LP_Index(2,1), summand.op_class.N)
                for lr in tabs:
                    print(f"{lr:fp}")
                    term = get_term_from_lr_tabs(summand, lr)
                    print(f"{term:c}")
                print("TEST")


    return single_terms