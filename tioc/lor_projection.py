import logging
import re
from fractions import Fraction
from pathlib import Path
from typing import List, Dict
from copy import copy

import sage.all
import sage.matrix as mx
from sage.rings.rational_field import QQ
from tioc.get_FORM_refactored.term import Term, TermType
from tioc.get_FORM_refactored.summand import Summand
from tioc.get_FORM_refactored.operators import Tensors
from tioc.get_FORM_refactored.indices import Indices_Operator
from tioc.get_FORM_refactored.index import Index
from tioc.get_FORM_refactored.read_write import get_terms
from tioc.sun_projection import get_type

from autoeft.io import load_basis
from autoeft.sun_projection import tensor_projection
from . import AUTOEFT_PATH, FORM_PATH, FORM_GENERAL_PATH, model, get_antisymEps, op_config, bosons, fermions, tensors, bosons_non_conj, fermions_non_conj, run_form, get_SUN_name
from .general import declaration_SL2C_sets

logger_autoeft = logging.getLogger("autoeft.projection")
logger = logger_autoeft.getChild(__name__)

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


def ibp_and_schouten_ids(single_terms):
    """
    Apply the integration by parts and Schouten identities to the lorentz structure.
    Parameters
    ----------
    single_terms

    Returns
    -------

    """
    # term_list = []  # flat list of Summands, which is later sorted by their types
    for type in single_terms.values():
        for term_mass_dim in type.values():
            for summand in term_mass_dim:
                pass

    return single_terms