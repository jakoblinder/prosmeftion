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
    TODO: Implement EOM of the SM-Lagrangian.
    Parameters
    ----------
    single_terms

    Returns
    -------

    """
    for type in single_terms.values():
        for term_mass_dim in type.values():
            if term_mass_dim.d != 6: continue
            for summand in term_mass_dim:
                check_derivatives = []
                if any([field.nD > 2 for field in summand.fields["H"]]): continue
                q = summand.fields[f"{op_config['fermionfields']['Q']['mathematica']['Q']}"]
                l = summand.fields[f"{op_config['fermionfields']['L']['mathematica']['L']}"]
                u = summand.fields[f"{op_config['fermionfields']['[u_C]']['mathematica']['uC']}"]
                d = summand.fields[f"{op_config['fermionfields']['[d_C]']['mathematica']['dC']}"]
                e = summand.fields[f"{op_config['fermionfields']['[e_C]']['mathematica']['eC']}"]
                bl = summand.fields[f"{op_config['bosonfields']['BL']['mathematica']['BL']}"]
                wl = summand.fields[f"{op_config['bosonfields']['WL']['mathematica']['WL']}"]

                fields_max_1_derivative = q + l + u + d + e + bl + wl
                if any([field.nD > 2 for field in fields_max_1_derivative]): continue

                a, b = bosons_non_conj, fermions_non_conj
                name_form = summand.name  # "".join([f"{name}{nD}" for name, nD in name_of_term])
                TERM_PATH = FORM_PATH / name_form
                TERM_PATH.mkdir(parents=True, exist_ok=True)  # Create directories if they don't exist.

                logger.info(f"Find EOMs in term {name_form}.")
                # The function which subsitute the EOMs should already be run here, so that all necessary indices are then declared afterwards.
                summand.form_higgsEOM(TERM_PATH)

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

                form += "#call antisymDerivative;\n"
                form += "#call spinorEOMidentification;\n"
                form += "#call fieldstrengthtensorEOMidentification;\n"
                form += "label 2;\n"

                # TODO: Here the function calls for all EOM substitutions should be

                # form += "* Bring indices of epsilons in order:\n"
                # epss = [list(tensor["mathematica"].values())[0] for tensor_name, tensor in
                #         op_config["tensors"].items() if
                #         "eps" in tensor_name]
                # for eps in epss:
                #     antisymeps = get_antisymEps(eps)
                #     form += f"Multiply replace_({eps},{antisymeps});\n"
                #     form += ".sort\n"
                #     form += f"Multiply replace_({antisymeps},{eps});\n"
                #     form += ".sort\n"
                # form += "\n"

                form += f"Bracket {', '.join(tensors + bosons + fermions)}, D2, EOM;\n\n"
                # form += f'#write <{TERM_PATH / "combined_term.h"}> "%E", expr\n'

                form += "Print +ss;\n"
                form += ".end\n"

                with open(TERM_PATH / f"{name_form}.frm", "w") as file:
                    file.write(form)

                # run_form(fp_cwd=TERM_PATH, filename=f"{name_form}.frm", fp_p=FORM_GENERAL_PATH)

                # terms = get_terms(TERM_PATH / "combined_term.h", as_one=True, name=name_form)
                op = op_config
                print("---------------")


    return single_terms
