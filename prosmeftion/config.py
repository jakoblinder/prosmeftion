import sys
import re
import subprocess
import logging
from fractions import Fraction
from pathlib import Path
from typing import List, Dict

from yaml import safe_load
from . import op_config

logger_autoeft = logging.getLogger("autoeft.projection")
logger = logger_autoeft.getChild(__name__)

# basis_path = AUTOEFT_PATH / Path("efts", "ssm-eft", "6", "basis")
# basis_file = io_basis.BasisFile(basis_path)
#
# basis = basis_file.get_basis()
# model = basis.model
# operator = basis[{"Q": 3, "L": 1}]
# print(model)
# print(operator)


def sort_model_fields(model):
    """
    Extract order from the helicity of the fields. Fields with identical helicity are ordered by their alphabetical order.
    Parameters
    ----------
    fields
        Dictionary containing fields.
    Returns
    -------
    Dict (ordered since Python 3.6) with the fields in the demanded order.
    """
    fields = list(model.fields.values())
    n = len(fields)
    for i in range(n - 1):
        for j in range(0, (n - 1) - i):
            if fields[j].helicity > fields[j + 1].helicity:
                fields[j], fields[j + 1] = fields[j + 1], fields[j]
            elif fields[j].helicity == fields[j + 1].helicity:
                # Note that due to the ASCII standard the letter 'A' stands before 'a' and so on
                if fields[j].name > fields[j + 1].name:
                    fields[j], fields[j + 1] = fields[j + 1], fields[j]
    field_dict = {}
    for field in fields:
        field_dict[field.name] = field

    model.fields = field_dict

    fields_sorted_ac = {}
    transl_autoeft_projection = {autoeft: proj for name, field in
                                 {**op_config["bosonfields"], **op_config["fermionfields"]}.items() if name != "D" for
                                 proj, autoeft in field["autoeft"].items()}


    ac = {name: field["ac"] for field in {**op_config["bosonfields"], **op_config["fermionfields"]}.values() for name in field["mathematica"].values() if "ac" in field.keys()}


    for name in model.fields.keys():
        proj_name = transl_autoeft_projection[name]
        fields_sorted_ac[proj_name] = ac[proj_name]
    del transl_autoeft_projection


    return fields_sorted_ac, model




# TODO: Move most parts of the __init__.py file here.

