import re
import sys
from pathlib import Path
from yaml import safe_load
from autoeft.model import Model

if Path.cwd().name == "projection":
    PROJECTION_PATH = Path.cwd()
elif Path.cwd().name == "tioc":
    PROJECTION_PATH = Path.cwd().parent
elif Path.cwd().name == "smeft-invariants":
    PROJECTION_PATH = Path.cwd() / "projection"

# Include the autoeft package in system path for easier import
AUTOEFT_PATH = PROJECTION_PATH.parent
sys.path.append(AUTOEFT_PATH/ "autoeft")


CONFIG_PATH = PROJECTION_PATH / "config"
FORM_PATH = PROJECTION_PATH / "form_files"
INPUT_PATH = PROJECTION_PATH / "BS"
LATEX_PATH = PROJECTION_PATH / "Latex"

def configurations(filename):
    filename = Path(filename)
    assert filename.suffix == ".yml", "The configurations file has to be in yml style."
    with open(CONFIG_PATH / filename, "r") as file:
        config = safe_load(file)
    return config

def get_values(yml_dict):
    """
    Get the all values of the nested dictionaries.
    """
    values = []
    for i in yml_dict.values():
        values += i.values()
    return list(values)

# Read in model file:
with open(AUTOEFT_PATH / "models/ssm.yml", "r") as infile:
    model = Model(**safe_load(infile))

def sort_model_fields(fields):
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
    fields = list(fields.values())
    n = len(fields)
    for i in range(n - 1):
        for j in range(0, (n - 1) - i):
            if fields[j].helicity > fields[j + 1].helicity:
                fields[j], fields[j + 1] = fields[j + 1], fields[j]
            elif fields[j].helicity == fields[j + 1].helicity:
                if fields[j].name.lower() > fields[j + 1].name.lower():
                    fields[j], fields[j + 1] = fields[j + 1], fields[j]
    field_dict = {}
    for field in fields:
        field_dict[field.name] = field

    return field_dict

model.fields = sort_model_fields(model.fields)

opname_sorted = configurations("opname.yml")
opname = {key: value for d in opname_sorted.values() for key, value in d.items()}
opvalues = get_values(opname_sorted)  # list(opname.values())

opnameSL2C_all = configurations("opnameSL2C.yml")
opnameSL2C = opnameSL2C_all["ordinary"]
opSL2Cvalues = list(opnameSL2C.values())
# Auxiliary commuting Weyl spinors:
spinorsSL2C_c = opnameSL2C_all["auxiliary"]
spSL2C_c_values = list(spinorsSL2C_c.values())

#TODO: Tex for Fields could be read from the modelfile
field_config = configurations("fields.yml")
coeff_config = configurations("coefficients.yml")

for field in field_config.keys():
    field_config[field]["tex"] = field_config[field]["tex"].split(" | ")
for field in field_config.keys():
    if field != "D":
        field_config[field]["autoeft"] = field_config[field]["autoeft"].split(" | ")

# one_dict = {key:item for i in coeff_config.values() for key,item in i.items()}
for dict_key, dict in coeff_config.items():
    for key in dict.keys():
        coeff_config[dict_key][key]["tex"] = coeff_config[dict_key][key]["tex"].split(" | ")

coeffname = coeff_config["constant"]
coeffvalues = [constant["FORM"] for key, constant in coeffname.items() if key != "I"]  # list of all FORM expressions that should be declared.
coeffvalues.append(coeff_config["abbreviation"]["lg"]["FORM"])
abbreviation = coeff_config["abbreviation"]
# coeffkeys = get_keys(coeffname)

def escape_regex(regex):
    """
    Escapes all common metacharacters of a given regular expression with a '\'.
    Parameters
    ----------
    regex: str
        Regular expression that should be excaped.
    Returns
    -------
        Escaped regular expression.
    """
    modyfied_regex = ""
    for a in regex:
        if a in r"^[]{}().$*\+|?<>=":
            modyfied_regex += f"\{a}"
        else:
            modyfied_regex += a
    return modyfied_regex


