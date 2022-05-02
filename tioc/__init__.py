import re
from pathlib import Path
from yaml import safe_load

if Path.cwd().name == "projection":
    PROJECTION_PATH = Path.cwd()
elif Path.cwd().name == "tioc":
    PROJECTION_PATH = Path.cwd().parent
elif Path.cwd().name == "smeft-invariants":
    PROJECTION_PATH = Path.cwd() / "projection"

CONFIG_PATH = PROJECTION_PATH / "config"

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
# def get_keys(yml_dict):
#     """
#     Get the all keys of the nested dictionaries.
#     """
#     keys = []
#     for i in yml_dict.values():
#         keys += i.keys()
#     return list(keys)

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


