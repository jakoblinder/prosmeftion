import sys
import re
import subprocess
import logging
from fractions import Fraction
from pathlib import Path
from typing import List, Dict

from yaml import safe_load

logger_autoeft = logging.getLogger("autoeft.projection")
logger = logger_autoeft.getChild(__name__)

PROJECTION_PATH = Path(__file__).parent.parent

# Include the autoeft package in system path for easier import
AUTOEFT_PATH = PROJECTION_PATH.parent
# sys.path.append(AUTOEFT_PATH/ "autoeft")
# sys.path.append(str(AUTOEFT_PATH))

import autoeft.io.basis as io_basis

basis_path = AUTOEFT_PATH / Path("efts", "ssm-eft", "6", "basis")
basis_file = io_basis.BasisFile(basis_path)

basis = basis_file.get_basis()
model = basis.model
operator = basis[{"Q": 3, "L": 1}]
# print(model)
# print(operator)


from autoeft.model import Model
from autoeft.io import load_basis

CONFIG_PATH = PROJECTION_PATH / "config"
CONFIG_PATH.mkdir(parents=True, exist_ok=True)  # Create directories if they don't exist.
FORM_PATH = PROJECTION_PATH / "form_files"
FORM_PATH.mkdir(parents=True, exist_ok=True)
INPUT_PATH = PROJECTION_PATH / "BS"
INPUT_PATH.mkdir(parents=True, exist_ok=True)
LATEX_PATH = PROJECTION_PATH / "Latex"
LATEX_PATH.mkdir(parents=True, exist_ok=True)
FORM_GENERAL_PATH = FORM_PATH / "general"
FORM_GENERAL_PATH.mkdir(parents=True, exist_ok=True)

def configurations(filename):
    filename = Path(filename)
    assert filename.suffix == ".yml", "The configuration file has to be written in yaml style."
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
# with open(AUTOEFT_PATH / "models/ssm.yml", "r") as infile:
#     model = Model(**safe_load(infile))

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
                # Note that due to the ASCII standard the letter 'A' stands before 'a' and so on
                if fields[j].name > fields[j + 1].name:
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

field_config = configurations("fields.yml")
coeff_config = configurations("coefficients.yml")

for field in field_config.keys():
    if field != "D":
        field_config[field]["helicity"] = str(field_config[field]["helicity"]).split(" | ")
        field_config[field]["helicity"] = list(map(float, map(Fraction, field_config[field]["helicity"])))
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
            modyfied_regex += fr"\{a}"
        else:
            modyfied_regex += a
    return modyfied_regex

###################
# FORM refactored #
###################

def n_der(mdim: int):
    """
    Calculates the maximum number of derivatives possible for a given mass dimension in the Standard Model (SM).
    At even mass dimensions, the maximum number of derivatives is d-2, since there are at least 2 Higgs bosons, with
    a mass dimension of 1 each, necessary to build a gauge invariant structure. Similar, there are at least 2 spinors,
    with a mass dimension of 3/2 necessary to build a gauge invariant structure. For odd mass dimensions, there are,
    thus, at most d-3 derivatives possible.

    Parameters
    ----------
    mdim
        Maximum mass dimension of the given Lagrangian.
    Returns
    -------
        Maximum number of derivatives for gauge invariant SM operators with mass dimension mdim.
    """
    if mdim%2:
        # Odd mass dimension
        return mdim - 3
    else:
        # Even mass dimension
        return mdim - 2

op_config = configurations("op_config.yml")

def save_set(dictionary, key, value):
    """Don't overwrite already set values."""
    if key not in dictionary.keys():
        dictionary[key] = value

def update_opname_and_modelfile(op_dict):
    """
    Specify default values in op_config and get helicities of the fields and other information from the model file
    and insert them into the op_config dictionary.
    Parameters
    ----------
    op_dict

    Returns
    -------
        Updated op_config dictionary.
    """
    # If autoeft name is not set, the form name without square brackets is taken.
    # def autoeft_transl(kind, skip = []):
    #     """
    #     Specify kind of fields (str) and which fields are to skip (List[str])
    #     and give autoeft translation if not specified.
    #     """
    #     for name, field in op_dict[kind].items():
    #         if name in skip:
    #             continue
    #         if "autoeft" not in field.keys():
    #             field["autoeft"] = {}
    #             field["autoeft"][name] = field["mathematica"][name]
    #             if f"conj[{name}]" in field["mathematica"].keys():
    #                 field["autoeft"][f"[{name}+]"] = field["mathematica"][f"conj[{name}]"][1:-1]  # "[]" are removed

    # namings like [su2eps] instead of su2eps doesn't allow the automatic filling of the autoeft dictionary.
    # autoeft_transl("tensors", ["T", "gamma"])
    # autoeft_transl("bosonfields")
    # autoeft_transl("fermionfields", ["D"])

    # get helicities values and so on
    # modelfields = model.fields

    # for model_field in modelfields.values():
    #     # write {} around daggers
    #     model_field.tex = re.sub(r"\^(?<!\{)(\\)+dagger(?!\})", "^{\\\\dagger}", model_field.tex)
    #     try:
    #         model_field.tex_hc = re.sub(r"\^(?<!\{)(\\)+dagger(?!\})", "^{\\\\dagger}", model_field.tex_hc)
    #     except AttributeError:
    #         pass
    #
    #     if "+" in model_field.name:
    #         # skip hermitian conjugated fields
    #         continue
    #     for name, field in {**op_dict["bosonfields"], **op_dict["fermionfields"]}.items():
    #         if model_field.name == name:
    #             save_set(field, "helicity", model_field.helicity)
    #             save_set(field, "ac", model_field.ac)
    #             # save_set(field, "conj", model_field.conj)
    #             save_set(field, "tex", model_field.tex)
    #             try:
    #                 save_set(field, "tex_hc", model_field.tex_hc)
    #             except KeyError:
    #                 pass
    for name, field in {**op_dict["bosonfields"], **op_dict["fermionfields"]}.items():
        try:
                field["helicity"] = Fraction(field["helicity"])
        except KeyError:
                pass
        save_set(field, "tex_hc", field["tex"] + "^{\\\\dagger}")

    # write index_structure as list of lists:
    for name, field in {**op_dict["tensors"], **op_dict["bosonfields"], **op_dict["fermionfields"]}.items():
        if type(field["index_structure"][0]) != list:
            field["index_structure"] = [field["index_structure"]]

    # Deduce massdimension of field from its helicity if helicity is defined: d = 1 + abs(helicity)
    for name, field in {**op_dict["bosonfields"], **op_dict["fermionfields"]}.items():
        if "helicity" not in field.keys():
            continue
        else:
            field["d"] = 1 + abs(field["helicity"])

    return op_dict

op_config = update_opname_and_modelfile(op_config)


mathematica = {name_typ: {key: value for fac in typ.values() for key, value in fac["mathematica"].items()} for name_typ, typ in op_config.items()}

# mathematica = {key: value for typ in op_config.values() for fac in typ.values() for key, value in fac["mathematica"].items()}

index_config = configurations("index.yml")
# If not explicitly given update SUN indices of index_config which are written in model file:
for name, index in index_config.items():
    if name == "gauge":
        ind = model.sun_groups["SU2_W"].indices
        save_set(index, "tex_indices", ind)
        save_set(index_config["gaugeadj"], "tex_indices", list(map(str.upper, ind)))
    elif name == "colf":
        ind = model.sun_groups["SU3_C"].indices
        save_set(index, "tex_indices", ind)
        save_set(index_config["cola"], "tex_indices", list(map(str.upper, ind)))

op_pattern = r"[a-zA-Z0-9,\(\)\[\]\+\_\?]+"
op_name_pattern = r"[a-zA-Z0-9\[\]\+\_]+"
index_number_pattern = r"[A-Za-z0-9]+"  # r"[A-Z0-9]+"
ambiguous_types = ["Lsl", "Usl", "gauge"]
negative_Assertion = [r"(?!dot)", r"(?!dot)", r"(?!adj)"]
index_pattern  = r"(?P<typ>"
index_pattern += r"|".join(indextyp for indextyp in index_config.keys() if indextyp not in ambiguous_types)
index_pattern += r"|" + r"|".join(map(lambda aType, negAssert: aType + negAssert, ambiguous_types, negative_Assertion))  # Lsl(?!dot)|Usl(?!dot)|gauge(?!adj)
index_pattern += r")"
# cannot use "index_number_pattern" in the following, because indices like Usldotgauge1234 occur.
index_pattern += r"(?P<id>" + index_number_pattern + r")"
dummy_index_pattern = r"N(?P<number>\d{1,2})\_\?"
del ambiguous_types
del negative_Assertion

coeff = [form_field for field in op_config["coefficients"].values() for form_field in field["mathematica"].values() if form_field != "i_"]
coeff += [form_field for field in op_config["abbreviation"].values() for form_field in field["mathematica"].values()]
bosons = [form_field for field in op_config["bosonfields"].values() for form_field in field["mathematica"].values()]
bosons_non_conj = [list(field["mathematica"].values())[0] for field in op_config["bosonfields"].values()]
fermions = [form_field for field in op_config["fermionfields"].values() for form_field in field["mathematica"].values()]
# the same like fermions, just without the derivative
fermionfields = [form_field for field in op_config["fermionfields"].values() for form_field in field["mathematica"].values() if "derivative" not in field["description"]]
fermions_non_conj = [list(field["mathematica"].values())[0] for field in op_config["fermionfields"].values()]
tensors = [form_field for field in op_config["tensors"].values() for form_field in field["mathematica"].values()]


fields_sorted = {}
transl_autoeft_projection = {autoeft: proj for name, field in {**op_config["bosonfields"], **op_config["fermionfields"]}.items() if name != "D" for proj, autoeft in field["autoeft"].items()}
for name, field in model.fields.items():
    proj_name = transl_autoeft_projection[name]
    field.form_name = proj_name
    fields_sorted[proj_name] = field
del transl_autoeft_projection

def get_basis(basispath: Path, max_dim: int):
    """Load basis from autoeft."""
    basispath = basispath.resolve()

    basis = {}  # dictionary with basis for each mass dimension from 4 to max_dim.
    for dim in range(4, max_dim + 1):
        # load_basis also returns some counters and the Hilbert series, which we don't need here...
        try:
            # basis[dim], _, _ = load_basis(basispath, dim)
            # basispath = AUTOEFT_PATH / Path("efts", "ssm-eft", "6", "basis")
            basispath = basispath / Path(f"{dim}", "basis")
            basisfile = io_basis.BasisFile(basispath)
            basis[dim] = basisfile.get_basis()
        except FileNotFoundError:
            logger.error(f"No model with the name {model.name} can be found in {AUTOEFT_PATH / Path('eft/')}.")
            sys.exit("STOP")

    return basis

def get_commuting_op(op):
    """eC -> eCc, [eC+] -> [eC+c], where eCc and [eC+c] are commuting functions."""
    if op[0] == "[" and op[-1] == "]":
        commuting_op = op[:-1] + "c]"
    else:
        commuting_op = op + "c"
    return commuting_op

def get_SUN_name(N):
    """Get Name of SU2_W out of model file."""
    for group_name, group_properties in model.sun_groups.items():
        if group_properties.N == N:
            return group_name

def get_antisymEps(eps:str):
    """
    Returns antisymmetric FORM expression of epsilon tensor.
    E.g.: su2eps -> su2epsA; [su2eps] -> [su2epsA].
    """
    if eps[-1] == ']':
        return eps[:-1] + "A]"
    else:
        return eps + "A"

def run_form(fp_cwd: Path, filename: Path, fp_p: Path = None, keep_backslash = False):
    """
    Run form in linux terminal.
    Parameters
    ----------
    fp_cwd
        Path of the working directory.
    filename
        name of the file which should be run.
    fp_p (optional)
        Path of additional input files for form.
    keep_backslash
        True in order to keep "\" which is important for tex outputs.
    Returns
    -------
        In console printed FORM output.
    """
    if fp_p:
        command = ["form", "-p", fp_p, f"{filename}"]
    else:
        command = ["form", f"{filename}"]
    try:
        formprocess = subprocess.run(
            command,
            cwd=fp_cwd,  # path of the working directory
            capture_output=True,  # If capture_output is true, stdout and stderr will be captured.
            text=True,  # output in stdout is now a string and not a byte sequence anymore
            check=True,  # If check is true, and the process exits with a non-zero exit code, a CalledProcessError
            # exception will be raised. Attributes of that exception hold the arguments, the exit code,
            # and stdout and stderr if they were captured.
        )
    except subprocess.CalledProcessError as exc:
        exc.cmd = list(map(str, list(exc.cmd)))
        logger.error(" ".join(exc.cmd) + "\n" + str(exc.stdout))
        logger.error(f"form returned non-zero exit status {exc.returncode}")
        sys.exit("STOP")
    else:
        # No Error occured
        logger.debug(" ".join(map(str, formprocess.args)))
        output = formprocess.stdout
        output = re.sub(r" *", "", output)  # Remove only all whitespaces
        # Sometimes FORM splits indices in long expressions with an backslash "\" which is discarded:
        if not keep_backslash:
            output = re.sub(r"\\", "", output)
        return output

def get_expression_from_FORM_output(output: str):
    """
    Return only the expression from the FORM output.
    Parameters
    ----------
    output

    Returns
    -------

    """
    pattern = r"Print(\+s{1,2})?;(.end)?\n{2}expr=\n{1,2}(?P<expression>(.|\n)*);"
    pattern_short = r"Print(\+s{1,2})?;(.end)?\n{2}expr=(?P<expression>(.|\n)*);"  # Pattern for extremely short expressions, i.e. fitting in one line.
    match = re.search(pattern, output)
    match_short = re.search(pattern_short, output)
    if match:
        # expression = re.sub(r"(\s)*", "", match.group("expression"))  # Replace all whitespaces and newlines: \s = [\t\n\r\f\v]
        return match.group("expression")
    elif match_short:
        # expression = re.sub(r"(\s)*", "", match_short.group("expression"))  # Replace all whitespaces and newlines: \s = [\t\n\r\f\v]
        return match_short.group("expression")
    else:
        logger.error(f"No output term has been found in {output}.")
        sys.exit("STOP")


# def retry(tries:int=5):
#     """
#     Decorator for running FORM via PIPE. Since the PIPE connection doesn't work always, the process sometimes needs
#     to be restarted. This is what this decorator establishes when it is used like in the following example:
#
#     @retry()
#     def form():
#         with PyFORM(form_path, 1, prompt="READY", input_dir=<FORM_Path>) as form
#
#             form.write(1, "STUFF")
#
#             return form.read_all(1)
#
#     res = form()
#
#     Parameters
#     ----------
#     tries:
#         Number of times the PIPE connection is started again.
#
#     Returns
#     -------
#
#     """
#     def decorator(func):
#         def wrapper(*args, **kwargs):
#             for i in range(tries):
#                 try:
#                     print(f"Try {i + 1:d}")
#                     return func(*args, **kwargs)
#                 except ConnectionError:
#                     continue
#         return wrapper
#     return decorator

