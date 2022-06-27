import sys
import re
import subprocess
import logging
from fractions import Fraction
from pathlib import Path

from yaml import safe_load

from autoeft.model import Model

logger_autoeft = logging.getLogger("autoeft.projection")
logger = logger_autoeft.getChild(__name__)

PROJECTION_PATH = Path(__file__).parent.parent

# Include the autoeft package in system path for easier import
AUTOEFT_PATH = PROJECTION_PATH.parent
sys.path.append(AUTOEFT_PATH/ "autoeft")

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
op_config = configurations("op_config.yml")

def save_set(dictionary, key, value):
    """Don't overwrite already set values. """
    if key not in dictionary.keys():
        dictionary[key] = value

def update_opname_and_modelfile(op_dict):
    """
    Specify default values in op_config and get helicities of the fields and other information from the model file and insert them into the op_config dictionary.
    Parameters
    ----------
    op_dict

    Returns
    -------
        Updated op_config dictionary.
    """
    # If autoeft name isn't set, the form name without square brackets is taken.
    def autoeft_transl(kind, skip = []):
        """
        Specify kind of fields (str) and which fields are to skip (List[str]) and give auteft translation if not specified.
        """
        for name, field in op_dict[kind].items():
            if name in skip:
                continue
            if "autoeft" not in field.keys():
                field["autoeft"] = {}
                field["autoeft"][name] = field["mathematica"][name]
                if f"conj[{name}]" in field["mathematica"].keys():
                    field["autoeft"][f"[{name}+]"] = field["mathematica"][f"conj[{name}]"][1:-1]  # "[]" are removed

    # namings like [su2eps] instead of su2eps doesn't allow the automatic filling of the autoeft dictionary.
    # autoeft_transl("tensors", ["T", "gamma"])
    # autoeft_transl("bosonfields")
    # autoeft_transl("fermionfields", ["D"])

    # get helicities values and so on
    modelfields = model.fields

    for model_field in modelfields.values():
        # write {} around daggers
        model_field.tex = re.sub(r"\^(?<!\{)(\\)+dagger(?!\})", "^{\\\\dagger}", model_field.tex)
        try:
            model_field.tex_hc = re.sub(r"\^(?<!\{)(\\)+dagger(?!\})", "^{\\\\dagger}", model_field.tex_hc)
        except:
            pass

        if "+" in model_field.name:
            # skip hermitian conjugated fields
            continue
        for name, field in {**op_dict["bosonfields"], **op_dict["fermionfields"]}.items():
            if model_field.name == name:
                save_set(field, "helicity", model_field.helicity)
                save_set(field, "ac", model_field.ac)
                save_set(field, "conj", model_field.conj)
                save_set(field, "tex", model_field.tex)
                try:
                    save_set(field, "tex_hc", model_field.tex_hc)
                except:
                    pass

    # write index_structure as list of lists:
    for name, field in {**op_dict["tensors"], **op_dict["bosonfields"], **op_dict["fermionfields"]}.items():
        if type(field["index_structure"][0]) != list:
            field["index_structure"] = [field["index_structure"]]

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
fermions = [form_field for field in op_config["fermionfields"].values() for form_field in field["mathematica"].values()]
tensors = [form_field for field in op_config["tensors"].values() for form_field in field["mathematica"].values()]


def get_SUN_name(N):
    """Get Name of SU2_W out of model file."""
    for group_name, group_properties in model.sun_groups.items():
        if group_properties.N == N:
            return group_name

def factorizeCoeff():
    """
    Write FORM function which replaces dimensional constants in the coefficient by dimensionless ones.
    Example
    -------
    #procedure factorizeCoeff
        id At = [At/Ms]*Ms;
        id mu = [mu/Ms]*Ms;
        id Mu = [Mu/Ms]*Ms;
        id muM = [muM/Ms]*Ms;
    #endprocedure
    Returns
    -------
    List of defined dimensionless constants.
    """
    ms = op_config["coefficients"]["Ms"]["mathematica"]["Ms"]  # FORM expression of EFT mass Ms

    dimlessConst = []
    idStatements = []
    for key, constant in op_config["coefficients"].items():
        if key == "I" or key == "Ms":
            continue
        if constant["massdim"] != 0:
            form_coeff = list(constant['mathematica'].values())[0]
            if constant["massdim"] == 1:
                dimlessC = f"[{form_coeff:s}/{ms:s}]"
                dimlessConst.append(dimlessC)
                idStatements.append(f"id {form_coeff:s} = {dimlessC:s}*{ms:s}")
            else:
                dimlessC = f"[{form_coeff:s}/{ms:s}^{constant['massdim']:d}]"
                dimlessConst.append(dimlessC)
                idStatements.append(f"id {form_coeff:s} = {dimlessC:s}*{ms:s}^{constant['massdim']:d}")

    form = "#procedure factorizeCoeff\n"
    for idStatement in idStatements:
        form += f"\t{idStatement:s};\n"
    form += "#endprocedure\n"
    with open(FORM_GENERAL_PATH / "factorizeCoeff.prc", "w") as file:
        file.write(form)

    return dimlessConst

def coefficient_handling():
    """
    Write FORM procedure necessary for formatting the coefficient of each Summand and Term.
    """
    ms = op_config["coefficients"]["Ms"]["mathematica"]["Ms"]  # FORM expression of EFT mass Ms
    # write fractions into frac(nominator, denominator) function
    form_fractorizeCoeff = """#procedure fractorizeCoeff
.sort
\tSymbols u, v, x, y;
* Get sign
\tPolyFun sign;
\t.sort
* Write positive numerical coefficient in coeff and sign in sign
\tPolyFun;
\tid sign(x?) = coeff(sig_(x)*x)*sign(sig_(x));
\t.sort
* write everything else in frac function
\trepeat;"""
    form_fractorizeCoeff += "\t\tid d?!{" + f"{ms}" + "}^n?pos_ = frac(d^n, 1);\n"
    form_fractorizeCoeff += "\t\tid d?!{" + f"{ms}" + "}^n?neg_ = frac(1, d^-n);\n"
    form_fractorizeCoeff += """\tendrepeat;
* combine the separate fractions
\trepeat;
\t\tid frac(u?, v?)*frac(x?, y?) = frac(u * x, v * y);
\tendrepeat;
\trepeat;
\t\tid frac(u?, 1) = u;
\t\tid sign(x?) = x;
\t\tid coeff(1) = 1;
\tendrepeat;
\t.sort\n"""
    # Write the expansion mass Ms outside of a bracket
    form_fractorizeCoeff += f"\tBracket {ms};\n"
    form_fractorizeCoeff += """\t.sort
\tCollect bbracket;
#endprocedure"""
    with open(FORM_GENERAL_PATH / "fractorizeCoeff.prc", "w") as file:
        file.write(form_fractorizeCoeff)
    #  Ensure that the first term in each bracket has a positive sign
    form_positiveTerm = """#procedure positiveTerm
.sort
* Ensure that first term in the bracket is positive
\tid bbracket(x?$inb) = bbracket(x);
* Save original expression in order to work only on the first term
\t$coeffic = coefficient;
\t.sort
\tCFunction bsign;
\tDrop coefficient;
* Get only the first term and its sign
\tLocal FirstTerm = firstterm_($inb);
\t.sort
\tPolyFun bsign;
\t.sort
\tPolyFun;
\tid bsign(x?$s) = bsign(x);
\t$sign = sig_($s);
*	Print "Sign: %$", $sign;
\t.sort
\tDrop FirstTerm;
* restore original expression
\tLocal coefficient = $coeffic;
* ensure that the first term in the bracket is positive
\tid bbracket(x?) = $sign*bbracket($sign*x);
#endprocedure"""
    with open(FORM_GENERAL_PATH / "positiveTerm.prc", "w") as file:
        file.write(form_positiveTerm)

def get_antisymEps(eps:str):
    """
    Returns antisymmetric FORM expression of epsilon tensor.
    E.g.: su2eps -> su2epsA; [su2eps] -> [su2epsA].
    """
    if eps[-1] == ']':
        return eps[:-1] + "A]"
    else:
        return eps + "A"

def form_declarations():
    """
    Contains all general declarations valid for any term, i.e. for example the declaration of all fields and indices.
    Returns
    -------
    form : str
        Content of the FORM file "declarations.h".
    """
    form = ""
    form += "*--#[ tensors :\n"
    form += "CFunction " + ", ".join(tensors) + ";\n"
    form += "* Indices and functions for derivatives in SL2C notation.\n"
    form += "CFunction sigma, sigmabar;\n"
    form += "CFunction sigma2, sigmabar2;\n"
    form += "* Auxiliary antisymmtric epsilons, used in combination with replace_.\n"
    eps = [list(tensor["mathematica"].values())[0] for tensor_name, tensor in op_config["tensors"].items() if
           "eps" in tensor_name]
    form += f"CFunction {', '.join(map(lambda text: get_antisymEps(text) + '(antisymmetric)', eps))};\n"  # sl2CepsA(antisymmetric), su2epsA(antisymmetric), su3epsA(antisymmetric)
    form += "\n"
    form += "* Declare Kronecker Delta symbol for Sl2C Indices, because built in can not handle upper and lower (un-)dottet indices.\n"
    form += "* Since two indices are also symmetric when they are cyclic and vice versa and pattern matching is not allowed for symmetric function but for cyclic it is, [sl2CdK] is declared as cyclic.\n"
    dK = [list(tensor["mathematica"].values())[0] for tensor_name, tensor in op_config["tensors"].items() if
          "dK" in tensor_name]
    form += f"CFunction {', '.join(map(lambda text: text + '(cyclic)', dK))};\n"  # sl2CdK(cyclic), su2dK(cyclic), su3dK(cyclic)
    form += "\n"
    form += "*--#] tensors :\n"
    form += "\n"
    form += "*--#[ coefficient :\n"
    # Coefficient
    form += f"Symbols {', '.join(coeff + ['n'])};\n" # Symbols d, eps, lambdah, At, g1, g2, g3, mu, lambdaphi, kappa, Ms, Mu, muM, [2L[Ms,muM]], n;

    # ms = op_config["coefficients"]["Ms"]["mathematica"]["Ms"]  # FORM expression of EFT mass Ms
    # dimlessConst = []
    # for key, constant in op_config["coefficients"].items():
    #     if key == "I" or key == "Ms":
    #         continue
    #     if constant["massdim"] != 0:
    #         form_coeff = list(constant['mathematica'].values())[0]
    #         if constant["massdim"] == 1:
    #             dimlessC = f"[{form_coeff:s}/{ms:s}]"
    #             dimlessConst.append(dimlessC)
    #         else:
    #             dimlessC = f"[{form_coeff:s}/{ms:s}^{constant['massdim']:d}]"
    #             dimlessConst.append(dimlessC)
    dimlessConst = factorizeCoeff()
    form += f"Symbols {', '.join(dimlessConst):s};\n"  # Symbols [At/Ms], [mu/Ms], [Mu/Ms], [muM/Ms];
    form += "*--#] coefficient :\n"
    form += "\n"
    form += "*--#[ operators :\n"
    form += "Off Statistics;\n"
    # write commuting and anti-commuting operators of each term in separate list for initialization. Thus
    # Duplicated operators are removed.
    def get_commuting_op(op):
        """eC -> eCc, [eC+] -> [eC+c], where eCc and [eC+c] are commuting functions."""
        if op[0] == "[" and op[-1] == "]":
            commuting_op = op[:-1] + "c]"
        else:
            commuting_op = op + "c"
        return commuting_op
    form += "Function " + ", ".join(bosons) + ";\n"
    # Auxiliary commuting bosons
    form += "CFunction " + ", ".join(map(get_commuting_op, bosons)) + ";\n"
    form += "\n"
    form += "Function " + ", ".join(fermions) + ";\n"
    # Auxiliary commuting fermions
    form += "CFunction " + ", ".join(map(get_commuting_op, fermions)) + ";\n"
    form += "\n"

    SL2C_fieldstrengths = [form_field for field in op_config["bosonfields"].values() for form_field in field["mathematica"].values() if "helicity" in field.keys() if field["helicity"] == -1]
    form += f"Set Fieldc: {', '.join(map(get_commuting_op, SL2C_fieldstrengths))};\n"
    form += "\n"
    form += "CFunction xi, [xi+], chi, [chi+];\n"
    form += "\n"
    spinors = [list(field["mathematica"].values())[0] for field_name, field in op_config["fermionfields"].items() if "helicity" not in field.keys() if field_name != "D"]
    adjspinors = [list(field["mathematica"].values())[1] for field_name, field in op_config["fermionfields"].items() if "helicity" not in field.keys() if field_name != "D"]
    form += f"Set spinors: {', '.join(spinors)};\n"
    form += f"Set spinorsAdj: {', '.join(adjspinors)};\n"
    form += f"Set spinorsAll: {', '.join(spinors + adjspinors)};\n"
    form += "\n"
    form += f"Set spinorsc: {', '.join(map(get_commuting_op, spinors))};\n"
    form += f"Set spinorsAdjc: {', '.join(map(get_commuting_op, adjspinors))};\n"
    form += f"Set spinorsAllc: {', '.join(map(get_commuting_op, spinors + adjspinors))};\n"
    form += "\n"
    form += "*--#] operators :\n"
    form += "\n"
    ind = index_config
    form += "*--#[ indices :\n"
    for index_name, index in index_config.items():
        form += f"AutoDeclare Indices {index_name:8s} = {index['dimension']}; * {index['description']}\n"
#     form += """AutoDeclare Indices lor      = 4; * 4d Lorentz index
# AutoDeclare Indices lorA     = 4; * Auxiliary 4d Lorentz index
# AutoDeclare Indices spin     = 4; * Index for Gamma matrices/ spinor index
# AutoDeclare Indices spinA    = 2; * Auxiliary index for Gamma matrices/ spinor index
# AutoDeclare Indices gauge    = 2; * SU(2)-index in fundamental
# AutoDeclare Indices gaugeA   = 2; * Auxiliary SU(2)-index in fundamental
# AutoDeclare Indices gaugeadj = 3; * SU(2)-index in adjoint
# AutoDeclare Indices colf     = 3; * SU(3)-index in fundamental
# AutoDeclare Indices colfA    = 3; * Auxiliary SU(3)-index in fundamental
# AutoDeclare Indices cola     = 8; * SU(3)-index in adjoint
# AutoDeclare Indices flav     = n; * flavor index
# AutoDeclare Indices Lsl      = 2; * SL2C Index
# AutoDeclare Indices Usl      = 2; * SL2C Index
# AutoDeclare Indices Lsldot   = 2; * SL2C Index
# AutoDeclare Indices Usldot   = 2; * SL2C Index\n"""
    form += "\n"
    form += "AutoDeclare Indices op; * auxiliary index for converting between commuting and noncommuting operators.\n\n"
    form += "* Declare some Symbols for pattern matching\n"
    form += "Symbols k,m;\n"
    form += "*--#] indices :\n"  # trailing "\n" important otherwise form will not find the "fold" declarations

    coefficient_handling()

    return form

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


