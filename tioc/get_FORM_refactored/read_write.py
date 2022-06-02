import re
import logging

from pathlib import Path
from yaml import safe_load
from typing import Dict, List, Tuple

from tioc import CONFIG_PATH, PROJECTION_PATH, FORM_GENERAL_PATH, FORM_PATH, opname, mathematica, escape_regex
from tioc import bosons, fermions, tensors, run_form, op_pattern, index_number_pattern, LATEX_PATH

from tioc.get_FORM_refactored.term import Term
from tioc import index_number_pattern as inp
INPUT_PATH = PROJECTION_PATH / "BS"

logger_autoeft = logging.getLogger("autoeft.projection")
logger = logger_autoeft.getChild(__name__)

def mathematica_to_form(inputfile: Path, outputfile: Path, header: int = 0):
    """

    Parameters
    ----------
    inputfile
        Input file in mathematica format.
    header
        Number of lines before the expression starts.
    header
        Path to the output file in FORM compatible format. If only a file name is given, the file will be saved at the
        same location as the inputfile.
    Returns
    -------
        Path to the output file.
    """
    logger.info("Convert mathematica input in FORM compatible expression.")
    inputfile = Path(inputfile)
    if inputfile.parent == Path("."):
        inputfile = INPUT_PATH / inputfile
    elif not inputfile.is_absolute():
        # Relative location w.r.t. to the working directory.
        inputfile = inputfile.resolve()

    outputfile = Path(outputfile)
    if outputfile.parent == Path("."):
        outputfile = INPUT_PATH / outputfile
    elif not outputfile.is_absolute():
        # Relative location w.r.t. to the working directory.
        outputfile = outputfile.resolve()

    with open(inputfile, "r") as infile:
        for i in range(header):
            infile.readline()
        line = infile.read()
    line = re.sub(r"(\s)*", "", line)  # Replace all whitespaces and newlines: \s = [\t\n\r\f\v]
    line = line[1:-1] # remove curly braces around expression

    # Replace complex conjugation of Yukawa matrices, by hermitiant conjugation, i.e. swap the indices:
    # conj[yd][{flav7784, flav7493}] -> conj[yd][{flav7493,flav7784}]
    yukawa = ["yu", "yd", "ye"]
    yukawa_conj = r'|'.join(escape_regex(list(opname["tensors"][y]["mathematica"].keys())[1]) for y in yukawa)
    line = re.sub(r"(?P<op>" + yukawa_conj + r")\[\{(?P<index1>flav" + index_number_pattern + r"),(?P<index2>flav" + index_number_pattern + r")\}\]",
                  "\g<op>[{\g<index2>,\g<index1>}]", line)

    sub_mathematica = {key: value for typ in mathematica.values() for key, value in typ.items()}
    placeholder = {}
    i = 0
    for ma, form in sub_mathematica.items():
        if r"[" in form or r"]" in form:
            sub = f"<{i}>"
            placeholder[sub] = form
            i += 1
        else:
            sub = form
        line = re.sub(escape_regex(ma), sub, line)

    # print(placeholder)
    # Remove curly braces around indices
    line = re.sub(r"\{(?P<index>[a-zA-Z0-9,]+)\}", "\g<index>", line)
    D = sub_mathematica["cov"]
    while re.search(D + r"\[(?P<squarebrackets>([a-zA-Z0-9,<>\[\]])+)\]", line):
        # cov[...] -> cov(...) by remembering that regex always tries to match the largest pattern
        line = re.sub(D + r"\[(?P<squarebrackets>([a-zA-Z0-9,<>\[\]])+)\]",
            D + "(\g<squarebrackets>)", line)
    # Replace [Index,Index,...] -> (Index,Index,...)
    line = re.sub(r"\[(?P<Index>([a-zA-Z0-9,])+)\]",
                  "(\g<Index>)", line)

    # replace double product symbol "**" by "*"
    line = re.sub(r"\*\*", r"*", line)

    # substitute placeholders again
    for p_holder, form in placeholder.items():
        line = re.sub(escape_regex(p_holder), form, line)

    with open(outputfile, "w") as outfile:
        outfile.write(line)

    return outputfile

# extract tensors fields, coefficients and so on from form file.

def read_form_1d_table(table_file: Path, table_label: str) -> Tuple[str]:
    """Return a map of table indices to content as given in the FORM file and table."""
    # pattern_2d_table = rf"Fill {table_label}\((?P<opPosition>(\d+))\)=(?P<op>([a-zA-Z0-9,\(\)\[\]\+]+);"
    op_list = []
    with open(table_file, "r") as infile:
        for line in infile.readlines():

            """v3.8
            if found := pattern_2d_table.fullmatch(line.strip()):
            """

            match = re.match(rf"Fill {table_label}\((?P<opPosition>(\d+))\) = \+ (?P<op>(" + op_pattern + "));", line.strip())
            if match:
                # print(f"Op {int(match.group('opPosition')):d}: {match.group('op'):s}")
                op = match.group("op")
                # Remove addtional bracket: D(lor1,(D(lor2,(H(gauge1234))))) -> D(lor1,D(lor2,H(gauge1234)))
                cov = escape_regex(opname["fermionfields"]["D"]["mathematica"]["cov"])
                while re.search(r"(" + cov + r"\(lor" + inp + r",\()+", op):
                    op = re.sub(cov + r"\((?P<index>lor" + inp + r"),\((?P<actedOnStuff>" + op_pattern + ")\)\)", cov + "(\g<index>,\g<actedOnStuff>)", op)
                op_list.append(op)
            elif line.strip():
                raise ValueError

    op_list = tuple(op_list)
    return op_list


def get_ops(expr: str, groupOps: List[List[str]], dir_name: Path, maxDimLagr:int=6, minDimOp:int=1):
    """
    Returns list of single operators.

    Parameters
    ----------
    expr
        FORM expression of only the operators.
    groupOps
        List of FORM compatible str representing operators which will be extracted from the term.
    dir_name
        Path of FORM files specific to the term.
    maxDimLagr
        Maximum mass dimension of the lagrangian operators.
    minDimOp
        Mass dimension of the operator with the minimal mass dimension.
    Returns
    -------
    """
    # maximum number of operators
    maxNOp = int(maxDimLagr// minDimOp)
    form = "Function " + ", ".join(tensors) + ";\n"
    form += "#include declarations.h # declarations\n"
    form += "\n"
    if type(groupOps) != list:
        groupOps = [groupOps]

    for i, ops in enumerate(groupOps):
        form += f"Set ops{i:d}: {', '.join(ops)};\n"
    form += "\n"
    form += f"Local expression = {expr:s};\n\n"
    form += "* Maximum numbers of terms in one operators\n"
    form += f'#define nterms "{3*maxNOp}"\n\n'
    form += "Format 255;\n"
    form += "\n"
    for i in range(len(groupOps)):
        form += f"#call getOps(ops{i:d}, `nterms')\n"
    form += "\n"
    form += "Print +ss;\n"
    form += ".end"
    TERM_PATH = FORM_PATH / dir_name
    TERM_PATH.mkdir(parents=True, exist_ok=True)  # Create directories if they don't exist.
    with open(TERM_PATH / "get_Ops.frm", "w") as file:
        file.write(form)
    run_form(fp_cwd=TERM_PATH, filename="get_Ops.frm", fp_p=FORM_GENERAL_PATH)
    ops = []
    for i in range(len(groupOps)):
        file = TERM_PATH / Path(f"ops{i:d}.t")
        ops.append(read_form_1d_table(file, f"ops{i:d}Tab"))
        # file.unlink()  # delete tab files

    return ops


def get_terms(filepath: Path, as_one=False, name:str=""):
    """
    The terms of the form output are extracted and written in individual Term and Summand objects. Since Terms may consist
    of multiple summands, each summand is stored as a Summand object and the coefficients are separated for each
    Summand object.
    For Example: A term like a*(c*A*B - d*C*D) is written in to two different Summand objects:
    1. Summand object: Coefficient = a*c, Fields = [A,B]
    2. Summand object: Coefficient = - a*d, Fields = [C,D].

    Parameters
    ----------
    filepath: Path
        Path to the inputfile which contains only the expression.
    as_one: bool
        If True, all founded terms are written as one Term object.
    name: str
        If everything is written in one Term, this name can be specified.
    Returns
    -------

    """
    if not as_one:
        assert not name, "The parameter name can only be set, when as_one is True."
    form = "#include declarations.h # tensors\n"
    form += "#include declarations.h # declarations\n"
    form += "\n"
    form += "Local expression = \n"
    form += f"#include {filepath.name}\n"
    form += ";\n"
    form += ".sort\n\n"
    form += "CFunction coeff;\n"
    form += f"Bracket {', '.join(tensors + bosons + fermions)};\n"
    form += ".sort\n"
    form += "collect coeff;\n"
    form += ".sort\n\n"
    form += "CFunction term;\n"
    form += "putinside term;\n"
    form += ".sort\n"
    form += "Format nospaces;\n\n"
    form += ""
    form += "Print +s;\n"
    form += ".end"

    with open(FORM_GENERAL_PATH / "term.frm", "w") as file:
        file.write(form)

    output = run_form(fp_cwd=FORM_GENERAL_PATH, filename="term.frm", fp_p=filepath.parent)
    pattern = r"Print(\+s{1,2})?;\n{2}expression=\n{1,2}(?P<expression>(.|\n)*);"
    pattern_short = r"Print(\+s{1,2})?;\n{2}expr=(?P<expression>(.|\n)*);"  # Pattern for extremely short expressions, i.e. fitting in one line.
    match = re.search(pattern, output)
    match_short = re.search(pattern_short, output)
    if match:
        expression = re.sub(r"(\s)*", "", match.group("expression"))  # Replace all whitespaces and newlines: \s = [\t\n\r\f\v]
        # expression = match.group("expression")
    elif match_short:
        expression = re.sub(r"(\s)*", "", match_short.group("expression"))  # Replace all whitespaces and newlines: \s = [\t\n\r\f\v]
        # expression = match_short.group("expression")
    else:
        logger.error(f"No output term has been found in {output}.")
        sys.exit("STOP")

    # Extract terms:
    matches = re.finditer(r"term", expression)  # Find Terms by identifying always "term".
    if matches:
        terms = []
        matches = list(matches)
        for i, match in enumerate(matches):
            if i < len(matches) - 1:
                term = expression[match.end():matches[i + 1].start()]
                term = term[1:-2]
            else:
                term = expression[match.end():]
                term = term[1:-1]
            terms.append(term)

    # Extract coefficients and operators:
    sorted_terms = []
    for i, term in enumerate(terms):
        match = re.search(r"coeff", term)
        coeff = term[match.end()+1:-1]
        ops = term[:match.start()-1]
        t, f = get_ops(expr=ops, groupOps=[tensors, bosons + fermions], dir_name=Path(f"term{i:d}"), maxDimLagr=6, minDimOp=1)
        sorted_terms.append({"tensors": t, "fields": f, "coefficient": coeff})

    if as_one:
        terms = Term(sorted_terms, name)
    else:
        # terms = [Term([term], f"term{i:d}") for i, term in enumerate(sorted_terms)]

        terms = list(map(Term, [[term] for term in sorted_terms], [f"term{i:d}" for i in range(len(sorted_terms))]))

    latex = ""
    for term in terms:
        latex += r"\paragraph{" + f"{term.name}" + "}\n"
        latex += r"\begin{dmath}" + "\n"
        latex += f"{term:tex} \n"
        latex += r"\end{dmath}" + "\n"
        # print(f"{term:tex}")
    with open(LATEX_PATH / "terms_all.tex", "w") as file:
        file.write(latex)
    # for term in terms:
    #     print(repr(term))
    return terms









