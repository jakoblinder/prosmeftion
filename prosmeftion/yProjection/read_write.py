import re
import logging
import sys

from pathlib import Path
from yaml import safe_load
from typing import Dict, List, Tuple

from prosmeftion import CONFIG_PATH, PROJECTION_PATH, FORM_GENERAL_PATH, FORM_PATH, op_config, mathematica, escape_regex
from prosmeftion import bosons, fermions, tensors, run_form, op_pattern, index_number_pattern, LATEX_PATH

from prosmeftion.pyform.pyform import PyFORM
from prosmeftion.pyform.newpyform import PyFORM as NPyFORM


from .term import Term
from prosmeftion import index_number_pattern as inp
INPUT_PATH = PROJECTION_PATH / "BS"

logger_autoeft = logging.getLogger("autoeft.projection")
logger = logger_autoeft.getChild(__name__)

def mathematica_to_form(ma_expr: str):
    """

    Parameters
    ----------
    ma_expr
        Expression in mathematica format as one string without any linebreaks and spaces.

    Returns
    -------
        Expression in FORM compatible format as one string without any linebreaks and spaces.
    """
    logger.info("Convert mathematica input in FORM compatible expression.")
    # Replace complex conjugation of Yukawa matrices, by hermitian conjugation, i.e. swap the indices:
    # conj[yd][{flav7784, flav7493}] -> conj[yd][{flav7493,flav7784}]
    yukawa = ["yu", "yd", "ye"]
    yukawa_conj = r'|'.join(escape_regex(list(op_config["tensors"][y]["mathematica"].keys())[1]) for y in yukawa)
    ma_expr = re.sub(r"(?P<op>" + yukawa_conj + r")\[\{(?P<index1>flav" + index_number_pattern + r"),(?P<index2>flav" + index_number_pattern + r")\}\]",
                  r"\g<op>[{\g<index2>,\g<index1>}]", ma_expr)

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
        ma_expr = re.sub(escape_regex(ma), sub, ma_expr)

    # print(placeholder)
    # Remove curly braces around indices
    ma_expr = re.sub(r"\{(?P<index>[a-zA-Z0-9,]+)\}", r"\g<index>", ma_expr)
    D = sub_mathematica["cov"]
    while re.search(D + r"\[(?P<squarebrackets>([a-zA-Z0-9,<>\[\]])+)\]", ma_expr):
        # cov[...] -> cov(...) by remembering that regex always tries to match the largest pattern
        ma_expr = re.sub(D + r"\[(?P<squarebrackets>([a-zA-Z0-9,<>\[\]])+)\]",
            D + r"(\g<squarebrackets>)", ma_expr)
    # Replace [Index,Index,...] -> (Index,Index,...)
    ma_expr = re.sub(r"\[(?P<Index>([a-zA-Z0-9,])+)\]",
                  r"(\g<Index>)", ma_expr)

    # replace double product symbol "**" by "*"
    ma_expr = re.sub(r"\*\*", r"*", ma_expr)

    # substitute placeholders again
    for p_holder, form in placeholder.items():
        ma_expr = re.sub(escape_regex(p_holder), form, ma_expr)

    return ma_expr

# extract tensors fields, coefficients and so on from form file.

def read_form_1d_table(table_file: Path, table_label: str) -> Tuple[str]:
    """Return a map of table indices to content as given in the FORM file and table."""
    # pattern_2d_table = rf"Fill {table_label}\((?P<opPosition>(\d+))\)=(?P<op>([a-zA-Z0-9,\(\)\[\]\+]+);"
    op_list = []
    with open(table_file, "r") as infile:
        for line in infile.readlines():
            match = re.match(rf"Fill {table_label}\((?P<opPosition>(\d+))\) = \+ (?P<op>(" + op_pattern + "));", line.strip())
            if match:
                # print(f"Op {int(match.group('opPosition')):d}: {match.group('op'):s}")
                op = match.group("op")
                # Remove additional brackets: D(lor1,(D(lor2,(H(gauge1234))))) -> D(lor1,D(lor2,H(gauge1234)))
                cov = escape_regex(op_config["fermionfields"]["D"]["mathematica"]["cov"])
                while re.search(r"(" + cov + r"\(lor" + inp + r",\()+", op):
                    op = re.sub(cov + r"\((?P<index>lor" + inp + r"),\((?P<actedOnStuff>" + op_pattern + r")\)\)", cov + r"(\g<index>,\g<actedOnStuff>)", op)
                # D(Lsl1, Usldot1, (D(Lsl2, Usldot2 (H(gauge1234))))) -> D(Lsl1, Usldot1, D(Lsl2, Usldot2, H(gauge1234)))
                while re.search(r"(" + cov + r"\(Lsl" + inp + r",Usl" + inp + r",\()+", op):
                    op = re.sub(cov + r"\((?P<index1>Lsl" + inp + r"),(?P<index2>Usl" + inp + r"),\((?P<actedOnStuff>" + op_pattern + r")\)\)", cov + r"(\g<index1>,\g<index2>,\g<actedOnStuff>)", op)
                op_list.append(op)
            elif line.strip():
                raise ValueError

    op_list = tuple(op_list)
    return op_list


def get_ops(expr: str, groupOps: List[List[str]]):
    """
    Returns list of single operators.

    Parameters
    ----------
    expr
        FORM expression of only the operators.
    groupOps
        List of FORM compatible str representing operators which will be extracted from the term.
    Returns
    -------
    """

    if type(groupOps) != list:
        groupOps = [groupOps]


    def retry(tries:int=5):
        def decorator(func):
            def wrapper(*args, **kwargs):
                for i in range(tries):
                    try:
                        print(f"Try {i + 1:d}")
                        return func(*args, **kwargs)
                    except ConnectionError:
                        continue
            return wrapper
        return decorator


    # Extract the operators via FORM:
    form_path = Path("getOps.frm")

    @retry()
    def form():
        with PyFORM(form_path, 1, prompt="READY", input_dir=FORM_GENERAL_PATH) as form:
            sets = ""
            for i, group in enumerate(groupOps):
                sets += f"Set ops{i:d}: {', '.join(group)};\n"

            form.write(1, sets)
            form.write(1, expr)
            form.write(1, len(groupOps) - 1)

            return form.read_all(1)

    res = form()

    # with NPyFORM(form_path, 1, FORM_GENERAL_PATH, "-Z") as form:
    #     sets = ""
    #     for i, group in enumerate(groupOps):
    #         sets += f"Set ops{i:d}: {', '.join(group)};\n"
    #
    #     form.write([sets], ["READY"])
    #     form.write([expr], ["READY"])
    #     form.write([len(groupOps) - 1], ["READY"])
    #     res = form.read()
    #     res = res[0]

    # Extract operator expressions from FORM output:
    matches = re.finditer(r"ops(\d{1,2}:)", res)
    if matches:
        op_group = []
        matches = list(matches)
        for i, match in enumerate(matches):
            if i < len(matches) - 1:
                ops = res[match.end():matches[i + 1].start()]
                ops = ops[1:-1]
            else:
                ops = res[match.end():]
                ops = ops[1:-1]
            if ops:
                op_group.append(ops.split("\n"))
            else:
                op_group.append([])
    else:
        logger.error("Should have found something.")
        sys.exit("STOP")

    for i, ops in enumerate(op_group):
        if not ops:
            continue
        else:
            for j, op in enumerate(ops):
                match = re.match(f"ops{i:d}" + r"-(\d{1,2}):\s", op)
                op_group[i][j] = op[match.end():]

    return op_group


def get_terms(expression: str, name:str=""):
    """
    The terms of the form output are extracted and written in individual Term and Summand objects. Since Terms may
    consist of multiple summands, each summand is stored as a Summand object and the coefficients are separated for each
    Summand object.
    For Example: A term like a*(c*A*B - d*C*D) is written in to two different Summand objects:
    1. Summand object: Coefficient = a*c, Fields = [A,B]
    2. Summand object: Coefficient = - a*d, Fields = [C,D].

    Parameters
    ----------
    expression: str
        Form compatible expression, without linebreaks and whitespaces.
    name: str
        If name is specified, everything is written in one Term object, with the specified name.
    Returns
    -------

    """
    with PyFORM(FORM_GENERAL_PATH / "getTerms.frm", 1) as form:
        form.write(1, expression)
        # Since FORM automatically adds a linebreak after 72 characters, one has to use read_all.
        # All linebreaks are removed afterwards.
        expression = form.read_all(1)
        expression = re.sub(r"(\s)*", "", expression)

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
        print("Identify term: ", i)
        if name:
            t, f = get_ops(expr=ops, groupOps=[tensors, bosons + fermions])
        else:
            t, f = get_ops(expr=ops, groupOps=[tensors, bosons + fermions])
        sorted_terms.append({"tensors": t, "fields": f, "coefficient": coeff})

    print(sorted_terms)

    if name:
        # write all written terms as one Term object.
        terms = Term(sorted_terms, name)
    else:
        # terms = [Term([term], f"term{i:d}") for i, term in enumerate(sorted_terms)]
        terms = []
        for i, term in enumerate(sorted_terms):
            terms.append(Term([term], f"term{i:d}"))
        # terms = list(map(Term, [[term] for term in sorted_terms], [f"term{i:d}" for i in range(len(sorted_terms))]))

    return terms









