import logging
import sys
from typing import List
from pathlib import Path
from typing import List, Dict
from itertools import permutations

from . import op_config, bosons, fermions, tensors, coeff, index_config, n_der
from . import get_antisymEps
from . import PROJECTION_PATH, CONFIG_PATH, FORM_PATH, FORM_GENERAL_PATH, INPUT_PATH, LATEX_PATH, AUTOEFT_PATH

logger_autoeft = logging.getLogger("autoeft.projection")
logger = logger_autoeft.getChild(__name__)

print("TEST")

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

def form_coefficient_handling():
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

def form_antisymDerivative(n_der: int):
    """
    Replace derivatives with 2 antisymmetric SL2C indices, e.g.
    eps_{adot,bdot}*D^{adot}_{a}*D^{bdot}_{b}
    by D^2 and a Fieldstrenghttensor.
    The D^2 could when acting on the Higgs be the EOM of the Higgs.
    Write the first with
    eps_{a,b}*D2(H(gauge1))
    and the second with
    FL(a,b,H(gauge1)).
    Note: FL denotes the total left handed part of the fieldstrengthtensor, which originates from the commutator of
    the covariant derivatives which acting on the Higgsfield. Therefore, only the part of the total Fieldstrength
    tensor for those groups under which the Higgsfield is charged are relavant, i.e. the SU(3) is for example not
    relevant for the Higgsfield.

    Parameters
    ----------
    n_der

    Returns
    -------

    """
    form = "#procedure antisymDerivative\n"

    assert n_der >= 2
    cov = op_config["fermionfields"]["D"]["mathematica"]["cov"]  # 'D'
    sl2Ceps = op_config["tensors"]["[sl2Ceps]"]["mathematica"]["sl2Ceps"]  # '[sl2Ceps]'
    # Replace dotted epsilon contracted with 2 derivatives by Fieldstrengthtensor and D^2:\n
    id_statements = []
    form += "* Replace dotted epsilon contracted with 2 derivatives by Fieldstrengthtensor and D^2:\n"
    for der in range(2, n_der + 1):
        # permutation of numbers 1,2,...,der:
        p = list(permutations(range(1, der + 1)))
        for per in p:
            if per.index(1) > per.index(2):
                minus_one_RHS = True
            else:
                minus_one_RHS = False

            def except_dotted_eps(fp_per, fp_minus_one_RHS):
                # Derivatives
                func = " * "
                for i in fp_per:
                    if i in (1,2):
                        func += f"{cov:s}(Lsl{i:d}?Lsl, Usldot{i:d}?ULsldot[i{i:d}sldot], "
                    else:
                        func += f"{cov:s}(Lsl{i:d}?Lsl, Usldot{i:d}?Usldot, "
                func += f"{'H'}?!" + "{" + f"{cov:s}" + "}(?a)"
                func += der * ")"
                func += " = "
                if fp_minus_one_RHS:
                    func += "(-1)*("
                # extract remaining indices, i.e. without 1 and 2 in the same order they occurred in the permutation.
                red_per = [i for i in fp_per if i not in (1,2)]
                red_der = [f"{cov:s}(Lsl{i:d}, Usldot{i:d}, " for i in red_per]
                func += f"{sl2Ceps:s}(Lsl1, Lsl2)"
                func += " * "
                func += f"{''.join(red_der)}D2(H(?a)){len(red_per)*')'}"
                func += f" - "
                func += f"i_ * FL(Lsl1, Lsl2, {''.join(red_der)}H(?a){len(red_per)*')'})"
                if fp_minus_one_RHS:
                    func += ")"
                return func

            id_statement1  = f"{sl2Ceps:s}(Lsldot1?LUsldot[i1sldot], Lsldot2?LUsldot[i2sldot])"
            id_statement1 += except_dotted_eps(per, minus_one_RHS)
            id_statement2 = f"{sl2Ceps:s}(Lsldot2?LUsldot[i2sldot], Lsldot1?LUsldot[i1sldot])"
            id_statement2 += except_dotted_eps(per, not minus_one_RHS)

            id_statements.append(f"id once ifmatch -> 2 {id_statement1:s};\n")
            id_statements.append(f"id once ifmatch -> 2 {id_statement2:s};\n")
            id_statements.append("\n")

    for id in id_statements:
        form += f"\t{id:s}"

    # Replace undotted epsilon contracted with 2 derivatives by Fieldstrengthtensor and D^2:\n
    id_statements = []
    form += "* Replace undotted epsilon contracted with 2 derivatives by Fieldstrengthtensor and D^2:\n"
    for der in range(2, n_der + 1):
        # permutation of numbers 1,2,...,der:
        p = list(permutations(range(1, der + 1)))
        for per in p:
            if per.index(1) > per.index(2):
                minus_one_RHS = True
            else:
                minus_one_RHS = False

            def except_unddotted_eps(fp_per, fp_minus_one_RHS):
                # Derivatives
                func = " * "
                for i in fp_per:
                    if i in (1, 2):
                        func += f"{cov:s}(Lsl{i:d}?LUsl[i{i:d}sl], Usldot{i:d}?Usldot, "
                    else:
                        func += f"{cov:s}(Lsl{i:d}?Lsl, Usldot{i:d}?Usldot, "
                func += f"{'H'}?!" + "{" + f"{cov:s}" + "}(?a)"
                func += der * ")"
                func += " = "
                if fp_minus_one_RHS:
                    func += "(-1)*("
                # extract remaining indices, i.e. without 1 and 2 in the same order they occurred in the permutation.
                red_per = [i for i in fp_per if i not in (1, 2)]
                red_der = [f"{cov:s}(Lsl{i:d}, Usldot{i:d}, " for i in red_per]
                func += f"{sl2Ceps:s}(Usldot1, Usldot2)"
                func += " * "
                func += f"{''.join(red_der)}D2(H(?a)){len(red_per) * ')'}"
                func += f" + "
                func += f"i_ * FR(Usldot1, Usldot2, {''.join(red_der)}H(?a){len(red_per) * ')'})"
                if fp_minus_one_RHS:
                    func += ")"
                return func

            id_statement1 = f"{sl2Ceps:s}(Usl1?ULsl[i1sl], Usl2?ULsl[i2sl])"
            id_statement1 += except_unddotted_eps(per, minus_one_RHS)
            id_statement2 = f"{sl2Ceps:s}(Usl2?ULsl[i2sl], Usl1?ULsl[i1sl])"
            id_statement2 += except_unddotted_eps(per, not minus_one_RHS)

            id_statements.append(f"id once ifmatch -> 2 {id_statement1:s};\n")
            id_statements.append(f"id once ifmatch -> 2 {id_statement2:s};\n")
            id_statements.append("\n")

    for id in id_statements:
        form += f"\t{id:s}"

    form += "label 2;\n"

    form += "\n"
    form += "* Discard fieldstrengthtensor which are contracted with an SL2C epsilontensor, because the SL2C indices are symmetric:\n"
    form += f"\tid FL(Lsl1?LUsl[i1], Lsl2?LUsl[i2], ?a) * {sl2Ceps:s}(Usl1?ULsl[i1], Usl2?ULsl[i2]) = 0;\n"
    form += f"\tid FL(Lsl2?LUsl[i2], Lsl1?LUsl[i1], ?a) * {sl2Ceps:s}(Usl1?ULsl[i1], Usl2?ULsl[i2]) = 0;\n"
    form += f"\tid FR(Usldot1?ULsldot[i1], Usldot2?ULsldot[i2], ?a) * {sl2Ceps:s}(Lsldot1?LUsldot[i1], Lsldot2?LUsldot[i2]) = 0;\n"
    form += f"\tid FR(Usldot2?ULsldot[i2], Usldot1?ULsldot[i1], ?a) * {sl2Ceps:s}(Lsldot1?LUsldot[i1], Lsldot2?LUsldot[i2]) = 0;\n"

    form += "#endprocedure"

    with open(FORM_GENERAL_PATH / "antisymDerivative.prc", "w") as file:
        file.write(form)

def form_simplifySigma2():
    """
    Assume that only sigma2 (i.e. sigma matrices with two lorentz indices) with two upper undotted and sigmabar2
    with two lower dotted indices exist, because this is the only relevant case for the fieldstrengthtensor replacement.
    Returns
    -------

    """
    form = ""
    # TODO: See method in class_term.py an rearrange for the now possible new generation of indices. -> Cannot be done in general folder.

def form_replaceSigmabyEps():
    """
    Replace contracted sigmas by SL2C epsilon tensors.
    Returns
    -------

    """
    sl2Ceps = op_config["tensors"]["[sl2Ceps]"]["mathematica"]["sl2Ceps"]
    sigma = "sigma"
    sigmabar = "sigmabar"
    form = ""
    form += "#procedure replaceSigmabyEps\n"
    form += "repeat;\n"
    form += "* Replace sigmabar by sigma.\n"
    form += "\t" + f"id {sigmabar}(?a, Lsldot1?Lsldot, Usl2?Usl, ?b) = {sigma}(?a, Usl2, Lsldot1, ?b);\n"
    form += "\t" + f"id {sigmabar}(?a, Lsldot1?Lsldot, Lsl2?Lsl, ?b) = {sigma}(?a, Lsl2, Lsldot1, ?b);\n"
    form += "\t" + f"id {sigmabar}(?a, Usldot1?Usldot, Usl2?Usl, ?b) = {sigma}(?a, Usl2, Usldot1, ?b);\n"
    form += "\t" + f"id {sigmabar}(?a, Usldot1?Usldot, Lsl2?Lsl, ?b) = {sigma}(?a, Lsl2, Usldot1, ?b);\n"
    form += "* Replace in lorentz indices contracted sigmas by SL2C epsilontensors.\n"
    form += "\t" + f"id {sigma}(lor1?lor, Usl1?Usl, Lsldot1?Lsldot) * {sigma}(lor1?lor, Usl2?Usl, Lsldot2?Lsldot) = - 2 * {sl2Ceps}(Usl1, Usl2) * {sl2Ceps}(Lsldot1, Lsldot2);\n"
    form += "\t" + f"id {sigma}(lor1?lor, Usl1?Usl, Usldot1?Usldot) * {sigma}(lor1?lor, Usl2?Usl, Usldot2?Usldot) = + 2 * {sl2Ceps}(Usl1, Usl2) * {sl2Ceps}(Usldot1, Usldot2);\n"
    form += "\t" + f"id {sigma}(lor1?lor, Lsl1?Lsl, Lsldot1?Lsldot) * {sigma}(lor1?lor, Lsl2?Lsl, Lsldot2?Lsldot) = + 2 * {sl2Ceps}(Lsl1, Lsl2) * {sl2Ceps}(Lsldot1, Lsldot2);\n"
    form += "\t" + f"id {sigma}(lor1?lor, Lsl1?Lsl, Usldot1?Usldot) * {sigma}(lor1?lor, Lsl2?Lsl, Usldot2?Usldot) = - 2 * {sl2Ceps}(Lsl1, Lsl2) * {sl2Ceps}(Usldot1, Usldot2);\n"
    form += "*\n"
    form += "\t" + f"id {sigma}(lor1?lor, Usl1?Usl, Usldot1?Usldot) * {sigma}(lor1?lor, Usl2?Usl, Lsldot2?Lsldot) = - 2 * {sl2Ceps}(Usl1, Usl2) * {sl2Ceps}(Usldot1, Lsldot2);\n"
    form += "\t" + f"id {sigma}(lor1?lor, Usl1?Usl, Lsldot1?Lsldot) * {sigma}(lor1?lor, Usl2?Usl, Usldot2?Usldot) = + 2 * {sl2Ceps}(Usl1, Usl2) * {sl2Ceps}(Lsldot1, Usldot2);\n"
    form += "\t" + f"id {sigma}(lor1?lor, Lsl1?Lsl, Usldot1?Usldot) * {sigma}(lor1?lor, Lsl2?Lsl, Lsldot2?Lsldot) = + 2 * {sl2Ceps}(Lsl1, Lsl2) * {sl2Ceps}(Usldot1, Lsldot2);\n"
    form += "\t" + f"id {sigma}(lor1?lor, Lsl1?Lsl, Lsldot1?Lsldot) * {sigma}(lor1?lor, Lsl2?Lsl, Usldot2?Usldot) = - 2 * {sl2Ceps}(Lsl1, Lsl2) * {sl2Ceps}(Lsldot1, Usldot2);\n"
    form += "*\n"
    form += "\t" + f"id {sigma}(lor1?lor, Usl1?Usl, Lsldot1?Lsldot) * {sigma}(lor1?lor, Lsl2?Lsl, Lsldot2?Lsldot) = + 2 * {sl2Ceps}(Usl1, Lsl2) * {sl2Ceps}(Lsldot1, Lsldot2);\n"
    form += "\t" + f"id {sigma}(lor1?lor, Usl1?Usl, Usldot1?Usldot) * {sigma}(lor1?lor, Lsl2?Lsl, Usldot2?Usldot) = - 2 * {sl2Ceps}(Usl1, Lsl2) * {sl2Ceps}(Usldot1, Usldot2);\n"
    form += "\t" + f"id {sigma}(lor1?lor, Usl1?Usl, Usldot1?Usldot) * {sigma}(lor1?lor, Lsl2?Lsl, Lsldot2?Lsldot) = + 2 * {sl2Ceps}(Usl1, Lsl2) * {sl2Ceps}(Usldot1, Lsldot2);\n"
    form += "\t" + f"id {sigma}(lor1?lor, Usl1?Usl, Lsldot1?Lsldot) * {sigma}(lor1?lor, Lsl2?Lsl, Usldot2?Usldot) = - 2 * {sl2Ceps}(Usl1, Lsl2) * {sl2Ceps}(Lsldot1, Usldot2);\n"
    form += "*\n"
    form += "\t" + f"id {sigma}(lor1?lor, Lsl1?Lsl, Lsldot1?Lsldot) * {sigma}(lor1?lor, Usl2?Usl, Lsldot2?Lsldot) = - 2 * {sl2Ceps}(Lsl1, Usl2) * {sl2Ceps}(Lsldot1, Lsldot2);\n"
    form += "\t" + f"id {sigma}(lor1?lor, Lsl1?Lsl, Usldot1?Usldot) * {sigma}(lor1?lor, Usl2?Usl, Usldot2?Usldot) = + 2 * {sl2Ceps}(Lsl1, Usl2) * {sl2Ceps}(Usldot1, Usldot2);\n"
    form += "\t" + f"id {sigma}(lor1?lor, Lsl1?Lsl, Usldot1?Usldot) * {sigma}(lor1?lor, Usl2?Usl, Lsldot2?Lsldot) = - 2 * {sl2Ceps}(Lsl1, Usl2) * {sl2Ceps}(Usldot1, Lsldot2);\n"
    form += "\t" + f"id {sigma}(lor1?lor, Lsl1?Lsl, Lsldot1?Lsldot) * {sigma}(lor1?lor, Usl2?Usl, Usldot2?Usldot) = + 2 * {sl2Ceps}(Lsl1, Usl2) * {sl2Ceps}(Lsldot1, Usldot2);\n"
    form += "endrepeat;\n"
    form += "#endprocedure\n"

    with open(FORM_GENERAL_PATH / "replaceSigmabyEps.prc", "w") as file:
        file.write(form)

def form_simplifySL2CEps():
    r"""
    Simplify contracted SL2C epsilon tensors:
    Examples:
    :math:`\epsilon^{\alpha \beta} \epsilon_{\beta \gamma} = \delta^{\alpha}_{\gamma}` and
    :math:`\delta^{\alpha}_{\alpha} = 2`
    Returns
    -------

    """
    sl2Ceps = op_config["tensors"]["[sl2Ceps]"]["mathematica"]["sl2Ceps"]
    sl2CepsA = get_antisymEps(sl2Ceps)
    sl2CdK = op_config["tensors"]["[sl2CdK]"]["mathematica"]["sl2CdK"]
    sigma = "sigma"
    sigmabar = "sigmabar"
    form = ""
    form += "#procedure simplifySL2CEps\n"
    form += "repeat;\n"
    form += "* Replace epsilons by Kronecker deltas [sl2CdK](,):\n"
    # TODO: Naming could be done more systematically like for example:
    # form += "\t" + f"id {sl2Ceps}(Usl1?Usl, Usl2?ULsl[k]) * {sl2Ceps}(Lsl3?LUsl[k], Lsl4?Lsl) = + {sl2CdK}(Usl1,Lsl4);\n"
    # Nevertheless, since this would again require a lot of testing this is postponed.

    form += "\t" + f"id {sl2Ceps}(Usl1?Usl, Usl2?ULsl[k]) * {sl2Ceps}(Lsl1?LUsl[k], Lsl3?Lsl) = + {sl2CdK}(Usl1,Lsl3);\n"
    form += "\t" + f"id {sl2Ceps}(Usl1?ULsl[k], Usl2?Usl) * {sl2Ceps}(Lsl1?LUsl[k], Lsl3?Lsl) = - {sl2CdK}(Usl2,Lsl3);\n"
    form += "\t" + f"id {sl2Ceps}(Lsl1?Lsl, Lsl2?LUsl[k]) * {sl2Ceps}(Usl2?ULsl[k], Usl3?Usl) = + {sl2CdK}(Usl3,Lsl1);\n"
    form += "\t" + f"id {sl2Ceps}(Lsl1?Lsl, Lsl2?LUsl[k]) * {sl2Ceps}(Usl3?Usl, Usl2?ULsl[k]) = - {sl2CdK}(Usl3,Lsl1);\n"

    form += "\t" + f"id {sl2Ceps}(Usldot1?Usldot, Usldot2?ULsldot[k]) * {sl2Ceps}(Lsldot1?LUsldot[k], Lsldot3?Lsldot) = + {sl2CdK}(Usldot1,Lsldot3);\n"
    form += "\t" + f"id {sl2Ceps}(Usldot1?ULsldot[k], Usldot2?Usldot) * {sl2Ceps}(Lsldot1?LUsldot[k], Lsldot3?Lsldot) = - {sl2CdK}(Usldot2,Lsldot3);\n"
    form += "\t" + f"id {sl2Ceps}(Lsldot1?Lsldot, Lsldot2?LUsldot[k]) * {sl2Ceps}(Usldot2?ULsldot[k], Usldot3?Usldot) = + {sl2CdK}(Usldot3,Lsldot1);\n"
    form += "\t" + f"id {sl2Ceps}(Lsldot1?Lsldot, Lsldot2?LUsldot[k]) * {sl2Ceps}(Usldot3?Usldot, Usldot2?ULsldot[k]) = - {sl2CdK}(Usldot3,Lsldot1);\n"

    form += "* Replace epsilons with one upper and one lower index by Kronecker-deltas:\n"
    form += "\t" + f"id {sl2Ceps}(Lsl1?Lsl,Usl2?Usl) = {sl2CdK}(Lsl1,Usl2);\n"
    form += "\t" + f"id {sl2Ceps}(Usl1?Usl,Lsl2?Lsl) = {sl2CdK}(Lsl2,Usl1);\n"
    form += "\t" + f"id {sl2Ceps}(Lsldot1?Lsldot,Usldot2?Usldot) = {sl2CdK}(Lsldot1,Usldot2);\n"
    form += "\t" + f"id {sl2Ceps}(Usldot1?Usldot,Lsldot2?Lsldot) = {sl2CdK}(Lsldot2,Usldot1);\n"
    form += "* Replace only contracted Kronecker-deltas which are contracted with eps and Kronecker-deltas themself, because contractions inside one building block are not wanted:\n"  # Replace contracted Kronecker-deltas (works for every function, not just eps):\n
    form += "\t" + f"id {sl2Ceps}?" + "{" + f"{sl2Ceps},{sl2CdK},{sigma},{sigmabar}" + "}" + f"(?a,Lsl1?LUsl[k],?b)*{sl2CdK}(?c,Usl1?ULsl[k],?d) = {sl2Ceps}(?a,?c,?d,?b);\n"  # [sl2Ceps]?
    form += "\t" + f"id {sl2Ceps}?" + "{" + f"{sl2Ceps},{sl2CdK},{sigma},{sigmabar}" + "}" + f"(?a,Lsldot1?LUsldot[k],?b)*{sl2CdK}(?c,Usldot1?ULsldot[k],?d) = {sl2Ceps}(?a,?c,?d,?b);\n"  # [sl2Ceps]?
    form += "* Replace self-contracted Kronecker-deltas by the dimension (=2):\n"
    form += "\t" + f"id {sl2CdK}(Usl1?ULsl[k],Lsl1?LUsl[k]) = d_(Lsl1,Lsl1);\n"
    form += "\t" + f"id {sl2CdK}(Usldot1?ULsldot[k],Lsldot1?LUsldot[k]) = d_(Lsldot1,Lsldot1);\n"
    form += "endrepeat;\n"
    form += "* Bring indices of epsilons in order:\n"
    form += f"Multiply replace_({sl2Ceps},{sl2CepsA});\n"
    form += ".sort\n"
    form += f"Multiply replace_({sl2CepsA},{sl2Ceps});\n"
    form += "#endprocedure\n"

    with open(FORM_GENERAL_PATH / "simplifySL2CEps.prc", "w") as file:
        file.write(form)

def form_replaceSUNGenerators(N:int):
    """
    Replace two in the adjoint representation contracted SUN Generators by Kronecker deltas in fundamental indices.
    Returns
    -------

    """
    indices_fund = {2: "gauge", 3: "colf"}
    indices_ad = {2: "gaugeadj", 3: "cola"}
    fund_index = indices_fund[N]
    ad_index = indices_ad[N]

    suNdK = op_config["tensors"][f"[su{N:d}dK]"]["mathematica"][f"su{N:d}dK"]

    form = f"#procedure replaceSU{N:d}Generators\n"
    form += f"* Replace two in the adjoint representation contracted SU{N:d} Generators by Kronecker deltas in fundamental indices.\n"
    form += "repeat;\n"
    form += f"\tid T({ad_index}1?{ad_index}, {fund_index}1?{fund_index}, {fund_index}2?{fund_index}) * T({ad_index}1?{ad_index}, {fund_index}3?{fund_index}, {fund_index}4?{fund_index})"
    form += f" = (1 / 2) * ({suNdK}({fund_index}1, {fund_index}4) * {suNdK}({fund_index}3, {fund_index}2) - (1 / {N:d}) * {suNdK}({fund_index}1, {fund_index}2) * {suNdK}({fund_index}3, {fund_index}4));\n"
    form += "endrepeat;\n"
    form += "#endprocedure\n"
    with open(FORM_GENERAL_PATH / f"replaceSU{N:d}Generators.prc", "w") as file:
        file.write(form)

# TODO: Generalise simplification of epsilon tensors for arbitrary SUN groups
def form_simplifyEpsSU2():
    r"""
    Simplify contracted SU2 epsilon tensors:
    Examples:
    :math:`\epsilon^{a b} \epsilon_{b c} = \delta^{a}_{c}` and
    :math:`\delta^{a}_{a} = 2`
    Returns
    -------

    """
    su2dK = op_config["tensors"][f"[su2dK]"]["mathematica"][f"su2dK"]
    su2eps = op_config["tensors"][f"[su2eps]"]["mathematica"][f"su2eps"]
    su2epsA = get_antisymEps(su2eps)
    # fundamental index:
    f_i = "gauge"

    form = ""
    form += "#procedure simplifyEpsSU2\n"
    form += "repeat;\n"
    form += "* Replace epsilons by Kronecker deltas [sl2CdK](,):\n"
    form += "\t" + f"id {su2eps}({f_i}1?{f_i}, {f_i}2?{f_i}) * {su2eps}({f_i}2?{f_i}, {f_i}3?{f_i}) = + {su2dK}({f_i}1,{f_i}3);\n"
    form += "\t" + f"id {su2eps}({f_i}1?{f_i}, {f_i}2?{f_i}) * {su2eps}({f_i}1?{f_i}, {f_i}3?{f_i}) = - {su2dK}({f_i}2,{f_i}3);\n"
    form += "\t" + f"id {su2eps}({f_i}1?{f_i}, {f_i}2?{f_i}) * {su2eps}({f_i}3?{f_i}, {f_i}1?{f_i}) = + {su2dK}({f_i}2,{f_i}3);\n"
    form += "\t" + f"id {su2eps}({f_i}1?{f_i}, {f_i}2?{f_i}) * {su2eps}({f_i}3?{f_i}, {f_i}2?{f_i}) = - {su2dK}({f_i}1,{f_i}3);\n"
    form += "* Replace only contracted Kronecker-deltas which are contracted with eps and Kronecker-deltas themself, because contractions inside one building block are not wanted:\n"
    form += "\t" + f"id {su2eps}?" + "{" + f"{su2eps},{su2dK}" + "}" + f"(?a,{f_i}1?{f_i},?b)*{su2dK}(?c,{f_i}1?{f_i},?d) = {su2eps}(?a,?c,?d,?b);\n"  # [su2eps]?
    form += "* Replace self-contracted Kronecker-deltas by the dimension (=2):\n"
    form += "\t" + f"id {su2dK}({f_i}1?{f_i},{f_i}1?{f_i}) = d_({f_i}1,{f_i}1);\n"
    form += "endrepeat;\n"
    form += "* Bring indices of epsilons in order:\n"
    form += f"Multiply replace_({su2eps},{su2epsA});\n"
    form += ".sort\n"
    form += f"Multiply replace_({su2epsA},{su2eps});\n"
    form += "#endprocedure\n"
    with open(FORM_GENERAL_PATH / f"simplifyEpsSU2.prc", "w") as file:
        file.write(form)

def form_simplifyEpsSU3():
    """
    Simplify SU3, i.e. 3 component epsilon tensor.
    #call simplifyEpsSU3
    Returns
    -------

    """
    su3dK = op_config["tensors"][f"[su3dK]"]["mathematica"][f"su3dK"]
    su3eps = op_config["tensors"][f"[su3eps]"]["mathematica"][f"su3eps"]
    su3epsA = get_antisymEps(su3eps)
    #fundamental index:
    f_i = "colf"

    form = ""
    form += "#procedure simplifyEpsSU3\n"
    form += "repeat;"
    form += "* Replace epsilons by Kronecker deltas [sl2CdK](,):\n"
    # 3 Cyclic permutations of first eps and first cyclic permutation of second eps.
    form +="* 3 Cyclic permutations of first eps and first cyclic permutation of second eps.\n"
    form += "\t" + f"id {su3eps}({f_i}1?{f_i}, {f_i}2?{f_i}, {f_i}3?{f_i}) * {su3eps}({f_i}1?{f_i}, {f_i}4?{f_i}, {f_i}5?{f_i}) = + ({su3dK}({f_i}2,{f_i}4)*{su3dK}({f_i}3,{f_i}5) - {su3dK}({f_i}2,{f_i}5)*{su3dK}({f_i}3,{f_i}4));\n"
    form += "\t" + f"id {su3eps}({f_i}3?{f_i}, {f_i}1?{f_i}, {f_i}2?{f_i}) * {su3eps}({f_i}1?{f_i}, {f_i}4?{f_i}, {f_i}5?{f_i}) = + ({su3dK}({f_i}2,{f_i}4)*{su3dK}({f_i}3,{f_i}5) - {su3dK}({f_i}2,{f_i}5)*{su3dK}({f_i}3,{f_i}4));\n"
    form += "\t" + f"id {su3eps}({f_i}2?{f_i}, {f_i}3?{f_i}, {f_i}1?{f_i}) * {su3eps}({f_i}1?{f_i}, {f_i}4?{f_i}, {f_i}5?{f_i}) = + ({su3dK}({f_i}2,{f_i}4)*{su3dK}({f_i}3,{f_i}5) - {su3dK}({f_i}2,{f_i}5)*{su3dK}({f_i}3,{f_i}4));\n"
    # 3 Antiyclic permutations of first eps and first cyclic permutation of second eps.
    form +="* 3 Antiyclic permutations of first eps and first cyclic permutation of second eps.\n"
    form += "\t" + f"id {su3eps}({f_i}1?{f_i}, {f_i}3?{f_i}, {f_i}2?{f_i}) * {su3eps}({f_i}1?{f_i}, {f_i}4?{f_i}, {f_i}5?{f_i}) = - ({su3dK}({f_i}2,{f_i}4)*{su3dK}({f_i}3,{f_i}5) - {su3dK}({f_i}2,{f_i}5)*{su3dK}({f_i}3,{f_i}4));\n"
    form += "\t" + f"id {su3eps}({f_i}2?{f_i}, {f_i}1?{f_i}, {f_i}3?{f_i}) * {su3eps}({f_i}1?{f_i}, {f_i}4?{f_i}, {f_i}5?{f_i}) = - ({su3dK}({f_i}2,{f_i}4)*{su3dK}({f_i}3,{f_i}5) - {su3dK}({f_i}2,{f_i}5)*{su3dK}({f_i}3,{f_i}4));\n"
    form += "\t" + f"id {su3eps}({f_i}3?{f_i}, {f_i}2?{f_i}, {f_i}1?{f_i}) * {su3eps}({f_i}1?{f_i}, {f_i}4?{f_i}, {f_i}5?{f_i}) = - ({su3dK}({f_i}2,{f_i}4)*{su3dK}({f_i}3,{f_i}5) - {su3dK}({f_i}2,{f_i}5)*{su3dK}({f_i}3,{f_i}4));\n"
    # 3 Cyclic permutations of first eps and second cyclic permutation of second eps.
    form +="* 3 Cyclic permutations of first eps and second cyclic permutation of second eps.\n"
    form += "\t" + f"id {su3eps}({f_i}1?{f_i}, {f_i}2?{f_i}, {f_i}3?{f_i}) * {su3eps}({f_i}5?{f_i}, {f_i}1?{f_i}, {f_i}4?{f_i}) = + ({su3dK}({f_i}2,{f_i}4)*{su3dK}({f_i}3,{f_i}5) - {su3dK}({f_i}2,{f_i}5)*{su3dK}({f_i}3,{f_i}4));\n"
    form += "\t" + f"id {su3eps}({f_i}3?{f_i}, {f_i}1?{f_i}, {f_i}2?{f_i}) * {su3eps}({f_i}5?{f_i}, {f_i}1?{f_i}, {f_i}4?{f_i}) = + ({su3dK}({f_i}2,{f_i}4)*{su3dK}({f_i}3,{f_i}5) - {su3dK}({f_i}2,{f_i}5)*{su3dK}({f_i}3,{f_i}4));\n"
    form += "\t" + f"id {su3eps}({f_i}2?{f_i}, {f_i}3?{f_i}, {f_i}1?{f_i}) * {su3eps}({f_i}5?{f_i}, {f_i}1?{f_i}, {f_i}4?{f_i}) = + ({su3dK}({f_i}2,{f_i}4)*{su3dK}({f_i}3,{f_i}5) - {su3dK}({f_i}2,{f_i}5)*{su3dK}({f_i}3,{f_i}4));\n"
    # 3 Antiyclic permutations of first eps and second cyclic permutation of second eps.
    form +="* 3 Antiyclic permutations of first eps and second cyclic permutation of second eps.\n"
    form += "\t" + f"id {su3eps}({f_i}1?{f_i}, {f_i}3?{f_i}, {f_i}2?{f_i}) * {su3eps}({f_i}5?{f_i}, {f_i}1?{f_i}, {f_i}4?{f_i}) = - ({su3dK}({f_i}2,{f_i}4)*{su3dK}({f_i}3,{f_i}5) - {su3dK}({f_i}2,{f_i}5)*{su3dK}({f_i}3,{f_i}4));\n"
    form += "\t" + f"id {su3eps}({f_i}2?{f_i}, {f_i}1?{f_i}, {f_i}3?{f_i}) * {su3eps}({f_i}5?{f_i}, {f_i}1?{f_i}, {f_i}4?{f_i}) = - ({su3dK}({f_i}2,{f_i}4)*{su3dK}({f_i}3,{f_i}5) - {su3dK}({f_i}2,{f_i}5)*{su3dK}({f_i}3,{f_i}4));\n"
    form += "\t" + f"id {su3eps}({f_i}3?{f_i}, {f_i}2?{f_i}, {f_i}1?{f_i}) * {su3eps}({f_i}5?{f_i}, {f_i}1?{f_i}, {f_i}4?{f_i}) = - ({su3dK}({f_i}2,{f_i}4)*{su3dK}({f_i}3,{f_i}5) - {su3dK}({f_i}2,{f_i}5)*{su3dK}({f_i}3,{f_i}4));\n"
    # 3 Cyclic permutations of first eps and third cyclic permutation of second eps.
    form +="* 3 Cyclic permutations of first eps and third cyclic permutation of second eps.\n"
    form += "\t" + f"id {su3eps}({f_i}1?{f_i}, {f_i}2?{f_i}, {f_i}3?{f_i}) * {su3eps}({f_i}4?{f_i}, {f_i}5?{f_i}, {f_i}1?{f_i}) = + ({su3dK}({f_i}2,{f_i}4)*{su3dK}({f_i}3,{f_i}5) - {su3dK}({f_i}2,{f_i}5)*{su3dK}({f_i}3,{f_i}4));\n"
    form += "\t" + f"id {su3eps}({f_i}3?{f_i}, {f_i}1?{f_i}, {f_i}2?{f_i}) * {su3eps}({f_i}4?{f_i}, {f_i}5?{f_i}, {f_i}1?{f_i}) = + ({su3dK}({f_i}2,{f_i}4)*{su3dK}({f_i}3,{f_i}5) - {su3dK}({f_i}2,{f_i}5)*{su3dK}({f_i}3,{f_i}4));\n"
    form += "\t" + f"id {su3eps}({f_i}2?{f_i}, {f_i}3?{f_i}, {f_i}1?{f_i}) * {su3eps}({f_i}4?{f_i}, {f_i}5?{f_i}, {f_i}1?{f_i}) = + ({su3dK}({f_i}2,{f_i}4)*{su3dK}({f_i}3,{f_i}5) - {su3dK}({f_i}2,{f_i}5)*{su3dK}({f_i}3,{f_i}4));\n"
    # 3 Antiyclic permutations of first eps and third cyclic permutation of second eps.
    form +="* 3 Antiyclic permutations of first eps and third cyclic permutation of second eps.\n"
    form += "\t" + f"id {su3eps}({f_i}1?{f_i}, {f_i}3?{f_i}, {f_i}2?{f_i}) * {su3eps}({f_i}4?{f_i}, {f_i}5?{f_i}, {f_i}1?{f_i}) = - ({su3dK}({f_i}2,{f_i}4)*{su3dK}({f_i}3,{f_i}5) - {su3dK}({f_i}2,{f_i}5)*{su3dK}({f_i}3,{f_i}4));\n"
    form += "\t" + f"id {su3eps}({f_i}2?{f_i}, {f_i}1?{f_i}, {f_i}3?{f_i}) * {su3eps}({f_i}4?{f_i}, {f_i}5?{f_i}, {f_i}1?{f_i}) = - ({su3dK}({f_i}2,{f_i}4)*{su3dK}({f_i}3,{f_i}5) - {su3dK}({f_i}2,{f_i}5)*{su3dK}({f_i}3,{f_i}4));\n"
    form += "\t" + f"id {su3eps}({f_i}3?{f_i}, {f_i}2?{f_i}, {f_i}1?{f_i}) * {su3eps}({f_i}4?{f_i}, {f_i}5?{f_i}, {f_i}1?{f_i}) = - ({su3dK}({f_i}2,{f_i}4)*{su3dK}({f_i}3,{f_i}5) - {su3dK}({f_i}2,{f_i}5)*{su3dK}({f_i}3,{f_i}4));\n"
    #
    # 3 Cyclic permutations of first eps and first anticyclic permutation of second eps.
    form +="* 3 Cyclic permutations of first eps and first anticyclic permutation of second eps.\n"
    form += "\t" + f"id {su3eps}({f_i}1?{f_i}, {f_i}2?{f_i}, {f_i}3?{f_i}) * {su3eps}({f_i}1?{f_i}, {f_i}5?{f_i}, {f_i}4?{f_i}) = - ({su3dK}({f_i}2,{f_i}4)*{su3dK}({f_i}3,{f_i}5) - {su3dK}({f_i}2,{f_i}5)*{su3dK}({f_i}3,{f_i}4));\n"
    form += "\t" + f"id {su3eps}({f_i}3?{f_i}, {f_i}1?{f_i}, {f_i}2?{f_i}) * {su3eps}({f_i}1?{f_i}, {f_i}5?{f_i}, {f_i}4?{f_i}) = - ({su3dK}({f_i}2,{f_i}4)*{su3dK}({f_i}3,{f_i}5) - {su3dK}({f_i}2,{f_i}5)*{su3dK}({f_i}3,{f_i}4));\n"
    form += "\t" + f"id {su3eps}({f_i}2?{f_i}, {f_i}3?{f_i}, {f_i}1?{f_i}) * {su3eps}({f_i}1?{f_i}, {f_i}5?{f_i}, {f_i}4?{f_i}) = - ({su3dK}({f_i}2,{f_i}4)*{su3dK}({f_i}3,{f_i}5) - {su3dK}({f_i}2,{f_i}5)*{su3dK}({f_i}3,{f_i}4));\n"
    # 3 Antiyclic permutations of first eps and first anticyclic permutation of second eps.
    form +="* 3 Antiyclic permutations of first eps and first anticyclic permutation of second eps.\n"
    form += "\t" + f"id {su3eps}({f_i}1?{f_i}, {f_i}3?{f_i}, {f_i}2?{f_i}) * {su3eps}({f_i}1?{f_i}, {f_i}5?{f_i}, {f_i}4?{f_i}) = + ({su3dK}({f_i}2,{f_i}4)*{su3dK}({f_i}3,{f_i}5) - {su3dK}({f_i}2,{f_i}5)*{su3dK}({f_i}3,{f_i}4));\n"
    form += "\t" + f"id {su3eps}({f_i}2?{f_i}, {f_i}1?{f_i}, {f_i}3?{f_i}) * {su3eps}({f_i}1?{f_i}, {f_i}5?{f_i}, {f_i}4?{f_i}) = + ({su3dK}({f_i}2,{f_i}4)*{su3dK}({f_i}3,{f_i}5) - {su3dK}({f_i}2,{f_i}5)*{su3dK}({f_i}3,{f_i}4));\n"
    form += "\t" + f"id {su3eps}({f_i}3?{f_i}, {f_i}2?{f_i}, {f_i}1?{f_i}) * {su3eps}({f_i}1?{f_i}, {f_i}5?{f_i}, {f_i}4?{f_i}) = + ({su3dK}({f_i}2,{f_i}4)*{su3dK}({f_i}3,{f_i}5) - {su3dK}({f_i}2,{f_i}5)*{su3dK}({f_i}3,{f_i}4));\n"
    # 3 Cyclic permutations of first eps and second anticyclic permutation of second eps.
    form +="* 3 Cyclic permutations of first eps and second anticyclic permutation of second eps.\n"
    form += "\t" + f"id {su3eps}({f_i}1?{f_i}, {f_i}2?{f_i}, {f_i}3?{f_i}) * {su3eps}({f_i}4?{f_i}, {f_i}1?{f_i}, {f_i}5?{f_i}) = - ({su3dK}({f_i}2,{f_i}4)*{su3dK}({f_i}3,{f_i}5) - {su3dK}({f_i}2,{f_i}5)*{su3dK}({f_i}3,{f_i}4));\n"
    form += "\t" + f"id {su3eps}({f_i}3?{f_i}, {f_i}1?{f_i}, {f_i}2?{f_i}) * {su3eps}({f_i}4?{f_i}, {f_i}1?{f_i}, {f_i}5?{f_i}) = - ({su3dK}({f_i}2,{f_i}4)*{su3dK}({f_i}3,{f_i}5) - {su3dK}({f_i}2,{f_i}5)*{su3dK}({f_i}3,{f_i}4));\n"
    form += "\t" + f"id {su3eps}({f_i}2?{f_i}, {f_i}3?{f_i}, {f_i}1?{f_i}) * {su3eps}({f_i}4?{f_i}, {f_i}1?{f_i}, {f_i}5?{f_i}) = - ({su3dK}({f_i}2,{f_i}4)*{su3dK}({f_i}3,{f_i}5) - {su3dK}({f_i}2,{f_i}5)*{su3dK}({f_i}3,{f_i}4));\n"
    # 3 Antiyclic permutations of first eps and second anticyclic permutation of second eps.
    form +="* 3 Antiyclic permutations of first eps and second anticyclic permutation of second eps.\n"
    form += "\t" + f"id {su3eps}({f_i}1?{f_i}, {f_i}3?{f_i}, {f_i}2?{f_i}) * {su3eps}({f_i}4?{f_i}, {f_i}1?{f_i}, {f_i}5?{f_i}) = + ({su3dK}({f_i}2,{f_i}4)*{su3dK}({f_i}3,{f_i}5) - {su3dK}({f_i}2,{f_i}5)*{su3dK}({f_i}3,{f_i}4));\n"
    form += "\t" + f"id {su3eps}({f_i}2?{f_i}, {f_i}1?{f_i}, {f_i}3?{f_i}) * {su3eps}({f_i}4?{f_i}, {f_i}1?{f_i}, {f_i}5?{f_i}) = + ({su3dK}({f_i}2,{f_i}4)*{su3dK}({f_i}3,{f_i}5) - {su3dK}({f_i}2,{f_i}5)*{su3dK}({f_i}3,{f_i}4));\n"
    form += "\t" + f"id {su3eps}({f_i}3?{f_i}, {f_i}2?{f_i}, {f_i}1?{f_i}) * {su3eps}({f_i}4?{f_i}, {f_i}1?{f_i}, {f_i}5?{f_i}) = + ({su3dK}({f_i}2,{f_i}4)*{su3dK}({f_i}3,{f_i}5) - {su3dK}({f_i}2,{f_i}5)*{su3dK}({f_i}3,{f_i}4));\n"
    # 3 Cyclic permutations of first eps and third anticyclic permutation of second eps.
    form +="* 3 Cyclic permutations of first eps and third anticyclic permutation of second eps.\n"
    form += "\t" + f"id {su3eps}({f_i}1?{f_i}, {f_i}2?{f_i}, {f_i}3?{f_i}) * {su3eps}({f_i}5?{f_i}, {f_i}4?{f_i}, {f_i}1?{f_i}) = - ({su3dK}({f_i}2,{f_i}4)*{su3dK}({f_i}3,{f_i}5) - {su3dK}({f_i}2,{f_i}5)*{su3dK}({f_i}3,{f_i}4));\n"
    form += "\t" + f"id {su3eps}({f_i}3?{f_i}, {f_i}1?{f_i}, {f_i}2?{f_i}) * {su3eps}({f_i}5?{f_i}, {f_i}4?{f_i}, {f_i}1?{f_i}) = - ({su3dK}({f_i}2,{f_i}4)*{su3dK}({f_i}3,{f_i}5) - {su3dK}({f_i}2,{f_i}5)*{su3dK}({f_i}3,{f_i}4));\n"
    form += "\t" + f"id {su3eps}({f_i}2?{f_i}, {f_i}3?{f_i}, {f_i}1?{f_i}) * {su3eps}({f_i}5?{f_i}, {f_i}4?{f_i}, {f_i}1?{f_i}) = - ({su3dK}({f_i}2,{f_i}4)*{su3dK}({f_i}3,{f_i}5) - {su3dK}({f_i}2,{f_i}5)*{su3dK}({f_i}3,{f_i}4));\n"
    # 3 Antiyclic permutations of first eps and third anticyclic permutation of second eps.
    form +="* 3 Antiyclic permutations of first eps and third anticyclic permutation of second eps.\n"
    form += "\t" + f"id {su3eps}({f_i}1?{f_i}, {f_i}3?{f_i}, {f_i}2?{f_i}) * {su3eps}({f_i}5?{f_i}, {f_i}4?{f_i}, {f_i}1?{f_i}) = + ({su3dK}({f_i}2,{f_i}4)*{su3dK}({f_i}3,{f_i}5) - {su3dK}({f_i}2,{f_i}5)*{su3dK}({f_i}3,{f_i}4));\n"
    form += "\t" + f"id {su3eps}({f_i}2?{f_i}, {f_i}1?{f_i}, {f_i}3?{f_i}) * {su3eps}({f_i}5?{f_i}, {f_i}4?{f_i}, {f_i}1?{f_i}) = + ({su3dK}({f_i}2,{f_i}4)*{su3dK}({f_i}3,{f_i}5) - {su3dK}({f_i}2,{f_i}5)*{su3dK}({f_i}3,{f_i}4));\n"
    form += "\t" + f"id {su3eps}({f_i}3?{f_i}, {f_i}2?{f_i}, {f_i}1?{f_i}) * {su3eps}({f_i}5?{f_i}, {f_i}4?{f_i}, {f_i}1?{f_i}) = + ({su3dK}({f_i}2,{f_i}4)*{su3dK}({f_i}3,{f_i}5) - {su3dK}({f_i}2,{f_i}5)*{su3dK}({f_i}3,{f_i}4));\n"

    form += "* Replace only contracted Kronecker-deltas which are contracted with eps and Kronecker-deltas themself, because contractions inside one building block are not wanted:\n"
    form += "\t" + f"id {su3eps}?" + "{" + f"{su3eps},{su3dK}" + "}" + f"(?a,{f_i}1?{f_i},?b)*{su3dK}(?c,{f_i}1?{f_i},?d) = {su3eps}(?a,?c,?d,?b);\n"  # [su3eps]?
    form += "* Replace self-contracted Kronecker-deltas by the dimension (=3):\n"
    form += "\t" + f"id {su3dK}({f_i}1?{f_i},{f_i}1?{f_i}) = d_({f_i}1,{f_i}1);\n"
    form += "endrepeat;\n"
    form += "* Bring indices of epsilons in order:\n"
    form += f"Multiply replace_({su3eps},{su3epsA});\n"
    form += ".sort\n"
    form += f"Multiply replace_({su3epsA},{su3eps});\n"
    # form += 1 * "\t" + f"id {su3eps}({f_i}1?{f_i}, {f_i}2?{f_i}, {f_i}3?{f_i}) = e_({f_i}1, {f_i}2, {f_i}3);\n"
    # form += "contract 0;\n"
    # form += 1 * "\t" + f"id e_({f_i}1?{f_i}, {f_i}2?{f_i}, {f_i}3?{f_i}) = {su3eps}({f_i}1, {f_i}2, {f_i}3);\n"
    # form += 1 * "\t" + f"id d_({f_i}1?{f_i}, {f_i}2?{f_i}) = {su3dK}({f_i}1, {f_i}2);\n"
    form += "#endprocedure"
    with open(FORM_GENERAL_PATH / f"simplifyEpsSU3.prc", "w") as file:
        file.write(form)

def form_declarations(n_der: int):
    """
    Contains all general declarations valid for any term, i.e. for example the declaration of all fields and indices.
    Creates also all other prewritten FORM files.
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
    form += f"Set Field: {', '.join(SL2C_fieldstrengths)};\n"
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
    form += "* D2 = D_mu * D^mu ;\n"
    form += "Function D2;\n"
    form += "* Total field strength tensor, necessary for EOM substitutions:\n"
    form += "Function FL, FR;\n"
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
    form += "Symbols k, m;\n"
    form += "Autodeclare Symbols i;\n"
    form += "*--#] indices :\n"  # trailing "\n" important otherwise form will not find the "fold" declarations

    form_coefficient_handling()

    # form_replaceSigmabyEps()
    form_simplifySL2CEps()
    # form_replaceSUNGenerators(N=2)
    # form_replaceSUNGenerators(N=3)
    # form_simplifyEpsSU2()
    # form_simplifyEpsSU3()
    form_antisymDerivative(n_der)

    return form

def declaration_SL2C_sets(indices) -> str:
    """
    Declare SL2C-indices of a Summand and write Sets for contraction of SL2C-indices.
    Write also Sets of other index types.
    Parameters
    ----------
    indices : Indices_Summand

    Returns
    -------

    """
    lsl, usl, lsldot, usldot = indices.get_sl2C_sets()

    max_ind = 5
    # Generate only indices if there aren't any, since the could be declared twice otherwise.
    def id_indices(typ: str, reference: List, max_index: int):
        return [f"{typ:s}{i:d}" for i in range(max_index) if f"{typ:s}{i:d}" not in [repr(index) for index in reference]]
    # if not (lsl and usl):
    lsl_aux    = id_indices("Lsl", lsl, max_ind)
    usl_aux    = id_indices("Usl", lsl, max_ind) #[f"Usl{i:d}" for i in range(max_ind)]
    # else:
    #     lsl_aux = []
    #     usl_aux = []
    # if not (lsldot and usldot):
    lsldot_aux = id_indices("Lsldot", lsldot, max_ind)  # [f"Lsldot{i:d}" for i in range(max_ind)]
    usldot_aux = id_indices("Usldot", usldot, max_ind)  # [f"Usldot{i:d}" for i in range(max_ind)]
    # else:
    #     lsldot_aux = []
    #     usldot_aux = []

    decl_lsl = ", ".join([f"{index}=2" for index in lsl_aux] + [f"{repr(index)}=2" for index in lsl])
    decl_usl = ", ".join([f"{index}=2" for index in usl_aux] + [f"{repr(index)}=2" for index in usl])
    decl_lsldot = ", ".join([f"{index}=2" for index in lsldot_aux] + [f"{repr(index)}=2" for index in lsldot])
    decl_usldot = ", ".join([f"{index}=2" for index in usldot_aux] + [f"{repr(index)}=2" for index in usldot])

    set_lusl = ", ".join([i for pair in list(zip(lsl_aux, usl_aux)) for i in pair] + [repr(i) for pair in list(zip(lsl, usl)) for i in pair])
    set_ulsl = ", ".join([i for pair in list(zip(usl_aux, lsl_aux)) for i in pair] + [repr(i) for pair in list(zip(usl, lsl)) for i in pair])
    set_lusldot = ", ".join([i for pair in list(zip(lsldot_aux, usldot_aux)) for i in pair] + [repr(i) for pair in list(zip(lsldot, usldot)) for i in pair])
    set_ulsldot = ", ".join([i for pair in list(zip(usldot_aux, lsldot_aux)) for i in pair] + [repr(i) for pair in list(zip(usldot, lsldot)) for i in pair])

    form = ""
    # Declaration SL2C_indices

    form += f"Indices {decl_lsl};\n"
    form += f"Indices {decl_usl};\n"
    form += f"Indices {decl_lsldot};\n"
    form += f"Indices {decl_usldot};\n"

    # Sets for contraction for SL2C-indices
    form += f"Set Lsl: {', '.join(map(repr, lsl))};\n"
    form += f"Set Usl: {', '.join(map(repr, usl))};\n"
    form += f"Set Lsldot: {', '.join(map(repr, lsldot))};\n"
    form += f"Set Usldot: {', '.join(map(repr, usldot))};\n"

    form += f"Set LUsl: {set_lusl};\n"
    form += f"Set ULsl: {set_ulsl};\n"
    form += f"Set LUsldot: {set_lusldot};\n"
    form += f"Set ULsldot: {set_ulsldot};\n"

    # Sets for contraction for all other indices except SL2C-indices:
    for index_name in index_config.keys():
        if index_name not in ["Lsl", "Usl", "Lsldot", "Usldot"]:
            indices_typ = indices[index_name]
            if indices_typ:
                ind_aux = id_indices(index_name, indices_typ, max_ind)
                form += f"Set {index_name}: {', '.join(ind_aux)}, {indices_typ};\n"

    return form

