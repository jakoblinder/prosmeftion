import logging
import sys
from typing import List
from pathlib import Path
from typing import List, Dict
from itertools import permutations
from fractions import Fraction

from . import op_config, bosons, fermions, fermionfields, tensors, coeff, index_config, n_der
from . import get_antisymEps, fields_sorted, get_commuting_op
from . import PROJECTION_PATH, CONFIG_PATH, FORM_PATH, FORM_GENERAL_PATH, INPUT_PATH, LATEX_PATH, AUTOEFT_PATH

logger_autoeft = logging.getLogger("autoeft.projection")
logger = logger_autoeft.getChild(__name__)

def create_procedure(function_path=Path("."), function_name=""):
    """
    Define Decorator to write and save the procedure by identifying the name of the function:
    form_testfunction
    wrapped by this decorator will generate a procedure of the name
    testfunction.
    The body of the procedure will be indented.
    Parameters of the decorator
    ----------
    function_path: str
        Path where the function should be stored.
    function_name
        Name of the function in FORM. If no name is given, the nam of the python function, i.e. testfunction
        in the above example is taken.
    Returns
    -------
    """
    function_path = Path(function_path)
    def decorator(func):
        def wrapper_with_func_args(*args, **kwargs):
            if not function_name:
                procedure_name = func.__name__[5:]
            else:
                procedure_name = function_name
            # Do something before the wrapped function is called.
            form = f"#procedure {procedure_name:s}\n"
            # Call wrapped function:
            func_output = func(*args, **kwargs)
            func_output_lines = func_output.split("\n")
            for i, line in enumerate(func_output_lines):
                if not line or line[0] == "*":
                    # Empty lines or lines which contain a comment are not indented.
                    continue
                else:
                    func_output_lines[i] = "\t" + line
            form += "\n".join(func_output_lines)
            # Do something with the output of the wrapped function.
            form += "\n#endprocedure"

            # Save the function in .prc file at the desired location.
            procedure_path = Path(function_path / f"{procedure_name}.prc")
            with open(procedure_path, "w") as file:
                file.write(form)
            # print(procedure_path)
            # return form
        return wrapper_with_func_args

    return decorator

def form_factorizeCoeff():
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

def form_getTerms():
    """
    Write FORM function which separates terms like a*(c*A*B - d*C*D) into two different terms:
        +term(a*c*A*B) -term(a*d*C*D)
    .

    Creates getTerms.frm which has to be called via a pipe proces.

    Returns
    -------

    """
    form = "#include declarations_general.h # coefficient\n"
    form += "#include declarations_general.h # indices\n"
    form += "#include declarations_general.h # tensors\n"
    form += "#include declarations_general.h # operators\n"
    form += "\n"
    form += "\n#setexternal `PIPE1_'\n"
    form += "Local expression = \n"
    form += f"#fromexternal\n"
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
    form += r'#toexternal "%E\n", expression'
    form += "\n"
    form += ".end"

    with open(FORM_GENERAL_PATH / "getTerms.frm", "w") as file:
        file.write(form)

def form_getOps():
    """
    Write FORM function which extracts from a term like a*A*c*B the coefficient tensors and fields:
        coefficient tensors: 1.) a
                             3.) c
        fields: 2.) A
                4.) B
    .

    Creates getOps.frm which has to be called via a pipe process.
    Creates also getOps.prc, a procedure called from within getOps.frm.

    Returns
    -------

    """
    def form_getOpsprc():
        """Creates getOps.prc, a procedure called from within getOps.frm."""
        form = """#procedure getOps(ops)
* ops: set of non commutable functions
*   specifies which fields are printed in the table
	.sort
	Function f;
* Max is a preprocessor variable, which tells the maximum possible number of fields and tensors in one operator.
* This number is allowed to be much larger then the number of ops which really occurs.
	#define Max "1000"
	#define matched "0"
	#define matchedsomething "0"
	Table `ops'Tab(1:`Max');

	multiply left f;
	.sort
	#write "Find ops from set `ops'."
	#do i= 1, `Max'
		redefine matchedsomething "0";
		if (match(f*H?`ops'(?b)));
			id once f*H?`ops'$fac(?b$arg) = H(?b)*f;
			redefine matched "1";
		else;
*			Print "STOP";
			redefine matched "0";
		endif;
		.sort
		#if `matched' == 1
			Fill `ops'Tab(`i') = `$fac'(`$arg');
			#write "Matched `i': %$(%$)", $fac,$arg
			#toexternal "`ops'-`i': %$(%$)\n", $fac,$arg
			redefine matchedsomething "1";
			goto 1;
		#endif
		if (match(f*H?!`ops'(?b)));
			id once f*H?!`ops'$fac(?b$arg) = H(?b)*f;
			redefine matchedsomething "1";
		endif;
		label 1;
		.sort
		#if ((`matched' == 0) && (`matchedsomething' == 1))
			#write "NoMatched `i': %$(%$)", $fac,$arg
		#endif

		#if `matchedsomething' == 0
			#breakdo
		#endif	
	#enddo
	id f = 1;
#endprocedure
"""
        return form

    with open(FORM_GENERAL_PATH / "getOps.prc", "w") as file:
        file.write(form_getOpsprc())

    form =  "#include declarations_general.h # NCtensors\n"  # Tensors as non-commuting functions.
    form += "#include declarations_general.h # coefficient\n"
    form += "#include declarations_general.h # indices\n"
    form += "#include declarations_general.h # operators\n"
    form += "\n"

    form += "#setexternal `PIPE1_'\n"
    form += "#prompt READY\n"

    # Define sets for extraction.
    form += "#fromexternal\n"
    form += "\n"

    form += "Local expression=\n"
    form += "#fromexternal\n"

    form += ";\n\n"

    form += ".sort\n"
    form += '#fromexternal "ngroups"\n'
    form += "Format 255;\n"
    form += "\n"

    form += "#do i=0, `ngroups'\n"
    form += "\t#toexternal"  + r' "ops`i' + r"'" + r':\n"' + "\n"
    form += "\t#call getOps(ops`i')\n"
    form += "#enddo\n"

    form += "\n"
    form += "Print +ss;\n"
    form += ".end"

    with open(FORM_GENERAL_PATH / "getOps.frm", "w") as file:
        file.write(form)

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

def form_converttoSL2C():
    """
    Function designed to convert the operators in the desired SL2C notation.

    Creates converttoSL2C.frm which has to be called via a pipe proces.

    Returns
    -------

    """
    form = "Off statistics;\n"
    form += "#setexternal `PIPE1_'\n"
    form += "#prompt READY\n"
    form += "#include declarations_general.h # coefficient\n"
    form += "#include declarations_general.h # indices\n"
    form += "#include declarations_general.h # tensors\n"
    form += "#include declarations_general.h # operators\n"
    form += ".sort\n"
    # Define sets for indices
    form += f"#fromexternal\n"
    form += "\n"
    form += "Local expression = \n"
    form += f"#fromexternal\n"
    form += ";\n"
    form += ".sort\n\n"
    form += "#call makecommutative\n"
    form += ".sort\n\n"

    # form += r'#toexternal "%E\n", expression'
    # form += "\n"
    # form += r'#toexternal "READY\n"'
    # form += "\n"

    form += "* SL2C Conversion\n"
    form += "* Convert Dirac spinors into Weyl spinors\n"
    form += f"#fromexternal\n"
    form += ".sort\n\n"
    form += "* Convert Derivatives into irreps\n"
    form += f"#fromexternal\n"
    form += ".sort\n\n"
    form += "* Convert fieldstrength tensors into irreps\n"
    form += f"#fromexternal\n"
    form += ".sort\n\n"

    form += "* Replace sigma2 and sigmabar2 which are contracted with an lorentz epsilon.\n"
    form += "#call lorepsandSigma2\n"
    form += "* Replace sigma2 and sigmabar2 by sigma and sigmabar.\n"
    form += f"#fromexternal\n"
    form += ".sort\n\n"

    form += "#call replaceSigmabyEps\n"
    form += "#call simplifySL2CEps\n"
    form += ".sort\n\n"

    form += "#call makenoncommutative\n"
    form += ".sort\n\n"

    form += "Format nospaces;\n"
    form += r'#toexternal "%E\n", expression'
    form += "\n\n"
    form += ".end"

    with open(FORM_GENERAL_PATH / "converttoSL2C.frm", "w") as file:
        file.write(form)

@create_procedure(FORM_GENERAL_PATH)
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
    Note: FL denotes the total left-handed part of the fieldstrengthtensor, which originates from the commutator of
    the covariant derivatives which acting on the Higgsfield. Therefore, only the part of the total Fieldstrength
    tensor for those groups under which the Higgsfield is charged are relavant, i.e. the SU(3) for example is not
    relevant for the Higgsfield.

    Call by
        #call antisymDerivative
    after including the SL2C-indices and groups:
        #include declaration_SL2C.h
    Also the label 2 has to be declared after all EOM finding functions:
        label 2;

    Parameters
    ----------
    n_der

    Returns
    -------

    """
    form = "$found = 0;\n"
    assert n_der >= 2

    n_der = 2  # TODO: At the moment only substitutions of at most 2 derivatives are possible, since otherwise derivatives would have to be permuted, i.e. this function has to be updated and the Leibniz (product) rule has to be implemented.

    cov = op_config["fermionfields"]["D"]["mathematica"]["cov"]  # 'D'
    sl2Ceps = op_config["tensors"]["[sl2Ceps]"]["mathematica"]["sl2Ceps"]  # '[sl2Ceps]'
    higgs = op_config["bosonfields"]["H"]["mathematica"]["H"]  # 'H'
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
                lHS = " * "
                func = ""
                for i in fp_per:
                    if i in (1,2):
                        lHS += f"{cov:s}(Lsl{i:d}?Lsl, Usldot{i:d}?ULsldot[i{i:d}sldot], "
                    else:
                        lHS += f"{cov:s}(Lsl{i:d}?Lsl, Usldot{i:d}?Usldot, "
                lHS += f"{higgs}?!" + "{" + f"{cov:s}" + "}(?a)"
                lHS += der * ")"
                func += lHS
                func += " = "
                if fp_minus_one_RHS:
                    func += "(-1)*("
                # extract remaining indices, i.e. without 1 and 2 in the same order they occurred in the permutation.
                red_per = [i for i in fp_per if i not in (1,2)]
                red_der = [f"{cov:s}(Lsl{i:d}, Usldot{i:d}, " for i in red_per]
                func += f"{sl2Ceps:s}(Lsl1, Lsl2)"
                func += " * "
                func += f"{''.join(red_der)}D2({higgs}(?a)){len(red_per)*')'}"
                func += f" - "
                func += f"i_ * FL(Lsl1, Lsl2, {''.join(red_der)}{higgs}(?a){len(red_per)*')'})"
                if fp_minus_one_RHS:
                    func += ")"

                return lHS, func

            id_statement1  = f"{sl2Ceps:s}(Lsldot1?LUsldot[i1sldot], Lsldot2?LUsldot[i2sldot])"
            lHS1, id_statement_tmp1 = except_dotted_eps(per, minus_one_RHS)
            lHS1 = id_statement1 + lHS1
            id_statement1 += id_statement_tmp1

            id_statement2 = f"{sl2Ceps:s}(Lsldot2?LUsldot[i2sldot], Lsldot1?LUsldot[i1sldot])"
            lHS2, id_statement_tmp2 = except_dotted_eps(per, not minus_one_RHS)
            lHS2 = id_statement2 + lHS2
            id_statement2 += id_statement_tmp2

            id_statements.append(f"if( match( {lHS1:s} ) ) $found = 1;\n")
            id_statements.append(f"id once ifmatch -> 1 {id_statement1:s};\n")
            id_statements.append(f"if( match( {lHS2:s} ) ) $found = 1;\n")
            id_statements.append(f"id once ifmatch -> 1 {id_statement2:s};\n")
            id_statements.append("\n")

    for id in id_statements:
        form += f"{id:s}"

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
                lHS = " * "
                func = ""
                for i in fp_per:
                    if i in (1, 2):
                        lHS += f"{cov:s}(Lsl{i:d}?LUsl[i{i:d}sl], Usldot{i:d}?Usldot, "
                    else:
                        lHS += f"{cov:s}(Lsl{i:d}?Lsl, Usldot{i:d}?Usldot, "
                lHS += f"{higgs}?!" + "{" + f"{cov:s}" + "}(?a)"
                lHS += der * ")"
                func += lHS
                func += " = "
                if fp_minus_one_RHS:
                    func += "(-1)*("
                # extract remaining indices, i.e. without 1 and 2 in the same order they occurred in the permutation.
                red_per = [i for i in fp_per if i not in (1, 2)]
                red_der = [f"{cov:s}(Lsl{i:d}, Usldot{i:d}, " for i in red_per]
                func += f"{sl2Ceps:s}(Usldot1, Usldot2)"
                func += " * "
                func += f"{''.join(red_der)}D2({higgs}(?a)){len(red_per) * ')'}"
                func += f" + "
                func += f"i_ * FR(Usldot1, Usldot2, {''.join(red_der)}{higgs}(?a){len(red_per) * ')'})"
                if fp_minus_one_RHS:
                   func += ")"

                return lHS, func

            id_statement1 = f"{sl2Ceps:s}(Usl1?ULsl[i1sl], Usl2?ULsl[i2sl])"
            lHS1, id_statement_tmp1 = except_unddotted_eps(per, minus_one_RHS)
            lHS1 = id_statement1 + lHS1
            id_statement1 += id_statement_tmp1

            id_statement2 = f"{sl2Ceps:s}(Usl2?ULsl[i2sl], Usl1?ULsl[i1sl])"
            lHS2, id_statement_tmp2 = except_unddotted_eps(per, not minus_one_RHS)
            lHS2 = id_statement2 + lHS2
            id_statement2 += id_statement_tmp2

            id_statements.append(f"if( match( {lHS1:s} ) ) $found = 1;\n")
            id_statements.append(f"id once ifmatch -> 1 {id_statement1:s};\n")
            id_statements.append(f"if( match( {lHS2:s} ) ) $found = 1;\n")
            id_statements.append(f"id once ifmatch -> 1 {id_statement2:s};\n")
            id_statements.append("\n")

    for id in id_statements:
        form += f"{id:s}"

    form += "label 1;\n"

    form += "\n"
    form += "* Discard fieldstrengthtensor which are contracted with an SL2C epsilontensor, because the SL2C indices are symmetric:\n"
    form += f"id FL(Lsl1?LUsl[i1], Lsl2?LUsl[i2], ?a) * {sl2Ceps:s}(Usl1?ULsl[i1], Usl2?ULsl[i2]) = 0;\n"
    form += f"id FL(Lsl2?LUsl[i2], Lsl1?LUsl[i1], ?a) * {sl2Ceps:s}(Usl1?ULsl[i1], Usl2?ULsl[i2]) = 0;\n"
    form += f"id FR(Usldot1?ULsldot[i1], Usldot2?ULsldot[i2], ?a) * {sl2Ceps:s}(Lsldot1?LUsldot[i1], Lsldot2?LUsldot[i2]) = 0;\n"
    form += f"id FR(Usldot2?ULsldot[i2], Usldot1?ULsldot[i1], ?a) * {sl2Ceps:s}(Lsldot1?LUsldot[i1], Lsldot2?LUsldot[i2]) = 0;\n"

    form += "if($found == 1) goto 2;\n"
    # form += "label 2;\n"

    return form

@create_procedure(FORM_GENERAL_PATH)
def form_spinorEOMidentification(n_der: int):
    """
    Replace derivatives acting on spinors, where one SL2C index of the derivative and one of the spinor
    are contracted by an epsilon tensor, e.g.
    eps^{a,b}*D^{adot}_{a}*L_b
    by the equation of motion of the lepton dublett.
    In a first step, abbreviate the EOM by
    EOM(L,adot,?x),
    where "?x" abbreviate all remaining indices of L, beginning first with super and then with subscript indices.
    Second example: Derivative acting on EOM
    eps^{a,b}*D^{cdot}_{c}*D^{adot}_{a}*D^{edot}_{e}*L_b
    would be abbreviated by
    D^{cdot}_{c}*D^{edot}_{e}*EOM(L,adot,?x)

    Note: The substituted expression EOM(L,adot,?x) ignores any signs that may occur. Like for example in
        \tensor{\eps}{^{\alpha \beta }} \tensor*{\cov}{_{\alpha} ^{\dot{\alpha}}} Q_{\beta, i, a, u}
		= - \tensor*{(\slashed{\cov} q_{\mathrm{L}})}{_{i, a, u} ^{\dot{\alpha}}}
		\equiv EOM(L,adot,?x)
	so that EOM(L,adot,?x) can directly substituted by the RHS
		+ \mi \tensor*{(y^{\mathrm{u}})}{_{u v}} H^{\dagger}_{i} u_{\mathbb{C} a, v}^{\dagger \dot{\alpha}}
		+ \mi \tensor*{(y^{\mathrm{d}})}{_{u v}} H_{i} d_{\mathbb{C} a, v}^{\dagger \dot{\alpha}}
	without introducing some kind of double signs.

    Permute all indices among the derivatives and take for the field just the next one. Further, contract with the
    epsilontensor in both permutations while considering the correct sign.

    Call by
        #call spinorEOMidentification
    after including the SL2C-indices and groups:
        #include declaration_SL2C.h
    Also the label 2 has to be declared after all EOM finding functions:
        label 2;

    Parameters
    ----------
    n_der
        Number of maximum derivatives before the spinor.
    Returns
    -------

    """
    form = ""

    assert n_der >= 1
    n_der = 1  # TODO: At the moment only substitutions of at most 1 derivative is possible, since otherwise derivatives would have to be permuted, i.e. this function has to be updated and the Leibniz (product) rule has to be implemented.

    cov = op_config["fermionfields"]["D"]["mathematica"]["cov"]  # 'D'
    sl2Ceps = op_config["tensors"]["[sl2Ceps]"]["mathematica"]["sl2Ceps"]  # '[sl2Ceps]'
    l = op_config["fermionfields"]["L"]["mathematica"]["L"]
    l_dagger = op_config["fermionfields"]["L"]["mathematica"]["conj[L]"]
    q = op_config["fermionfields"]["Q"]["mathematica"]["Q"]
    q_dagger = op_config["fermionfields"]["Q"]["mathematica"]["conj[Q]"]

    d = op_config["fermionfields"]["[d_C]"]["mathematica"]["dC"]
    d_dagger = op_config["fermionfields"]["[d_C]"]["mathematica"]["conj[dC]"]
    u = op_config["fermionfields"]["[u_C]"]["mathematica"]["uC"]
    u_dagger = op_config["fermionfields"]["[u_C]"]["mathematica"]["conj[uC]"]
    e = op_config["fermionfields"]["[e_C]"]["mathematica"]["eC"]
    e_dagger = op_config["fermionfields"]["[e_C]"]["mathematica"]["conj[eC]"]

    # Replace epsilon contracted with derivative and Q, L, d, u or e by EOM:
    id_statements = []
    form += "* Replace epsilon contracted with derivative and Q, L, d, u or e by EOM:\n"
    form += "#do field = {" + f"{l}, {q}, {d}, {u}, {e}" + "}\n"
    # L(Lslspin64762,gaugeF1I1,flav6474)
    # Q(Lslspin647912,gaugeF1I1,colfF1I1,flav6477)
    # [d_C](Lslspin27682,colfF2I1,colfF2I2,flav2767)
    # [u_C](Lslspin648111,colfF1I1,colfF1I2,flav7031)
    # [e_C](Lslspin27662,flav2765)
    for der in range(1, n_der + 1):
        field_ind = der + 1
        # permutation of numbers 1,2,...,der:
        p = list(permutations(range(1, der + 1)))
        for per in p:
            def except_undotted_eps(fp_per: List, fp_minus_one_RHS: bool):
                """
                
                Parameters
                ----------
                fp_per
                    specific permutation
                fp_minus_one_RHS
                    RHS get a minus or not.
                Returns
                    List of id statements
                -------

                """
                # Derivatives
                func = " * "
                der_ind = 1
                for i in fp_per:
                    if i == der_ind:
                        #  Index on contracted derivative
                        func += f"{cov:s}(Lsl{i:d}?LUsl[i{i:d}sl], Usldot{i:d}?Usldot, "
                    else:
                        func += f"{cov:s}(Lsl{i:d}?Lsl, Usldot{i:d}?Usldot, "
                func += f"`field'(Lsl{field_ind:d}?LUsl[i{field_ind:d}sl], ?a)"
                func += der * ")"
                func += " = "
                if fp_minus_one_RHS:
                    func += "(-1)*("
                # extract remaining indices, i.e. without 1 and 2 in the same order they occurred in the permutation.
                red_per = [i for i in fp_per if i not in (1,)]
                red_der = [f"{cov:s}(Lsl{i:d}, Usldot{i:d}, " for i in red_per]
                func += f"{''.join(red_der)}EOM(`field', Usldot{der_ind:d}, ?a){len(red_per)*')'}"
                if fp_minus_one_RHS:
                    func += ")"
                return func

            id_statement1  = f"{sl2Ceps:s}(Usl1?ULsl[i1sl], Usl{field_ind:d}?ULsl[i{field_ind:d}sl])"
            id_statement1 += except_undotted_eps(per, fp_minus_one_RHS=False)
            id_statement2 = f"{sl2Ceps:s}(Usl{field_ind:d}?ULsl[i{field_ind:d}sl], Usl1?ULsl[i1sl])"
            id_statement2 += except_undotted_eps(per, fp_minus_one_RHS=True)

            id_statements.append(f"\tid once ifmatch -> 2 {id_statement1:s};\n")
            id_statements.append(f"\tid once ifmatch -> 2 {id_statement2:s};\n")
            id_statements.append("\n")

    for id in id_statements:
        form += f"{id:s}"

    form += "#enddo\n\n"

    # Replace dotted epsilon contracted with derivative and Q+, L+, d+, u+ or e+ by EOM:
    id_statements = []
    form += "* Replace dotted epsilon contracted with derivative and Q+, L+, d+, u+ or e+ by EOM:\n"
    form += "#do field = {" + f"{l_dagger}, {q_dagger}, {d_dagger}, {u_dagger}, {e_dagger}" + "}\n"
    # [Q+](Usldotspin64781,gaugeF4I1,colfF4I1,colfF4I2,flav7010)
    # [L+](Usldotspin64751,gaugeF4I1,flav6978)
    # [d_C+](Usldotspin60452,colfF4I1,flav6044)
    # [u_C+](Usldotspin64822,colfF4I1,flav6480)
    # [e_C+](Usldotspin649512,flav6493)
    for der in range(1, n_der + 1):
        field_ind = der + 1
        # permutation of numbers 1,2,...,der:
        p = list(permutations(range(1, der + 1)))
        for per in p:
            def except_dotted_eps(fp_per: List, fp_minus_one_RHS: bool):
                """

                Parameters
                ----------
                fp_per
                    specific permutation
                fp_minus_one_RHS
                    RHS get a minus or not.
                Returns
                    List of id statements
                -------

                """
                # Derivatives
                func = " * "
                der_ind = 1
                for i in fp_per:
                    if i == der_ind:
                        #  Index on contracted derivative
                        func += f"{cov:s}(Lsl{i:d}?Lsl, Usldot{i:d}?ULsldot[i{i:d}sldot], "
                    else:
                        func += f"{cov:s}(Lsl{i:d}?Lsl, Usldot{i:d}?Usldot, "
                func += f"`field'(Usldot{field_ind:d}?ULsldot[i{field_ind:d}sldot], ?a)"
                func += der * ")"
                func += " = "
                if fp_minus_one_RHS:
                    func += "(-1)*("
                # extract remaining indices, i.e. without 1 and 2 in the same order they occurred in the permutation.
                red_per = [i for i in fp_per if i not in (1,)]
                red_der = [f"{cov:s}(Lsl{i:d}, Usldot{i:d}, " for i in red_per]
                func += f"{''.join(red_der)}EOM(`field', Lsl{der_ind:d}, ?a){len(red_per) * ')'}"
                if fp_minus_one_RHS:
                    func += ")"
                return func

            id_statement1 = f"{sl2Ceps:s}(Lsldot1?LUsldot[i1sldot], Lsldot{field_ind:d}?LUsldot[i{field_ind:d}sldot])"
            id_statement1 += except_dotted_eps(per, fp_minus_one_RHS=False)
            id_statement2 = f"{sl2Ceps:s}(Lsldot{field_ind:d}?LUsldot[i{field_ind:d}sldot], Lsldot1?LUsldot[i1sldot])"
            id_statement2 += except_dotted_eps(per, fp_minus_one_RHS=True)

            id_statements.append(f"\tid once ifmatch -> 2 {id_statement1:s};\n")
            id_statements.append(f"\tid once ifmatch -> 2 {id_statement2:s};\n")
            id_statements.append("\n")

    for id in id_statements:
        form += f"{id:s}"

    form += "#enddo\n\n"

    # form += "label 2;\n"
    form += "\n"

    return form

@create_procedure(FORM_GENERAL_PATH)
def form_fieldstrengthtensorEOMidentification(n_der: int):
    """
    Replace derivatives acting on field strength tensors, where one SL2C index of the derivative and one of
    the field strength tensor are contracted by an epsilon tensor, e.g.
    eps^{a,b}*D^{adot}_{a}*BL_bc
    by the equation of motion of the field strength tensor.
    In a first step, abbreviate the EOM by
    EOM(BL,b,adot,?x),
    where "?x" abbreviate all remaining indices of BL, beginning first with super and then with subscript indices.
    Second example: Derivative acting on EOM
    eps^{a,b}*D^{cdot}_{c}*D^{adot}_{a}*D^{edot}_{e}*BL_bf
    would be abbreviated by
    D^{cdot}_{c}*D^{edot}_{e}*EOM(L,f, adot,?x)

    Permute all indices among the derivatives and take for the field just the next two. Further, contract with the
    epsilon tensor in both permutations while considering the correct sign.

    Call by
        #call fieldstrengthtensorEOMidentification
    after including the SL2C-indices and groups:
        #include declaration_SL2C.h
    Also the label 2 has to be declared after all EOM finding functions:
        label 2;

    Parameters
    ----------
    n_der
        Number of maximum derivatives before the spinor.
    Returns
    -------

    """
    form = ""

    assert n_der >= 1
    n_der = 1  # TODO: At the moment only substitutions of at most 1 derivative is possible, since otherwise derivatives would have to be permuted, i.e. this function has to be updated and the Leibniz (product) rule has to be implemented.

    cov = op_config["fermionfields"]["D"]["mathematica"]["cov"]  # 'D'
    sl2Ceps = op_config["tensors"]["[sl2Ceps]"]["mathematica"]["sl2Ceps"]  # '[sl2Ceps]'
    bl = op_config["bosonfields"]["BL"]["mathematica"]["BL"]
    br = op_config["bosonfields"]["BL"]["mathematica"]["conj[BL]"]
    wl = op_config["bosonfields"]["WL"]["mathematica"]["WL"]
    wr = op_config["bosonfields"]["WL"]["mathematica"]["conj[WL]"]

    # Replace epsilon contracted with derivative and BL or WL by EOM:
    id_statements = []
    form += "* Replace epsilon contracted with derivative and BL or WL by EOM:\n"
    form += "#do field = {" + f"{bl}, {wl}" + "}\n"
    # BL(Lsl1,Lsl2)
    # WL(gaugeF2I1,gaugeF2I2,Lsl10,Lsl11)
    for der in range(1, n_der + 1):
        field_ind1 = der + 1
        field_ind2 = der + 2
        # permutation of numbers 1,2,...,der:
        p = list(permutations(range(1, der + 1)))
        for per in p:
            def except_undotted_eps(fp_per: List, fp_minus_one_RHS: bool, fst_permutation: int):
                """

                Parameters
                ----------
                fp_per
                    specific permutation
                fp_minus_one_RHS
                    RHS get a minus or not.
                fst_permutation
                    Since the field strength tensor SL2C-indices are symmetric, those need to be written in all permutations.
                Returns
                    List of id statements
                -------

                """
                # Derivatives
                func = " * "
                der_ind = 1
                for i in fp_per:
                    if i == der_ind:
                        #  Index on contracted derivative
                        func += f"{cov:s}(Lsl{i:d}?LUsl[i{i:d}sl], Usldot{i:d}?Usldot, "
                    else:
                        func += f"{cov:s}(Lsl{i:d}?Lsl, Usldot{i:d}?Usldot, "
                if fst_permutation == 1:
                    func += f"`field'(?a, Lsl{field_ind1:d}?LUsl[i{field_ind1:d}sl], Lsl{field_ind2:d}?Lsl, ?b)"
                elif fst_permutation == 2:
                    func += f"`field'(?a, Lsl{field_ind2:d}?Lsl, Lsl{field_ind1:d}?LUsl[i{field_ind1:d}sl], ?b)"
                func += der * ")"
                func += " = "
                if fp_minus_one_RHS:
                    func += "(-1)*("
                # extract remaining indices, i.e. without 1 in the same order they occurred in the permutation.
                red_per = [i for i in fp_per if i not in (1,)]
                red_der = [f"{cov:s}(Lsl{i:d}, Usldot{i:d}, " for i in red_per]
                # TODO: Change index order, so that SL2C-indices are always the first (change in sl2c folder first.
                func += f"{''.join(red_der)}EOM(`field', Lsl{field_ind2:d}, Usldot{der_ind:d}, ?a, ?b){len(red_per) * ')'}"
                if fp_minus_one_RHS:
                    func += ")"
                return func

            id_statement1 = f"{sl2Ceps:s}(Usl1?ULsl[i1sl], Usl{field_ind1:d}?ULsl[i{field_ind1:d}sl])"
            id_statement1 += except_undotted_eps(per, fp_minus_one_RHS=False, fst_permutation=1)
            id_statement2 = f"{sl2Ceps:s}(Usl1?ULsl[i1sl], Usl{field_ind1:d}?ULsl[i{field_ind1:d}sl])"
            id_statement2 += except_undotted_eps(per, fp_minus_one_RHS=False, fst_permutation=2)
            id_statement3 = f"{sl2Ceps:s}(Usl{field_ind1:d}?ULsl[i{field_ind1:d}sl], Usl1?ULsl[i1sl])"
            id_statement3 += except_undotted_eps(per, fp_minus_one_RHS=True, fst_permutation=1)
            id_statement4 = f"{sl2Ceps:s}(Usl{field_ind1:d}?ULsl[i{field_ind1:d}sl], Usl1?ULsl[i1sl])"
            id_statement4 += except_undotted_eps(per, fp_minus_one_RHS=True, fst_permutation=2)

            id_statements.append(f"\tid once ifmatch -> 2 {id_statement1:s};\n")
            id_statements.append(f"\tid once ifmatch -> 2 {id_statement2:s};\n")
            id_statements.append(f"\tid once ifmatch -> 2 {id_statement3:s};\n")
            id_statements.append(f"\tid once ifmatch -> 2 {id_statement4:s};\n")
            id_statements.append("\n")

    for id in id_statements:
        form += f"{id:s}"

    form += "#enddo\n\n"

    # Replace dotted epsilon contracted with derivative and Q+, L+, d+, u+ or e+ by EOM:
    id_statements = []
    form += "* Replace dotted epsilon contracted with derivative and BR or WR by EOM:\n"
    form += "#do field = {" + f"{br}, {wr}" + "}\n"
    # BR(Usldot1,Usldot2)
    # WR(gaugeA10,gaugeA12,Usldot10,Usldot11)
    for der in range(1, n_der + 1):
        field_ind1 = der + 1
        field_ind2 = der + 2
        # permutation of numbers 1,2,...,der:
        p = list(permutations(range(1, der + 1)))
        for per in p:
            def except_dotted_eps(fp_per: List, fp_minus_one_RHS: bool, fst_permutation: int):
                """

                Parameters
                ----------
                fp_per
                    specific permutation
                fp_minus_one_RHS
                    RHS get a minus or not.
                fst_permutation
                    Since the field strength tensor SL2C-indices are symmetric, those need to be written in all permutations.
                Returns
                    List of id statements
                -------

                """
                # Derivatives
                func = " * "
                der_ind = 1
                for i in fp_per:
                    if i == der_ind:
                        #  Index on contracted derivative
                        func += f"{cov:s}(Lsl{i:d}?Lsl, Usldot{i:d}?ULsldot[i{i:d}sldot], "
                    else:
                        func += f"{cov:s}(Lsl{i:d}?Lsl, Usldot{i:d}?Usldot, "
                if fst_permutation == 1:
                    func += f"`field'(?a, Usldot{field_ind1:d}?ULsldot[i{field_ind1:d}sldot], Usldot{field_ind2:d}?Usldot, ?b)"
                elif fst_permutation == 2:
                    func += f"`field'(?a, Usldot{field_ind2:d}?Usldot, Usldot{field_ind1:d}?ULsldot[i{field_ind1:d}sldot], ?b)"
                func += der * ")"
                func += " = "
                if fp_minus_one_RHS:
                    func += "(-1)*("
                # extract remaining indices, i.e. without 1 and 2 in the same order they occurred in the permutation.
                red_per = [i for i in fp_per if i not in (1,)]
                red_der = [f"{cov:s}(Lsl{i:d}, Usldot{i:d}, " for i in red_per]
                func += f"{''.join(red_der)}EOM(`field', Lsl{der_ind:d}, Usldot{field_ind2:d}, ?a, ?b){len(red_per) * ')'}"
                if fp_minus_one_RHS:
                    func += ")"
                return func

            id_statement1 = f"{sl2Ceps:s}(Lsldot1?LUsldot[i1sldot], Lsldot{field_ind1:d}?LUsldot[i{field_ind1:d}sldot])"
            id_statement1 += except_dotted_eps(per, fp_minus_one_RHS=False, fst_permutation=1)
            id_statement2 = f"{sl2Ceps:s}(Lsldot1?LUsldot[i1sldot], Lsldot{field_ind1:d}?LUsldot[i{field_ind1:d}sldot])"
            id_statement2 += except_dotted_eps(per, fp_minus_one_RHS=False, fst_permutation=2)
            id_statement3 = f"{sl2Ceps:s}(Lsldot{field_ind1:d}?LUsldot[i{field_ind1:d}sldot], Lsldot1?LUsldot[i1sldot])"
            id_statement3 += except_dotted_eps(per, fp_minus_one_RHS=True, fst_permutation=1)
            id_statement4 = f"{sl2Ceps:s}(Lsldot{field_ind1:d}?LUsldot[i{field_ind1:d}sldot], Lsldot1?LUsldot[i1sldot])"
            id_statement4 += except_dotted_eps(per, fp_minus_one_RHS=True, fst_permutation=2)

            id_statements.append(f"\tid once ifmatch -> 2 {id_statement1:s};\n")
            id_statements.append(f"\tid once ifmatch -> 2 {id_statement2:s};\n")
            id_statements.append(f"\tid once ifmatch -> 2 {id_statement3:s};\n")
            id_statements.append(f"\tid once ifmatch -> 2 {id_statement4:s};\n")
            id_statements.append("\n")

    for id in id_statements:
        form += f"{id:s}"

    form += "#enddo\n\n"

    # form += "label 2;\n"
    form += "\n"

    return form

@create_procedure(FORM_GENERAL_PATH)
def form_simplifySymmetricFieldstrengthtensor():
    """
    Since the SL2C-indices of a field strength tensor are symmetric, a contraction with an epsilon tensor yields a zero.
    Returns
    -------

    """
    sl2Ceps = op_config["tensors"]["[sl2Ceps]"]["mathematica"]["sl2Ceps"]  # '[sl2Ceps]'
    bl = op_config["bosonfields"]["BL"]["mathematica"]["BL"]
    br = op_config["bosonfields"]["BL"]["mathematica"]["conj[BL]"]
    # wl = op_config["bosonfields"]["WL"]["mathematica"]["WL"]
    # wr = op_config["bosonfields"]["WL"]["mathematica"]["conj[WL]"]

    form = ""
    form += "* Discard fieldstrengthtensor which are contracted with an SL2C epsilontensor, because the SL2C indices are symmetric:\n"
    form += f"id {bl}?Field(?a, Lsl1?LUsl[i1], Lsl2?LUsl[i2], ?b) * {sl2Ceps:s}(Usl1?ULsl[i1], Usl2?ULsl[i2]) = 0;\n"  # The set Field contains all Field strength tensors, LH and RH.
    form += f"id {bl}?Field(?a, Lsl2?LUsl[i2], Lsl1?LUsl[i1], ?b) * {sl2Ceps:s}(Usl1?ULsl[i1], Usl2?ULsl[i2]) = 0;\n"
    form += f"id {br}?Field(?a, Usldot1?ULsldot[i1], Usldot2?ULsldot[i2], ?b) * {sl2Ceps:s}(Lsldot1?LUsldot[i1], Lsldot2?LUsldot[i2]) = 0;\n"
    form += f"id {br}?Field(?a, Usldot2?ULsldot[i2], Usldot1?ULsldot[i1], ?b) * {sl2Ceps:s}(Lsldot1?LUsldot[i1], Lsldot2?LUsldot[i2]) = 0;"

    return form

@create_procedure(FORM_GENERAL_PATH)
def form_lorepsandSigma2():
    # tensors
    loreps = op_config["tensors"]["[loreps]"]["mathematica"]["loreps"]
    sigma2 = op_config["tensors"]["sigma2"]["mathematica"]["sigma2"]
    sigmabar2 = op_config["tensors"]["sigmabar2"]["mathematica"]["sigmabar2"]
    # imaginary unit
    im = op_config["coefficients"]["I"]["mathematica"]["I"]

    permutations = [([1, 2, 3, 4], 1),
                    ([1, 2, 4, 3], -1),
                    ([1, 3, 2, 4], -1),
                    ([1, 3, 4, 2], 1),
                    ([1, 4, 2, 3], 1),
                    ([1, 4, 3, 2], -1),
                    ([2, 1, 3, 4], -1),
                    ([2, 1, 4, 3], 1),
                    ([2, 3, 1, 4], 1),
                    ([2, 3, 4, 1], -1),
                    ([2, 4, 1, 3], -1),
                    ([2, 4, 3, 1], 1),
                    ([3, 1, 2, 4], 1),
                    ([3, 1, 4, 2], -1),
                    ([3, 2, 1, 4], -1),
                    ([3, 2, 4, 1], 1),
                    ([3, 4, 1, 2], 1),
                    ([3, 4, 2, 1], -1),
                    ([4, 1, 2, 3], -1),
                    ([4, 1, 3, 2], 1),
                    ([4, 2, 1, 3], 1),
                    ([4, 2, 3, 1], -1),
                    ([4, 3, 1, 2], -1),
                    ([4, 3, 2, 1], 1)]


    form = "repeat;\n"
    id_statements = []
    for per, sign in permutations:
        id = f"{loreps}(lor{per[0]}?lor,lor{per[1]}?lor,lor{per[2]}?lor,lor{per[3]}?lor)"
        id += f"*{sigma2}(lor3?lor,lor4?lor,Usl1?Usl,Usl2?Usl)"
        id += " = "
        id += f"({sign:d})*(-2*{im})*{sigma2}(lor1,lor2,Usl1,Usl2)"
        id_statements.append(id)
        id = f"{loreps}(lor{per[0]}?lor,lor{per[1]}?lor,lor{per[2]}?lor,lor{per[3]}?lor)"
        id += f"*{sigmabar2}(lor3?lor,lor4?lor,Lsldot1?Lsldot,Lsldot2?Lsldot)"
        id += " = "
        id += f"({sign:d})*(2*{im})*{sigmabar2}(lor1,lor2,Lsldot1,Lsldot2)"
        id_statements.append(id)

    for id in id_statements:
        form += f"\tid once {id};\n"
    form += "endrepeat;"

    return form


@create_procedure(FORM_GENERAL_PATH)
def form_replaceSigmabyEps():
    """
    Replace contracted sigmas by SL2C epsilon tensors.
    Returns
    -------

    """
    sl2Ceps = op_config["tensors"]["[sl2Ceps]"]["mathematica"]["sl2Ceps"]
    sl2CdK = op_config["tensors"]["[sl2CdK]"]["mathematica"]["sl2CdK"]
    sigma = op_config["tensors"]["sigma"]["mathematica"]["sigma"]
    sigmabar = op_config["tensors"]["sigmabar"]["mathematica"]["sigmabar"]
    form = ""
    form += "repeat;\n"
    form += "* Replace sigmabar by sigma.\n"
    form += "\t" + f"id {sigmabar}(?a, Lsldot1?Lsldot, Usl2?Usl, ?b) = {sigma}(?a, Usl2, Lsldot1, ?b);\n"
    form += "\t" + f"id {sigmabar}(?a, Lsldot1?Lsldot, Lsl2?Lsl, ?b) = {sigma}(?a, Lsl2, Lsldot1, ?b);\n"
    form += "\t" + f"id {sigmabar}(?a, Usldot1?Usldot, Usl2?Usl, ?b) = {sigma}(?a, Usl2, Usldot1, ?b);\n"
    form += "\t" + f"id {sigmabar}(?a, Usldot1?Usldot, Lsl2?Lsl, ?b) = {sigma}(?a, Lsl2, Usldot1, ?b);\n"

    form += "\t" + f"id {sigma}(?a, Lsldot1?Lsldot, Usl2?Usl, ?b) = {sigma}(?a, Usl2, Lsldot1, ?b);\n"
    form += "\t" + f"id {sigma}(?a, Usldot1?Usldot, Lsl2?Lsl, ?b) = {sigma}(?a, Lsl2, Usldot1, ?b);\n"

    form += "* Replace in lorentz indices contracted sigmas by SL2C epsilontensors.\n"
    form += "\t" + f"id {sigma}(lor1?lor, Usl1?Usl, Lsldot1?Lsldot) * {sigma}(lor1?lor, Usl2?Usl, Lsldot2?Lsldot) = - 2 * {sl2Ceps}(Usl1, Usl2) * {sl2Ceps}(Lsldot1, Lsldot2);\n"
    form += "\t" + f"id {sigma}(lor1?lor, Usl1?Usl, Usldot1?Usldot) * {sigma}(lor1?lor, Usl2?Usl, Usldot2?Usldot) = + 2 * {sl2Ceps}(Usl1, Usl2) * {sl2Ceps}(Usldot1, Usldot2);\n"
    form += "\t" + f"id {sigma}(lor1?lor, Lsl1?Lsl, Lsldot1?Lsldot) * {sigma}(lor1?lor, Lsl2?Lsl, Lsldot2?Lsldot) = + 2 * {sl2Ceps}(Lsl1, Lsl2) * {sl2Ceps}(Lsldot1, Lsldot2);\n"
    form += "\t" + f"id {sigma}(lor1?lor, Lsl1?Lsl, Usldot1?Usldot) * {sigma}(lor1?lor, Lsl2?Lsl, Usldot2?Usldot) = - 2 * {sl2Ceps}(Lsl1, Lsl2) * {sl2Ceps}(Usldot1, Usldot2);\n"
    form += "*\n"
    form += "\t" + f"id {sigma}(lor1?lor, Usl1?Usl, Usldot1?Usldot) * {sigma}(lor1?lor, Usl2?Usl, Lsldot2?Lsldot) = - 2 * {sl2Ceps}(Usl1, Usl2) * {sl2CdK}(Usldot1, Lsldot2);\n"
    form += "\t" + f"id {sigma}(lor1?lor, Usl1?Usl, Lsldot1?Lsldot) * {sigma}(lor1?lor, Usl2?Usl, Usldot2?Usldot) = + 2 * {sl2Ceps}(Usl1, Usl2) * {sl2CdK}(Lsldot1, Usldot2);\n"
    form += "\t" + f"id {sigma}(lor1?lor, Lsl1?Lsl, Usldot1?Usldot) * {sigma}(lor1?lor, Lsl2?Lsl, Lsldot2?Lsldot) = + 2 * {sl2Ceps}(Lsl1, Lsl2) * {sl2CdK}(Usldot1, Lsldot2);\n"
    form += "\t" + f"id {sigma}(lor1?lor, Lsl1?Lsl, Lsldot1?Lsldot) * {sigma}(lor1?lor, Lsl2?Lsl, Usldot2?Usldot) = - 2 * {sl2Ceps}(Lsl1, Lsl2) * {sl2CdK}(Lsldot1, Usldot2);\n"
    form += "*\n"
    form += "\t" + f"id {sigma}(lor1?lor, Usl1?Usl, Lsldot1?Lsldot) * {sigma}(lor1?lor, Lsl2?Lsl, Lsldot2?Lsldot) = + 2 * {sl2CdK}(Usl1, Lsl2) * {sl2Ceps}(Lsldot1, Lsldot2);\n"
    form += "\t" + f"id {sigma}(lor1?lor, Usl1?Usl, Usldot1?Usldot) * {sigma}(lor1?lor, Lsl2?Lsl, Usldot2?Usldot) = - 2 * {sl2CdK}(Usl1, Lsl2) * {sl2Ceps}(Usldot1, Usldot2);\n"
    form += "\t" + f"id {sigma}(lor1?lor, Usl1?Usl, Usldot1?Usldot) * {sigma}(lor1?lor, Lsl2?Lsl, Lsldot2?Lsldot) = + 2 * {sl2CdK}(Usl1, Lsl2) * {sl2CdK}(Usldot1, Lsldot2);\n"
    form += "\t" + f"id {sigma}(lor1?lor, Usl1?Usl, Lsldot1?Lsldot) * {sigma}(lor1?lor, Lsl2?Lsl, Usldot2?Usldot) = - 2 * {sl2CdK}(Usl1, Lsl2) * {sl2CdK}(Lsldot1, Usldot2);\n"
    form += "*\n"
    form += "\t" + f"id {sigma}(lor1?lor, Lsl1?Lsl, Lsldot1?Lsldot) * {sigma}(lor1?lor, Usl2?Usl, Lsldot2?Lsldot) = - 2 * {sl2CdK}(Lsl1, Usl2) * {sl2Ceps}(Lsldot1, Lsldot2);\n"
    form += "\t" + f"id {sigma}(lor1?lor, Lsl1?Lsl, Usldot1?Usldot) * {sigma}(lor1?lor, Usl2?Usl, Usldot2?Usldot) = + 2 * {sl2CdK}(Lsl1, Usl2) * {sl2Ceps}(Usldot1, Usldot2);\n"
    form += "\t" + f"id {sigma}(lor1?lor, Lsl1?Lsl, Usldot1?Usldot) * {sigma}(lor1?lor, Usl2?Usl, Lsldot2?Lsldot) = - 2 * {sl2CdK}(Lsl1, Usl2) * {sl2CdK}(Usldot1, Lsldot2);\n"
    form += "\t" + f"id {sigma}(lor1?lor, Lsl1?Lsl, Lsldot1?Lsldot) * {sigma}(lor1?lor, Usl2?Usl, Usldot2?Usldot) = + 2 * {sl2CdK}(Lsl1, Usl2) * {sl2CdK}(Lsldot1, Usldot2);\n"
    form += "endrepeat;\n"

    return form

@create_procedure(FORM_GENERAL_PATH)
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

    # form += "* Replace epsilons with one upper and one lower index by Kronecker-deltas:\n"
    # form += "\t" + f"id {sl2Ceps}(Lsl1?Lsl,Usl2?Usl) = {sl2CdK}(Lsl1,Usl2);\n"
    # form += "\t" + f"id {sl2Ceps}(Usl1?Usl,Lsl2?Lsl) = {sl2CdK}(Lsl2,Usl1);\n"
    # form += "\t" + f"id {sl2Ceps}(Lsldot1?Lsldot,Usldot2?Usldot) = {sl2CdK}(Lsldot1,Usldot2);\n"
    # form += "\t" + f"id {sl2Ceps}(Usldot1?Usldot,Lsldot2?Lsldot) = {sl2CdK}(Lsldot2,Usldot1);\n"
    form += "* Simplify only Kronecker-deltas which are contracted with eps and Kronecker-deltas, since contractions inside one building block are not wanted:\n"  # Replace contracted Kronecker-deltas (works for every function, not just eps):\n
    form += "\t" + f"id {sl2Ceps}?" + "{" + f"{sl2Ceps},{sl2CdK},{sigma}" + "}" + f"(?a,Lsl1?LUsl[k],?b)*{sl2CdK}(?c,Usl1?ULsl[k],?d) = {sl2Ceps}(?a,?c,?d,?b);\n"  # [sl2Ceps]?
    form += "\t" + f"id {sl2Ceps}?" + "{" + f"{sl2Ceps},{sl2CdK},{sigma}" + "}" + f"(?a,Lsldot1?LUsldot[k],?b)*{sl2CdK}(?c,Usldot1?ULsldot[k],?d) = {sl2Ceps}(?a,?c,?d,?b);\n"  # [sl2Ceps]?
    form += "* Replace self-contracted Kronecker-deltas by the dimension (=2):\n"
    form += "\t" + f"id {sl2CdK}(Usl1?ULsl[k],Lsl1?LUsl[k]) = d_(Lsl1,Lsl1);\n"
    form += "\t" + f"id {sl2CdK}(Usldot1?ULsldot[k],Lsldot1?LUsldot[k]) = d_(Lsldot1,Lsldot1);\n"
    form += "endrepeat;\n"
    form += "* Bring indices of epsilons in order:\n"
    form += f"Multiply replace_({sl2Ceps},{sl2CepsA});\n"
    form += ".sort\n"
    form += f"Multiply replace_({sl2CepsA},{sl2Ceps});\n"

    return form

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
@create_procedure(FORM_GENERAL_PATH)
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
    # form += "#procedure simplifyEpsSU2\n"
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

    return form

    # form += "#endprocedure\n"
    # with open(FORM_GENERAL_PATH / f"simplifyEpsSU2.prc", "w") as file:
    #     file.write(form)

@create_procedure(FORM_GENERAL_PATH)
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
    # form += "#procedure simplifyEpsSU3\n"
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

    return form

    # form += "#endprocedure"
    # with open(FORM_GENERAL_PATH / f"simplifyEpsSU3.prc", "w") as file:
    #     file.write(form)

@create_procedure(FORM_GENERAL_PATH)
def form_makecommutative(n_der: int=4):
    """
    FORM function which converts every non commuting Funtion and derivatives of it into a Placeholder and a commuting
    Function, like:
        D(lor1,lbar(spin1,gauge1,flav1)) -> D(lor1,lbar(op1))*lbarc(op1,spin1,gauge1,flav1)
    Parameters
    ----------
    n_der
        Maximum number of derivatives occuring.
    Returns
    -------

    """
    form = ".sort\n"
    form += "Function tmp, Field;\n"
    form += "Multiply left tmp;\n\n"
    form += "#$i = 1;\n"
    def derivative(nD, withset=True):
        return ''.join([f"D(lor{i}{'?lor' if withset else ''}," for i in range(1, nD + 1)])

    while_conditions = [f"(match(tmp*{derivative(nD)}Field?AllFields(?a){nD * ')'}))" for nD in range(n_der + 1)]
    id_left = [f"tmp*{derivative(nD)}Field?AllFields[k?](?a){nD * ')'}" for nD in range(n_der + 1)]
    id_right = [f"{derivative(nD, False)}Field(op[$i]){nD * ')'}*AllcFields[k](op[$i],?a)*tmp" for nD in range(n_der + 1)]
    form += f"while ( {' || '.join(while_conditions)} );\n"
    # (match(tmp*Field?AllFields(?a))) || (match( tmp*D(lor1?lor,Field?AllFields(?a)) ))
    for nD, while_condition in enumerate(while_conditions):
        if nD == 0:
            form += "\tif "
        else:
            form += "\telseif "
        form += f"{while_condition};\n"
        form += f"\t\tid once {id_left[nD]} = {id_right[nD]};\n"
        form += "\t\t$i = $i + 1;\n"

    form += "\tendif;\n"
    form += "endwhile;\n\n"
    form += "id tmp = 1;"

    return form

@create_procedure(FORM_GENERAL_PATH)
def form_makenoncommutative(sets:int, n_der:int=4, sl2c=True):
    """
    FORM function which reastablishes the non-commuting Funtions (s. form_makecommutative()):
        D(lor1,lbar(op1))*[L+c](op1,Lsldot1,gauge1,flav1) -> D(lor1,[L+](Lsldot1,gauge1,flav1))
    Parameters
    ----------
    n_der
        Maximum number of derivatives occuring.
    Returns
    -------

    """
    form = ".sort\n"
    form += "Function Field;\n"
    form += "CFunction CField;\n\n"

    def derivative(nD, sl2c, withset=True):
        if sl2c:
            return ''.join([f"D(Lsl{i}{'?Lsl' if withset else ''},Usldot{i}{'?Usldot' if withset else ''}," for i in range(1, nD + 1)])
        else:
            return ''.join([f"D(lor{i}{'?lor' if withset else ''}," for i in range(1, nD + 1)])

    form += "repeat;\n"
    for i in range(sets):
        for nD in range(n_der + 1):
            form += "\tid once "
            form += f"{derivative(nD, sl2c)}Field?AllFields{i}[k](op1?op[m]){nD * ')'}"
            form += "*"
            form += f"CField?AllcconFields{i}[k](op1?op[m],?a)"
            form += " = "
            form += f"{derivative(nD,sl2c,False)}AllconFields{i}[k](?a){nD*')'};\n"
    form += "endrepeat;"

    return form

@create_procedure(FORM_GENERAL_PATH)
def form_derivativeasIndex(n_der: int=4):
    """
    Write derivative which act on any field as indices of that field.
    E.g.: D(a, adot, H(gauge1)) -> H(a,adot,gauge)

    This is done to sort the fields easily by their helicity.

    Parameters
    ----------
    n_der

    Returns
    -------

    """
    cov = op_config["fermionfields"]["D"]["mathematica"]["cov"]  # 'D'
    l = op_config["fermionfields"]["L"]["mathematica"]["L"]  # 'L'
    def derivativestoInd(argumentLeft, n):
        """
        Write for example: id D(?a, D(?b, L?(?x))) = L(?a, ?b, ?x);
        Note: 'x', 'y', 'z' are free to use them for internal field indices.
        Parameters
        ----------
        argumentLeft
        argumentRight
        n

        Returns
        -------
        lHS
        rHS_indices
        """
        abc = "abcdefghijklmnopqrstuvw"
        lHS = ""
        rHS_indices = ""
        for a in abc[:n]:
            lHS += f"{cov}(?{a:s}, "
        lHS += argumentLeft
        lHS += n * ")"

        for a in abc[:n]:
            rHS_indices += f"?{a:s}, "
        return lHS, rHS_indices

    form = "repeat;\n"
    for i in range(1, n_der + 1):
        lHS, rHS_Indices = derivativestoInd(f"{l}?(?x)", n=i)
        form += f"\tid {lHS} = {l}({rHS_Indices}?x);\n"
    form += "endrepeat;"
    return form

@create_procedure(FORM_GENERAL_PATH)
def form_indexasDerivative(n_der: int=4):
    """
    Since, derivatives are the only objects that have always the index structure [Lsl, Usldot], it is always uniquely
    possible to extract the derivatives of a fields by the indices:
    E.g.: H(a,adot,gauge) -> D(a, adot, H(gauge1))
    """
    cov = op_config["fermionfields"]["D"]["mathematica"]["cov"]  # 'D'
    l = op_config["fermionfields"]["L"]["mathematica"]["L"]  # 'L'
    def indtoDerivatives(n):
        """
        'x', 'y', 'z' are free to use them for internal field indices.
        Parameters
        ----------
        n

        Returns
        -------
        lHS_indices
        rHS
        brackets
        """
        abc = "abcdefghijklmnopqrstuvw"
        lHS_indices = []
        rHS = []
        for i in range(1, n+1):
            lHS_indices.append(f"Lsl{i:d}?Lsl, Usldot{i:d}?Usldot")
            rHS.append(f"{cov}(Lsl{i:d}, Usldot{i:d}")
        rHS = ", ".join(rHS)
        rHS += ", "
        brackets = n * ")"

        lHS_indices = ", ".join(lHS_indices)

        return lHS_indices, rHS, brackets

    form = "repeat;\n"
    for i in reversed(range(1, n_der + 1)):
        lHS_Indices, rHS, brackets = indtoDerivatives(n=i)
        form += f"\tid {l}?!" + "{" + f"{cov}" + "}" + f"({lHS_Indices}, ?x) = {rHS}{l}(?x){brackets};\n"
    form += "endrepeat;"

    return form

@create_procedure(FORM_GENERAL_PATH)
def form_sortfields(order: Dict):
    """

    Parameters
    ----------
    order : Dict[Field]

    Returns
    -------

    """
    form = ""
    form += "repeat;\n"
    ordered_fields = list(order.values())[::-1]  # Reversed list
    for i in range(len(ordered_fields)):
        for j in range(len(ordered_fields)):
            if j > i:
                if ordered_fields[i].ac and ordered_fields[j].ac:
                    sign = "-"
                else:
                    sign = "+"
                form += f"\tid {ordered_fields[i].form_name}(?a)*{ordered_fields[j].form_name}(?b) = {sign}{ordered_fields[j].form_name}(?b)*{ordered_fields[i].form_name}(?a);\n"
        form += "\n"
    form += "endrepeat;"
    return form

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
    tensors_without_Kronecker = [tensor for tensor in tensors if "dK" not in tensor]
    form += "CFunction " + ", ".join(tensors_without_Kronecker) + ";\n"

    form += "* Auxiliary antisymmetric epsilons, used in combination with replace_.\n"
    eps = [list(tensor["mathematica"].values())[0] for tensor_name, tensor in op_config["tensors"].items() if
           "eps" in tensor_name]
    form += f"CFunction {', '.join(map(lambda text: get_antisymEps(text) + '(antisymmetric)', eps))};\n"  # sl2CepsA(antisymmetric), su2epsA(antisymmetric), su3epsA(antisymmetric)
    form += "\n"
    form += "* Declare Kronecker Delta symbol for Sl2C Indices, because built in can not handle upper and lower (un-)dottet indices.\n"
    form += "* Since two indices are also symmetric when they are cyclic and vice versa, and pattern matching is not allowed for symmetric function but for cyclic, [sl2CdK] is declared as cyclic.\n"
    dK = [list(tensor["mathematica"].values())[0] for tensor_name, tensor in op_config["tensors"].items() if
          "dK" in tensor_name]
    form += f"CFunction {', '.join(map(lambda text: text + '(cyclic)', dK))};\n"  # sl2CdK(cyclic), su2dK(cyclic), su3dK(cyclic)
    form += "*--#] tensors :\n"

    #
    form += "\n"
    form += "*--#[ NCtensors :\n"
    form += "Function " + ", ".join(tensors) + ";\n"
    form += "* Declare Kronecker Delta symbol for Sl2C Indices, because built in can not handle upper and lower (un-)dottet indices.\n"
    form += "* Since two indices are also symmetric when they are cyclic and vice versa, and pattern matching is not allowed for symmetric function but for cyclic, [sl2CdK] is declared as cyclic.\n"
    dK = [list(tensor["mathematica"].values())[0] for tensor_name, tensor in op_config["tensors"].items() if
          "dK" in tensor_name]
    form += f"Function {', '.join(map(lambda text: text + '(cyclic)', dK))};\n"  # sl2CdK(cyclic), su2dK(cyclic), su3dK(cyclic)
    form += "*--#] NCtensors :\n"

    form += "\n"
    form += "*--#[ coefficient :\n"
    # Coefficient
    form += f"Symbols {', '.join(coeff + ['n', 'N'])};\n" # Symbols d, eps, lambdah, At, g1, g2, g3, mu, lambdaphi, kappa, Ms, Mu, muM, [2L[Ms,muM]], n;

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
    dimlessConst = form_factorizeCoeff()
    form += f"Symbols {', '.join(dimlessConst):s};\n"  # Symbols [At/Ms], [mu/Ms], [Mu/Ms], [muM/Ms];
    form += "*--#] coefficient :\n"
    form += "\n"
    form += "*--#[ operators :\n"
    # write commuting and anti-commuting operators of each term in separate list for initialization. Thus
    # duplicated operators are removed.

    def commutative_fields(ops):
        non_singlets = {}
        singlets = {}
        for name, field in {**ops["bosonfields"], **ops["fermionfields"]}.items():
            if "helicity" not in field.keys():
                if "derivative" not in field["description"]:
                    non_singlets[name] = field
            else:
                if field["helicity"] == 0:
                    non_singlets[name] = field
                    singlets[name] = field
                else:
                    singlets[name] = field

        allFields = []
        allconFields = []
        for singlet_field in singlets.values():
            if singlet_field["helicity"] == 0:
                allFields += singlet_field["mathematica"].values()
                allconFields += singlet_field["mathematica"].values()

            elif singlet_field["helicity"] == Fraction(-1,2):
                non_singlet_field = non_singlets[singlet_field["dirac"]]
                autoeft_form = {autoeft: form for form, autoeft in singlet_field["autoeft"].items()}
                for dirac, weyl in non_singlet_field["autoeft"].items():
                    allFields.append(dirac)
                    allconFields.append(autoeft_form[weyl])

            elif singlet_field["helicity"] == -1:
                non_singlet_field = list(non_singlets[singlet_field["non_singlet"]]["mathematica"].values())[0]
                for sfield in singlet_field["mathematica"].values():
                    allFields.append(non_singlet_field)
                    allconFields.append(sfield)

        allcFields = list(map(get_commuting_op,allFields))

        allcconFields = list(map(get_commuting_op, allconFields))

        return allFields, allcFields, allcconFields, allconFields

    allFields, allcFields, allcconFields, allconFields = commutative_fields(op_config)
    # Since FORM is confused when something appears twice in a Set, multiple Sets are needed and the id-statements
    # referring to those just need to be copied.
    allFieldsSet, allcFieldsSet, allcconFieldsSet, allconFieldsSet = [[]], [[]], [[]], [[]]
    for aF, acF, acconF, aconF in zip(allFields, allcFields, allcconFields, allconFields):
        if (aF not in allFieldsSet[0]) and (acF not in allcFieldsSet[0]) and (acconF not in allcconFieldsSet[0]) and (aconF not in allconFieldsSet[0]):
            allFieldsSet[0].append(aF)
            allcFieldsSet[0].append(acF)
            allcconFieldsSet[0].append(acconF)
            allconFieldsSet[0].append(aconF)
        else:
            n = 0
            while (aF in allFieldsSet[n]) or (acF in allcFieldsSet[n]) or (acconF in allcconFieldsSet[n]) or (aconF in allconFieldsSet[n]):
                n += 1
                if len(allFieldsSet) < n+1 or len(allcFieldsSet) < n+1 or len(allcconFieldsSet) < n+1 or len(allconFieldsSet) < n+1:
                    allFieldsSet.append([])
                    allcFieldsSet.append([])
                    allcconFieldsSet.append([])
                    allconFieldsSet.append([])
            assert ((aF not in allFieldsSet[n]) and (acF not in allcFieldsSet[n]) and (acconF not in allcconFieldsSet[n]) and (aconF not in allconFieldsSet[n]))
            allFieldsSet[n].append(aF)
            allcFieldsSet[n].append(acF)
            allcconFieldsSet[n].append(acconF)
            allconFieldsSet[n].append(aconF)


    form += "* Define all fields and auxiliary fields:\n"
    form += f"Function {', '.join(set(allFields))};\n"
    form += f"CFunction {', '.join(set(allcFields))};\n"
    form += f"CFunction {', '.join(set(allcconFields))};\n"
    form += f"Function {', '.join(set(allconFields))};\n\n"
    form += "* Derivative:\n"
    form += f"Function {op_config['fermionfields']['D']['mathematica']['cov']};\n"

    # form += "Function " + ", ".join(bosons) + ";\n"
    # Auxiliary commuting bosons
    # form += "CFunction " + ", ".join(map(get_commuting_op, bosons)) + ";\n"
    # form += "\n"
    # form += "Function " + ", ".join(fermions) + ";\n"

    # Auxiliary commuting fermions
    # form += "CFunction " + ", ".join(map(get_commuting_op, fermionfields)) + ";\n"
    # form += "\n"

    SL2C_fieldstrengths = [form_field for field in op_config["bosonfields"].values() for form_field in field["mathematica"].values() if "helicity" in field.keys() if field["helicity"] == -1]
    form += f"Set fieldstrengthsc: {', '.join(map(get_commuting_op, SL2C_fieldstrengths))};\n"
    form += f"Set fieldstrengths: {', '.join(SL2C_fieldstrengths)};\n"
    form += "\n"
    form += "CFunction xi1, [xi1+], chi1, [chi1+], xi2, [xi2+], chi2, [chi2+];\n"
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
    form += "* List of all possibly occurring fields bevor they are converted into Lorentz irreps\n"
    form += "* Note: If a field (like e.g. a fieldstrength tensor B) is converted into two different Lorentz irreps (BL, BR) is has to occur twice in this list\n"
    form += f"Set AllFields: {', '.join(allFields)};\n"
    for i, allFs in enumerate(allFieldsSet):
        form += f"Set AllFields{i}: {', '.join(allFs)};\n"
    form += "\n"
    form += "* List of the same fields as in AllFields (in the same order!) but defined as a commuting Function\n"
    form += f"Set AllcFields: {', '.join(allcFields)};\n"
    for i, allcFs in enumerate(allcFieldsSet):
        form += f"Set AllcFields{i}: {', '.join(allcFs)};\n"
    form += "\n"
    form += "* List of the same fields (first as commutative fields) as in AllFields (in the same order!) but with possible replacements like Dirac to Weyl spinors and so on\n"
    form += f"Set AllcconFields: {', '.join(allcconFields)};\n"
    for i, allcconFs in enumerate(allcconFieldsSet):
        form += f"Set AllcconFields{i}: {', '.join(allcconFs)};\n"
    form += "\n"
    form += "* Now as noncommutative fields\n"
    form += f"Set AllconFields: {', '.join(allconFields)};\n"
    for i, allconFs in enumerate(allconFieldsSet):
        form += f"Set AllconFields{i}: {', '.join(allconFs)};\n"
    form += "\n"
    form += "* Set for convenient insertion of op1, op2, ... indices\n"
    form += "Set op: op1,...,op100;\n"
    form += "\n"
    form += "* D2 = D_mu * D^mu:\n"
    form += "Function D2;\n"
    form += "* Intern abbreviation for equation of motion:\n"
    form += "Function EOM;\n"
    form += "* Total field strength tensor, necessary for EOM substitutions:\n"
    form += "Function FL, FR;\n"
    form += "\n"
    form += "*--#] operators :\n"
    form += "\n"
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
    form += "Autodeclare Symbols i, j;\n"
    form += "*--#] indices :\n"  # trailing "\n" important otherwise form will not find the "fold" declarations

    form_getTerms()
    form_getOps()
    form_coefficient_handling()
    form_converttoSL2C()

    form_makecommutative(n_der)
    form_makenoncommutative(len(allFieldsSet), n_der)

    form_lorepsandSigma2()
    form_replaceSigmabyEps()
    form_simplifySL2CEps()
    # form_replaceSUNGenerators(N=2)
    # form_replaceSUNGenerators(N=3)
    form_simplifyEpsSU2()
    form_simplifyEpsSU3()
    form_antisymDerivative(n_der)
    form_spinorEOMidentification(n_der)
    form_fieldstrengthtensorEOMidentification(n_der)
    form_simplifySymmetricFieldstrengthtensor()
    form_derivativeasIndex(n_der)
    form_indexasDerivative(n_der)

    # TODO: Create sorted field new from own model file.
    form_sortfields(fields_sorted)


    return form

def declaration_SL2C_sets(indices, n_der) -> str:
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

    max_ind = n_der + 4
    # Generate only indices if there aren't any, since they could be declared twice otherwise.
    def id_indices(typ: str, reference: List, max_index: int):
        return [f"{typ:s}{i:d}" for i in range(max_index) if f"{typ:s}{i:d}" not in [repr(index) for index in reference]]
    # if not (lsl and usl):
    lsl_aux    = id_indices("Lsl", lsl, max_ind)
    usl_aux    = id_indices("Usl", usl, max_ind) #[f"Usl{i:d}" for i in range(max_ind)]
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

    set_lsl = ", ".join([f"{index}" for index in lsl_aux] + [f"{repr(index)}" for index in lsl])
    set_usl = ", ".join([f"{index}" for index in usl_aux] + [f"{repr(index)}" for index in usl])
    set_lsldot = ", ".join([f"{index}" for index in lsldot_aux] + [f"{repr(index)}" for index in lsldot])
    set_usldot = ", ".join([f"{index}" for index in usldot_aux] + [f"{repr(index)}" for index in usldot])

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
    form += f"Set Lsl: {set_lsl};\n"
    form += f"Set Usl: {set_usl};\n"
    form += f"Set Lsldot: {set_lsldot};\n"
    form += f"Set Usldot: {set_usldot};\n"

    form += f"Set LUsl: {set_lusl};\n"
    form += f"Set ULsl: {set_ulsl};\n"
    form += f"Set LUsldot: {set_lusldot};\n"
    form += f"Set ULsldot: {set_ulsldot};\n"
    form += "\n"

    form += "* Auxiliary sets\n"
    # Sets for contraction for all other indices except SL2C-indices:
    for index_name in index_config.keys():
        if index_name not in ["Lsl", "Usl", "Lsldot", "Usldot"]:
            indices_typ = indices[index_name]
            ind_aux = id_indices(index_name, indices_typ, max_ind)
            form += f"Set {index_name}: {', '.join(ind_aux + list(map(str, indices_typ)))};\n"

    return form

