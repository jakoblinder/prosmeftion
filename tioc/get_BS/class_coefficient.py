import re
import logging
import sys
import subprocess

from pathlib import Path

from tioc import coeffname, abbreviation, coeffvalues, escape_regex, PROJECTION_PATH, FORM_PATH
logger_autoeft = logging.getLogger("autoeft.projection")
logger = logger_autoeft.getChild("coefficient")

class Coefficient():
    """
    A class describing the coefficient in a term.
    Example coefficient (of term 18):
        +At*(-72*(4+3*eps)*kappa*Ms^4*mu+36*At^2*mu*(-3*eps*Ms^2+eps^2*Ms^2-8*Mu^2-8*eps*Mu^2)
         -36*At^3*((-4+5*eps+3*eps^2)*Ms^2-8*(2+eps)*Mu^2)
         +At*(36*Ms^2*((4+9*eps)*lambdah*Ms^2-2*(2+eps)*lambdaphi*Ms^2+4*(1+eps)*mu^2)
         -g2^2*Ms^2*((168+58*eps+d*(12+11*eps))*Ms^2+72*(2+eps)*Mu^2)
         -36*kappa*((-16-5*eps+eps^2)*Ms^4+2*(-4+eps)*Ms^2*Mu^2-8*(2+eps)*Mu^4))
         +6*eps*(24*kappa*Ms^4*mu-3*At^2*mu*(eps*Ms^2-8*Mu^2)+12*At^3*((-1+eps)*Ms^2-4*Mu^2)
         +At*(g2^2*((14+d)*Ms^4+12*Ms^2*Mu^2)+3*(-4*Ms^2*(lambdah*Ms^2-lambdaphi*Ms^2+mu^2)
         +kappa*((-16+eps)*Ms^4-8*Ms^2*Mu^2-16*Mu^4))))*L[Ms^2,\[Mu]^2])/(72*eps*Ms^8)
         with L[M1,M2] := Log[M1/M2].
    Attributes
    ----------
    expression: str
        Explicit expression of the operator.
    Magicmethods
    ------------
    __str__()
        Specify the format for printing an Operator object with the print command or converting a Operator object in
        string with str() function.
    __format__()
        Specify the format for "format" function in print statement: Here the same as the print statement itself.
    Methods
    -------
    replace()
        Rewrites all constants and the logarithm occurring in the coefficient.
    """
    expression : str
    def __init__(self, expression, name=None, first_read_in = True):
        self.expression = expression
        if first_read_in:
            self.name = name
            self.replace() # rewrites expression in FORM compatible way
            self.dim_operator = self.get_dim()
        self.tex_expression = self.get_tex()

    def __str__(self):
        """Specify the format for printing an Operator object with the print command or converting a Operator object in
        string with str() function. """
        return f"{self.expression:s}"

    def __repr__(self):
        return self.__str__()

    def __format__(self, key):
        """Specify the format for "format" function in print statement: Here the same as the print statement itself."""
        return self.__str__()

    def __mul__(self, other):
        """
        Multiply a string from the right to the expression of the Coefficient.
        Parameters
        ----------
        other: str
            Example: other = '(1/2)'
        Returns
        -------
        """
        self.expression = self.expression + f"*{other:s}"
        return self

    def replace(self):
        """
        Rewrite coefficient of each term in such a way that it is compatible with FORM. This means essentially that all
        constants occurring in the term are replaced by the expressions written in "coeffname".
        Returns
        -------
        """
        # Since both 'Mu' and '\[Mu]' occur in the output, the latter is substituted first to remove any ambiguities
        # in the substitution.
        for key in coeffname.keys():
            if key != "Mu":
                self.expression = re.sub(escape_regex(key), coeffname[key]["FORM"], self.expression)
        self.expression = re.sub(escape_regex("Mu"), coeffname["Mu"]["FORM"], self.expression)

        # Rewrite Logarithms L[Ms^2,\[Mu]^2] as specified by lg = "[2L[Ms,muM]]" in coefficients.yml:
        muM = coeffname['\\[Mu]']["FORM"]
        Ms = coeffname['Ms']["FORM"]
        self.expression = re.sub(escape_regex(f"L[{Ms:s}^2,{muM:s}^2]"), abbreviation["lg"]["FORM"], self.expression)

    def get_dim(self):
        """
        Determines the dimension of the operator by investigation of the coefficient.
        Parameters
        ----------
        dimensionalConstants: Dict{str: int}
            Dictionary of
        Returns
        -------

        """
        TERM_PATH = FORM_PATH / self.name
        TERM_PATH.mkdir(parents=True, exist_ok=True)  # Create directories if they don't exist.
        form = ""
        symbols = "Symbols " + ", ".join(coeffvalues + ["n"]) + ";\n"
        form += symbols  # Symbols d, lambdah, At, g1, g2, g3, [c_mu], lambdaphi, kappa, I, [m_{Ms}], [m_{Mu}], [Mu], eps, n, [2L[Ms,Mu]];
        form += "\n"
        form += f"Local coefficient = {self.expression:s};\n"
        form += ".sort\n"
        dimlessConst = []
        idStatements = []
        for key, constant in coeffname.items():
            if key == "I" or key == "Ms":
                continue
            if constant["massdim"] != 0:
                if constant["massdim"] == 1:
                    dimlessC = f"[{constant['FORM']:s}/{coeffname['Ms']['FORM']:s}]"
                    dimlessConst.append(dimlessC)
                    idStatements.append(f"id {constant['FORM']:s} = {dimlessC:s}*{coeffname['Ms']['FORM']:s}")
                else:
                    dimlessC = f"[{constant['FORM']:s}/{coeffname['Ms']['FORM']:s}^{constant['massdim']:d}]"
                    dimlessConst.append(dimlessC)
                    idStatements.append(f"id {constant['FORM']:s} = {dimlessC:s}*{coeffname['Ms']['FORM']:s}^{constant['massdim']:d}")
        form += f"Symbols {', '.join(dimlessConst):s};\n"
        form_factorizeCoeff = "#procedure factorizeCoeff\n"
        for idStatement in idStatements:
            form_factorizeCoeff += f"\t{idStatement:s};\n"
        form_factorizeCoeff += "#endprocedure\n"
        with open(TERM_PATH / "factorizeCoeff.prc", "w") as file:
            file.write(form_factorizeCoeff)

        form += "\n"
        form += "#call factorizeCoeff\n"
        form += "\n"
        form += "id Ms^n?$dim = Ms^n;\n"
        form += "$dim = 4 - $dim;\n"
        form += "\n"
        form += ".sort\n"
        form += f"Bracket {coeffname['Ms']['FORM']:s};\n"
        FILE_PATH = TERM_PATH / Path("dimension.aux")
        form += f'#write <{str(FILE_PATH):s}> "dim = %$", $dim\n'
        form += "\n"
        form += "Print +s;\n"
        form += ".end"

        with open(TERM_PATH / "get_dimension.frm", "w") as file:
            file.write(form)

        try:
            formprocess = subprocess.run(["form", "-p", TERM_PATH, TERM_PATH / "get_dimension.frm"],
                                          capture_output=True, text=True, check=True,)
        except subprocess.CalledProcessError as exc:
            exc.cmd = list(map(str, list(exc.cmd)))
            logger.error(" ".join(exc.cmd) + "\n" + str(exc.stdout))
            logger.error(f"form returned non-zero exit status {exc.returncode}")
            sys.exit("STOP")
        else:
            # No Error occured
            logger.debug(" ".join(map(str, formprocess.args)))
            output = formprocess.stdout
            with open(FILE_PATH, "r") as file:
                match = re.match(r"dim = (?P<dim>\d{1,2})",file.read())
                if match:
                    dim = int(match.group("dim"))
                    if dim < 4 and dim > 0:
                        logger.warning(f"{self.name:s} has a mass dimension of only {dim}.")
                    elif dim == 0:
                        logger.info(f"{self.name:s} has a mass dimension zero.")
                    elif dim < 0:
                        logger.error(f"{self.name:s} has a negativ mass dimension: {dim}.")
                        sys.exit("STOP")
                else:
                    logger.error("Format of the dimension.aux not as expected.")
                    sys.exit("STOP")
            # print(f"{self.name}: {dim}")
            # write the newly formatted expression in the expression attribute.
            output = re.sub(r"(\s)*", "", output)  # Remove all whitespaces and newlines: \s = [\t\n\r\f\v]
            output = re.sub(r"\\", "",
                            output)  # Sometimes FORM splits indices in long expressions with an backslash "\" which is discarded.
            pattern = r"(?<!Local)coefficient=(?P<coefficient>.*);"
            match = re.search(pattern, output)
            if match:
                self.expression = match.group("coefficient")
            else:
                logger.error(f"No coefficient has been found in {output}.")
                sys.exit("STOP")
            return dim

    def update_coefficient(self, yukawa = False):
        """
        Write expression again inside of Form to include appended Coefficients and so on.
        Does not include Yukawa matrices (TODO)
        Parameters
        ----------
        yukawa

        Returns
        -------

        """
        TERM_PATH = FORM_PATH / self.name
        TERM_PATH.mkdir(parents=True, exist_ok=True)  # Create directories if they don't exist.
        form = ""
        symbols = "Symbols " + ", ".join(coeffvalues + ["n"]) + ";\n"
        form += symbols  # Symbols d, lambdah, At, g1, g2, g3, [c_mu], lambdaphi, kappa, I, [m_{Ms}], [m_{Mu}], [Mu], eps, n, [2L[Ms,Mu]];

        dimlessConst = []
        for key, constant in coeffname.items():
            if key == "I" or key == "Ms":
                continue
            if constant["massdim"] != 0:
                if constant["massdim"] == 1:
                    dimlessConst.append(f"[{constant['FORM']:s}/{coeffname['Ms']['FORM']:s}]")
                else:
                    dimlessConst.append(f"[{constant['FORM']:s}/{coeffname['Ms']['FORM']:s}^{constant['massdim']:d}]")
        form += f"Symbols {', '.join(dimlessConst):s};\n"

        form += "\n"
        form += f"Local coefficient = {self.expression:s};\n"
        form += ".sort\n"
        form += "\n"
        form += "#call factorizeCoeff\n"
        form += "\n"
        form += ".sort\n"
        form += f"Bracket {coeffname['Ms']['FORM']:s};\n"
        form += "\n"
        form += "Print +s;\n"
        form += ".end"
        with open(TERM_PATH / "update_coefficient.frm", "w") as file:
            file.write(form)

        try:
            formprocess = subprocess.run(["form", "-p", TERM_PATH, TERM_PATH / "update_coefficient.frm"],
                                         capture_output=True, text=True, check=True)
        except subprocess.CalledProcessError as exc:
            exc.cmd = list(map(str, list(exc.cmd)))
            logger.error(" ".join(exc.cmd) + "\n" + str(exc.stdout))
            logger.error(f"form returned non-zero exit status {exc.returncode}")
            sys.exit("STOP")
        else:
            # No Error occured
            logger.debug(" ".join(map(str, formprocess.args)))
            output = formprocess.stdout
            # write the newly formatted expression in the expression attribute.
            output = re.sub(r"(\s)*", "", output)  # Remove all whitespaces and newlines: \s = [\t\n\r\f\v]
            output = re.sub(r"\\", "",
                            output)  # Sometimes FORM splits indices in long expressions with an backslash "\" which is discarded.
            pattern = r"(?<!Local)coefficient=(?P<coefficient>.*);"
            match = re.search(pattern, output)
            if match:
                self.expression = match.group("coefficient")
            else:
                logger.error(f"No coefficient has been found in {output}.")
                sys.exit("STOP")

        self.tex_expression = self.get_tex()

    def get_tex(self):
        """
        Write Tex expression of coefficient.
        Returns
        -------

        """
        TERM_PATH = FORM_PATH / self.name
        TERM_PATH.mkdir(parents=True, exist_ok=True)  # Create directories if they don't exist.
        form = "Off statistics;\n"
        symbols = "Symbols " + ", ".join(coeffvalues + ["n"]) + ";\n"
        form += symbols  # Symbols d, lambdah, At, g1, g2, g3, [c_mu], lambdaphi, kappa, I, [m_{Ms}], [m_{Mu}], [Mu], eps, n, [2L[Ms,Mu]];

        dimlessConst = []
        for key, constant in coeffname.items():
            if key == "I" or key == "Ms":
                continue
            if constant["massdim"] != 0:
                if constant["massdim"] == 1:
                    dimlessConst.append(f"[{constant['FORM']:s}/{coeffname['Ms']['FORM']:s}]")
                else:
                    dimlessConst.append(f"[{constant['FORM']:s}/{coeffname['Ms']['FORM']:s}^{constant['massdim']:d}]")
        form += f"Symbols {', '.join(dimlessConst):s};\n"

        # form += "\n"
        #         form += f"Local coefficient = {self.expression:s};\n"
        #         form += ".sort\n"
        #         form += "\n"
        #         form += "#call factorizeCoeff\n"
        #         form += "\n"
        #         form += ".sort\n"
        #         form += "CFunction sign, coeff, frac;\n"
        #         form += "Symbols u, v, x, y;\n"
        #         form += "* Get sign\n"
        #         form += "PolyFun sign;\n"
        #         form += ".sort\n"
        #         form += "* Write positive numerical coefficient in coeff\n"
        #         form += "PolyFun coeff;\n"
        #         form += "id sign(x?) = sig_(x)*x*sign(sig_(x));\n"
        #         form += ".sort\n"
        #         form += "PolyFun;\n"
        #         form += ".sort\n"
        #         form += "* write numerical fraction in frac function\n"
        #         form += "PolyRatFun frac;\n"
        #         form += ".sort\n"
        #         form += "PolyRatFun;\n"
        #         form += "* write everything else in frac function\n"
        #         form += "repeat;\n"
        #         form += "\tid d?!{Ms}^n?pos_ = frac(d^n, 1);\n"
        #         form += "\tid d?!{Ms}^n?neg_ = frac(1, d^-n);\n"
        #         form += "endrepeat;\n"
        #         form += "* combine the separate fractions\n"
        #         form += "repeat;\n"
        #         form += "\tid frac(u?, v?)*frac(x?, y?) = frac(u * x, v * y);\n"
        #         form += "endrepeat;\n"
        #         form += "repeat;\n"
        #         form += "\tid frac(u?, 1) = u;\n"
        #         form += "\tid sign(x?) = x;\n"
        #         form += "\tid coeff(1) = 1;\n"
        #         form += "endrepeat;\n"
        #         form += ".sort\n"
        #         form += "* Rewrite everything in Latex notation\n"
        #         form += "#OpenDictionary constants\n"
        form += "\n"
        form += f"Local coefficient = {self.expression:s};\n"
        form += ".sort\n"
        form += "\n"
        form += "#call factorizeCoeff\n"
        form += "\n"
        form += ".sort\n"
        form += "CFunction sign, coeff, frac;\n"
        form += "Symbols u, v, x, y;\n"
        form += "* Get sign\n"
        form += "PolyFun sign;\n"
        form += ".sort\n"
        form += "* Write positive numerical coefficient in coeff and sign in sign\n"
        form += "PolyFun;\n"
        form += "id sign(x?) = coeff(sig_(x)*x)*sign(sig_(x));\n"
        form += ".sort\n"
        form += "* write everything else in frac function\n"
        form += "repeat;\n"
        form += "\tid d?!{Ms}^n?pos_ = frac(d^n, 1);\n"
        form += "\tid d?!{Ms}^n?neg_ = frac(1, d^-n);\n"
        form += "endrepeat;\n"
        form += "* combine the separate fractions\n"
        form += "repeat;\n"
        form += "\tid frac(u?, v?)*frac(x?, y?) = frac(u * x, v * y);\n"
        form += "endrepeat;\n"
        form += "repeat;\n"
        form += "\tid frac(u?, 1) = u;\n"
        form += "\tid sign(x?) = x;\n"
        form += "\tid coeff(1) = 1;\n"
        form += "endrepeat;\n"
        form += ".sort\n"
        # Write the expansion mass Ms outside of a bracket
        form += "CFunction bbracket;\n"
        form += "Bracket Ms;\n"
        form += ".sort\n"
        form += "Collect bbracket;\n"
        #  Ensure that the first term in each bracket has a positive sign
        form += ".sort\n"
        form += """#procedure positiveTerm
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
#endprocedure
#call positiveTerm"""
        # Rewrite everything in Latex notation
        form += "* Rewrite everything in Latex notation\n"
        form += "#OpenDictionary constants\n"
        tex_statements = {}

        for key, constant in coeffname.items():
            if key not in ["I"]:
                if constant["massdim"] == 0 or key == "Ms":
                    tex_statements[constant["FORM"]] = constant["tex"][0]
                elif constant["massdim"] == 1 and key != "Ms":
                    tex_statements[f"[{constant['FORM']:s}/{coeffname['Ms']['FORM']:s}]"] = constant["tex_dimless"]
                elif constant["massdim"] > 1 and key != "Ms":
                    tex_statements[f"[{constant['FORM']:s}/{coeffname['Ms']['FORM']:s}^{constant['massdim']:d}]"] = constant["tex_dimless"]

        for key, constant in abbreviation.items():
            tex_statements[constant["FORM"]] = constant["tex"][0]

        for form_expr, tex_expr in tex_statements.items():
            form += f'\t#add {form_expr:s}: "{tex_expr:s}"\n'
        form += "#CloseDictionary\n"
        form += "\n"
        # form += f"Bracket {coeffname['Ms']['FORM']:s};\n"
        # form += "\n"
        form += "#UseDictionary constants\n"
        form += "Print +s;\n"
        form += ".end"
        with open(TERM_PATH / "tex_coefficient.frm", "w") as file:
            file.write(form)

        try:
            formprocess = subprocess.run(["form", "-p", TERM_PATH, TERM_PATH / "tex_coefficient.frm"],
                                         capture_output=True, text=True, check=True)
        except subprocess.CalledProcessError as exc:
            exc.cmd = list(map(str, list(exc.cmd)))
            logger.error(" ".join(exc.cmd) + "\n" + str(exc.stdout))
            logger.error(f"form returned non-zero exit status {exc.returncode}")
            sys.exit("STOP")
        else:
            # No Error occured
            logger.debug(" ".join(map(str, formprocess.args)))
            output = formprocess.stdout
            # write the newly formatted expression in the expression attribute.
            output = re.sub(r"(\s)*", "", output)  # Remove all whitespaces and newlines: \s = [\t\n\r\f\v]
            pattern = r"(?<!Local)coefficient=(?P<coefficient>.*);"
            match = re.search(pattern, output)
            if match:
                tex = match.group("coefficient")
                # substitute coeff(x/y) by \frac{x}{y}
                tex = re.sub(r"coeff\((?P<nominator>\d+)/(?P<denominator>\d+)\)",
                             r"\\" + r"frac{\g<nominator>}{\g<denominator>}", tex)
                # substitute coeff(x) by x
                tex = re.sub(r"coeff\((?P<number>\d+)\)",
                             r"\g<number>", tex)

                # First, substitute fractions standing at the end of the term
                tex = re.sub(r"frac\((?P<nominator>[^+-]+),(?P<denominator>[^+-]+)\)" + r"\)", #  + r"(?=\)|\*" + escape_regex(coeffname["Ms"]["tex"]) + r")",
                             r"\\" + r"frac{\g<nominator>}{\g<denominator>}" + ")", tex)
                # substitute all remaining fractions
                tex = re.sub(r"frac\((?P<nominator>[^+-]+),(?P<denominator>[^+-]+)\)",
                             r"\\" + r"frac{\g<nominator>}{\g<denominator>}", tex)

                # The Expression should look like: bbracket(<a bunch of terms here>)*Ms^-n -> \left(<a bunch of terms here>\right)*M_{s}^{-n}
                tex = re.sub(r"bbracket\((?P<summands>.+)\)\*" + escape_regex(coeffname["Ms"]["tex"]) + r"\^(?P<exponent>-?\d{1,2})",
                             r"\\" + r"left(" + r"\g<summands>" + r"\\" + r"right)*" + escape_regex(coeffname["Ms"]["tex"]) + r"^{\g<exponent>}", tex)

                # Subistute terms where n = 0, i.e. Ms^n = 1:
                tex = re.sub(r"bbracket\((?P<summands>.+)\)\Z",
                             r"\\" + r"left(" + r"\g<summands>" + r"\\" + r"right)",
                             tex)
            else:
                logger.error(f"No coefficient has been found in {output}.")
                sys.exit("STOP")
        return tex
