import re
import logging
import sys
import subprocess

from pathlib import Path
from abc import ABC, abstractmethod

from tioc import run_form, op_config, coeff
from tioc import escape_regex, PROJECTION_PATH, FORM_PATH, FORM_GENERAL_PATH
logger_autoeft = logging.getLogger("autoeft.projection")
logger = logger_autoeft.getChild("coefficient")

class Coefficient_Model(ABC):
    r"""
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
    """
    expr : str
    output_dimless = True # Specifies whether tex output contains only dimensionless variables or not
    @abstractmethod
    def __init__(self, expr):
        self.expr = expr

    @abstractmethod
    def __repr__(self):
        """Spedify general string representation of the coefficient."""
        return f"{self.expr:s}"

    def __str__(self):
        """Specify the format for printing an Operator object with the print command or converting a Operator object in
        string with str() function. """
        return self.__repr__()

    def __format__(self, key):
        """Specify the format for "format" function in print statement: Here the same as the print statement itself."""
        if key == "tex":
            return self.tex
        else:
            return self.__repr__()

    @abstractmethod
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
        if self.name and other.name:
            new_name = self.name
        elif self.name and not other.name:
            new_name = self.name
        elif not self.name and other.name:
            new_name = self.name
        else:
            login.error("At least one coefficient has to be named.")
            sys.exit("STOP")

        return Coefficient(f"({self.expr})*({other.expr})", new_name)

    @abstractmethod
    def __add__(self, other):
        """
        Add a string from the right to the expression of the Coefficient.
        Parameters
        ----------
        other: str
            Example: other = '(1/2)'
        Returns
        -------
        """
        if self.name and other.name:
            new_name = self.name
        elif self.name and not other.name:
            new_name = self.name
        elif not self.name and other.name:
            new_name = self.name
        else:
            login.error("At least one coefficient has to be named.")
            sys.exit("STOP")

        return Coefficient(f"{self.expr}+{other.expr}", new_name)

    @property
    def tex(self):
        """
        Write Tex expression of coefficient.
        Returns
        -------
        """
        TERM_PATH = FORM_PATH / self.name
        TERM_PATH.mkdir(parents=True, exist_ok=True)  # Create directories if they don't exist.
        ms = op_config["coefficients"]["Ms"]["mathematica"]["Ms"]  # FORM expression of EFT mass Ms

        form = "Off statistics;\n"
        form += "#include declarations_general.h # coefficient\n"
        form += "\n"
        form += f"Local coefficient = {self.expr:s};\n"
        form += ".sort\n"
        form += "\n"

        if self.output_dimless:
            form += "#call factorizeCoeff\n"
            form += ".sort\n"

        form += "CFunction sign, coeff, frac;\n"
        form += "CFunction bbracket;\n"
        # write fractions into frac(nominator, denominator) function
        form += "#call fractorizeCoeff\n"
        #  Ensure that the first term in each bracket has a positive sign
        form += "#call positiveTerm\n"
        # Rewrite everything in Latex notation
        form += "* Rewrite everything in Latex notation\n"
        form += "#OpenDictionary constants\n"
        tex_statements = {}

        if self.output_dimless:
            # only dimless constants appear.
            for key, constant in op_config["coefficients"].items():
                form_coeff = list(constant['mathematica'].values())[0]
                if constant["massdim"] == 0 or key == "Ms":
                    tex_statements[form_coeff] = constant["tex"]
                elif constant["massdim"] == 1 and key != "Ms":
                    tex_statements[f"[{form_coeff:s}/{ms:s}]"] = constant["tex_dimless"]
                elif constant["massdim"] > 1 and key != "Ms":
                    tex_statements[f"[{form_coeff:s}/{ms:s}^{constant['massdim']:d}]"] = constant["tex_dimless"]
        else:
            for key, constant in op_config["coefficients"].items():
                form_coeff = list(constant['mathematica'].values())[0]
                tex_statements[form_coeff] = constant["tex"]

        for key, constant in op_config["abbreviation"].items():
            form_coeff = list(constant['mathematica'].values())[0]
            tex_statements[form_coeff] = constant["tex"]

        for form_expr, tex_expr in tex_statements.items():
            form += f'\t#add {form_expr:s}: "{tex_expr:s}"\n'
        form += "#CloseDictionary\n"
        form += "\n"

        form += "#UseDictionary constants\n"
        form += "Print +s;\n"
        form += ".end"
        with open(TERM_PATH / "tex_coefficient.frm", "w") as file:
            file.write(form)

        output = run_form(fp_cwd=TERM_PATH, filename="tex_coefficient.frm", fp_p=FORM_GENERAL_PATH, keep_backslash=True)
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
            tex = re.sub(r"bbracket\((?P<summands>.+)\)\*" + escape_regex(op_config["coefficients"]["Ms"]["tex"]) + r"\^(?P<exponent>-?\d{1,2})",
                         r"\\" + r"left(" + r"\g<summands>" + r"\\" + rf"right)*" + op_config["coefficients"]["Ms"]["tex"] + r"^{\g<exponent>}", tex)

            # Subistute terms where n = 0, i.e. Ms^n = 1:
            tex = re.sub(r"bbracket\((?P<summands>.+)\)\Z",
                         r"\\" + r"left(" + r"\g<summands>" + r"\\" + r"right)",
                         tex)
            multSign = ""
            tex = re.sub(r"\*", f"{multSign} ", tex)
        else:
            logger.error(f"No coefficient has been found in {output}.")
            sys.exit("STOP")
        return tex

class Coefficient(Coefficient_Model):
    def __init__(self, expr, name):
        super().__init__(expr)
        self.name = name
        # self.mdim
        # self.tex

    def __repr__(self):
        """Spedify general string representation of the coefficient."""
        return super().__repr__()

    def __mul__(self, other):
        return super().__mul__(other)

    def __add__(self, other):
        return super().__add__(other)

    @property
    def mdim(self):
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

        ms = op_config["coefficients"]["Ms"]["mathematica"]["Ms"]  # FORM expression of EFT mass Ms

        form = "Off statistics;\n"
        form += "#include declarations_general.h # coefficient\n"
        form += "\n"

        form += f"Local coefficient = {self.expr:s};\n"
        form += ".sort\n"
        form += "\n"
        if self.output_dimless:
            form += "#call factorizeCoeff\n"
            form += ".sort\n"
        form += "\n"
        form += f"id {ms}^n?$dim = {ms}^n;\n"
        form += "$dim = 4 - $dim;\n"
        form += "\n"
        form += ".sort\n"
        form += f"Bracket {ms:s};\n"
        FILE_PATH = TERM_PATH / Path("dimension.aux")
        form += f'#write <{str(FILE_PATH):s}> "dim = %$", $dim\n'
        form += "\n"
        form += "Print +s;\n"
        form += ".end"

        with open(TERM_PATH / "get_dimension.frm", "w") as file:
            file.write(form)

        run_form(fp_cwd=TERM_PATH, fp_p=FORM_GENERAL_PATH, filename="get_dimension.frm")

        with open(FILE_PATH, "r") as file:
            match = re.match(r"dim = (?P<dim>\d{1,2})", file.read())
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

        return dim

class Factor(Coefficient_Model):
    """
    Class for small factors which are multiplied to a Coefficient of a Summand.
    """
    def __init__(self, expr):
        super().__init__(expr)
        self.name = None
    def __repr__(self):
        """Spedify general string representation of the coefficient."""
        return super().__repr__()
    def __mul__(self, other):
        return Factor(f"{self.expr}*{other.expr}")
    def __add__(self, other):
        return Factor(f"{self.expr}+{other.expr}")
