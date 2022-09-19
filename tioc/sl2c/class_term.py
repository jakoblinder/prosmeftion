import logging
import re
import sys
from typing import List, Tuple

from tioc import opname, opnameSL2C, spinorsSL2C_c, field_config

from .class_coefficient import Coefficient
from .class_operator import Operator, OperatorModel

logger_autoeft = logging.getLogger("autoeft.projection")
logger = logger_autoeft.getChild("term")
# logger = logger_autoeft.getChild(__name__)

class Term(Operator, Coefficient):
    """
    A class describing one term of Benjamin Summs output, which consist of approximately 100 terms each containing
    some operators and a matching coefficient.

    ...

    Attributes
    ----------
    coeff : str
        Coefficient of the operator.
    cops : str
        Expression of the fully contracted operators themself matching the coefficient. This expression is rewritten
        such that is compatible with FORM.
    operators : [Operator, Operator, ...]
        List of the operators in the expression as objects of the Operator class.
    indices : [str, str, ...]
        Indices occurring in the expression, without duplicates. Important for the sum statement in FORM.
    form : str
        Formatted string that should be written in a file name.frm and executed in FORM. Converts term in the y-Basis.
    Magic Methods
    -------
    __str__()
        Specifies the str() function.
    __format__()
        Specifies the format() function used in print().
    Methods
    -------
    get_operators()
    get_formCops()
    """
    manual: bool
    coeff: Coefficient
    cops: str
    operators: Tuple[Operator]
    indices: List[str]
    cops_original: str
    sorted_indices: List[str]
    possible_indices: List[str]

    def __init__(self, *args, **kwargs):
        """
        Constructor of Class Term.

        Parameters
        ----------
        args : tuple
            Should contain the following three variables if the Terms are generated automatically:
                    args[0]: str
                        Coefficient of the operator.
                    args[1]: str
                        Expression of the fully contracted operators themselves matching the coefficient. This expression
                        is rewritten such that is compatible with FORM.
                    args[2]: str
                        Name of the term, e.g. term12.
        kwargs["operators"] : tuple
            If this argument exists, the expression is written manually, i.e. each Operator has to be written in a FORM
            compatible way, containing all Indices, in the correct order in a tuple.
            Note: The operators don't have to be fully contracted.
            Tuple containing the operators in the correct order.
        kwargs["operators"] : list containing all indices which occur in the operator. Indices which are contracted should occur twice
            in the declaration.
        Returns
        -------
        None.

        """
        if "operators" not in kwargs.keys():
            self.manual = False
            self.coeff = Coefficient(args[0], name=args[2])
            self.cops = args[1]  # Fully contracted operators from Mathematica output.
            self.name = args[2]
            self.operators = self.get_operators()  # Separate operators from the fully contracted expression.
        else:
            # manual = True
            self.manual = True
            self.operators = ()
            try:
                self.operators = []
                assert len(kwargs["operators"]) == len(
                    kwargs["indices"]), "Each Operator musst have a list containing its Index, even if it is empty."
                indices = kwargs["indices"]
                for i, v in enumerate(kwargs["operators"]):
                    v = re.sub(r"(\s)*", "", v)  # Replace all whitespaces and newlines: \s = [\t\n\r\f\v]
                    try:
                        self.operators.append(Operator(expression=v, id=i, manual=True, indices=indices[i], name=kwargs["name"][i]))
                    except KeyError:
                        self.operators.append(Operator(expression=v, id=i, manual=True, indices=indices[i]))
                    except IndexError:
                        self.operators.append(Operator(expression=v, id=i, manual=True, indices=indices[i]))
                self.operators = tuple(self.operators)
                self.coeff = None
            except KeyError:
                raise ValueError("Operators or Indices are not found in the kwargs.")
            try:
                self.name = kwargs["name"]
            except KeyError:
                logger.error("Name of the term has to be given.")
                sys.exit("STOP")

        self.indices = None  # List with occurring indices over which one should sum.
        # In the end every index should appear only ones, but in order to have an additional check
        # let them stay in the beginning to check that every index occurs exactly twice. Only then
        # remove the doubles.

        # Replace fully contracted operators:
        self.cops, self.indices = self.get_formCops()
        self.cops_original = self.cops
        logger.debug(f"New Term: {self.cops}")
        # sort indices
        self.sorted_indices = self.sort_indices()
        # In order to declare all necessary indices and the corresponding sets, define an over-defined set of indices:
        self.possible_indices = self.poss_indices(1, 100)

    def __str__(self):
        """Specify the format for printing a Term object with the print command or converting a Term object in
        string with str() function. """
        if self.coeff == None:
            return "{0:s}".format(self.cops)
        else:
            return "{0:s}*{1:s}".format(self.coeff, self.cops)
    def __repr__(self):
        return self.__str__()
    def __format__(self, key):
        """Specify the format for "format" function in print statement: Here the same as the print statement itself."""
        return self.__str__()  # "{0:s} * {1:s}".format(self.coeff, self.cops)

    def get_operators(self):
        """Utility rsync
        Extract each operator from the fully contracted operator product in each term. Note: An operators denotes in
        this case everything which carries an index, e.g. the yukawa matrix or the su2epsilon are listed as
        non fermionic operators.
        Returns
        -------
        Tuple of operators.
        """
        match = re.finditer(r"\*\*",
                            self.cops)  # Find operators by identifying always non-commuting multiplication "**".
        if match:
            matches = list(match)
            operators = []
            for i, v in enumerate(matches):
                if i > 0:
                    operators.append(self.cops[matches[i - 1].end():v.start()])
                else:
                    operators.append(self.cops[:v.start()])  # Append first operator as the FIRST element in the list.
            operators.append(self.cops[matches[-1].end():])  # Append also last operator.
        else:
            raise ValueError("Only one operator was found which is not supported until now.")
        operatorobj = []
        for i, v in enumerate(operators):
            operatorobj.append(Operator(expression=v,id=i))  # Operators are realised as objects from the Operator Class.
        return tuple(operatorobj)  # Rewrite operator list as tuple to ensure order.

    def sort_indices(self):
        """
        Sort indices of different kind, each in a list.
        Returns
        -------

        """
        # dictionary for indices:
        indDic = {}

        # Lorentz indices
        indDic["lor"] = self.extract_indices(r"lor\d{1,5}")
        indDic["lor"] += self.extract_indices(r"lorA\d{1,5}")
        # spin indices - for gamma matrices and so on
        indDic["spin"] = self.extract_indices(r"spin\d{1,5}")
        # Gauge indices - SU2 fundamental
        indDic["gauge"] = self.extract_indices(r"gauge\d{1,5}")
        indDic["gauge"] += self.extract_indices(r"gaugeA\d{1,5}")
        # gauge adjunct indices - SU2 adjoint
        indDic["gaugeadj"] = self.extract_indices(r"gaugeadj\d{1,5}")
        # color indices - SU3 fundamental
        indDic["colf"] = self.extract_indices(r"colf\d{1,5}")
        indDic["colf"] += self.extract_indices(r"colfA\d{1,5}")
        # color indices - SU3 adjoint
        indDic["cola"] = self.extract_indices(r"cola\d{1,5}")
        # flavor indices
        indDic["flav"] = self.extract_indices(r"flav\d{1,5}")
        # if indDic["flav"] == []: del indDic["flav"]
        # SL2C-Indices:
        indDic["Usl"] = self.extract_indices(r"Usl\d{1,3}")
        indDic["Lsl"] = self.extract_indices(r"Lsl\d{1,3}")
        indDic["Usldot"] = self.extract_indices(r"Usldot\d{1,3}")
        indDic["Lsldot"] = self.extract_indices(r"Lsldot\d{1,3}")
        indDic["Usl"] += self.extract_indices(r"Uslop\d{1,2}\d{1,3}")
        indDic["Lsl"] += self.extract_indices(r"Lslop\d{1,2}\d{1,3}")
        indDic["Usldot"] += self.extract_indices(r"Usldotop\d{1,2}\d{1,3}")
        indDic["Lsldot"] += self.extract_indices(r"Lsldotop\d{1,2}\d{1,3}")
        return indDic

    def poss_indices(self, min=1, max=100):
        # dictionary for indices:
        indDic = {}
        # SL2C-Indices:
        indDic["Usl"] = ["Usl{0:d}".format(i) for i in range(min, max + 1)]
        indDic["Lsl"] = ["Lsl{0:d}".format(i) for i in range(min, max + 1)]
        indDic["Usldot"] = ["Usldot{0:d}".format(i) for i in range(min, max + 1)]
        indDic["Lsldot"] = ["Lsldot{0:d}".format(i) for i in range(min, max + 1)]
        for i in self.sorted_indices["Usl"]:
            if i not in indDic["Usl"]:
                match = re.match(r"Usl(?P<usl>(op\d{1,2})?\d{1,3})", i)
                if match:
                    slN = match.group("usl")
                    indDic["Usl"].append(f"Usl{slN:s}")
                    # Need to append also the regarding Lsl index in order to get a correct matching in the ULsl and LUsl sets.
                    indDic["Lsl"].append(f"Lsl{slN:s}")
        for i in self.sorted_indices["Lsl"]:
            if i not in indDic["Lsl"]:
                match = re.match(r"Lsl(?P<lsl>(op\d{1,2})?\d{1,3})", i)
                if match:
                    slN = match.group("lsl")
                    indDic["Lsl"].append(f"Lsl{slN:s}")
                    indDic["Usl"].append(f"Usl{slN:s}")
        for i in self.sorted_indices["Usldot"]:
            if i not in indDic["Usldot"]:
                match = re.match(r"Usldot(?P<usldot>(op\d{1,2})?\d{1,3})", i)
                if match:
                    sldotN = match.group("usldot")
                    indDic["Usldot"].append(f"Usldot{sldotN:s}")
                    indDic["Lsldot"].append(f"Lsldot{sldotN:s}")
        for i in self.sorted_indices["Lsldot"]:
            if i not in indDic["Lsldot"]:
                match = re.match(r"Lsldot(?P<lsldot>(op\d{1,2})?\d{1,3})", i)
                if match:
                    sldotN = match.group("lsldot")
                    indDic["Lsldot"].append(f"Lsldot{sldotN:s}")
                    indDic["Usldot"].append(f"Usldot{sldotN:s}")

        # Lorentz indices
        indDic["lor"] = self.sorted_indices["lor"].copy()
        # spin indices - for gamma matrices and so on
        indDic["spin"] = self.sorted_indices["spin"].copy()
        # Auxiliary spin indices for substituteSpinors method
        spinA = []
        for i in self.sorted_indices["spin"]:
            matchSpin = re.match(r"spin(\d{1,5})", i)
            if matchSpin:
                spinA.append("spinA" + matchSpin.group(1))
        indDic["spinA"] = tuple(spinA)
        assert len(indDic["spinA"]) == len(indDic["spin"])
        # Gauge indices - SU2 fundamental
        indDic["gauge"] = self.sorted_indices["gauge"].copy()
        temp_gauge = self.sorted_indices["gauge"].copy()
        for ind in temp_gauge:
            indDic["gauge"].append(ind + "a")
            indDic["gauge"].append(ind + "b")

        # gauge adjoint indices - SU2 adjoint
        indDic["gaugeadj"] = self.sorted_indices["gaugeadj"].copy()
        # color indices - SU3 fundamental
        indDic["colf"] = self.sorted_indices["colf"].copy()
        temp_colf = self.sorted_indices["colf"].copy()
        for ind in temp_colf:
            indDic["colf"].append(ind + "a")
            indDic["colf"].append(ind + "b")
            indDic["colf"].append(ind + "c")
            indDic["colf"].append(ind + "d")

        # color indices - SU3 adjoint
        indDic["cola"] = self.sorted_indices["cola"].copy()
        # flavor indices
        indDic["flav"] = self.sorted_indices["flav"].copy()

        return indDic

    def update_sortIndices(self):
        """
        This method should be called if indices are attached from functions inside the form function of the projection file.
        This is especially important for the fieldstrengthtensors.
        """
        cindices = []
        for i in self.operators:
            cindices += i.indices
        if not self.manual:
            # Ensure that every index occurs exactly twice.
            con, uncon = self.contracted(cindices)
            assert con, f"All Indices should be contracted, but {str(uncon):s} are not contracted."
        # Remove all duplicates of indices
        self.indices = list(set(cindices))

        # sort indices
        self.sorted_indices = self.sort_indices()
        # In order to declare all necessary indices and the corresponding sets, define an over-defined set of indices:
        self.possible_indices = self.poss_indices(1, 100)

    def contracted(self, cindices):
        contract = []
        uncontractedInd = []
        cindices = list(map(str,cindices))
        for i in cindices:
            if cindices.count(i) == 2:
                contract.append(True)
            else:
                uncontractedInd.append(i)
        # Copy list for iteration in such away that only elements from "uncontractedInd" are removed in an "uncontractedInd.remove("something")" order.
        uncontractedInd_tmp = uncontractedInd.copy()
        # Consider now the possible uncontracted indices which can only be SL2C-indices:
        for i in uncontractedInd_tmp:
            # try:
            matchUsl = re.match(r"Usl(?P<usl>(op\d{1,2})?\d{1,3})", i)
            # except TypeError:
            #     print(i)
            #     sys.exit("STOP")
            if matchUsl:
                if "Lsl{0:s}".format(matchUsl.group("usl")) in uncontractedInd:
                    contract.append(True)
                    uncontractedInd.remove("Usl{0:s}".format(matchUsl.group("usl")))
                    uncontractedInd.remove("Lsl{0:s}".format(matchUsl.group("usl")))
            matchUsldot = re.match(r"Usldot(?P<usldot>(op\d{1,2})?\d{1,3})", i)
            if matchUsldot:
                if "Lsldot{0:s}".format(matchUsldot.group("usldot")) in uncontractedInd:
                    contract.append(True)
                    uncontractedInd.remove("Usldot{0:s}".format(matchUsldot.group("usldot")))
                    uncontractedInd.remove("Lsldot{0:s}".format(matchUsldot.group("usldot")))
        if len(uncontractedInd) == 0 and all(contract):
            return True, []
        else:
            return False, uncontractedInd

    def get_formCops(self):
        """
        Construct whole FORM expression by joining the Operators from the Operator class. Further, also the indices are
        collected and it is checked that each index occurs exactly twice. Afterwards duplicates are deleted.
        Returns
        -------
        conOp : str
            Fully contracted FORM compatible Operator expression of the term.
        cindices : [str, str, ...]
            List with indices over which is summed in the fully contracted Operator.
        """
        conOp = ""
        cindices = []
        for i in self.operators:
            conOp += i.expression
            conOp += "*"
            cindices += i.indices
        # Remove last redundant "*"
        conOp = conOp[:-1]
        if not self.manual:
            # Ensure that every index occurs exactly twice.
            con, uncon = self.contracted(cindices)
            assert con, "All Indices should be contracted, but {0:s} are not contracted.".format(str(uncon))
            # for i in cindices:
            #     assert cindices.count(i)==2, "All Indices should be contracted."
        # Remove all duplicates of indices
        cindices = list(set(cindices))
        return conOp, cindices

    def get_form(self, cout=False):
        """
        Writes the 'Local expr = D(lor1,..)...;' statement for each term in a form file termI.frm which can be executed
        and further manipulated in FORM.

        Returns
        -------

        """
        form = ""
        coperator = "Local expr = {0:s};\n".format(
            self.cops)  # self.cops gives only the fully contracted operators, but rather
        # only self gives the whole expression including the coefficient.
        form += coperator
        # summ = "sum " + ", ".join(v.indices) + ";\n"
        # form += summ
        if cout: form += "Print;\n"
        form += ".sort\n\n"
        return form

    def form_higgsDerivativetoCommutative(self, nDer=4):
        """
        Rewrite (if there exists a gauge index in the term) Higgfields inside of a derivative as:
            D(lor1?, D(lor2?, H('gaugeindex'))) = D(lor1, D(lor2, [H('gaugeindex')])) * HC('gaugeindex'),
        where [H('gaugeindex')] is some placeholder to remember the position of the Higgsfield and HC('gaugeindex')
        is a commuting field. This makes it easier to implement id-statement for contracted Higgsfields.
        Note: Also a "lonely" Higgsfield is replaced,
            H('gaugeindex') = [H('gaugeindex')]*HC('gaugeindex'),
        in order to ensure also the order of the "lonely" Higgsfields which is of course much more important for
        non-commuting, e.g. fermionic objects, but also does not harm anything.

        Parameters
        ----------
        nDer : int
            Maximum number of derivatives.
        Returns
        -------
        initialize : str
            Contains essentially the declaration of the dummy function [H('gaugeindex')] and [[H+]('gaugeindex')] with
            the correct indices.
        form : str
            Formatted string of the FORM procedure file "higgsDerivativetoCommutative.prc" with the correct indices.
            in the do-loop index list.

        """
        initialize = ""
        form = ""
        if len(self.sorted_indices["gauge"]) != 0:
            # Initialize auxiliary Higgsfields:
            initialize += "Function "
            for v in self.sorted_indices["gauge"]:
                initialize += "[H({0:s})], [[H+]({0:s})], ".format(v)
            # Replace last comma by semicolon.
            initialize = initialize[:-2] + ";\n"
            # Form expression of the function which should be stored in "higgsDerivativetoCommutative.prc":
            form_list = []
            form_list_conj = []
            for i in range(0, nDer + 1):
                form_temp = "id " + Term.derivativesLorentz("H('gaugeindex')", n=i)
                form_temp += " = " + Term.derivativesLorentz("[H('gaugeindex')]", n=i, questionmark=False)
                form_temp += "*HC('gaugeindex');"
                form_list.append(form_temp)
            for i in range(0, nDer + 1):
                form_temp = "id " + Term.derivativesLorentz("[H+]('gaugeindex')", n=i)
                form_temp += " = " + Term.derivativesLorentz("[[H+]('gaugeindex')]", n=i, questionmark=False)
                form_temp += "*[HC+]('gaugeindex');"
                form_list_conj.append(form_temp)
            form += "#procedure higgsDerivativetoCommutative\n"
            form += "#do gaugeindex={"
            for v in self.sorted_indices["gauge"]:
                form += v + ","
            # Remove last comma.
            form = form[:-1]
            form += "}\n"
            form += "* Replacements of ordinary Higgs field by commuting Higgs field and dummy Higgs.\n"
            # form += "\tid H('gaugeindex') = [H('gaugeindex')]*HC('gaugeindex');\n"
            for i in form_list:
                form += 1 * "\t" + i + "\n"
            form += "* Replacements of Hermitian conjugated field by commuting Higgs field and dummy Higgs.\n"
            # form += "\tid [H+]('gaugeindex') = [[H+]('gaugeindex')]*[HC+]('gaugeindex');\n"
            for i in form_list_conj:
                form += 1 * "\t" + i + "\n"
            #             form += """\tid H('gaugeindex') = [H('gaugeindex')]*HC('gaugeindex');
            # \tid D(lor1?, H('gaugeindex')) = D(lor1, [H('gaugeindex')]) * HC('gaugeindex');
            # \tid D(lor1?, D(lor2?, H('gaugeindex'))) = D(lor1, D(lor2, [H('gaugeindex')])) * HC('gaugeindex');
            # \tid D(lor1?, D(lor2?, D(lor3?, H('gaugeindex')))) = D(lor1, D(lor2, D(lor3, [H('gaugeindex')]))) * HC('gaugeindex');
            # \tid D(lor1?, D(lor2?, D(lor3?, D(lor4?, H('gaugeindex'))))) = D(lor1, D(lor2, D(lor3, D(lor4, [H('gaugeindex')])))) * HC('gaugeindex');
            # \tid [H+]('gaugeindex') = [[H+]('gaugeindex')]*[HC+]('gaugeindex');
            # \tid D(lor1?, [H+]('gaugeindex')) = D(lor1, [[H+]('gaugeindex')]) * [HC+]('gaugeindex');
            # \tid D(lor1?, D(lor2?, [H+]('gaugeindex'))) = D(lor1, D(lor2, [[H+]('gaugeindex')])) * [HC+]('gaugeindex');
            # \tid D(lor1?, D(lor2?, D(lor3?, [H+]('gaugeindex')))) = D(lor1, D(lor2, D(lor3, [[H+]('gaugeindex')]))) * [HC+]('gaugeindex');
            # \tid D(lor1?, D(lor2?, D(lor3?, D(lor4?, [H+]('gaugeindex'))))) = D(lor1, D(lor2, D(lor3, D(lor4, [[H+]('gaugeindex')])))) * [HC+]('gaugeindex');"""
            #             form += "\n"
            form += "#enddo\n"
            form += "#endprocedure"
            if initialize == "" and form == "":
                return
            else:
                return initialize, form

    def form_higgsCommutativetoDerivative(self, nDer=4):
        """
        Rewrite (if there exists a gauge index in the term) Higgfields again in inside of the derivative where it
        came from as:
            D(lor1?,[H('gaugeindex1')])*HC('gaugeindex2') = D(lor1,H('gaugeindex2')),
        where [H('gaugeindex1')] is some placeholder to remember the position of the Higgsfield and HC('gaugeindex2')
        is a commuting field which is no contracted with an SU2 Epsilon Tensor. This makes it easier to implement
        id-statement for contracted Higgsfields.
        Note: Also a "lonely" Higgsfield is replaced backwards,
            [H('gaugeindex1')]*HC('gaugeindex2') = H('gaugeindex2'),
        in order to ensure also the order of the "lonely" Higgsfields which is of course much mor important for
        non-commuting, e.g. fermionic objects, but also does not harm anything.

        Returns
        -------
        form : str
            Formatted string of the FORM procedure file "higgsDerivativetoCommutative.prc" with the correct indices.
            in the do-loop indice list.

        """
        form = ""
        if len(self.sorted_indices["gauge"]) != 0:
            form += "#procedure higgsCommutativetoDerivative\n"
            form += "#do gaugeindex1={" + ",".join([gauge for gauge in self.sorted_indices["gauge"]]) + "}\n"
            form += "\t#do gaugeindex2={,a,b}\n"
            form_list = []
            form_list_conj = []
            for i in range(nDer + 1):
                lHS, rHS = Term.derivativesSLbackwards("[H('gaugeindex1')]", "H('gaugeindex1''gaugeindex2')", n=i)
                form_temp = "id " + lHS
                form_temp += "*HC('gaugeindex1''gaugeindex2')"
                form_temp += " = " + rHS + ";"
                form_list.append(form_temp)
            for i in range(nDer + 1):
                lHS, rHS = Term.derivativesSLbackwards("[[H+]('gaugeindex1')]", "[H+]('gaugeindex1''gaugeindex2')", n=i)
                form_temp = "id " + lHS
                form_temp += "*[HC+]('gaugeindex1''gaugeindex2')"
                form_temp += " = " + rHS + ";"
                form_list_conj.append(form_temp)

            form += "* Replacements of ordinary Higgs field.\n"
            for i in form_list:
                form += 2 * "\t" + i + "\n"
            form += "* Replacements of Hermitian conjugated field.\n"
            for i in form_list_conj:
                form += 2 * "\t" + i + "\n"
            form += "\t#enddo\n"
            form += "#enddo\n"
            form += "#endprocedure"

            return form

    def form_higgsReplacebyEps(self):
        """
        Replace contracted Higgsfields by Higgsfields and an SU2 epsilon tensor.
        Returns
        -------

        """
        form = ""
        if len(self.sorted_indices["gauge"]) != 0:
            form += "* Substitute Higgsfields and other anti-fundamental SU2 field (Q^dagger and L^dagger) with Higgsfields and SU2 epsilon tensors.\n"
            form += "#do gaugeindex={"
            for v in self.sorted_indices["gauge"]:
                form += v + ","
            # remove last comma
            form = form[:-1]
            form += "}\n"
            form += "\tid once [HC+]?{[HC+], " + f"{spinorsSL2C_c['[L+]']}, {spinorsSL2C_c['[Q+]']}" + "}(?c, 'gaugeindex', ?d) = [su2eps]('gaugeindex'a, 'gaugeindex') * [HC+](?c, 'gaugeindex'a, ?d);\n"
            form += "\tid once [HC+]?{[HC+], " + f"{spinorsSL2C_c['[L+]']}, {spinorsSL2C_c['[Q+]']}" + "}(?c, 'gaugeindex', ?d) = [su2eps]('gaugeindex'b, 'gaugeindex') * [HC+](?c, 'gaugeindex'b, ?d);\n"
            # form += "\tid HC('gaugeindex')*[HC+]('gaugeindex') = [su2eps]('gaugeindex'a, 'gaugeindex'b) * [HC+]('gaugeindex'b) * HC('gaugeindex'a);\n"
            form += "#enddo\n"

            return form
    #
    # Fieldstrengthtensor
    #
    def form_fieldstrengthtensorDerivativeHandling(self, nDer=4):  # fieldstrengthtensorDerivativetoCommutative
        """
        Rewrite fieldstrengthtensors inside a derivative as:
            D(lorA1?lor, B(lor1,lor2))
            = D(lorA1, [B(op0,lor1,lor2)])*(-i_/4)*(
              BCR(op0, Usldot1, Usldot2)*sigmabar2(lor1, lor2, Lsldot1, Lsldot2)
              - BCL(op0, Lsl1, Lsl2)*sigma2(lor1, lor2, Usl1, Usl2)
              ),
        and
            D(lorA1?lor, W(gaugeadj1,lor1,lor2))
            = D(lorA1, [W(op3,gaugeadj1,lor1,lor2)])*T(gaugeadj1,gaugeA10,gaugeA11)*[su2eps](gaugeA11,gaugeA12)*(-i_/4)*(
              WCR(op3, gaugeA10, gaugeA12, Usldot10, Usldot11)*sigmabar2(lor1, lor2, Lsldot10, Lsldot11)
              - WCL(op3, gaugeA10, gaugeA12, Lsl10, Lsl11)*sigma2(lor1, lor2, Usl10, Usl11)
              ).
        Note 1: The [B(op0,lor1,lor2)] are some placeholder to remember the position of the fieldstrengthtensor and
                the BCL, BCR, WCL, WCR, GCL and GCR are a commuting fields. This makes it easier to implement
                id-statement for contracted fields.
        Note 2: The additional index 'op0' specifies the position of the fields, which is important for the replacement
                of commuting fields again with not commuting fields if for example two identical fields are involved.
                TODO: This form of replacement should be done for all fields not only the fieldstrengthtensors.
                TODO: The index should corresponce to the field position used also later. Thus, the first index should have the position index 'op1' and not 'op0'.
        Note 3: Also a "lonely" fieldstrengthtensor is replaced, not only to ensure also the order of the "lonely"
                fieldstrengthtensors, but for the SU(2) and SU(3) fieldstrengthtensors also to replace the adjoint
                gauge group indices by fundamental ones:
                W(gaugeadj1,lor1,lor2)
                = [W(op3,gaugeadj1,lor1,lor2)]*T(gaugeadj1,gaugeA10,gaugeA11)*[su2eps](gaugeA11,gaugeA12)*(-i_/4)*(
                  WCR(op3, gaugeA10, gaugeA12, Usldot10, Usldot11)*sigmabar2(lor1, lor2, Lsldot10, Lsldot11)
                  - WCL(op3, gaugeA10, gaugeA12, Lsl10, Lsl11)*sigma2(lor1, lor2, Usl10, Usl11)
                  ).
        Note 4: The fundamental representation generator T(gaugeadj1,gauge2,gauge1) has one adjoint, one fundamental
                and one anti-fundamental index. Since the BSUOLEA-output contains terms like
                    [H+](gauge2)*T(gaugeadj1,gauge2,gauge1)*H(gauge1),
                and an fundamental index always has to be contracted with an anti-fundamental one, it is clear that
                the first index of T (here 'gauge2') is fundamental and second (here 'gauge1') is anti-fundamental.
                -> \tensor{(T^{A})}{_{a}^{b}}, with 'a' fundamental and 'b' anti-fundamental.
                The SU(2) and SU(3) fieldstrengthtensors are therefore replaced as follows:
                W^{I} = \eps^{j k} \tensor{(T^{I})}{_{k} ^{i}} W_{(i j)},
                <=> W(gaugeadj1,lor1,lor2) = [su2eps](gauge2,gauge3)*WT(gaugeadj1,gauge3,gauge1)*(gauge1,gauge2,lor1,lor2)
                G^{A} = \frac{\tensor{(T^{A})}{_{d} ^{b}}}{2} \eps^{a c d} G_{a b c}.
                <=> G(cola1,lor1,lor2) = (1/2)*T(cola1,colf4,colf2)*[su3eps](colf1,colf3,colf4)*G(colf1,colf2,colf3,lor1,lor2)

        Parameters
        ----------
        nDer : int
            Maximum number of derivatives.
        Returns: 4 FORM strings
        -------
        initialize_DtC : str
            Contains essentially the declaration of the dummy function [q('spin', gauge1, colf2, flav3)] and [[qbar]('spin', gauge1, colf2, flav3)] with
            the correct indices.
        form_DtC : str
            Formatted string of the FORM procedure file "spinorDerivativetoCommutative.prc" with the correct indices.
            in the do-loop index list.
            -> Replacements of fieldstrengthtensor by commuting fieldstrengthtensor and dummy fieldstrengthtensor.
        form_DtSL : str
            Formatted string of the FORM procedure file "fieldstrengthtensorDerivativetoSL2C.prc" with the correct indices.
            in the do-loop index list.
            -> Replacements of derivative acting on a fieldstrengthtensor by SL2C derivative and sigma matrix.
        form_CtD : str
            Formatted string of the FORM procedure file "fieldstrengthtensorDerivativetoSL2C.prc" with the correct indices.
            in the do-loop index list.
            -> Replacements of commuting fieldstrengthtensor and dummy fieldstrengthtensor by fieldstrengthtensor.

        """
        initialize_DtC = ""
        form_DtC = ""
        form_DtSL = ""
        form_CtD = ""
        if any([True if op.name in [opname["F"], opname["V"], opname["G"]] else False for op in self.operators]):
            form_list_DtC = []
            form_list_DtSL = []
            form_list_CtD = []
            # Initialize auxiliary spinors:
            initialize_DtC += "Function "
            for j, operator in enumerate(self.operators):
                ind2_1 = 2 * j + 1  # ind2_1 = 1,3,5,7,9,...
                ind2_2 = 2 * j + 2  # ind2_2 = 2,4,6,8,10,...

                ind3_1 = 3 * j + 1  # ind3_1 = 1,4,7,10,13,...
                ind3_2 = 3 * j + 2  # ind3_2 = 2,5,8,11,14,...
                ind3_3 = 3 * j + 3  # ind3_3 = 3,6,9,12,15,...

                ind4_1 = 4 * j + 1  # ind4_1 = 1,5, 9,13,17,...
                ind4_2 = 4 * j + 2  # ind4_2 = 2,6,10,14,18,...
                ind4_3 = 4 * j + 3  # ind4_3 = 3,7,11,15,19,...
                ind4_4 = 4 * j + 4  # ind4_4 = 4,8,12,16,20,...

                fieldstrengthtensor = ""
                if operator.name == opname["F"]:
                    # F(lor, lor) - Fieldstrengthtensor of U(1)
                    fieldstrengthtensor = opname["F"]
                    patternF = fieldstrengthtensor + r"\((?P<lor1>lor\d{1,5}),(?P<lor2>lor\d{1,5})\)"
                    matchF = re.search(patternF, operator.expression)
                    if matchF:
                        lor = [matchF.group("lor1"), matchF.group("lor2")]
                        auxfield = f"[{fieldstrengthtensor:s}({','.join([operator.id] + lor):s})]"
                    else:
                        logger.error(f"Operator {operator.expression} couldn't be matched as an fieldstrengthtensor {fieldstrengthtensor}.")
                        sys.exit("STOP")
                elif operator.name == opname["V"]:
                    # V(gaugeadj, lor, lor) - Fieldstrengthtensor of SU(2)
                    fieldstrengthtensor = opname["V"]
                    patternV = fieldstrengthtensor + r"\((?P<gaugeadjindex>gaugeadj\d{1,5}),(?P<lor1>lor\d{1,5}),(?P<lor2>lor\d{1,5})\)"
                    matchV = re.search(patternV, operator.expression)
                    if matchV:
                        lor = [matchV.group("lor1"), matchV.group("lor2")]
                        gaugeadj = matchV.group("gaugeadjindex")
                        auxfield = f"[{fieldstrengthtensor:s}({','.join([operator.id] + [gaugeadj] + lor):s})]"
                    else:
                        logger.error(f"Operator {operator.expression} couldn't be matched as an fieldstrengthtensor {fieldstrengthtensor}.")
                        sys.exit("STOP")
                elif operator.name == opname["G"]:
                    # G(cola, lor, lor) - Fieldstrengthtensor of SU(3)
                    fieldstrengthtensor = opname["G"]
                    patternG = fieldstrengthtensor + r"\((?P<colaindex>cola\d{1,5}),(?P<lor1>lor\d{1,5}),(?P<lor2>lor\d{1,5})\)"
                    matchG = re.search(patternG, operator.expression)
                    if matchG:
                        lor = [matchG.group("lor1"), matchG.group("lor2")]
                        cola = matchG.group("colaindex")
                        auxfield = f"[{fieldstrengthtensor:s}({','.join([operator.id] + [cola] + lor):s})]"
                    else:
                        logger.error(f"Operator {operator.expression} couldn't be matched as an fieldstrengthtensor {fieldstrengthtensor}.")
                        sys.exit("STOP")
                #
                # Replacements of fieldstrengthtensor by commuting fieldstrengthtensor and dummy fieldstrengthtensor:
                #
                if fieldstrengthtensor:
                    initialize_DtC += f"{auxfield}, "
                    if fieldstrengthtensor == opname["F"]:
                        for j in range(0, nDer + 1):
                            form_temp = "id once " + Term.derivativesLorentz(f"{fieldstrengthtensor:s}({','.join(lor):s})", n=j)
                            form_temp += " = " + Term.derivativesLorentz(auxfield, n=j, questionmark=False)
                            form_temp += "*" + f"(-i_/4)*({opname['F']:s}CR({operator.id:s}, Usldot{ind2_1:d}, Usldot{ind2_2:d})*sigmabar2({lor[0]:s}, {lor[1]:s}, Lsldot{ind2_1:d}, Lsldot{ind2_2:d}) " \
                                                     f"- {opname['F']:s}CL({operator.id:s}, Lsl{ind2_1:d}, Lsl{ind2_2:d})*sigma2({lor[0]:s}, {lor[1]:s}, Usl{ind2_1:d}, Usl{ind2_2:d}));"
                            form_list_DtC.append(form_temp)
                        operator.indices += [f"Usldot{ind2_1:d}",
                                             f"Usldot{ind2_2:d}",
                                             f"Lsldot{ind2_1:d}",
                                             f"Lsldot{ind2_2:d}",
                                             f"Lsl{ind2_1:d}",
                                             f"Lsl{ind2_2:d}",
                                             f"Usl{ind2_1:d}",
                                             f"Usl{ind2_2:d}"]
                    elif fieldstrengthtensor == opname["V"]:
                        # TODO: Change index order, so that SL2C-indices are always the first.
                        for j in range(0, nDer + 1):
                            # Adjoint gauge index substitution is of the following form:
                            # W(gaugeadj1,lor1,lor2) = T(gaugeadj1,gauge2,gauge1)*[su2eps](gauge2,gauge3)*W(gauge1,gauge3,lor1,lor2)
                            form_temp = "id once " + Term.derivativesLorentz(f"{fieldstrengthtensor:s}({','.join([gaugeadj] + lor):s})", n=j)
                            form_temp += " = " + Term.derivativesLorentz(auxfield, n=j, questionmark=False)
                            form_temp += "*" + f"[su2eps](gaugeA{ind3_3:d},gaugeA{ind3_2:d})*T({gaugeadj:s},gaugeA{ind3_2:d},gaugeA{ind3_1:d})*(-i_/4)" \
                                         f"*({opname['V']:s}CR({operator.id:s}, gaugeA{ind3_1:d}, gaugeA{ind3_3:d}, Usldot{ind3_1:d}, Usldot{ind3_2:d})*sigmabar2({lor[0]:s}, {lor[1]:s}, Lsldot{ind3_1:d}, Lsldot{ind3_2:d})" \
                                         f"- {opname['V']:s}CL({operator.id:s}, gaugeA{ind3_1:d}, gaugeA{ind3_3:d}, Lsl{ind3_1:d}, Lsl{ind3_2:d})*sigma2({lor[0]:s}, {lor[1]:s}, Usl{ind3_1:d}, Usl{ind3_2:d}));"
                            form_list_DtC.append(form_temp)
                        operator.indices += [f"Usldot{ind3_1:d}",
                                             f"Usldot{ind3_2:d}",
                                             f"Lsldot{ind3_1:d}",
                                             f"Lsldot{ind3_2:d}",
                                             f"Lsl{ind3_1:d}",
                                             f"Lsl{ind3_2:d}",
                                             f"Usl{ind3_1:d}",
                                             f"Usl{ind3_2:d}",
                                             f"gaugeA{ind3_1:d}",
                                             f"gaugeA{ind3_1:d}",
                                             f"gaugeA{ind3_2:d}",
                                             f"gaugeA{ind3_2:d}",
                                             f"gaugeA{ind3_3:d}",
                                             f"gaugeA{ind3_3:d}"]
                    elif fieldstrengthtensor == opname["G"]:
                        for j in range(0, nDer + 1):
                            # Adjoint gauge index substitution is of the following form:
                            # G(cola1,lor1,lor2) = (1/2)*T(cola1,colf4,colf2)*[su3eps](colf1,colf3,colf4)*G(colf1,colf2,colf3,lor1,lor2)
                            form_temp = "id once " + Term.derivativesLorentz(f"{fieldstrengthtensor:s}({','.join([cola] + lor):s})", n=j)
                            form_temp += " = " + Term.derivativesLorentz(auxfield, n=j, questionmark=False)
                            form_temp += "*" + f"(1/2)*T({cola:s},colfA{ind4_4:d},colfA{ind4_2:d})*[su3eps](colfA{ind4_1:d},colfA{ind4_3:d},colfA{ind4_4:d})*(-i_/4)" \
                                         f"*({opname['G']:s}CR({operator.id:s}, colfA{ind4_1:d}, colfA{ind4_2:d}, colfA{ind4_3:d}, Usldot{ind4_1:d}, Usldot{ind4_2:d})*sigmabar2({lor[0]:s}, {lor[1]:s}, Lsldot{ind4_1:d}, Lsldot{ind4_2:d})" \
                                         f"- {opname['G']:s}CL({operator.id:s}, colfA{ind4_1:d}, colfA{ind4_2:d}, colfA{ind4_3:d}, Lsl{ind4_1:d}, Lsl{ind4_2:d})*sigma2({lor[0]:s}, {lor[1]:s}, Usl{ind4_1:d}, Usl{ind4_2:d}));"
                            form_list_DtC.append(form_temp)
                        operator.indices += [f"Usldot{ind4_1:d}",
                                             f"Usldot{ind4_2:d}",
                                             f"Lsldot{ind4_1:d}",
                                             f"Lsldot{ind4_2:d}",
                                             f"Lsl{ind4_1:d}",
                                             f"Lsl{ind4_2:d}",
                                             f"Usl{ind4_1:d}",
                                             f"Usl{ind4_2:d}",
                                             f"colfA{ind4_1:d}",
                                             f"colfA{ind4_1:d}",
                                             f"colfA{ind4_2:d}",
                                             f"colfA{ind4_2:d}",
                                             f"colfA{ind4_3:d}",
                                             f"colfA{ind4_3:d}",
                                             f"colfA{ind4_4:d}",
                                             f"colfA{ind4_4:d}"]

                    #
                    # Replacements of derivative acting on a fieldstrengthtensor by SL2C derivative and sigma matrix:
                    # form_fieldstrengthtensorDerivativetoSL2C -> DtSL
                    #
                    if fieldstrengthtensor == opname["F"]:
                        for j in range(1, nDer + 1):
                            form_temp_DtSL = "id once " + Term.derivativesLorentz(auxfield, n=j)
                            form_temp_DtSL += f"*{opname['F']:s}CR({operator.id:s}, Usldot{ind2_1:d}, Usldot{ind2_2:d})"
                            form_temp_DtSL += " = "
                            form_temp_DtSL += Term.derivativesSL(auxfield, n=j, loop=operator.id)
                            form_temp_DtSL += f"*{opname['F']:s}CR({operator.id:s}, Usldot{ind2_1:d}, Usldot{ind2_2:d});"
                            form_list_DtSL.append(form_temp_DtSL)
                            form_temp_DtSL = "id once " + Term.derivativesLorentz(auxfield, n=j)
                            form_temp_DtSL += f"*{opname['F']:s}CL({operator.id:s}, Lsl{ind2_1:d}, Lsl{ind2_2:d})"
                            form_temp_DtSL += " = "
                            form_temp_DtSL += Term.derivativesSL(auxfield, n=j, loop=operator.id)
                            form_temp_DtSL += f"*{opname['F']:s}CL({operator.id:s}, Lsl{ind2_1:d}, Lsl{ind2_2:d});"
                            form_list_DtSL.append(form_temp_DtSL)
                            operator.indices += [f"Lsl{operator.id}{j:d}",
                                                 f"Usl{operator.id}{j:d}",
                                                 f"Lsldot{operator.id}{j:d}",
                                                 f"Usldot{operator.id}{j:d}"]
                    elif fieldstrengthtensor == opname["V"]:
                        for j in range(1, nDer + 1):
                            form_temp_DtSL = "id once " + Term.derivativesLorentz(auxfield, n=j)
                            form_temp_DtSL += f"*{opname['V']:s}CR({operator.id:s}, gaugeA{ind3_1:d}, gaugeA{ind3_3:d}, Usldot{ind3_1:d}, Usldot{ind3_2:d})"
                            form_temp_DtSL += " = "
                            form_temp_DtSL += Term.derivativesSL(auxfield, n=j, loop=operator.id)
                            form_temp_DtSL += f"*{opname['V']:s}CR({operator.id:s}, gaugeA{ind3_1:d}, gaugeA{ind3_3:d}, Usldot{ind3_1:d}, Usldot{ind3_2:d});"
                            form_list_DtSL.append(form_temp_DtSL)
                            form_temp_DtSL = "id once " + Term.derivativesLorentz(auxfield, n=j)
                            form_temp_DtSL += f"*{opname['V']:s}CL({operator.id:s}, gaugeA{ind3_1:d}, gaugeA{ind3_3:d}, Lsl{ind3_1:d}, Lsl{ind3_2:d})"
                            form_temp_DtSL += " = "
                            form_temp_DtSL += Term.derivativesSL(auxfield, n=j, loop=operator.id)
                            form_temp_DtSL += f"*{opname['V']:s}CL({operator.id:s}, gaugeA{ind3_1:d}, gaugeA{ind3_3:d}, Lsl{ind3_1:d}, Lsl{ind3_2:d});"
                            form_list_DtSL.append(form_temp_DtSL)
                            operator.indices += [f"Lsl{operator.id}{j:d}",
                                                 f"Usl{operator.id}{j:d}",
                                                 f"Lsldot{operator.id}{j:d}",
                                                 f"Usldot{operator.id}{j:d}"]
                    elif fieldstrengthtensor == opname["G"]:
                        for j in range(1, nDer + 1):
                            form_temp_DtSL = "id once " + Term.derivativesLorentz(auxfield, n=j)
                            form_temp_DtSL += f"*{opname['G']:s}CR({operator.id:s}, colfA{ind4_1:d}, colfA{ind4_2:d}, colfA{ind4_3:d}, Usldot{ind4_1:d}, Usldot{ind4_2:d})"
                            form_temp_DtSL += " = "
                            form_temp_DtSL += Term.derivativesSL(auxfield, n=j, loop=operator.id)
                            form_temp_DtSL += f"*{opname['G']:s}CR({operator.id:s}, colfA{ind4_1:d}, colfA{ind4_2:d}, colfA{ind4_3:d}, Usldot{ind4_1:d}, Usldot{ind4_2:d});"
                            form_list_DtSL.append(form_temp_DtSL)
                            form_temp_DtSL = "id once " + Term.derivativesLorentz(auxfield, n=j)
                            form_temp_DtSL += f"*{opname['G']:s}CL({operator.id:s}, colfA{ind4_1:d}, colfA{ind4_2:d}, colfA{ind4_3:d}, Lsl{ind4_1:d}, Lsl{ind4_2:d})"
                            form_temp_DtSL += " = "
                            form_temp_DtSL += Term.derivativesSL(auxfield, n=j, loop=operator.id)
                            form_temp_DtSL += f"*{opname['G']:s}CL({operator.id:s}, colfA{ind4_1:d}, colfA{ind4_2:d}, colfA{ind4_3:d}, Lsl{ind4_1:d}, Lsl{ind4_2:d});"
                            form_list_DtSL.append(form_temp_DtSL)
                            operator.indices += [f"Lsl{operator.id}{j:d}",
                                                 f"Usl{operator.id}{j:d}",
                                                 f"Lsldot{operator.id}{j:d}",
                                                 f"Usldot{operator.id}{j:d}"]
                    #
                    # Replacements of commuting fieldstrengthtensor and dummy fieldstrengthtensor by fieldstrengthtensor:
                    # form_fieldstrengthtensorCommutativetoDerivative -> CtD
                    #
                    form_temp_CtD = ""
                    if fieldstrengthtensor == opname["F"]:
                        for i in range(nDer + 1):
                            lHS, rHS = Term.derivativesSLbackwards(auxfield,
                                                                   f"{opname['F']:s}L(Lsl{ind2_1:d}, Lsl{ind2_2:d})",
                                                                   n=i)
                            form_temp_CtD = "\t" + "id once " + lHS
                            form_temp_CtD += f"*{opname['F']:s}CL({operator.id:s}, Lsl{ind2_1:d}, Lsl{ind2_2:d})"
                            form_temp_CtD += " = " + rHS + ";"
                            form_list_CtD.append(form_temp_CtD)
                            lHS, rHS = Term.derivativesSLbackwards(auxfield,
                                                                   f"{opname['F']:s}R(Usldot{ind2_1:d}, Usldot{ind2_2:d})",
                                                                   n=i)
                            form_temp_CtD = "\t" + "id once " + lHS
                            form_temp_CtD += f"*{opname['F']:s}CR({operator.id:s}, Usldot{ind2_1:d}, Usldot{ind2_2:d})"
                            form_temp_CtD += " = " + rHS + ";"
                            form_list_CtD.append(form_temp_CtD)
                    elif fieldstrengthtensor == opname["V"]:
                        for i in range(nDer + 1):
                            lHS, rHS = Term.derivativesSLbackwards(auxfield,
                                                                   f"{opname['V']:s}L(gaugeA{ind3_1:d}, gaugeA{ind3_3:d}, Lsl{ind3_1:d}, Lsl{ind3_2:d})",
                                                                   n=i)
                            form_temp_CtD = "\t" + "id once " + lHS
                            form_temp_CtD += f"*{opname['V']:s}CL({operator.id:s}, gaugeA{ind3_1:d}, gaugeA{ind3_3:d}, Lsl{ind3_1:d}, Lsl{ind3_2:d})"
                            form_temp_CtD += " = " + rHS + ";"
                            form_list_CtD.append(form_temp_CtD)
                            lHS, rHS = Term.derivativesSLbackwards(auxfield,
                                                                   f"{opname['V']:s}R(gaugeA{ind3_1:d}, gaugeA{ind3_3:d}, Usldot{ind3_1:d}, Usldot{ind3_2:d})",
                                                                   n=i)
                            form_temp_CtD = "\t" + "id once " + lHS
                            form_temp_CtD += f"*{opname['V']:s}CR({operator.id:s}, gaugeA{ind3_1:d}, gaugeA{ind3_3:d}, Usldot{ind3_1:d}, Usldot{ind3_2:d})"
                            form_temp_CtD += " = " + rHS + ";"
                            form_list_CtD.append(form_temp_CtD)
                    elif fieldstrengthtensor == opname["G"]:
                        for i in range(nDer + 1):
                            lHS, rHS = Term.derivativesSLbackwards(auxfield,
                                                                   f"{opname['G']:s}L(colfA{ind4_1:d}, colfA{ind4_2:d}, colfA{ind4_3:d}, Lsl{ind4_1:d}, Lsl{ind4_2:d})",
                                                                    n = i)
                            form_temp_CtD = "\t" + "id once " + lHS
                            form_temp_CtD += f"*{opname['G']:s}CL({operator.id:s}, colfA{ind4_1:d}, colfA{ind4_2:d}, colfA{ind4_3:d}, Lsl{ind4_1:d}, Lsl{ind4_2:d})"
                            form_temp_CtD += " = " + rHS + ";"
                            form_list_CtD.append(form_temp_CtD)
                            lHS, rHS = Term.derivativesSLbackwards(auxfield,
                                                                   f"{opname['G']:s}R(colfA{ind4_1:d}, colfA{ind4_2:d}, colfA{ind4_3:d}, Usldot{ind4_1:d}, Usldot{ind4_2:d})",
                                                                    n = i)
                            form_temp_CtD = "\t" + "id once " + lHS
                            form_temp_CtD += f"*{opname['G']:s}CR({operator.id:s}, colfA{ind4_1:d}, colfA{ind4_2:d}, colfA{ind4_3:d}, Usldot{ind4_1:d}, Usldot{ind4_2:d})"
                            form_temp_CtD += " = " + rHS + ";"
                            form_list_CtD.append(form_temp_CtD)

            # Replace last comma by semicolon.
            initialize_DtC = initialize_DtC[:-2] + ";\n"
            # Form expression of the function which should be stored in "fieldstrengthtensorDerivativetoCommutative.prc":
            form_DtC += "#procedure fieldstrengthtensorDerivativetoCommutative\n"
            form_DtC += "* Replacements of fieldstrengthtensor by commuting fieldstrengthtensor and dummy fieldstrengthtensor.\n"
            for i in form_list_DtC:
                form_DtC += 1 * "\t" + i + "\n"
            form_DtC += "#endprocedure"

            # Form expression of the function which should be stored in "fieldstrengthtensorDerivativetoSL2C.prc":
            form_DtSL += "#procedure fieldstrengthtensorDerivativetoSL2C\n"
            form_DtSL += "* Replacements of derivative acting on a fieldstrengthtensor by SL2C derivative and sigma matrix.\n"
            for i in form_list_DtSL:
                form_DtSL += 1 * "\t" + i + "\n"
            form_DtSL += "#endprocedure"

            # Form expression of the function which should be stored in "fieldstrengthtensorDerivativetoSL2C.prc":
            form_CtD += "#procedure fieldstrengthtensorCommutativetoDerivative\n"
            form_CtD += "* Replacements of commuting fieldstrengthtensor and dummy fieldstrengthtensor by fieldstrengthtensor.\n"
            for i in form_list_CtD:
                form_CtD += 1 * "\t" + i + "\n"
            form_CtD += "#endprocedure"
            return initialize_DtC, form_DtC, form_DtSL, form_CtD

    @staticmethod
    def derivativesLorentz(argument, n=4, indOffset=0, questionmark=True):
        d = ""
        for i in range(n + indOffset, 0 + indOffset, -1):  # Inverted range
            if questionmark:
                d += f"D(lorA{i:d}?lor, "
            else:
                d += f"D(lorA{i:d}, "
        d += argument
        d += n * ")"
        return d

    @staticmethod
    def derivativesSL(argument, n=4, indOffset=0, loop=""):
        if n == 0:
            d = ""
        else:
            d = f"(-(1/2))^{n:d}*"
        for i in range(n + indOffset, 0 + indOffset, -1):  # Inverted range
            d += f"D(Lsl{loop:s}{i:d}, Usldot{loop:s}{i:d}, "
        d += argument
        d += n * ")"
        for i in range(n + indOffset, 0 + indOffset, -1):  # Inverted range
            d += f"*sigma(lorA{i:d}, Usl{loop:s}{i:d}, Lsldot{loop:s}{i:d})"
        return d

    @staticmethod
    def derivativesSLbackwards(argumentLeft, argumentRight, n=4):
        abc = "abcdefghijklmnopqrstuvwxyz"
        lHS = ""
        rHS = ""
        for a in abc[:n]:
            lHS += f"D(?{a:s}, "
        lHS += argumentLeft
        lHS += n * ")"

        for a in abc[:n]:
            rHS += f"D(?{a:s}, "
        rHS += argumentRight
        rHS += n * ")"

        return lHS, rHS

    @staticmethod
    def derivativestoInd(argumentLeft, n=4):
        """
        'x', 'y', 'z' are free to use them for internal field indices.
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
            lHS += f"D(?{a:s}, "
        lHS += argumentLeft
        lHS += n * ")"

        for a in abc[:n]:
            rHS_indices += f"?{a:s}, "
        return lHS, rHS_indices

    @staticmethod
    def indtoDerivatives(n=4):
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
            j1 = 2 * i - 1
            j2 = 2*i
            lHS_indices.append(f"Lsl{j1:d}?Lsl, Usldot{j2:d}?Usldot")
            rHS.append(f"D(Lsl{j1:d}, Usldot{j2:d}")
        rHS = ", ".join(rHS)
        rHS += ", "
        brackets = n * ")"

        lHS_indices = ", ".join(lHS_indices)

        return lHS_indices, rHS, brackets

    def form_higgsDerivativetoSL2C(self, nDer=4):
        """
        Replace derivatives acting on a Higgsfields SL2C notated derivative and sigma matrices.
        Returns
        -------

        """
        if len(self.sorted_indices["gauge"]) != 0:
            form_list = []
            form_list_conj = []
            for i in range(1, nDer + 1):
                form_temp = "id " + Term.derivativesLorentz("[H('gaugeindex1')]", n=i)
                form_temp += "*HC('gaugeindex2')"
                form_temp += " = "
                form_temp += Term.derivativesSL("[H('gaugeindex1')]", n=i, loop="'gaugeindex1'")
                form_temp += "*HC('gaugeindex2');"
                form_list.append(form_temp)
            for i in range(1, nDer + 1):
                form_temp = "id " + Term.derivativesLorentz("[[H+]('gaugeindex1')]", n=i, indOffset=nDer)
                form_temp += "*[HC+]('gaugeindex2')"
                form_temp += " = "
                form_temp += Term.derivativesSL("[[H+]('gaugeindex1')]", n=i, indOffset=nDer, loop="'gaugeindex1'")
                form_temp += "*[HC+]('gaugeindex2');"
                form_list_conj.append(form_temp)
            form = ""
            form += "#procedure higgsDerivativetoSL2C\n"
            form += "#do gaugeindex1={"
            for v in self.sorted_indices["gauge"]:
                form += v + ","
            # remove last comma
            form = form[:-1]
            form += "}\n"
            form += "\t#do gaugeindex2={"
            for v in self.sorted_indices["gauge"]:
                form += v + ","
                form += v + "a" + ","
                # form += v + "b" + ","
            # remove last comma
            form = form[:-1]
            form += "}\n"
            form += "* Replacements of ordinary Higgs field - Indices range form 1 to 'maximum number of derivatives'.\n"
            for i in form_list:
                form += 2 * "\t" + i + "\n"
            form += "* Replacements of Hermitian conjugated field - Indices range form 'maximum number of derivatives' + 1 to 2*'maximum number of derivatives'.\n"
            for i in form_list_conj:
                form += 2 * "\t" + i + "\n"
            form += "\t#enddo\n"
            form += "#enddo\n"
            form += "#endprocedure"
            return form
        else:
            return

    def form_higgsDeclarations(self, nDer=4):
        """
        Write needed indices in corresponding lists in order to automatically generate all declarations for indices
        needed for derivative replacement of Higgs fields.
        Returns
        -------

        """
        if len(self.sorted_indices["gauge"]) != 0:
            for gauge in self.sorted_indices["gauge"]:
                for nindex in range(1, 2 * nDer + 1):
                    self.possible_indices["Usl"].append("Usl{gaugeindex:s}{i:d}".format(gaugeindex=gauge, i=nindex))
                    self.possible_indices["Lsl"].append("Lsl{gaugeindex:s}{i:d}".format(gaugeindex=gauge, i=nindex))
                    self.possible_indices["Usldot"].append(
                        "Usldot{gaugeindex:s}{i:d}".format(gaugeindex=gauge, i=nindex))
                    self.possible_indices["Lsldot"].append(
                        "Lsldot{gaugeindex:s}{i:d}".format(gaugeindex=gauge, i=nindex))

    def form_replaceSigmabyEps(self):
        """
        Replace contracted sigmas by SL2C epsilon tensors.
        Returns
        -------

        """
        form = ""
        form += "#procedure replaceSigmabyEps\n"
        form += "repeat;\n"
        form += "* Replace sigmabar by sigma.\n"
        form += "\t" + "id sigmabar(?a, Lsldot1?Lsldot, Usl2?Usl, ?b) = sigma(?a, Usl2, Lsldot1, ?b);\n"
        form += "\t" + "id sigmabar(?a, Lsldot1?Lsldot, Lsl2?Lsl, ?b) = sigma(?a, Lsl2, Lsldot1, ?b);\n"
        form += "\t" + "id sigmabar(?a, Usldot1?Usldot, Usl2?Usl, ?b) = sigma(?a, Usl2, Usldot1, ?b);\n"
        form += "\t" + "id sigmabar(?a, Usldot1?Usldot, Lsl2?Lsl, ?b) = sigma(?a, Lsl2, Usldot1, ?b);\n"
        form += "* Replace in lorentz indices contracted sigmas by SL2C epsilontensors.\n"
        form += "\t" + "id sigma(lor1?, Usl1?Usl, Lsldot1?Lsldot) * sigma(lor1?, Usl2?Usl, Lsldot2?Lsldot) = - 2 * [sl2Ceps](Usl1, Usl2) * [sl2Ceps](Lsldot1, Lsldot2);\n"
        form += "\t" + "id sigma(lor1?, Usl1?Usl, Usldot1?Usldot) * sigma(lor1?, Usl2?Usl, Usldot2?Usldot) = + 2 * [sl2Ceps](Usl1, Usl2) * [sl2Ceps](Usldot1, Usldot2);\n"
        form += "\t" + "id sigma(lor1?, Lsl1?Lsl, Lsldot1?Lsldot) * sigma(lor1?, Lsl2?Lsl, Lsldot2?Lsldot) = + 2 * [sl2Ceps](Lsl1, Lsl2) * [sl2Ceps](Lsldot1, Lsldot2);\n"
        form += "\t" + "id sigma(lor1?, Lsl1?Lsl, Usldot1?Usldot) * sigma(lor1?, Lsl2?Lsl, Usldot2?Usldot) = - 2 * [sl2Ceps](Lsl1, Lsl2) * [sl2Ceps](Usldot1, Usldot2);\n"
        form += "*\n"
        form += "\t" + "id sigma(lor1?, Usl1?Usl, Usldot1?Usldot) * sigma(lor1?, Usl2?Usl, Lsldot2?Lsldot) = - 2 * [sl2Ceps](Usl1, Usl2) * [sl2CdK](Usldot1, Lsldot2);\n"
        form += "\t" + "id sigma(lor1?, Usl1?Usl, Lsldot1?Lsldot) * sigma(lor1?, Usl2?Usl, Usldot2?Usldot) = + 2 * [sl2Ceps](Usl1, Usl2) * [sl2CdK](Lsldot1, Usldot2);\n"
        form += "\t" + "id sigma(lor1?, Lsl1?Lsl, Usldot1?Usldot) * sigma(lor1?, Lsl2?Lsl, Lsldot2?Lsldot) = + 2 * [sl2Ceps](Lsl1, Lsl2) * [sl2CdK](Usldot1, Lsldot2);\n"
        form += "\t" + "id sigma(lor1?, Lsl1?Lsl, Lsldot1?Lsldot) * sigma(lor1?, Lsl2?Lsl, Usldot2?Usldot) = - 2 * [sl2Ceps](Lsl1, Lsl2) * [sl2CdK](Lsldot1, Usldot2);\n"
        form += "*\n"
        form += "\t" + "id sigma(lor1?, Usl1?Usl, Lsldot1?Lsldot) * sigma(lor1?, Lsl2?Lsl, Lsldot2?Lsldot) = + 2 * [sl2CdK](Usl1, Lsl2) * [sl2Ceps](Lsldot1, Lsldot2);\n"
        form += "\t" + "id sigma(lor1?, Usl1?Usl, Usldot1?Usldot) * sigma(lor1?, Lsl2?Lsl, Usldot2?Usldot) = - 2 * [sl2CdK](Usl1, Lsl2) * [sl2Ceps](Usldot1, Usldot2);\n"
        form += "\t" + "id sigma(lor1?, Usl1?Usl, Usldot1?Usldot) * sigma(lor1?, Lsl2?Lsl, Lsldot2?Lsldot) = + 2 * [sl2CdK](Usl1, Lsl2) * [sl2CdK](Usldot1, Lsldot2);\n"
        form += "\t" + "id sigma(lor1?, Usl1?Usl, Lsldot1?Lsldot) * sigma(lor1?, Lsl2?Lsl, Usldot2?Usldot) = - 2 * [sl2CdK](Usl1, Lsl2) * [sl2CdK](Lsldot1, Usldot2);\n"
        form += "*\n"
        form += "\t" + "id sigma(lor1?, Lsl1?Lsl, Lsldot1?Lsldot) * sigma(lor1?, Usl2?Usl, Lsldot2?Lsldot) = - 2 * [sl2CdK](Lsl1, Usl2) * [sl2Ceps](Lsldot1, Lsldot2);\n"
        form += "\t" + "id sigma(lor1?, Lsl1?Lsl, Usldot1?Usldot) * sigma(lor1?, Usl2?Usl, Usldot2?Usldot) = + 2 * [sl2CdK](Lsl1, Usl2) * [sl2Ceps](Usldot1, Usldot2);\n"
        form += "\t" + "id sigma(lor1?, Lsl1?Lsl, Usldot1?Usldot) * sigma(lor1?, Usl2?Usl, Lsldot2?Lsldot) = - 2 * [sl2CdK](Lsl1, Usl2) * [sl2CdK](Usldot1, Lsldot2);\n"
        form += "\t" + "id sigma(lor1?, Lsl1?Lsl, Lsldot1?Lsldot) * sigma(lor1?, Usl2?Usl, Usldot2?Usldot) = + 2 * [sl2CdK](Lsl1, Usl2) * [sl2CdK](Lsldot1, Usldot2);\n"
        form += "endrepeat;\n"
        form += "#endprocedure\n"

        return form

    def form_simplifyEps(self):
        r"""
        Simplify contracted SL2C epsilon tensors:
        Examples:
        :math:`\epsilon^{\alpha \beta} \epsilon_{\beta \gamma} = \delta^{\alpha}_{\gamma}` and
        :math:`\delta^{\alpha}_{\alpha} = 2`
        Returns
        -------

        """
        form = ""
        form += "#procedure simplifyEps\n"
        form += "repeat;\n"
        form += "* Replace epsilons by Kronecker deltas [sl2CdK](,):\n"
        form += "\t" + "id [sl2Ceps](Usl1?Usl, Usl2?ULsl[k]) * [sl2Ceps](Lsl1?LUsl[k], Lsl3?Lsl) = + [sl2CdK](Usl1,Lsl3);\n"
        form += "\t" + "id [sl2Ceps](Usl1?ULsl[k], Usl2?Usl) * [sl2Ceps](Lsl1?LUsl[k], Lsl3?Lsl) = - [sl2CdK](Usl2,Lsl3);\n"
        form += "\t" + "id [sl2Ceps](Lsl1?Lsl, Lsl2?LUsl[k]) * [sl2Ceps](Usl2?ULsl[k], Usl3?Usl) = + [sl2CdK](Usl3,Lsl1);\n"
        form += "\t" + "id [sl2Ceps](Lsl1?Lsl, Lsl2?LUsl[k]) * [sl2Ceps](Usl3?Usl, Usl2?ULsl[k]) = - [sl2CdK](Usl3,Lsl1);\n"
        form += "\t" + "id [sl2Ceps](Usldot1?Usldot, Usldot2?ULsldot[k]) * [sl2Ceps](Lsldot1?LUsldot[k], Lsldot3?Lsldot) = + [sl2CdK](Usldot1,Lsldot3);\n"
        form += "\t" + "id [sl2Ceps](Usldot1?ULsldot[k], Usldot2?Usldot) * [sl2Ceps](Lsldot1?LUsldot[k], Lsldot3?Lsldot) = - [sl2CdK](Usldot2,Lsldot3);\n"
        form += "\t" + "id [sl2Ceps](Lsldot1?Lsldot, Lsldot2?LUsldot[k]) * [sl2Ceps](Usldot2?ULsldot[k], Usldot3?Usldot) = + [sl2CdK](Usldot3,Lsldot1);\n"
        form += "\t" + "id [sl2Ceps](Lsldot1?Lsldot, Lsldot2?LUsldot[k]) * [sl2Ceps](Usldot3?Usldot, Usldot2?ULsldot[k]) = - [sl2CdK](Usldot3,Lsldot1);\n"
        # form += "* Replace epsilons with one upper and one lower index by Kronecker-deltas:\n"
        # form += "\t" + "id [sl2Ceps](Lsl1?Lsl,Usl2?Usl) = [sl2CdK](Lsl1,Usl2);\n"
        # form += "\t" + "id [sl2Ceps](Usl1?Usl,Lsl2?Lsl) = [sl2CdK](Lsl2,Usl1);\n"
        # form += "\t" + "id [sl2Ceps](Lsldot1?Lsldot,Usldot2?Usldot) = [sl2CdK](Lsldot1,Usldot2);\n"
        # form += "\t" + "id [sl2Ceps](Usldot1?Usldot,Lsldot2?Lsldot) = [sl2CdK](Lsldot2,Usldot1);\n"
        form += "* Replace only contracted Kronecker-deltas which are contracted with eps and Kronecker-deltas themself, because contractions inside one building block are not wanted:\n"  # Replace contracted Kronecker-deltas (works for every function, not just eps):\n
        form += "\t" + "id [sl2Ceps]?{[sl2Ceps],[sl2CdK],sigma,sigmabar}(?a,Lsl1?LUsl[k],?b)*[sl2CdK](?c,Usl1?ULsl[k],?d) = [sl2Ceps](?a,?c,?d,?b);\n"  # [sl2Ceps]?
        form += "\t" + "id [sl2Ceps]?{[sl2Ceps],[sl2CdK],sigma,sigmabar}(?a,Lsldot1?LUsldot[k],?b)*[sl2CdK](?c,Usldot1?ULsldot[k],?d) = [sl2Ceps](?a,?c,?d,?b);\n"  # [sl2Ceps]?
        form += "* Replace self-contracted Kronecker-deltas by the dimension (=2):\n"
        form += "\t" + "id [sl2CdK](Usl1?ULsl[k],Lsl1?LUsl[k]) = d_(Lsl1,Lsl1);\n"
        form += "\t" + "id [sl2CdK](Usldot1?ULsldot[k],Lsldot1?LUsldot[k]) = d_(Lsldot1,Lsldot1);\n"
        form += "endrepeat;\n"
        form += "* Bring indices of epsilons in order:\n"
        form += "Multiply replace_([sl2Ceps],[sl2CepsA]);\n"
        form += ".sort\n"
        form += "Multiply replace_([sl2CepsA],[sl2Ceps]);\n"
        form += "#endprocedure\n"

        return form

    def form_substituteSU2dKbySU2Eps(self):
        pass
        # TODO: Substitute Kronecker-deltas by 2 contracted epsilon tensors. In order to get the simplest form of
        #  epsilon tensors, the substitution is done after the simplification of epsilon tensors. If turns out that the
        #  form is not important, the lines in form_simplifyEps which introduce the epsilontensors just could be
        #  removed.

    def form_simplifySigma2(self):
        """
        Assume that only sigma2 (i.e. sigma matrices with two lorentz indices) with two upper undotted and sigmabar2
        with two lower dotted indices exist, because this is the only relevant case for the fieldstrengthtensor replacement.
        Returns
        -------

        """
        form = ""
        if len(self.sorted_indices["Usl"]) == 0 and len(self.sorted_indices["Lsldot"]) == 0:
            return
        # .pop(0)
        form += "#procedure simplifySigma2\n"
        uppersl = self.sorted_indices["Usl"]
        if len(uppersl) != 0:
            ind = []
            for usl1 in uppersl:
                match1 = re.match(r"Usl(\d{1,3})", usl1)
                if match1:
                    sl1N = int(match1.group(1))
                    ind.append(sl1N)
            ind.sort()
            # Ind contains now a sorted list with integer numbers, which denote the index number of the SL2C indices. E.g. Usl5 has number 5 and so on.
            for sl1N in ind:
                self.possible_indices["Usl"].append("Usl{0:d}a".format(sl1N))
                self.possible_indices["Lsl"].append("Lsl{0:d}a".format(sl1N))
                self.possible_indices["Usldot"].append("Usldot{0:d}a".format(sl1N))
                self.possible_indices["Lsldot"].append("Lsldot{0:d}a".format(sl1N))
            for usl1 in uppersl:
                for usl2 in uppersl:
                    if usl1 != usl2:
                        match1 = re.match(r"Usl(\d{1,3})", usl1)
                        if match1:
                            sl1N = int(match1.group(1))
                            form += "id once sigma2(lor1?,lor2?,{0:s},{1:s}) = [sl2Ceps]({0:s},Usl{slj:d}a)*(i_/2)*(sigma(lor1,Lsl{slj:d}a,Lsldot{slj:d}a)*sigmabar(lor2,Usldot{slj:d}a,{1:s}) - sigma(lor2, Lsl{slj:d}a, Lsldot{slj:d}a)*sigmabar(lor1,Usldot{slj:d}a,{1:s}));\n".format(
                                usl1, usl2, slj=sl1N)

        lowersldot = self.sorted_indices["Lsldot"]
        if len(lowersldot) != 0:
            ind = []
            for lsldot1 in lowersldot:
                match1 = re.match(r"Lsldot(\d{1,3})", lsldot1)
                if match1:
                    sldot1N = int(match1.group(1))
                    ind.append(sldot1N)
            ind.sort()
            # Ind contains now a sorted list with integer numbers, which denote the index number of the SL2C indices. E.g. Lsldot5 has number 5 and so on.
            for sldot1N in ind:
                self.possible_indices["Usl"].append("Usl{0:d}a".format(sldot1N))
                self.possible_indices["Lsl"].append("Lsl{0:d}a".format(sldot1N))
                self.possible_indices["Usldot"].append("Usldot{0:d}a".format(sldot1N))
                self.possible_indices["Lsldot"].append("Lsldot{0:d}a".format(sldot1N))
            for lsldot1 in lowersldot:
                for lsldot2 in lowersldot:
                    if lsldot1 != lsldot2:
                        match1 = re.match(r"Lsldot(\d{1,3})", lsldot1)
                        if match1:
                            sldot1N = int(match1.group(1))
                            form += "id once sigmabar2(lor1?,lor2?,{0:s},{1:s}) = [sl2Ceps]({0:s},Lsldot{sldotj:d}a)*(i_/2)*(sigmabar(lor1,Usldot{sldotj:d}a,Usl{sldotj:d}a)*sigma(lor2,Lsl{sldotj:d}a,{1:s}) - sigmabar(lor2, Usldot{sldotj:d}a, Usl{sldotj:d}a)*sigma(lor1,Lsl{sldotj:d}a,{1:s}));\n".format(
                                lsldot1, lsldot2, sldotj=sldot1N)

        form += "#endprocedure\n"
        return form

    def form_simplifyEpsSU2(self):
        r"""
        Simplify contracted SU2 epsilon tensors:
        Examples:
        :math:`\epsilon^{a b} \epsilon_{b c} = \delta^{a}_{c}` and
        :math:`\delta^{a}_{a} = 2`
        Returns
        -------

        """
        form = ""
        if len(self.sorted_indices["gauge"]) != 0:
            form += "#procedure simplifyEpsSU2\n"
            form += "repeat;\n"
            form += "* Replace epsilons by Kronecker deltas [sl2CdK](,):\n"
            form += "\t" + "id [su2eps](gauge1?gauge, gauge2?gauge) * [su2eps](gauge2?gauge, gauge3?gauge) = + [su2dK](gauge1,gauge3);\n"
            form += "\t" + "id [su2eps](gauge1?gauge, gauge2?gauge) * [su2eps](gauge1?gauge, gauge3?gauge) = - [su2dK](gauge2,gauge3);\n"
            form += "\t" + "id [su2eps](gauge1?gauge, gauge2?gauge) * [su2eps](gauge3?gauge, gauge1?gauge) = + [su2dK](gauge2,gauge3);\n"
            form += "\t" + "id [su2eps](gauge1?gauge, gauge2?gauge) * [su2eps](gauge3?gauge, gauge2?gauge) = - [su2dK](gauge1,gauge3);\n"
            form += "* Replace only contracted Kronecker-deltas which are contracted with eps and Kronecker-deltas themself, because contractions inside one building block are not wanted:\n"
            form += "\t" + "id [su2eps]?{[su2eps],[su2dK]}(?a,gauge1?gauge,?b)*[su2dK](?c,gauge1?gauge,?d) = [su2eps](?a,?c,?d,?b);\n"  # [su2eps]?
            form += "* Replace self-contracted Kronecker-deltas by the dimension (=2):\n"
            form += "\t" + "id [su2dK](gauge1?gauge,gauge1?gauge) = d_(gauge1,gauge1);\n"
            form += "endrepeat;\n"
            form += "* Bring indices of epsilons in order:\n"
            form += "Multiply replace_([su2eps],[su2epsA]);\n"
            form += ".sort\n"
            form += "Multiply replace_([su2epsA],[su2eps]);\n"
            form += "#endprocedure\n"
            return form

    def form_replaceSU2Generators(self):
        """
        Replace two in the adjoint representation contracted SU2 Generators by Kronecker deltas in fundamental indices.
        Returns
        -------

        """
        form = ""
        if len(self.sorted_indices["gaugeadj"]) != 0:
            form += "* Replace two in the adjoint representation contracted SU2 Generators by Kronecker deltas in fundamental indices.\n"
            form += "repeat;\n"
            form += "\tid T(gaugeadj1?gaugeadj, gauge1?gauge, gauge2?gauge) * T(gaugeadj1?gaugeadj, gauge3?gauge, gauge4?gauge)"
            form += " = (1 / 2) * ([su2dK](gauge1, gauge4) * [su2dK](gauge3, gauge2) - (1 / 2) * [su2dK](gauge1, gauge2) * [su2dK](gauge3, gauge4));\n"
            form += "endrepeat;\n"
            return form

    def form_replaceSU3Generators(self):
        """
        Replace two in the adjoint representation contracted SU2 Generators by Kronecker deltas in fundamental indices.
        Returns
        -------

        """
        form = ""
        if len(self.sorted_indices["cola"]) != 0:
            form += "* Replace two in the adjoint representation contracted SU3 Generators by Kronecker deltas in fundamental indices.\n"
            form += "repeat;\n"
            form += "\tid T(cola1?cola, colf1?colf, colf2?colf) * T(cola1?cola, colf3?colf, colf4?colf)"
            form += " = (1 / 2) * ([su3dK](colf1, colf4) * [su3dK](colf3, colf2) - (1 / 3) * [su3dK](colf1, colf2) * [su3dK](colf3, colf4));\n"
            form += "endrepeat;\n"
            return form

    #
    # Spinors
    #
    def form_spinorDeclarations(self, nDer):
        """
        Write needed indices in corresponding lists in order to automatically generate all declarations for indices
        needed for derivative replacement of spinor fields.
        Define further sets of spinorfunction for id-statements.
        Returns
        -------

        """
        if len(self.sorted_indices["spin"]) != 0:
            for spin in self.sorted_indices["spin"]:
                for nindex in range(1, 2 * nDer + 1):
                    usl = "Usl{spinindex:s}{i:d}".format(spinindex=spin, i=nindex)
                    if usl not in self.possible_indices["Usl"]:
                        self.possible_indices["Usl"].append(usl)
                    lsl = "Lsl{spinindex:s}{i:d}".format(spinindex=spin, i=nindex)
                    if lsl not in self.possible_indices["Lsl"]:
                        self.possible_indices["Lsl"].append(lsl)
                    usldot = "Usldot{spinindex:s}{i:d}".format(spinindex=spin, i=nindex)
                    if usldot not in self.possible_indices["Usldot"]:
                        self.possible_indices["Usldot"].append(usldot)
                    lsldot = "Lsldot{spinindex:s}{i:d}".format(spinindex=spin, i=nindex)
                    if lsldot not in self.possible_indices["Lsldot"]:
                        self.possible_indices["Lsldot"].append(lsldot)

                usl = "Usl{spinindex:s}".format(spinindex=spin)
                if usl not in self.possible_indices["Usl"]:
                    self.possible_indices["Usl"].append(usl)
                lsl = "Lsl{spinindex:s}".format(spinindex=spin)
                if lsl not in self.possible_indices["Lsl"]:
                    self.possible_indices["Lsl"].append(lsl)
                usldot = "Usldot{spinindex:s}".format(spinindex=spin)
                if usldot not in self.possible_indices["Usldot"]:
                    self.possible_indices["Usldot"].append(usldot)
                lsldot = "Lsldot{spinindex:s}".format(spinindex=spin)
                if lsldot not in self.possible_indices["Lsldot"]:
                    self.possible_indices["Lsldot"].append(lsldot)

                for nindex in range(1, len(self.sorted_indices["spin"])):
                    if nindex == 0:
                        nindex = ""
                    else:
                        nindex = str(nindex)
                    for n in range(1, 3):
                        usl = "Usl{spinindex:s}{i:s}{fs:d}".format(spinindex=spin, i=nindex, fs=n)
                        if usl not in self.possible_indices["Usl"]:
                            self.possible_indices["Usl"].append(usl)
                        lsl = "Lsl{spinindex:s}{i:s}{fs:d}".format(spinindex=spin, i=nindex, fs=n)
                        if lsl not in self.possible_indices["Lsl"]:
                            self.possible_indices["Lsl"].append(lsl)
                        usldot = "Usldot{spinindex:s}{i:s}{fs:d}".format(spinindex=spin, i=nindex, fs=n)
                        if usldot not in self.possible_indices["Usldot"]:
                            self.possible_indices["Usldot"].append(usldot)
                        lsldot = "Lsldot{spinindex:s}{i:s}{fs:d}".format(spinindex=spin, i=nindex, fs=n)
                        if lsldot not in self.possible_indices["Lsldot"]:
                            self.possible_indices["Lsldot"].append(lsldot)

    def form_spinorDerivativetoCommutative(self, nDer=4):
        """
        Rewrite (if there exists a spin index in the term) spinor fields inside of a derivative as:
            D(lor1?, D(lor2?, [qbar]('spin', gauge1, colf2, flav3))) = D(lor1, D(lor2, [[qbar]('spin', gauge1, colf2, flav3)])) * [qCbar]('spin', gauge1, colf2, flav3),
        where [[qbar]('spin', gauge1, colf2, flav3)] is some placeholder to remember the position of the spinor field and [qCbar]('spin', gauge1, colf2, flav3)
        is a commuting field. This makes it easier to implement id-statement for contracted spinor fields.
        Note: Also a "lonely" spinor field is replaced,
            [qbar]('spin', gauge1, colf2, flav3) = [[qbar]('spin', gauge1, colf2, flav3)]*[qCbar]('spin', gauge1, colf2, flav3),
        in order to ensure also the order of the "lonely" spinor fields.

        Parameters
        ----------
        nDer : int
            Maximum number of derivatives.
        Returns
        -------
        initialize : str
            Contains essentially the declaration of the dummy function [q('spin', gauge1, colf2, flav3)] and [[qbar]('spin', gauge1, colf2, flav3)] with
            the correct indices.
        form : str
            Formatted string of the FORM procedure file "spinorDerivativetoCommutative.prc" with the correct indices.
            in the do-loop index list.

        """
        initialize = ""
        form = ""
        if len(self.sorted_indices["spin"]) != 0:
            form_list = []
            form_list_bar = []
            # Initialize auxiliary spinors:
            initialize += "Function "
            for operator in self.operators:
                spinor = ""
                if operator.name == opname["e"] or operator.name == opname["bar[e]"]:
                    # e(spin, flav) - RH leptons
                    spinor = "e"
                elif operator.name == opname["u"] or operator.name == opname["bar[u]"]:
                    # u(spin, colf, flav) - RH up-type quarks
                    spinor = "u"
                elif operator.name == opname["b"] or operator.name == opname["bar[b]"]:
                    # b(spin, colf, flav) - RH down-type quarks
                    spinor = "b"
                elif operator.name == opname["l"] or operator.name == opname["bar[l]"]:
                    # l(spin, gauge, flav) - SU2 lepton dublett
                    spinor = "l"
                elif operator.name == opname["q"] or operator.name == opname["bar[q]"]:
                    # q(spin, gauge, colf, flav) - SU2 quark dublett
                    spinor = "q"
                if spinor:
                    indices = ",".join(operator.spinorIndices)
                    initialize += "[{op:s}({i:s})], [{opbar:s}({i:s})], ".format(op=opname[spinor],
                                                                                 opbar=opname["bar[" + spinor + "]"],
                                                                                 i=indices)

                    for j in range(0, nDer + 1):
                        form_temp = "id " + Term.derivativesLorentz(
                            "{op:s}({i:s})".format(op=opname[spinor], i=indices), n=j)
                        form_temp += " = " + Term.derivativesLorentz(
                            "[{op:s}({i:s})]".format(op=opname[spinor], i=indices), n=j, questionmark=False)
                        form_temp += "*{op:s}C({i:s});".format(op=spinor, i=indices)
                        form_list.append(form_temp)
                    for j in range(0, nDer + 1):
                        form_temp = "id " + Term.derivativesLorentz(
                            "{opbar:s}({i:s})".format(opbar=opname["bar[" + spinor + "]"], i=indices), n=j)
                        form_temp += " = " + Term.derivativesLorentz(
                            "[{opbar:s}({i:s})]".format(opbar=opname["bar[" + spinor + "]"], i=indices), n=j,
                            questionmark=False)
                        form_temp += "*[{op:s}Cbar]({i:s});".format(op=spinor, i=indices)
                        form_list_bar.append(form_temp)

            # Replace last comma by semicolon.
            initialize = initialize[:-2] + ";\n"

            # Form expression of the function which should be stored in "spinorDerivativetoCommutative.prc":
            form += "#procedure spinorDerivativetoCommutative\n"
            form += "* Replacements of ordinary spinor field by commuting spinor field and dummy spinor.\n"
            for i in form_list:
                form += 1 * "\t" + i + "\n"
            form += "* Replacements of Hermitian conjugated field by commuting spinor field and dummy spinor.\n"
            for i in form_list_bar:
                form += 1 * "\t" + i + "\n"
            form += "#endprocedure"
            if initialize != "" and form != "":
                return initialize, form

    def form_spinorDerivativetoSL2C(self, nDer=4):
        """
        Replace derivatives acting on a spinorfield by SL2C notated derivative and sigma matrices.
        Place directly after the spinorDerivativetoCommutative call.
        Returns
        -------

        """
        if len(self.sorted_indices["spin"]) != 0:
            form_list = []
            form_list_bar = []
            for operator in self.operators:
                spinor = ""
                if operator.name == opname["e"] or operator.name == opname["bar[e]"]:
                    # e(spin, flav) - RH leptons
                    spinor = "e"
                elif operator.name == opname["u"] or operator.name == opname["bar[u]"]:
                    # u(spin, colf, flav) - RH up-type quarks
                    spinor = "u"
                elif operator.name == opname["b"] or operator.name == opname["bar[b]"]:
                    # b(spin, colf, flav) - RH down-type quarks
                    spinor = "b"
                elif operator.name == opname["l"] or operator.name == opname["bar[l]"]:
                    # l(spin, gauge, flav) - SU2 lepton dublett
                    spinor = "l"
                elif operator.name == opname["q"] or operator.name == opname["bar[q]"]:
                    # q(spin, gauge, colf, flav) - SU2 quark dublett
                    spinor = "q"

                if spinor:
                    indices = ",".join(operator.spinorIndices)
                    spin = operator.spinorIndices[0]
                    for j in range(0, nDer + 1):
                        form_temp = "id " + Term.derivativesLorentz(
                            "[{op:s}({i:s})]".format(op=opname[spinor], i=indices), n=j)
                        form_temp += "*{op:s}C({i:s})".format(op=spinor, i=indices)
                        form_temp += " = "
                        form_temp += Term.derivativesSL("[{op:s}({i:s})]".format(op=opname[spinor], i=indices), n=j,
                                                        loop=spin)
                        form_temp += "*{op:s}C({i:s});".format(op=spinor, i=indices)
                        form_list.append(form_temp)

                    for j in range(0, nDer + 1):
                        form_temp = "id " + Term.derivativesLorentz(
                            "[{opbar:s}({i:s})]".format(opbar=opname["bar[" + spinor + "]"], i=indices),
                            n=j, indOffset=nDer)
                        form_temp += "*[{op:s}Cbar]({i:s})".format(op=spinor, i=indices)
                        form_temp += " = "
                        form_temp += Term.derivativesSL(
                            "[{opbar:s}({i:s})]".format(opbar=opname["bar[" + spinor + "]"], i=indices), n=j,
                            indOffset=nDer, loop=spin)
                        form_temp += "*[{op:s}Cbar]({i:s});".format(op=spinor, i=indices)
                        form_list_bar.append(form_temp)

            form = ""
            form += "#procedure spinorDerivativetoSL2C\n"
            form += "* Replacements of ordinary spinor field - Indices range form 1 to 'maximum number of derivatives'.\n"
            for i in form_list:
                form += "\t" + i + "\n"
            form += "* Replacements of Hermitian conjugated field - Indices range form 'maximum number of derivatives' + 1 to 2*'maximum number of derivatives'.\n"
            for i in form_list_bar:
                form += "\t" + i + "\n"
            form += "#endprocedure"
            return form

    def form_substituteSpinors(self):
        """
        Substitute 4 Dirac spinors by SL2C-/ Weyl-spinors.
        Returns
        -------

        """
        form = ""
        if len(self.sorted_indices["spin"]) != 0:
            form += "#procedure substituteSpinors\n"
            form += "* Add auxiliary spin index\n"
            form += "id eC?spinorsAll(spin1?spin[k?],?a) = eC(spinA[k],spin1,?a);\n"
            form += "id gamma(lor1?lor,spin1?spin[k?],spin2?spin[m?]) = gamma(lor1,spinA[k],spinA[m],spin1,spin2);\n"
            form += "\n"
            form += "* Sum over all occurring spin indices.\n"
            nSpinAInd = len(self.possible_indices["spinA"])
            for i, spinA in enumerate(self.possible_indices["spinA"]):
                if i == 0:
                    form += "sum {spA:s} 1, 2;\n".format(spA=spinA)
                else:
                    form += "sum {spA:s} {ind:d}1, {ind:d}2;\n".format(ind=i, spA=spinA)
            form += "\n"
            form += """repeat;
\t$sXi= 0;
\t$sChiD = 0;
\t$sbarChi = 0;
\t$sbarXiD = 0;
\t$sXiArg= 0;
\t$sChiDArg = 0;
\t$sbarChiArg = 0;
\t$sbarXiDArg = 0;
\t$s = 0;
\t$sbar = 0;"""
            form += "\n"
            form += "* First, substitute directly contracted spinors:\n"
            form += 1 * "\t" + "#do i={," + ",".join(
                [str(i) for i in range(1, nSpinAInd)]) + "}\n"  # '#do i={,1,2,3,4}'
            form += 1 * "\t" + "#do spin={" + ",".join([spin for spin in self.sorted_indices[
                "spin"]]) + "}\n"  # '#do spin={spin6954,spin6953,spin695,spin694,spin2768}'
            form += 2 * "\t" + "id [eCbar]?spinorsAdj('i'1, 'spin', ?a) * eC?spinors('i'2, 'spin', ?b) = 0;\n"
            form += 2 * "\t" + "id [eCbar]?spinorsAdj('i'2, 'spin', ?a) * eC?spinors('i'1, 'spin', ?b) = 0;\n"
            form += "\n"
            form += """\t\tif (match([eCbar]?spinorsAdj[k?$sbar]('i'1, 'spin', ?a) * eC?spinors[m?$s]('i'1, 'spin', ?b)));
\t\t\tid [eCbar]?spinorsAdj[k?$sbar]('i'1, 'spin', ?a$b) * eC?spinors[m?$s]('i'1, 'spin', ?b$c) = [sl2Ceps](Usl'spin'1,Usl'spin'2)*chi(Lsl'spin'2,?a)*xi(Lsl'spin'1,?b);
\t\t\t$sbarChi = chi(Lsl'spin'2,$b);
\t\t\t$sXi = xi(Lsl'spin'1,$c);
\t\t\tinside $sbarChi,$sXi;
\t\t\t\tid chi(?b$sbarChiArg)=chi(?b);
\t\t\t\tid xi(?a$sXiArg)=xi(?a);
\t\t\tendinside;
\t\t\tgoto 1;
\t\tendif;
\t\tif (match([eCbar]?spinorsAdj[k?$sbar]('i'2, 'spin', ?a) * eC?spinors[m?$s]('i'2, 'spin', ?b)));
\t\t\tid [eCbar]?spinorsAdj[k?$sbar]('i'2, 'spin', ?a$b) * eC?spinors[m?$s]('i'2, 'spin', ?b$c) = [sl2Ceps](Lsldot'spin'1,Lsldot'spin'2)*[xi+](Usldot'spin'2,?a)*[chi+](Usldot'spin'1,?b);
\t\t\t$sbarXiD = [xi+](Usldot'spin'2,$b);
\t\t\t$sChiD = [chi+](Usldot'spin'1,$c);
\t\t\tinside $sbarXiD,$sChiD;
\t\t\t\tid [xi+](?a$sbarXiDArg)=xi(?a);
\t\t\t\tid [chi+](?b$sChiDArg)=chi(?b);
\t\t\tendinside;
\t\t\tgoto 1;
\t\t\tendif;"""
            form += "\n"
            form += 2 * "\t#enddo\n"
            form += "\n"
            form += "* Insert Weyl representation of gamma matrix\n"
            form += 1 * "\t" + "#do i={," + ",".join(
                [str(i) for i in range(1, nSpinAInd)]) + "}\n"  # '#do i={,1,2,3,4}'
            form += 1 * "\t" + "#do j={," + ",".join(
                [str(i) for i in range(1, nSpinAInd)]) + "}\n"  # '#do j={,1,2,3,4}'
            form += 1 * "\t" + "#do spin1={" + ",".join([spin for spin in self.sorted_indices[
                "spin"]]) + "}\n"  # '#do spin1={spin6954,spin6953,spin695,spin694,spin2768}'
            form += 1 * "\t" + "#do spin2={" + ",".join([spin for spin in self.sorted_indices[
                "spin"]]) + "}\n"  # '#do spin2={spin6954,spin6953,spin695,spin694,spin2768}'
            form += """\t\tid gamma(lor1?lor, 'i'1, 'j'1, 'spin1', 'spin2') = 0;
\t\tid gamma(lor1?lor, 'i'2, 'j'2, 'spin1', 'spin2') = 0;
\t\tid gamma(lor1?lor, 'i'1, 'j'2, 'spin1', 'spin2') = sigma(lor1, Lsl'spin1''i'1, Lsldot'spin2''j'2);
\t\tid gamma(lor1?lor, 'i'2, 'j'1, 'spin1', 'spin2') = sigmabar(lor1, Usldot'spin1''i'1, Usl'spin2''j'2);"""
            form += "\n"
            form += 4 * "\t#enddo\n"
            form += "\n"
            form += 1 * "\t" + "#do i={," + ",".join(
                [str(i) for i in range(1, nSpinAInd)]) + "}\n"  # '#do i={,1,2,3,4}'
            form += 1 * "\t" + "#do spin={" + ",".join([spin for spin in self.sorted_indices["spin"]]) + "}\n"  # '#do spin={spin6954,spin6953,spin695,spin694,spin2768}'
            form += "* Assign ordinary spinor\n"
            form += """\t\tif (match(eC?spinors('i'1,'spin',?a)));
\t\t\tid eC?spinors[k?$s]('i'1,'spin',?a$b) = xi(Lsl'spin''i'2,?a);
\t\t\t$sXi = xi(Lsl'spin''i'2, $b);
\t\t\tinside $sXi;
\t\t\t\tid xi(?a$sXiArg)=xi(?a);
\t\t\tendinside;
\t\t\tgoto 1;
\t\tendif;
\t\tif (match(eC?spinors[k?$s]('i'2,'spin',?a)));
\t\t\tid eC?spinors[k?$s]('i'2,'spin',?a$b) = [chi+](Usldot'spin''i'2,?a);
\t\t\t$sChiD = [chi+](Usldot'spin''i'2,$b);
\t\t\tinside $sChiD;
\t\t\t\tid [chi+](?a$sChiDArg)=[chi+](?a);
\t\t\tendinside;
\t\t\tgoto 1;
\t\tendif;"""
            form += "\n"
            form += "* Assign adjoint spinor\n"
            form += """\t\tif (match([eCbar]?spinorsAdj('i'1,'spin',?a)));
\t\t\tid [eCbar]?spinorsAdj[m?$sbar]('i'1,'spin',?a$b) = chi(Usl'spin''i'1,?a);
\t\t\t$sbarChi = chi(Usl'spin''i'1,$b);
\t\t\tinside $sbarChi;
\t\t\t\tid chi(?a$sbarChiArg)=chi(?a);
\t\t\tendinside;
\t\t\tgoto 1;
\t\tendif;
\t\tif (match([eCbar]?spinorsAdj('i'2,'spin',?a)));
\t\t\tid [eCbar]?spinorsAdj[m?$sbar]('i'2,'spin',?a$b) = [xi+](Lsldot'spin''i'1,?a);
\t\t\t$sbarXiD = [xi+](Lsldot'spin''i'1,$b);
\t\t\tinside $sbarXiD;
\t\t\t\tid [xi+](?a$sbarXiDArg)=[xi+](?a);
\t\t\tendinside;
\t\t\tgoto 1;
\t\tendif;"""
            form += "\n"
            form += 2 * "\t#enddo\n"
            form += "\n"
            form += "* Identify spinors and substitute by corresponding SL2C (in first place commuting auxiliary) expression:\n"
            form += "*[e_C], [u_C], [d_C], L, Q, [e_C+], [u_C+], [d_C+], [L+], [Q+];\n"
            form += "* eC = 1, uC = 2, bC = 3, lC = 4, qC = 5;\n"
            form += 1 * "\t" + "label 1;\n"
            form += """\tif (($s == 1)&&($sXiArg)); * eC
\t\tid xi($sXiArg) = 0;
\telseif (($s == 2)&&($sXiArg)); * uC
\t\tid xi($sXiArg) = 0;
\telseif (($s == 3)&&($sXiArg)); * bC
\t\tid xi($sXiArg) = 0;
\telseif (($s == 4)&&($sXiArg)); * lC
\t\tid xi($sXiArg) = {L:s}($sXiArg);
\telseif (($s == 5)&&($sXiArg)); * qC
\t\tid xi($sXiArg) = {Q:s}($sXiArg);
\tendif;""".format(L=spinorsSL2C_c["L"], Q=spinorsSL2C_c["Q"])
            form += "\n"

            form += """\tif (($s == 1)&&($sChiDArg)); * eC
\t\tid [chi+]($sChiDArg) = {eD:s}($sChiDArg);
\telseif (($s == 2)&&($sChiDArg)); * uC
\t\tid [chi+]($sChiDArg) = {uD:s}($sChiDArg);
\telseif (($s == 3)&&($sChiDArg)); * bC
\t\tid [chi+]($sChiDArg) = {dD:s}($sChiDArg);
\telseif (($s == 4)&&($sChiDArg)); * lC
\t\tid [chi+]($sChiDArg) = 0;
\telseif (($s == 5)&&($sChiDArg)); * qC
\t\tid [chi+]($sChiDArg) = 0;
\tendif;""".format(eD=spinorsSL2C_c["[e_C+]"], uD=spinorsSL2C_c["[u_C+]"], dD=spinorsSL2C_c["[d_C+]"])
            form += "\n"
            form += "\n"
            form += "* [eCbar] = 1, [uCbar] = 2, [bCbar] = 3, [lCbar] = 4, [qCbar] = 5;\n"
            form += """\tif (($sbar == 1)&&($sbarChiArg)); * [eCbar]
\t\tid chi($sbarChiArg) = {e:s}($sbarChiArg);
\telseif (($sbar == 2)&&($sbarChiArg)); * [uCbar]
\t\tid chi($sbarChiArg) = {u:s}($sbarChiArg);
\telseif (($sbar == 3)&&($sbarChiArg)); * [bCbar]
\t\tid chi($sbarChiArg) = {d:s}($sbarChiArg);
\telseif (($sbar == 4)&&($sbarChiArg)); * [lCbar]
\t\tid chi($sbarChiArg) = 0;
\telseif (($sbar == 5)&&($sbarChiArg)); * [qCbar]
\t\tid chi($sbarChiArg) = 0;
\tendif;""".format(e=spinorsSL2C_c["[e_C]"], u=spinorsSL2C_c["[u_C]"], d=spinorsSL2C_c["[d_C]"])
            form += "\n"

            form += """\tif (($sbar == 1)&&($sbarXiDArg)); * [eCbar]
\t\tid [xi+]($sbarXiDArg) = 0;
\telseif (($sbar == 2)&&($sbarXiDArg)); * [uCbar]
\t\tid [xi+]($sbarXiDArg) = 0;
\telseif (($sbar == 3)&&($sbarXiDArg)); * [bCbar]
\t\tid [xi+]($sbarXiDArg) = 0;
\telseif (($sbar == 4)&&($sbarXiDArg)); * [lCbar]
\t\tid [xi+]($sbarXiDArg) = {LD:s}($sbarXiDArg);
\telseif (($sbar == 5)&&($sbarXiDArg)); * [qCbar]
\t\tid [xi+]($sbarXiDArg) = {QD:s}($sbarXiDArg);
\tendif;""".format(LD=spinorsSL2C_c["[L+]"], QD=spinorsSL2C_c["[Q+]"])
            form += "\n"
            form += "* Raise and lower indices in the correct way.\n"
            form += 1 * "\t" + "#do i={," + ",".join(
                [str(i) for i in range(1, nSpinAInd)]) + "}\n"  # '#do i={,1,2,3,4}'
            form += "\t#do j={1,2}\n"
            form += 1 * "\t" + "#do spin={" + ",".join([spin for spin in self.sorted_indices[
                "spin"]]) + "}\n"  # '#do spin={spin6954,spin6953,spin695,spin694,spin2768}'
            form += "\t\tid {e:s}?spinorsAllN(Usl'spin''i''j',?a)".format(e=spinorsSL2C_c["[e_C]"]) + "*sigma?{sigma,sigmabar,[sl2Ceps],sigma2,sigmabar2}(?b,Lsl'spin''i''j',?c) =" + " -{e:s}(Lsl'spin''i''j',?a)*sigma(?b,Usl'spin''i''j',?c);\n".format(e=spinorsSL2C_c["[e_C]"])
            form += "\t\tid {e:s}?spinorsAllN(Lsldot'spin''i''j',?a)".format(e=spinorsSL2C_c["[e_C]"]) + "*sigma?{sigma,sigmabar,[sl2Ceps],sigma2,sigmabar2}(?b,Usldot'spin''i''j',?c) =" + " -{e:s}(Usldot'spin''i''j',?a)*sigma(?b,Lsldot'spin''i''j',?c);\n".format(e=spinorsSL2C_c["[e_C]"])
            form += 3 * "\t#enddo\n"

            form += "endrepeat;\n"
            form += "#endprocedure"
            return form

    def form_spinorCommutativetoDerivative(self, nDer):
        """
        Rewrite (if there exists a spin index in the term) spinor fields again in inside the derivative where it
        came from as:
            D(lor1?,[q(spin7786, gauge7124, colf7492, flav7371)])*Q(Lslspin778612,gauge7124,colf7492,flav7371) = D(lor1,Q(Lslspin778612,gauge7124,colf7492,flav7371)),
        where [q(spin7786, gauge7124, colf7492, flav7371)] is some placeholder to remember the position of the
        spinor field and Q(Lslspin778612,gauge7124,colf7492,flav7371) is an weyl spinor field which coming from the
        contracted Dirac spinors.
        Note: Also a "lonely" spinor field is replaced backwards,
            [q(spin7786, gauge7124, colf7492, flav7371)]*Q(Lslspin778612,gauge7124,colf7492,flav7371) = Q(Lslspin778612,gauge7124,colf7492,flav7371),
        in order to ensure also the order of the "lonely" spinor fields.

        Returns
        -------
        form : str
            Formatted string of the FORM procedure file "spinorCommutativetoDerivative.prc".

        """
        if len(self.sorted_indices["spin"]) != 0:
            form_list = []
            form_list_bar = []
            for operator in self.operators:
                spinor = ""
                if operator.name == opname["e"] or operator.name == opname["bar[e]"]:
                    spinor = "e"  # e(spin, flav) - RH leptons
                    sL2Cspinor = "[e_C+]"
                    sL2Cadjspinor = "[e_C]"
                elif operator.name == opname["u"] or operator.name == opname["bar[u]"]:
                    spinor = "u"  # u(spin, colf, flav) - RH up-type quarks
                    sL2Cspinor = "[u_C+]"
                    sL2Cadjspinor = "[u_C]"
                elif operator.name == opname["b"] or operator.name == opname["bar[b]"]:
                    spinor = "b"  # b(spin, colf, flav) - RH down-type quarks
                    sL2Cspinor = "[d_C+]"
                    sL2Cadjspinor = "[d_C]"
                elif operator.name == opname["l"] or operator.name == opname["bar[l]"]:
                    spinor = "l"  # l(spin, gauge, flav) - SU2 lepton dublett
                    sL2Cspinor = "L"
                    sL2Cadjspinor = "[L+]"
                elif operator.name == opname["q"] or operator.name == opname["bar[q]"]:
                    spinor = "q"  # q(spin, gauge, colf, flav) - SU2 quark dublett
                    sL2Cspinor = "Q"
                    sL2Cadjspinor = "[Q+]"

                if spinor:
                    indices = ",".join(operator.spinorIndices)
                    # Replacements of ordinary spinor field.
                    if spinor in ["u","b"]:
                        form_temp = "\t#do colfindex={" + ",".join(
                            ["", "a", "b", "c", "d"]) + "}\n"  # colf1234a...colf1234d possible.
                        indices_without_spin = operator.spinorIndices[1] + "'colfindex'" + "," + \
                                               ",".join(operator.spinorIndices[2:])
                        for i in range(nDer + 1):
                            lHS, rHS = Term.derivativesSLbackwards(f"[{opname[spinor]:s}({indices:s})]",
                                                                   f"{opnameSL2C[sL2Cspinor]:s}(Usldot1, {indices_without_spin:s})",
                                                                   n=i)
                            form_temp += 2 * "\t" + "id " + lHS
                            form_temp += f"*{spinorsSL2C_c[sL2Cspinor]:s}(Usldot1?Usldot, {indices_without_spin:s})"
                            form_temp += " = " + rHS + ";\n"
                        form_temp += "\t" + "#enddo"
                    elif spinor == "e":
                        form_temp = ""
                        for i in range(nDer + 1):
                            indices_without_spin = ",".join(operator.spinorIndices[1:])
                            lHS, rHS = Term.derivativesSLbackwards(f"[{opname[spinor]:s}({indices:s})]",
                                                                   f"{opnameSL2C[sL2Cspinor]:s}(Usldot1, {indices_without_spin:s})",
                                                                   n=i)
                            form_temp += "\t" + "id " + lHS
                            form_temp += f"*{spinorsSL2C_c[sL2Cspinor]:s}(Usldot1?Usldot, {indices_without_spin:s})"
                            form_temp += " = " + rHS + ";"
                    elif spinor == "l":
                        form_temp = "\t#do gaugeindex={" + ",".join(["","a", "b"]) + "}\n"  # gauge1234a and gauge1234b possible.
                        # The operator in [] with indices "indices" is the same since the beginning. Therefore,
                        # auxiliary indices "a" and "b" may only occur in "indices_without_spin".
                        indices_without_spin = operator.spinorIndices[1] + "'gaugeindex'" + "," + ",".join(
                            operator.spinorIndices[2:])
                        for i in range(nDer + 1):
                            lHS, rHS = Term.derivativesSLbackwards(f"[{opname[spinor]:s}({indices:s})]",
                                                                   f"{opnameSL2C[sL2Cspinor]:s}(Lsl1, {indices_without_spin:s})",
                                                                   n=i)
                            form_temp += 2*"\t" + "id " + lHS
                            form_temp += f"*{spinorsSL2C_c[sL2Cspinor]:s}(Lsl1?Lsl, {indices_without_spin:s})"
                            form_temp += " = " + rHS + ";\n"
                        form_temp += "\t" + "#enddo"
                    elif spinor == "q":
                        form_temp = "\t#do gaugeindex={" + ",".join(["","a", "b"]) + "}\n"  # gauge1234a and gauge1234b possible.
                        form_temp += "\t#do colfindex={" + ",".join(["", "a", "b", "c", "d"]) + "}\n"  # colf1234a...colf1234d possible.
                        indices_without_spin = operator.spinorIndices[1] + "'gaugeindex'" + ","
                        indices_without_spin += operator.spinorIndices[2] + "'colfindex'" + "," + ",".join(operator.spinorIndices[3:])
                        for i in range(nDer + 1):
                            lHS, rHS = Term.derivativesSLbackwards(f"[{opname[spinor]:s}({indices:s})]",
                                                                   f"{opnameSL2C[sL2Cspinor]:s}(Lsl1, {indices_without_spin:s})",
                                                                   n=i)
                            form_temp += 2*"\t" + "id " + lHS
                            form_temp += f"*{spinorsSL2C_c[sL2Cspinor]:s}(Lsl1?Lsl, {indices_without_spin:s})"
                            form_temp += " = " + rHS + ";\n"
                        form_temp += "\t" + "#enddo\n"
                        form_temp += "\t" + "#enddo"

                    form_list.append(form_temp)


                    # Replacements of Hermitian conjugated spinor fields.
                    if spinor in ["u","b"]:
                        form_temp = "\t#do colfindex1={" + ",".join(
                            ["", "a", "b", "c", "d"]) + "}\n"  # colf1234a...colf1234d possible.
                        form_temp += "\t#do colfindex2={" + ",".join(
                            ["", "a", "b", "c", "d"]) + "}\n"  # colf1234a...colf1234d possible.
                        indices_without_spin = operator.spinorIndices[1] + "'colfindex1'" + "," + \
                                               operator.spinorIndices[1] + "'colfindex2'" + "," + ",".join(operator.spinorIndices[2:])
                        for i in range(nDer + 1):
                            barspinor = "bar[" + spinor + "]"
                            lHS, rHS = Term.derivativesSLbackwards(f"[{opname[barspinor]:s}({indices:s})]",
                                                                   f"{opnameSL2C[sL2Cadjspinor]:s}(Lsl1, {indices_without_spin:s})",
                                                                   n=i)
                            form_temp += 2 * "\t" + "id " + lHS
                            form_temp += f"*{spinorsSL2C_c[sL2Cadjspinor]:s}(Lsl1?Lsl, {indices_without_spin:s})"
                            form_temp += " = " + rHS + ";\n"
                        form_temp += "\t" + "#enddo\n"
                        form_temp += "\t" + "#enddo\n"

                        form_temp += "\t#do colfindex={" + ",".join(["", "a", "b", "c", "d"]) + "}\n"  # colf1234a...colf1234d possible.
                        indices_without_spin = operator.spinorIndices[1] + "'colfindex'" + "," + \
                                               ",".join(operator.spinorIndices[2:])
                        for i in range(nDer + 1):
                            barspinor = "bar[" + spinor + "]"
                            lHS, rHS = Term.derivativesSLbackwards(f"[{opname[barspinor]:s}({indices:s})]",
                                                                   f"{opnameSL2C[sL2Cadjspinor]:s}(Lsl1, {indices_without_spin:s})",
                                                                   n=i)
                            form_temp += 2 * "\t" + "id " + lHS
                            form_temp += f"*{spinorsSL2C_c[sL2Cadjspinor]:s}(Lsl1?Lsl, {indices_without_spin:s})"
                            form_temp += " = " + rHS + ";\n"
                        form_temp += "\t" + "#enddo"
                    elif spinor == "e":
                        form_temp = ""
                        for i in range(nDer + 1):
                            indices_without_spin = ",".join(operator.spinorIndices[1:])
                            barspinor = "bar[" + spinor + "]"
                            lHS, rHS = Term.derivativesSLbackwards(f"[{opname[barspinor]:s}({indices:s})]",
                                                                   f"{opnameSL2C[sL2Cadjspinor]:s}(Lsl1, {indices_without_spin:s})",
                                                                   n=i)
                            form_temp += "\t" + "id " + lHS
                            form_temp += f"*{spinorsSL2C_c[sL2Cadjspinor]:s}(Lsl1?Lsl, {indices_without_spin:s})"
                            form_temp += " = " + rHS + ";"
                    elif spinor == "l":
                        form_temp = "\t#do gaugeindex={" + ",".join(["", "a", "b"]) + "}\n"  # gauge1234a and gauge1234b possible.
                        indices_without_spin = operator.spinorIndices[1] + "'gaugeindex'" + "," + ",".join(
                            operator.spinorIndices[2:])
                        for i in range(nDer + 1):
                            barspinor = "bar[" + spinor + "]"
                            lHS, rHS = Term.derivativesSLbackwards(f"[{opname[barspinor]:s}({indices:s})]",
                                                                   f"{opnameSL2C[sL2Cadjspinor]:s}(Usldot1, {indices_without_spin:s})",
                                                                   n=i)
                            form_temp += 2*"\t" + "id " + lHS
                            form_temp += f"*{spinorsSL2C_c[sL2Cadjspinor]:s}(Usldot1?Usldot, {indices_without_spin:s})"
                            form_temp += " = " + rHS + ";\n"
                        form_temp += "\t" + "#enddo"
                    elif spinor == "q":
                        form_temp = "\t#do gaugeindex={" + ",".join(["", "a", "b"]) + "}\n"  # gauge1234a and gauge1234b possible.
                        form_temp += "\t#do colfindex1={" + ",".join(["", "a", "b", "c", "d"]) + "}\n"  # colf1234a...colf1234d possible.
                        form_temp += "\t#do colfindex2={" + ",".join(["", "a", "b", "c", "d"]) + "}\n"  # colf1234a...colf1234d possible.
                        indices_without_spin = operator.spinorIndices[1] + "'gaugeindex'" + ","
                        indices_without_spin += operator.spinorIndices[2] + "'colfindex1'" + "," + operator.spinorIndices[2] + "'colfindex2'" + "," + ",".join(operator.spinorIndices[3:])
                        for i in range(nDer + 1):
                            barspinor = "bar[" + spinor + "]"
                            lHS, rHS = Term.derivativesSLbackwards(f"[{opname[barspinor]:s}({indices:s})]",
                                                                   f"{opnameSL2C[sL2Cadjspinor]:s}(Usldot1, {indices_without_spin:s})",
                                                                   n=i)
                            form_temp += 2*"\t" + "id " + lHS
                            form_temp += f"*{spinorsSL2C_c[sL2Cadjspinor]:s}(Usldot1?Usldot, {indices_without_spin:s})"
                            form_temp += " = " + rHS + ";\n"
                        form_temp += "\t" + "#enddo\n"
                        form_temp += "\t" + "#enddo\n"
                        form_temp += "\t" + "#enddo\n"

                        form_temp += "\t#do gaugeindex={" + ",".join(
                            ["", "a", "b"]) + "}\n"  # gauge1234a and gauge1234b possible.
                        form_temp += "\t#do colfindex1={" + ",".join(
                            ["", "a", "b", "c", "d"]) + "}\n"  # colf1234a...colf1234d possible.
                        indices_without_spin = operator.spinorIndices[1] + "'gaugeindex'" + ","
                        indices_without_spin += operator.spinorIndices[2] + "'colfindex1'" + "," + \
                                                ",".join(operator.spinorIndices[3:])
                        for i in range(nDer + 1):
                            barspinor = "bar[" + spinor + "]"
                            lHS, rHS = Term.derivativesSLbackwards(f"[{opname[barspinor]:s}({indices:s})]",
                                                                   f"{opnameSL2C[sL2Cadjspinor]:s}(Usldot1, {indices_without_spin:s})",
                                                                   n=i)
                            form_temp += 2 * "\t" + "id " + lHS
                            form_temp += f"*{spinorsSL2C_c[sL2Cadjspinor]:s}(Usldot1?Usldot, {indices_without_spin:s})"
                            form_temp += " = " + rHS + ";\n"
                        form_temp += "\t" + "#enddo\n"
                        form_temp += "\t" + "#enddo"


                    form_list_bar.append(form_temp)

            form = ""
            form += "#procedure spinorCommutativetoDerivative\n"
            form += "* Replacements of ordinary spinor field.\n"
            for i in form_list:
                form += i + "\n"
            form += "* Replacements of Hermitian conjugated spinor fields.\n"
            for i in form_list_bar:
                form += i + "\n"
            form += "#endprocedure"

            return form

    def form_substituteSU3Spinors(self):
        """
        Substitute hermitian conjugated fields wich have colour indices and transform in anti-fundamental representations
        by SU3 epsilon tensor and field with 2 color indices of the fundamental representation.
        #call substituteSU3Spinors
        Returns
        -------

        """
        form = ""
        if len(self.sorted_indices["colf"]) != 0:
            form += "#procedure substituteSU3Spinors\n"
            form += "#do colf={" + ",".join([colf for colf in self.sorted_indices[
                "colf"]]) + "}\n"  # '#do colf={colf5554, colf6674}'
            form += 1 * "\t" + f"id once {spinorsSL2C_c['[Q+]']}?" + "{" + f"{spinorsSL2C_c['[Q+]']},{spinorsSL2C_c['[u_C]']},{spinorsSL2C_c['[d_C]']}" + "}" + f"(?a, 'colf', ?b) = (1/2)*[su3eps]('colf', 'colf'a, 'colf'b)*{spinorsSL2C_c['[Q+]']}(?a, 'colf'a, 'colf'b,?b);\n"
            form += 1 * "\t" + f"id once {spinorsSL2C_c['[Q+]']}?" + "{" + f"{spinorsSL2C_c['[Q+]']},{spinorsSL2C_c['[u_C]']},{spinorsSL2C_c['[d_C]']}" + "}" + f"(?a, 'colf', ?b) = (1/2)*[su3eps]('colf', 'colf'c, 'colf'd)*{spinorsSL2C_c['[Q+]']}(?a, 'colf'c, 'colf'd,?b);\n"
            form += "#enddo\n"
            form += "#endprocedure"
            return form

    def form_simplifyEpsSU3(self):
        """
        Simplify SU3, i.e. 3 component epsilon tensor.
        #call simplifyEpsSU3
        Returns
        -------

        """
        form = ""
        if len(self.sorted_indices["colf"]) != 0:
            form += "#procedure simplifyEpsSU3\n"
            form += "repeat;"
            form += "* Replace epsilons by Kronecker deltas [sl2CdK](,):\n"
            # 3 Cyclic permutations of first eps and first cyclic permutation of second eps.
            form +="* 3 Cyclic permutations of first eps and first cyclic permutation of second eps.\n"
            form += "\t" + "id [su3eps](colf1?colf, colf2?colf, colf3?colf) * [su3eps](colf1?colf, colf4?colf, colf5?colf) = + ([su3dK](colf2,colf4)*[su3dK](colf3,colf5) - [su3dK](colf2,colf5)*[su3dK](colf3,colf4));\n"
            form += "\t" + "id [su3eps](colf3?colf, colf1?colf, colf2?colf) * [su3eps](colf1?colf, colf4?colf, colf5?colf) = + ([su3dK](colf2,colf4)*[su3dK](colf3,colf5) - [su3dK](colf2,colf5)*[su3dK](colf3,colf4));\n"
            form += "\t" + "id [su3eps](colf2?colf, colf3?colf, colf1?colf) * [su3eps](colf1?colf, colf4?colf, colf5?colf) = + ([su3dK](colf2,colf4)*[su3dK](colf3,colf5) - [su3dK](colf2,colf5)*[su3dK](colf3,colf4));\n"
            # 3 Antiyclic permutations of first eps and first cyclic permutation of second eps.
            form +="* 3 Antiyclic permutations of first eps and first cyclic permutation of second eps.\n"
            form += "\t" + "id [su3eps](colf1?colf, colf3?colf, colf2?colf) * [su3eps](colf1?colf, colf4?colf, colf5?colf) = - ([su3dK](colf2,colf4)*[su3dK](colf3,colf5) - [su3dK](colf2,colf5)*[su3dK](colf3,colf4));\n"
            form += "\t" + "id [su3eps](colf2?colf, colf1?colf, colf3?colf) * [su3eps](colf1?colf, colf4?colf, colf5?colf) = - ([su3dK](colf2,colf4)*[su3dK](colf3,colf5) - [su3dK](colf2,colf5)*[su3dK](colf3,colf4));\n"
            form += "\t" + "id [su3eps](colf3?colf, colf2?colf, colf1?colf) * [su3eps](colf1?colf, colf4?colf, colf5?colf) = - ([su3dK](colf2,colf4)*[su3dK](colf3,colf5) - [su3dK](colf2,colf5)*[su3dK](colf3,colf4));\n"
            # 3 Cyclic permutations of first eps and second cyclic permutation of second eps.
            form +="* 3 Cyclic permutations of first eps and second cyclic permutation of second eps.\n"
            form += "\t" + "id [su3eps](colf1?colf, colf2?colf, colf3?colf) * [su3eps](colf5?colf, colf1?colf, colf4?colf) = + ([su3dK](colf2,colf4)*[su3dK](colf3,colf5) - [su3dK](colf2,colf5)*[su3dK](colf3,colf4));\n"
            form += "\t" + "id [su3eps](colf3?colf, colf1?colf, colf2?colf) * [su3eps](colf5?colf, colf1?colf, colf4?colf) = + ([su3dK](colf2,colf4)*[su3dK](colf3,colf5) - [su3dK](colf2,colf5)*[su3dK](colf3,colf4));\n"
            form += "\t" + "id [su3eps](colf2?colf, colf3?colf, colf1?colf) * [su3eps](colf5?colf, colf1?colf, colf4?colf) = + ([su3dK](colf2,colf4)*[su3dK](colf3,colf5) - [su3dK](colf2,colf5)*[su3dK](colf3,colf4));\n"
            # 3 Antiyclic permutations of first eps and second cyclic permutation of second eps.
            form +="* 3 Antiyclic permutations of first eps and second cyclic permutation of second eps.\n"
            form += "\t" + "id [su3eps](colf1?colf, colf3?colf, colf2?colf) * [su3eps](colf5?colf, colf1?colf, colf4?colf) = - ([su3dK](colf2,colf4)*[su3dK](colf3,colf5) - [su3dK](colf2,colf5)*[su3dK](colf3,colf4));\n"
            form += "\t" + "id [su3eps](colf2?colf, colf1?colf, colf3?colf) * [su3eps](colf5?colf, colf1?colf, colf4?colf) = - ([su3dK](colf2,colf4)*[su3dK](colf3,colf5) - [su3dK](colf2,colf5)*[su3dK](colf3,colf4));\n"
            form += "\t" + "id [su3eps](colf3?colf, colf2?colf, colf1?colf) * [su3eps](colf5?colf, colf1?colf, colf4?colf) = - ([su3dK](colf2,colf4)*[su3dK](colf3,colf5) - [su3dK](colf2,colf5)*[su3dK](colf3,colf4));\n"
            # 3 Cyclic permutations of first eps and third cyclic permutation of second eps.
            form +="* 3 Cyclic permutations of first eps and third cyclic permutation of second eps.\n"
            form += "\t" + "id [su3eps](colf1?colf, colf2?colf, colf3?colf) * [su3eps](colf4?colf, colf5?colf, colf1?colf) = + ([su3dK](colf2,colf4)*[su3dK](colf3,colf5) - [su3dK](colf2,colf5)*[su3dK](colf3,colf4));\n"
            form += "\t" + "id [su3eps](colf3?colf, colf1?colf, colf2?colf) * [su3eps](colf4?colf, colf5?colf, colf1?colf) = + ([su3dK](colf2,colf4)*[su3dK](colf3,colf5) - [su3dK](colf2,colf5)*[su3dK](colf3,colf4));\n"
            form += "\t" + "id [su3eps](colf2?colf, colf3?colf, colf1?colf) * [su3eps](colf4?colf, colf5?colf, colf1?colf) = + ([su3dK](colf2,colf4)*[su3dK](colf3,colf5) - [su3dK](colf2,colf5)*[su3dK](colf3,colf4));\n"
            # 3 Antiyclic permutations of first eps and third cyclic permutation of second eps.
            form +="* 3 Antiyclic permutations of first eps and third cyclic permutation of second eps.\n"
            form += "\t" + "id [su3eps](colf1?colf, colf3?colf, colf2?colf) * [su3eps](colf4?colf, colf5?colf, colf1?colf) = - ([su3dK](colf2,colf4)*[su3dK](colf3,colf5) - [su3dK](colf2,colf5)*[su3dK](colf3,colf4));\n"
            form += "\t" + "id [su3eps](colf2?colf, colf1?colf, colf3?colf) * [su3eps](colf4?colf, colf5?colf, colf1?colf) = - ([su3dK](colf2,colf4)*[su3dK](colf3,colf5) - [su3dK](colf2,colf5)*[su3dK](colf3,colf4));\n"
            form += "\t" + "id [su3eps](colf3?colf, colf2?colf, colf1?colf) * [su3eps](colf4?colf, colf5?colf, colf1?colf) = - ([su3dK](colf2,colf4)*[su3dK](colf3,colf5) - [su3dK](colf2,colf5)*[su3dK](colf3,colf4));\n"
            #
            # 3 Cyclic permutations of first eps and first anticyclic permutation of second eps.
            form +="* 3 Cyclic permutations of first eps and first anticyclic permutation of second eps.\n"
            form += "\t" + "id [su3eps](colf1?colf, colf2?colf, colf3?colf) * [su3eps](colf1?colf, colf5?colf, colf4?colf) = - ([su3dK](colf2,colf4)*[su3dK](colf3,colf5) - [su3dK](colf2,colf5)*[su3dK](colf3,colf4));\n"
            form += "\t" + "id [su3eps](colf3?colf, colf1?colf, colf2?colf) * [su3eps](colf1?colf, colf5?colf, colf4?colf) = - ([su3dK](colf2,colf4)*[su3dK](colf3,colf5) - [su3dK](colf2,colf5)*[su3dK](colf3,colf4));\n"
            form += "\t" + "id [su3eps](colf2?colf, colf3?colf, colf1?colf) * [su3eps](colf1?colf, colf5?colf, colf4?colf) = - ([su3dK](colf2,colf4)*[su3dK](colf3,colf5) - [su3dK](colf2,colf5)*[su3dK](colf3,colf4));\n"
            # 3 Antiyclic permutations of first eps and first anticyclic permutation of second eps.
            form +="* 3 Antiyclic permutations of first eps and first anticyclic permutation of second eps.\n"
            form += "\t" + "id [su3eps](colf1?colf, colf3?colf, colf2?colf) * [su3eps](colf1?colf, colf5?colf, colf4?colf) = + ([su3dK](colf2,colf4)*[su3dK](colf3,colf5) - [su3dK](colf2,colf5)*[su3dK](colf3,colf4));\n"
            form += "\t" + "id [su3eps](colf2?colf, colf1?colf, colf3?colf) * [su3eps](colf1?colf, colf5?colf, colf4?colf) = + ([su3dK](colf2,colf4)*[su3dK](colf3,colf5) - [su3dK](colf2,colf5)*[su3dK](colf3,colf4));\n"
            form += "\t" + "id [su3eps](colf3?colf, colf2?colf, colf1?colf) * [su3eps](colf1?colf, colf5?colf, colf4?colf) = + ([su3dK](colf2,colf4)*[su3dK](colf3,colf5) - [su3dK](colf2,colf5)*[su3dK](colf3,colf4));\n"
            # 3 Cyclic permutations of first eps and second anticyclic permutation of second eps.
            form +="* 3 Cyclic permutations of first eps and second anticyclic permutation of second eps.\n"
            form += "\t" + "id [su3eps](colf1?colf, colf2?colf, colf3?colf) * [su3eps](colf4?colf, colf1?colf, colf5?colf) = - ([su3dK](colf2,colf4)*[su3dK](colf3,colf5) - [su3dK](colf2,colf5)*[su3dK](colf3,colf4));\n"
            form += "\t" + "id [su3eps](colf3?colf, colf1?colf, colf2?colf) * [su3eps](colf4?colf, colf1?colf, colf5?colf) = - ([su3dK](colf2,colf4)*[su3dK](colf3,colf5) - [su3dK](colf2,colf5)*[su3dK](colf3,colf4));\n"
            form += "\t" + "id [su3eps](colf2?colf, colf3?colf, colf1?colf) * [su3eps](colf4?colf, colf1?colf, colf5?colf) = - ([su3dK](colf2,colf4)*[su3dK](colf3,colf5) - [su3dK](colf2,colf5)*[su3dK](colf3,colf4));\n"
            # 3 Antiyclic permutations of first eps and second anticyclic permutation of second eps.
            form +="* 3 Antiyclic permutations of first eps and second anticyclic permutation of second eps.\n"
            form += "\t" + "id [su3eps](colf1?colf, colf3?colf, colf2?colf) * [su3eps](colf4?colf, colf1?colf, colf5?colf) = + ([su3dK](colf2,colf4)*[su3dK](colf3,colf5) - [su3dK](colf2,colf5)*[su3dK](colf3,colf4));\n"
            form += "\t" + "id [su3eps](colf2?colf, colf1?colf, colf3?colf) * [su3eps](colf4?colf, colf1?colf, colf5?colf) = + ([su3dK](colf2,colf4)*[su3dK](colf3,colf5) - [su3dK](colf2,colf5)*[su3dK](colf3,colf4));\n"
            form += "\t" + "id [su3eps](colf3?colf, colf2?colf, colf1?colf) * [su3eps](colf4?colf, colf1?colf, colf5?colf) = + ([su3dK](colf2,colf4)*[su3dK](colf3,colf5) - [su3dK](colf2,colf5)*[su3dK](colf3,colf4));\n"
            # 3 Cyclic permutations of first eps and third anticyclic permutation of second eps.
            form +="* 3 Cyclic permutations of first eps and third anticyclic permutation of second eps.\n"
            form += "\t" + "id [su3eps](colf1?colf, colf2?colf, colf3?colf) * [su3eps](colf5?colf, colf4?colf, colf1?colf) = - ([su3dK](colf2,colf4)*[su3dK](colf3,colf5) - [su3dK](colf2,colf5)*[su3dK](colf3,colf4));\n"
            form += "\t" + "id [su3eps](colf3?colf, colf1?colf, colf2?colf) * [su3eps](colf5?colf, colf4?colf, colf1?colf) = - ([su3dK](colf2,colf4)*[su3dK](colf3,colf5) - [su3dK](colf2,colf5)*[su3dK](colf3,colf4));\n"
            form += "\t" + "id [su3eps](colf2?colf, colf3?colf, colf1?colf) * [su3eps](colf5?colf, colf4?colf, colf1?colf) = - ([su3dK](colf2,colf4)*[su3dK](colf3,colf5) - [su3dK](colf2,colf5)*[su3dK](colf3,colf4));\n"
            # 3 Antiyclic permutations of first eps and third anticyclic permutation of second eps.
            form +="* 3 Antiyclic permutations of first eps and third anticyclic permutation of second eps.\n"
            form += "\t" + "id [su3eps](colf1?colf, colf3?colf, colf2?colf) * [su3eps](colf5?colf, colf4?colf, colf1?colf) = + ([su3dK](colf2,colf4)*[su3dK](colf3,colf5) - [su3dK](colf2,colf5)*[su3dK](colf3,colf4));\n"
            form += "\t" + "id [su3eps](colf2?colf, colf1?colf, colf3?colf) * [su3eps](colf5?colf, colf4?colf, colf1?colf) = + ([su3dK](colf2,colf4)*[su3dK](colf3,colf5) - [su3dK](colf2,colf5)*[su3dK](colf3,colf4));\n"
            form += "\t" + "id [su3eps](colf3?colf, colf2?colf, colf1?colf) * [su3eps](colf5?colf, colf4?colf, colf1?colf) = + ([su3dK](colf2,colf4)*[su3dK](colf3,colf5) - [su3dK](colf2,colf5)*[su3dK](colf3,colf4));\n"

            form += "* Replace only contracted Kronecker-deltas which are contracted with eps and Kronecker-deltas themself, because contractions inside one building block are not wanted:\n"
            form += "\t" + "id [su3eps]?{[su3eps],[su3dK]}(?a,colf1?colf,?b)*[su3dK](?c,colf1?colf,?d) = [su3eps](?a,?c,?d,?b);\n"  # [su3eps]?
            form += "* Replace self-contracted Kronecker-deltas by the dimension (=3):\n"
            form += "\t" + "id [su3dK](colf1?colf,colf1?colf) = d_(colf1,colf1);\n"
            form += "endrepeat;\n"
            form += "* Bring indices of epsilons in order:\n"
            form += "Multiply replace_([su3eps],[su3epsA]);\n"
            form += ".sort\n"
            form += "Multiply replace_([su3epsA],[su3eps]);\n"
            # form += 1 * "\t" + "id [su3eps](colf1?colf, colf2?colf, colf3?colf) = e_(colf1, colf2, colf3);\n"
            # form += "contract 0;\n"
            # form += 1 * "\t" + "id e_(colf1?colf, colf2?colf, colf3?colf) = [su3eps](colf1, colf2, colf3);\n"
            # form += 1 * "\t" + "id d_(colf1?colf, colf2?colf) = [su3dK](colf1, colf2);\n"
            form += "#endprocedure"
            return form

    #
    # Flavor
    #
    def form_simplifyFlavorMatrices(self):
        """
        It can be noticed that the flavormatrices often occur in terms like:
            [yd+](spin1, spin2)*yd(spin2, spin3).
        For no further reason than readability terms like this are substituted by
            [[yd+]*yd](spin1, spin3).
        Returns
        -------

        """
        form = ""
        initialize = ""
        if len(self.sorted_indices["flav"]) != 0:
            initialize += "CFunction "
            form += "*\n* Simplify flavormatrices\n*\n"
            for flM in ["yu", "yd", "ye"]:
                form += "id {yDagger:s}(flav1?flav,flav2?flav)*{y:s}(flav2?flav,flav3?flav) = [{yDagger:s}*{y:s}](flav1,flav3) ;\n".format(
                    y=opname[flM], yDagger=opname["conj[" + flM + "]"])
                form += "id {y:s}(flav1?flav,flav2?flav)*{yDagger:s}(flav2?flav,flav3?flav) = [{y:s}*{yDagger:s}](flav1,flav3) ;\n".format(
                    y=opname[flM], yDagger=opname["conj[" + flM + "]"])
                initialize += "[{yDagger:s}*{y:s}], ".format(y=opname[flM], yDagger=opname["conj[" + flM + "]"])
                initialize += "[{y:s}*{yDagger:s}], ".format(y=opname[flM], yDagger=opname["conj[" + flM + "]"])
            initialize = initialize[:-2] + ";\n"
            return form, initialize
    #
    # Sort fields
    #
    def form_derivativeasIndex(self, nDer=4):
        form = "#procedure derivativeasIndex\n"
        form += "repeat;\n"
        for i in range(1, nDer + 1):
            lHS, rHS_Indices = Term.derivativestoInd(f"{opnameSL2C['[e_C]']}?(?x)", n=i)
            form += f"\tid {lHS} = {opnameSL2C['[e_C]']}({rHS_Indices}?x);\n"
        form += "endrepeat;\n"
        form += "#endprocedure"
        return form

    def form_indexasDerivative(nDer=4):
        """
        Derivative are the only object that have always the index structure [Lsl, Usldot].
        """
        form = "#procedure indexasDerivative\n"
        form += "repeat;\n"
        for i in reversed(range(1, nDer + 1)):
            lHS_Indices, rHS, brackets = Term.indtoDerivatives(n=i)
            form += f"\tid {opnameSL2C['[e_C]']}?!" + "{D}" + f"({lHS_Indices}, ?x) = {rHS}{opnameSL2C['[e_C]']}(?x){brackets};\n"
        form += "endrepeat;\n"
        form += "#endprocedure"
        return form

    #     repeat;
    # 	id [e_C]?!{D}(Lsl1?Lsl, Usldot2?Usldot, Lsl3?Lsl, Usldot4?Usldot, Lsl5?Lsl, Usldot6?Usldot, Lsl7?Lsl, Usldot8?Usldot, gauge1?) = D(Lsl1, Usldot2, D(Lsl3, Usldot4, D(Lsl5, Usldot6, D(Lsl7, Usldot8, [e_C](gauge1)))));
    # 	id [e_C]?!{D}(Lsl1?Lsl, Usldot2?Usldot, Lsl3?Lsl, Usldot4?Usldot, Lsl5?Lsl, Usldot6?Usldot, ?x) = D(Lsl1, Usldot2, D(Lsl3, Usldot4, D(Lsl5, Usldot6, [e_C](?x))));
    # 	id [e_C]?!{D}(Lsl1?Lsl, Usldot2?Usldot, Lsl3?Lsl, Usldot4?Usldot, ?x) = D(Lsl1, Usldot2, D(Lsl3, Usldot4, [e_C](?x)));
    # 	id [e_C]?!{D}(Lsl1?Lsl, Usldot2?Usldot, ?x) = D(Lsl1, Usldot2, [e_C](?x));
    # endrepeat;

    @staticmethod
    def extractOrder():
        """
        Extract order from the helicity of the fields specified in config/fields.yml.
        Returns
        -------
        Tuple of str, i.e. names of the fields in the demanded order.
        """
        form = {'GL' : f"{opname['G']:s}L", 'GR' : f"{opname['G']:s}R",
               'WL' : f"{opname['V']:s}L", 'WR' : f"{opname['V']:s}R",
               'BL' : f"{opname['F']:s}L", 'BR' : f"{opname['F']:s}R",
               'L' : [opnameSL2C['L'], opnameSL2C['[L+]']],
               'eC' : [opnameSL2C['[e_C]'], opnameSL2C['[e_C+]']],
               'Q' : [opnameSL2C['Q'], opnameSL2C['[Q+]']],
               'uC' : [opnameSL2C['[u_C]'], opnameSL2C['[u_C+]']],
               'dC' : [opnameSL2C['[d_C]'], opnameSL2C['[d_C+]']],
               'H' : [opname['H'], opname['conj[H]']],
               'D' : opname['cov']}

        fields = []
        for field in field_config.keys():
            if field != "D":
                if "ac" in field_config[field].keys():
                    if field_config[field]["ac"]:
                        fields.append(OperatorModel(field, form[field][0], field_config[field]["helicity"][0], True))
                        fields.append(OperatorModel(field, form[field][1], field_config[field]["helicity"][1], True))
                else:
                    if field == "H":
                        fields.append(OperatorModel(field, form[field][0], field_config[field]["helicity"][0], False))
                        fields.append(OperatorModel(field, form[field][1], field_config[field]["helicity"][0], False))
                    else:
                        fields.append(OperatorModel(field, form[field], field_config[field]["helicity"][0], False))

        # Sort listed field, first by their helicity and then by the alphabetical order of the names. Use bubble sort:
        n = len(fields)
        # Traverse through all array elements
        for i in range(n - 1):
            # range(n) also work but outer loop will repeat one time more than needed.
            # Last i elements are already in place
            for j in range(0, (n - 1) - i):
                # traverse the array from 0 to n-i-1
                # Swap if the element found is greater
                # than the next element
                if fields[j].helicity > fields[j + 1].helicity:
                    fields[j], fields[j + 1] = fields[j + 1], fields[j]
                elif fields[j].helicity == fields[j + 1].helicity:
                    # Note that due to the ASCII standard the letter 'A' stands before 'a' and so on
                    if fields[j].name > fields[j + 1].name:
                        fields[j], fields[j + 1] = fields[j + 1], fields[j]
        logger.debug(f"Fields are ordered by helicity like: {' '.join([f.form_name for f in fields])}")
        return fields

    @staticmethod
    def form_sortfields(order):
        """

        Parameters
        ----------
        order : List[OperatorModel]

        Returns
        -------

        """

        form = ""
        form += "#procedure sortfields\n"
        form += "repeat;\n"
        ordered_fields = order[::-1]  # Reversed list
        for i in range(len(ordered_fields)):
            for j in range(len(ordered_fields)):
                if j > i:
                    if ordered_fields[i].fermion and ordered_fields[j].fermion:
                        sign = "-"
                    else:
                        sign = "+"
                    form += f"\tid {ordered_fields[i].form_name}(?a)*{ordered_fields[j].form_name}(?b) = {sign}{ordered_fields[j].form_name}(?b)*{ordered_fields[i].form_name}(?a);\n"
            form += "\n"
        form += "endrepeat;\n"
        form += "#endprocedure"
        return form


