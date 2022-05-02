import re
import logging
import sys
from yaml import safe_load
from typing import Dict, List, Tuple

from tioc import opname, CONFIG_PATH
# from tioc.class_index import Index

logger_autoeft = logging.getLogger("autoeft.projection")
logger = logger_autoeft.getChild("operator")

class OperatorModel:
    """
    Class necessary for ordering of the fields.
    """
    def __init__(self, name, form_name, helicity, fermion):
        self.name = name
        self.form_name = form_name
        self.helicity = helicity
        self.fermion = fermion


class Operator():
    """
    A class describing the operators in a term.
    Possible Operators:
        {H, e, u, b, l, q} & {conj[H], bar[e], bar[u], bar[b], bar[l], bar[q]}
        field strength tensors of U(1),SU(2) and SU(3) named F, V and G
        covariant derivative cov[{lor},Field] of any field called Field
        completely antisymmetric symbol su2eps[{gauge1,gauge2}]
        group generators TT[{gaugeadj1},{gauge1, gauge2}]
        Gammamatrices, e.g.: gamma[{lor1}, {spin6497, spin6498}]
        couplings that carry indices (only Yukawa couplings): {yu,yd,ye,conj[yu],conj[yd],conj[ye]}


    Attributes
    ----------
    name : str
        Name of the Operator.
    indices : [str, str, ...]
        All type of indices in an operator. Occurring always at most one time in this least even if there would be
        two contracted indices. Note: All indices in a term are always fully contracted.
        Example for indices: [lor1, lor34532, spin34,...].
    fermion : bool
        Determines whether an operator is an fermion, i.e. a anti-commuting object.
    expression : str
        Explicit expression of the operator. If manual is true. The expression has to contain an valid FORM expression,
        containing only functions which are defined by the declaration functions.
    manual : str
        default: False.
    Magicmethods
    ------------
    __str__()
        Specify the format for printing an Operator object with the print command or converting a Operator object in
        string with str() function.
    __format__()
        Specify the format for "format" function in print statement: Here the same as the print statement itself.
    Methods
    -------
    rewriteCov()
        Rewrites all covariant derivatives.
    identifyHiggs()
        Identifies all occurring Higgs and rewrites them according to the "opname" dictionary.
    identifyFlavorMatrices()
        Identifies all occurring flavor matrices and rewrites them according to the "opname" dictionary.
    testIndex()
        Assert number of digits in each Index <= 5.
    identifySU2eps()
        Identifies all occurring SU2 epsilon tensor and rewrites them according to the "opname" dictionary.
    identifyGroupGenerators()
        Identifies all occurring group generators and rewrites them according to the "opname" dictionary.
    identifyGammaMatrices()
        Identifies all occurring gamma matrices and rewrites them according to the "opname" dictionary.
    identifyFieldstrengthtensors()
        Identifies all occurring fieldstrengthtensors and rewrites them according to the "opname" dictionary.
    identifySpinors()
        Identifies all occurring spinors and rewrites them according to the "opname" dictionary.

    """
    expression : str
    id : str  # op_numID
    numID: int  # number of the id
    indices : List[str]
    fermion : bool
    name : str
    description : str

    def __init__(self, expression, id:int, manual=False, **kwargs):
        self.expression = expression
        # self.numID = id
        self.id = f"op{id}"
        self.indices = []  # list with str of used indices in operator which are always fully contracted, i.e. if indices
        # are contracted inside of one operators their corresponding index should occur twice.
        # e.g. [lor1, lor34532, spin34,...]
        self.fermion = None
        self.name = None
        self.description = None
        if not manual:
            self.rewriteCov(derivativeasmultiplication=False)
            self.identifyHiggs()
            self.identifyFlavorMatrices()
            self.testIndex()
            self.identifySU2eps()
            self.identifyGroupGenerators()
            self.identifyGammaMatrices()
            self.identifyFieldstrengthtensors()
            self.identifySpinors()
        else:
            # manual case:
            try:
                self.indices = kwargs["indices"]
            except KeyError:
                pass
            try:
                self.fermion = kwargs["fermion"]
            except KeyError:
                pass
            try:
                self.name = kwargs["name"]
            except KeyError:
                pass
            try:
                self.description = kwargs["description"]
            except KeyError:
                pass

    @property
    def numID(self):
        return self._numID

    @numID.setter
    def numID(self, value):
        try:
            self._numID = int(value)
        except ValueError:
            logger.error(f"The numerical ID has to be an integer and not of the type {type(value)}.")
            sys.exit("STOP")

    @property
    def id(self):
        return f"op{self.numID:d}"
    @id.setter
    def id(self, value):
        match = re.match(r"op(?P<numID>\d+)", value)
        if match:
            self._numID = int(match.group("numID"))
        else:
            logger.error(f"The id has to be in the format 'op1234'.")
            sys.exit("STOP")


    def __str__(self):
        """Specify the format for printing an Operator object with the print command or converting a Operator object in
        string with str() function. """
        return "{0:s}".format(self.expression)
    def __repr__(self):
        return self.__str__()
    def __format__(self, key):
        """Specify the format for "format" function in print statement: Here the same as the print statement itself."""
        return self.__str__()

    def extract_indices(self, regex, *indices):
        """
        Return list of indices contained in indices or self.indices and matching the given pattern regex.
        Parameters
        ----------
        regex : str
            Regular expression which matches the index.
        indices : list(optional)
            If given, this list of indices is used instead of self.indices
        Returns
        -------
            List of extracted indices.
        """
        indexlist = []
        if indices:
            ind = list(indices)
        else:
            ind = self.indices
        for v in ind:
            match = re.match(regex, v)
            if match:
                indexlist.append(v)
        return indexlist

    def rewriteCov(self, derivativeasmultiplication=True):
        """

        Parameters
        ----------
        derivativeasmultiplication

        Returns
        -------

        """
        """
        Convert cov[{lor1234},cov[{lor5678}], ...]] into FORM compatible form and extract Lorentzindices of
        covariant derivatives.
        Parameters
        ----------
        derivativeasmultiplication : bool
            If true (default): Convert cov[{lor1234},cov[{lor5678}], ...]] into cov(lor1234)*cov(5678) ....
            If False: Convert cov[{lor1234},cov[{lor5678}], ...]] into cov(lor1234,cov(lor5678, ...)).
        Returns -------
            None.
        """
        indices_temp = []
        expression_temp = ""
        pattern = r"cov\[\{(lor\d{1,5})\},"
        match = re.finditer(pattern, self.expression)
        matches = list(match)
        if any(matches):
            for v in matches:
                index = v.group(1)
                indices_temp.append(index)
                if derivativeasmultiplication:
                    expression_temp += "{0:s}({1:s})*".format(opname["cov"], index)
                else:
                    expression_temp += "{0:s}({1:s},".format(opname["cov"], index)
            op = self.expression[matches[-1].end():-1*len(matches)]
            expression_temp += op
            if derivativeasmultiplication:
                self.expression = expression_temp
            else:
                self.expression = expression_temp + ")"*len(matches)
            self.indices += indices_temp

    def identifyHiggs(self):
        """
        Search for Higgs (Higgsboson) H and rewrite:
            H[{gauge1234}] -> H(gauge1234).
        Search for hermitian conjugated Higgs and rewrite:
            conj[H][{gauge1234}] -> [H+](gauge1234).
        Returns
        -------

        """
        # Search for normal Higgs first.
        pattern = r"H\[\{(gauge\d{1,5})\}\](?P<brackets>\)*)$"
        match = re.search(pattern, self.expression)
        if match:
            self.fermion = False
            index = match.group(1)
            self.indices.append(index)
            self.expression = self.expression[:match.start()] + "{0:s}({1:s})".format(opname["H"], index) + match.group("brackets")
            self.description = "Higgsfield"
            self.name = opname["H"]
        else:
            # Search for hermitian conjugated Higgs.
            pattern = r"conj\[H\]\[\{(gauge\d{1,5})\}\](?P<brackets>\)*)$"
            match = re.search(pattern, self.expression)
            if match:
                self.fermion = False
                index = match.group(1)
                self.indices.append(index)
                self.expression = self.expression[:match.start()] + "{0:s}({1:s})".format(opname["conj[H]"], index) + match.group("brackets")
                self.description = "conjugated Higgsfield"
                self.name = opname["conj[H]"]

    def testIndex(self):
        """
        Assert number of digits in each Index <= 5.
        Returns
        -------
            None.
        """
        pattern = r"\d{6,}"
        match = re.search(pattern, self.expression)
        if match:
            raise AssertionError(
                "Only indices with at most 5 digits are allow, but it was found an index with 6 or more digits.")

    def identifyFlavorMatrices(self):
        """
        Search for y[ude] (FlavorMatrices, i.e. Yukawa matrices) {yu, yd, ye} and rewrite:
            yu[{flav1234,flav5678}] -> yu(flav1234,flav5678).
        Search for complex conjugate Yukawa matrix {conj[yu], conj[yd], conj[ye]} and rewrite:
            conj[yu][{flav1234,flav5678}] -> [yu+](flav5678, flav1234).
            Note that they are rewritten as hermitian conjugates. Thus the indices are transposed.
        Returns
        -------

        """
        pattern = r"(y[ude])\[\{(flav\d{1,5}),(flav\d{1,5})\}\]"
        match = re.match(pattern, self.expression)  # Use "match", because the Yukawa coupling matrix should be all
        # in the expression of the operator.
        # Search for normal Yukawa matrix:
        if match:
            self.fermion = False
            index1 = match.group(2)
            index2 = match.group(3)
            self.indices.append(index1)
            self.indices.append(index2)
            self.expression = "{0:s}({1:s},{2:s})".format(opname[match.group(1)], index1, index2)
            self.description = "Flavor Matrix"
            self.name = opname[match.group(1)]
        else:
            # Search for hermitian conjugated Yukawa matrix:
            pattern2 = r"(conj\[y[ude]\])\[\{(flav\d{1,5}),(flav\d{1,5})\}\]"
            match2 = re.match(pattern2, self.expression)
            if match2:
                self.fermion = False
                index1 = match2.group(2)
                index2 = match2.group(3)
                self.indices.append(index1)
                self.indices.append(index2)
                self.expression = "{0:s}({2:s},{1:s})".format(opname[match2.group(1)], index1, index2)
                self.description = "conjugated Flavor Matrix"
                self.name = opname[match2.group(1)]

    def identifySU2eps(self):
        """
        Search for su2eps (Epsilontensor of SU2) and rewrite:
            su2eps[{gauge6153,gauge5345}] -> [su2eps](gauge6153,gauge5345)
        Search also for su3eps (Epsilontensor of SU3) and rewrite:
            su3eps[{colf1234,colf5678,colf2468}] -> [su3eps](colf1234,colf5678,colf2468)
        Returns
        -------

        """
        patternSU2 = r"su2eps\[\{(?P<gauge1>gauge\d{1,5}),(?P<gauge2>gauge\d{1,5})\}\]"
        patternSU3 = r"su3eps\[\{(?P<colf1>colf\d{1,5}),(?P<colf2>colf\d{1,5}),(?P<colf3>colf\d{1,5})\}\]"
        matchSU2 = re.match(patternSU2, self.expression)  # Use "matchSU2", because the su2eps should be all in the expression
        # of the operator.
        if matchSU2:
            self.fermion = False
            index1 = matchSU2.group("gauge1")
            index2 = matchSU2.group("gauge2")
            self.indices.append(index1)
            self.indices.append(index2)
            self.description = "SU2 Epsilon"
            self.expression = f"{opname['su2eps']:s}({index1:s},{index2:s})"
            self.name = opname["su2eps"]
        else:
            matchSU3 = re.match(patternSU3, self.expression)
            if matchSU3:
                self.fermion = False
                index1 = matchSU3.group("colf1")
                index2 = matchSU3.group("colf2")
                index3 = matchSU3.group("colf3")
                self.indices.append(index1)
                self.indices.append(index2)
                self.indices.append(index3)
                self.description = "SU3 Epsilon"
                self.expression = f"{opname['su3eps']:s}({index1:s},{index2:s},{index3:s})"
                self.name = opname["su3eps"]


    def identifyGroupGenerators(self):
        """
        Search for TT (GroupGenerators) and rewrite for SU2:
            TT[{gaugeadj8013},{gauge8019,gauge7973}] -> T(gaugeadj8013,gauge8019,gauge7973),
        for SU3:
            TT[{cola1234},{colf2345,colf8765}] -> T(cola1234,colf2345,colf8765),
        Return
        -------

        """
        patternSU2 = r"TT\[\{(?P<gaugeadj>gaugeadj\d{1,5})\},\{(?P<gauge1>gauge\d{1,5}),(?P<gauge2>gauge\d{1,5})\}\]"
        matchSU2 = re.match(patternSU2, self.expression)  # Use "match", because the TT should be all in the expression
        patternSU3 = r"TT\[\{(?P<adjcolorindex>cola\d{1,5})\},\{(?P<colf1>colf\d{1,5}),(?P<colf2>colf\d{1,5})\}\]"
        matchSU3 = re.match(patternSU3, self.expression)
        # of the operator.
        if matchSU2:
            self.fermion = False
            index1 = matchSU2.group("gaugeadj")
            index2 = matchSU2.group("gauge1")
            index3 = matchSU2.group("gauge2")
            self.indices.append(index1)
            self.indices.append(index2)
            self.indices.append(index3)
            self.description = "Group Generator of SU2."
            self.expression = "{0:s}({1:s},{2:s},{3:s})".format(opname["TT"], index1, index2, index3)
            self.name = opname["TT"]
        if matchSU3:
            self.fermion = False
            index1 = matchSU3.group("adjcolorindex")
            index2 = matchSU3.group("colf1")
            index3 = matchSU3.group("colf2")
            self.indices.append(index1)
            self.indices.append(index2)
            self.indices.append(index3)
            self.description = "Group Generator of SU3."
            self.expression = "{0:s}({1:s},{2:s},{3:s})".format(opname["TT"], index1, index2, index3)
            self.name = opname["TT"]

    def identifyGammaMatrices(self):
        """
        Search for gamma matrices and rewrite:
            gamma[{lor1}, {spin6497, spin6498}] -> gamma(lor1, spin6497,spin6498).
        Returns
        -------

        """
        pattern = r"gamma\[\{(lor\d{1,5})\},\{(spin\d{1,5}),(spin\d{1,5})\}\]"
        match = re.match(pattern, self.expression)  # Use "match", because the gamma matrix should be all in the
        # expression of the operator.
        if match:
            self.fermion = False
            index1 = match.group(1)
            index2 = match.group(2)
            index3 = match.group(3)
            self.indices.append(index1)
            self.indices.append(index2)
            self.indices.append(index3)
            self.description = "Gamma Matrix"
            self.expression = "{0:s}({1:s},{2:s},{3:s})".format(opname["gamma"], index1, index2, index3)
            self.name = opname["gamma"]

    def identifyFieldstrengthtensors(self):
        """
        Search for Fieldstrengthtensors {V, F, G} and rewrite:
        For SU(2):
            V[{gaugeadj8012},{lor1,lor2}] -> V(gaugeadj8012,lor1,lor2).
        For U(1):
            F[{lor1,lor2}] -> F(lor1,lor2).
        For SU(3):
            G[{cola?},{lor1,lor2}] -> G(cola?,lor1,lor2) but doesn't seem to occur.
        Returns
        -------

        """
        pattern = r"(?P<name>[VFG])\[(?P<gaugeadjindex>\{gaugeadj\d{1,5}\},)?(?P<colaindex>\{cola\d{1,5}\},)?\{(?P<lor1>lor\d{1,5}),(?P<lor2>lor\d{1,5})\}\](?P<brackets>\)*)$"  #
        match = re.search(pattern, self.expression)
        if match:
            self.fermion = False
            name = opname[match.group("name")]
            indexadj = match.group("gaugeadjindex")  # indexadj = None, for field F[{lor1,lor2}] and G[{cola?},{lor1,lor2}]
            indexcola = match.group("colaindex") # indexcola = None, for field F[{lor1,lor2}] and V[{gaugeadj8012},{lor1,lor2}]
            indices = []
            if indexadj != None and indexcola == None:
                # Additional index for V
                indexadj = indexadj[1:-2]
                indices.append(indexadj)
            elif indexcola != None and indexadj == None:
                # Additional index for G
                indexcola = indexcola[1:-2]
                indices.append(indexcola)
            elif indexcola != None and indexadj != None:
                raise AssertionError("A fieldstrengthtensor with both an color and an SU2 adjoint index is matched which shouldnt exist.")
            indexlor1 = match.group("lor1")
            indexlor2 = match.group("lor2")
            indices.append(indexlor1)
            indices.append(indexlor2)
            self.expression = self.expression[:match.start()] + f"{name:s}({','.join(indices):s})" + match.group("brackets")
            self.indices = self.indices + indices
            group = {opname["F"]: "U(1)", opname["V"]: "SU(2)", opname["G"]: "SU(3)"}
            self.description = f"Fieldstrengthtensor of {group[name]:s}"
            self.name = name

    def identifySpinors(self):
        """
        Search for spinors {e, u, b, l, q} & {bar[e], bar[u], bar[b], bar[l], bar[q]} and rewrite:

        right-handed leptons ("gauge" and "colf" is None):
          e[{spin5610},{flav5608}] -> e(spin5610,flav5608).
        right-handed up-type quarks ("gauge" is None):
          u[{spin3806},{colf3079},{flav2775}] -> u(spin3806,colf3079,flav2775).
        right-handed down-type quarks ("gauge" is None):
          b[{spin3787},{colf3065},{flav2791}] -> b(spin3787,colf3065,flav2791).
        SU2 Lepton Dublett ("colf" is None):
          l[{spin5952},{gauge6153},{flav5951}] -> l(spin5952,gauge6153,flav5951).
        SU2 Quark Dublett (Nothing is None):
          bar[q][{spin2763}, {gauge2764} , {colf2891} ,{flav2762}] -> [q+](spin2763,gauge2764,colf2891,flav2762).

        Returns
        -------

        """

        pattern = r"(?P<name>[eublq])\[\{(?P<spin>spin\d{1,5})\},(?P<gauge>\{gauge\d{1,5}\},)?(?P<colf>\{colf\d{1,5}\},)?\{(?P<flav>flav\d{1,5})\}\](?P<brackets>\)*)$"
        match = re.search(pattern, self.expression)
        if match:
            self.fermion = True
            name = opname[match.group("name")]
            indexSpin = match.group("spin")
            indexGauge = match.group("gauge")  # None
            indexColf = match.group("colf")  # maybe None
            indexFlav = match.group("flav")
            self.spinorIndices=[]
            self.indices.append(indexSpin)
            self.spinorIndices.append(indexSpin)
            if indexGauge:
                indexGauge = indexGauge[1:-2]
                self.indices.append(indexGauge)
                self.spinorIndices.append(indexGauge)
            if indexColf:
                indexColf = indexColf[1:-2]
                self.indices.append(indexColf)
                self.spinorIndices.append(indexColf)
            self.indices.append(indexFlav)
            self.spinorIndices.append(indexFlav)
            self.spinorIndices = tuple(self.spinorIndices)
            spinorName = "unidentified"
            if name == opname["e"]:
                spinorName = "of right-handed lepton"
            elif name == opname["u"]:
                spinorName = "of right-handed up-type quark"
            elif name == opname["b"]:
                spinorName = "of right-handed down-type quark"
            elif name == opname["l"]:
                spinorName = "of SU2 Lepton Dublett"
            elif name == opname["q"]:
                spinorName = "of SU2 Quark Dublett"
            self.description = "Spinor " + spinorName
            self.expression = self.expression[:match.start()] + "{0:s}({1:s}{2:s},{3:s}{4:s})".format(name,
                                                                                                      indexSpin,
                                                                                                      "" if indexGauge == None else "," + indexGauge,
                                                                                                      "" if indexColf == None else indexColf + ",",
                                                                                                      indexFlav)
            self.expression = self.expression + match.group("brackets")
            self.name = name
        else:
            pattern2 = r"(?P<name>bar\[[eublq]\])\[\{(?P<spin>spin\d{1,5})\},(?P<gauge>\{gauge\d{1,5}\},)?(?P<colf>\{colf\d{1,5}\},)?\{(?P<flav>flav\d{1,5})\}\](?P<brackets>\)*)$"
            match2 = re.search(pattern2, self.expression)
            if match2:
                self.fermion = True
                name = opname[match2.group("name")]
                indexSpin = match2.group("spin")
                indexGauge = match2.group("gauge")  # None
                indexColf = match2.group("colf")  # maybe None
                indexFlav = match2.group("flav")
                self.spinorIndices = []
                self.indices.append(indexSpin)
                self.spinorIndices.append(indexSpin)
                if indexGauge:
                    indexGauge = indexGauge[1:-2]
                    self.indices.append(indexGauge)
                    self.spinorIndices.append(indexGauge)
                if indexColf:
                    indexColf = indexColf[1:-2]
                    self.indices.append(indexColf)
                    self.spinorIndices.append(indexColf)
                self.indices.append(indexFlav)
                self.spinorIndices.append(indexFlav)
                self.spinorIndices = tuple(self.spinorIndices)
                spinorName = "unidentified"
                if name == opname["bar[e]"]:
                    spinorName = "of conjugated right-handed lepton"
                elif name == opname["bar[u]"]:
                    spinorName = "of conjugated right-handed up-type quark"
                elif name == opname["bar[b]"]:
                    spinorName = "of conjugated right-handed down-type quark"
                elif name == opname["bar[l]"]:
                    spinorName = "of conjugated SU2 Lepton Dublett"
                elif name == opname["bar[q]"]:
                    spinorName = "of conjugated SU2 Quark Dublett"
                self.description = "Spinor " + spinorName
                self.expression = self.expression[:match2.start()] + "{0:s}({1:s}{2:s},{3:s}{4:s})".format(name,
                                                                                                           indexSpin,
                                                                                                           "" if indexGauge == None else "," + indexGauge,
                                                                                                           "" if indexColf == None else indexColf + ",",
                                                                                                           indexFlav)
                self.expression = self.expression + match2.group("brackets")
                self.name = name
