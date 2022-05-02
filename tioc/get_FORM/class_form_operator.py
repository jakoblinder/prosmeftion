import re
import logging
import sys
from yaml import safe_load
from typing import Dict, List, Tuple

from tioc import CONFIG_PATH, opname_sorted, opname, opvalues, opSL2Cvalues, escape_regex, field_config, coeff_config
from .class_index import Index

logger_autoeft = logging.getLogger("autoeft.projection")
logger = logger_autoeft.getChild("form_operator")

class Operator():
    """
    A class describing the operators in a term which is previously transformed via form.
    Possible Operators:
        {H, e, u, b, l, q} & {[H+], [e_C+], [u_C+], [d_C+], [L+], [Q+]}
        field strength tensors of U(1),SU(2) and SU(3) named B, W and G
        completely antisymmetric symbol su2eps(gauge1,gauge2), su3eps, ...
        completely symmetric symbol su2dK(gauge1,gauge2), su3dK, ...
        couplings that carry indices (only Yukawa couplings): {yu, yd, ye, [yu+], [yd+], [ye+]}

    Attributes
    ----------
    tex : str
    description : str

    name : str
        Name of the Operator.
    indices : [Index, Index, ...]
        All type of indices in an operator. Occurring always at most one time in this list even, BECAUSE there shouldn't exist. Note: All indices in a term are always fully contracted.
        Example for indices: [lor1, lor34532, spin34,...].
    fermion : bool
        Determines whether an operator is an fermion, i.e. a anti-commuting object.
    expression : str
        Explicit expression of the operator. If manual is true. The expression has to contain an valid FORM expression,
        containing only functions which are defined by the declaration functions.
    tex : str
        LaTex expression of the operator.
    description : str
        Describes the type of operator.

    Magicmethods
    ------------
    __str__()
        Specify the format for printing an Operator object with the print command or converting a Operator object in
        string with str() function.
    __repr__()
        Specify the string representation of the object: Here the same as the str function itself.
    __format__()
        Specify the format for "format" function in print statement: Here the same as the print statement itself.
    Methods
    -------
    identifyHiggs()
        Identifies all occurring Higgs and rewrites them according to the "opname" dictionary.
    identifyFieldstrengthtensors()
        Identifies all occurring fieldstrengthtensors and rewrites them according to the "opname" dictionary.
    identifySpinors()
        Identifies all occurring spinors and rewrites them according to the "opname" dictionary.

    identifyFlavorMatrices()
        Identifies all occurring flavor matrices and rewrites them according to the "opname" dictionary.
    identifySUNepsandKroneckerDelta()
        Search for suNeps (Epsilontensor of SUN) and suNdK (Kronecker-Delta of SUN) and sl2C tensors and identify indices:
            [su2eps](gauge6153,gauge5345), [su2dK](gauge6153,gauge5345), ...
    """
    expression : str
    indices : List[Index]
    tex : str
    fermion : bool
    name : str
    description : str

    def __init__(self, expression, **kwargs):
        self.expression = expression
        self.indices = []  # list with Index-typed objects of used indices in operator. Note: At this stage an
        # operator shouldn't have itself contracted indices.
        # e.g. [lor1, lor34532, spin34,...]
        self.fermion = None
        self.name = None
        self.description = None

        if kwargs["type"] == "field":
            self.numID = kwargs["numID"]  # Denotes the number of the field, beginning at 1. Tensors don't have such a numID.
            self.identifyHiggs()
            self.identifyFieldstrengthtensors()
            self.identifySpinors()
            self.gaugeIndicesforProjection()
            self.n_D = self.get_NDerivative()
        elif kwargs["type"] == "tensor":
            pass
            self.identifyFlavorMatrices()
            self.identifySUNepsandKroneckerDelta()

    def __str__(self):
        """Specify the format for printing with str() or print() statement function. """
        return f"{self.expression:s}"
    def __repr__(self):
        return self.__str__()
    def __format__(self, key):
        """Specify the format for "format" function in print statement: Here the same as the print statement itself."""
        return self.__str__()

    def identifyHiggs(self):
        """
        Identify Higgs (Higgsboson) H and indices.
        Returns
        -------
        """
        # Search for normal Higgs first.
        pattern = r"(?P<Higgs>" + f"{opname['H']}" + r"|" + f"{escape_regex(opname['conj[H]'])}" + r")\((?P<Indices>[a-zA-Z\d,]+)\)"
        match = re.match(pattern, self.expression)
        if match:
            indices = match.group('Indices').split(',')
            name = match.group('Higgs')
            # print(f"{name}{indices}")
            self.fermion = False
            for ind in indices:
                self.indices.append(Index(ind))
            self.indices = tuple(self.indices)
            self.description = field_config["H"]["description"]
            self.name = name
            if name == "H":
                self.tex = field_config["H"]["tex"][0]
            elif name == "[H+]":
                self.tex = field_config["H"]["tex"][1]

    def identifyFieldstrengthtensors(self):
        """
        Identify for Fieldstrengthtensors {V, F, G} and Indices:
        For SU(2):
            VL & VR(gaugeA13,gaugeA15,Usldot13,Usldot14).
        For U(1):
            FL(Lsl5,Lsl6) & FR
        For SU(3):
            GL & GR(cola?,Usldot13,Usldot14) but doesn't seem to occur.
        Returns
        -------
        """
        pattern = r"(?P<name>[" + f"{opname['V']}{opname['F']}{opname['G']}" + r"][LR])\((?P<Indices>[a-zA-Z\d,]+)\)"  #
        match = re.match(pattern, self.expression)
        if match:
            self.fermion = False
            name = match.group("name")
            indices = match.group("Indices").split(',')
            for ind in indices:
                self.indices.append(Index(ind))
            self.indices = tuple(self.indices)
            # print(f"{name}{indices}")
            self.tex = field_config[name]["tex"][0]
            groupName = "unidentified Group"
            if name == f"{opname['F']}L" or name == f"{opname['F']}R":
                groupName = "U(1)"
            elif name == f"{opname['V']}L" or name == f"{opname['V']}R":
                groupName = "SU(2)"
            elif name == f"{opname['G']}L" or name == f"{opname['G']}R":
                groupName = "SU(3)"
            self.description = "Fieldstrengthtensor for {0:s}.".format(groupName)
            self.name = name

    def identifySpinors(self):
        """
        Search for spinors {[e_C], [u_C], [d_C], L, Q} & {[e_C+], [u_C+], [d_C+], [L+], [Q+]} and identify indices:

        right-handed leptons: [e_C](spin5610,flav5608).
        right-handed up-type quarks: [u_C](spin3806,colf3079,flav2775).
        right-handed down-type quarks: [d_C](spin3787,colf3065,flav2791).
        SU2 Lepton Dublett: L(spin5952,gauge6153,flav5951).
        SU2 Quark Dublett: [Q+](spin2763,gauge2764,colf2891,flav2762).

        Returns
        -------

        """
        pattern = r"(?P<name>(" + f"{'|'.join(map(escape_regex,opSL2Cvalues))}" + r"))\((?P<Indices>[a-zA-Z\d,]+)\)"
        match = re.match(pattern, self.expression)
        if match:
            self.fermion = True
            name = match.group("name")
            indices = match.group("Indices").split(',')
            for ind in indices:
                self.indices.append(Index(ind))
            self.indices = tuple(self.indices)
            # print(f"{name}{indices}")
            spinorName = "unidentified"
            if name == "L":
                self.tex = field_config["L"]["tex"][0]
                spinorName = "of SU2 Lepton Dublett"
            elif name == "[L+]":
                self.tex = field_config["L"]["tex"][1]
                spinorName = "of conjugated SU2 Lepton Dublett"
            elif name == "Q":
                self.tex = field_config["Q"]["tex"][0]
                spinorName = "of SU2 Quark Dublett"
            elif name == "[Q+]":
                self.tex = field_config["Q"]["tex"][1]
                spinorName = "of conjugated SU2 Quark Dublett"
            elif name == "[e_C]":
                self.tex = field_config["eC"]["tex"][0]
                spinorName = "of conjugated right-handed lepton"
            elif name == "[e_C+]":
                self.tex = field_config["eC"]["tex"][1]
                spinorName = "of right-handed lepton"
            elif name == "[u_C]":
                self.tex = field_config["uC"]["tex"][0]
                spinorName = "of conjugated right-handed up-type quark"
            elif name == "[u_C+]":
                self.tex = field_config["uC"]["tex"][1]
                spinorName = "of right-handed up-type quark"
            elif name == "[d_C]":
                self.tex = field_config["dC"]["tex"][0]
                spinorName = "of conjugated right-handed down-type quark"
            elif name == "[d_C+]":
                self.tex = field_config["dC"]["tex"][1]
                spinorName = "of right-handed down-type quark"
            self.description = "Spinor " + spinorName
            self.name = name

    def identifyFlavorMatrices(self):
        """
        Search for y[ude] (FlavorMatrices, i.e. Yukawa matrices) {yu, yd, ye} and identify indices:
            yu(flav1234,flav5678).
        Search for complex conjugate Yukawa matrix {[yu+], [yd+], [ye+]} and identify indices:
            [yu+](flav5678, flav1234).
        Returns
        -------

        """
        flavormatrices = list(opname_sorted["coefficients"].values())[:6]
        pattern = r"(?P<name>(" + f"{'|'.join(map(escape_regex, flavormatrices))}" + "))\((?P<Indices>[a-zA-Z\d,]+)\)"
        match = re.match(pattern, self.expression)  # Use "match", because the Yukawa coupling matrix should be all
        if match:
            self.fermion = False
            name = match.group("name")
            indices = match.group("Indices").split(',')
            for ind in indices:
                self.indices.append(Index(ind))
            self.indices = tuple(self.indices)
            self.name = name
            # print(f"{name}{indices}")
            if "yu" in name:
                self.description = coeff_config["flavormatrix"]["yu"]["description"]
                if name == "yu":
                    self.tex = coeff_config["flavormatrix"]["yu"]["tex"][0]
                else:
                    self.tex = coeff_config["flavormatrix"]["yu"]["tex"][1]
            elif "ye" in name:
                self.description = coeff_config["flavormatrix"]["ye"]["description"]
                if name == "ye":
                    self.tex = coeff_config["flavormatrix"]["ye"]["tex"][0]
                else:
                    self.tex = coeff_config["flavormatrix"]["ye"]["tex"][1]
            elif "yd" in name:
                self.description = coeff_config["flavormatrix"]["yd"]["description"]
                if name == "yd":
                    self.tex = coeff_config["flavormatrix"]["yd"]["tex"][0]
                else:
                    self.tex = coeff_config["flavormatrix"]["yd"]["tex"][1]

    def identifySUNepsandKroneckerDelta(self):
        """
        Search for suNeps (Epsilontensor of SUN) and suNdK (Kronecker-Delta of SUN) and sl2C tensors and identify indices:
            [su2eps](gauge6153,gauge5345), [su2dK](gauge6153,gauge5345), ...
        Returns
        -------
        """
        tensors = [tensor["FORM"] for tensor in coeff_config["tensor"].values()]
        pattern = r"(?P<name>(" + f"{'|'.join(map(escape_regex, tensors))}" + "))\((?P<Indices>[a-zA-Z\d,]+)\)"
        match = re.match(pattern, self.expression)  # Use "match", because the su2eps should be all in the expression of the operator.
        if match:
            self.fermion = False
            name = match.group("name")
            indices = match.group("Indices").split(',')
            for ind in indices:
                self.indices.append(Index(ind))
            self.indices = tuple(self.indices)
            self.name = name[1:-1]
            self.description = coeff_config["tensor"][name[1:-1]]["description"]
            self.tex = coeff_config["tensor"][name[1:-1]]["tex"][0]

    def gaugeIndicesforProjection(self):
        """
        The gauge indices of a field should be denoted by the pattern idxF2I1, where A denotes the second field
        (remember the unique ordering by helicity) and one the first index of this field. I.e. idxF2I1 denotes the first
        gauge index of the second field in this term.
        Note that the indices of the different gauge group are not distinguished and thus should always be part of the
        index object.

        Returns
        -------

        """
        pre = "idx"
        counter = {"gauge": 1, "colf": 1}
        # indtype = "gauge"
        for indtype in ["gauge", "colf"]:
            for index in self.indices:
                if index.typ == indtype:
                    index.projection = f"{pre:s}F{self.numID:d}I{counter[indtype]:d}"
                    counter[indtype] += 1

    def get_NDerivative(self):
        """
        Determine number of derivatives acting on a field by remembering that each derivative must have an usldot and an lsl index.
        By taking the minimal number of usldot and lsl indices the number of derivatives can be unambiguously determined.
        Returns
        -------
        Number of derivatives acting on a field.
        """
        n_usldot = 0  # number of usldot indices
        n_lsl = 0  # number of lsl indices
        for ind in self.indices:
            if ind.typ == "Usldot":
                n_usldot += 1
            elif ind.typ == "Lsl":
                n_lsl += 1

        return min(n_usldot, n_lsl)