# -*- coding: utf-8 -*-
"""
Created on Mon Nov  29 8:22:04 2021.

@author: Jakob Linder
"""
import os
import subprocess
import sys
import re
import logging.config
import multiprocessing as mp

from itertools import chain
from pathlib import Path
from typing import List

from .sl2c.class_term import Term
from . import coeffvalues, opname_sorted, opname, opnameSL2C, opvalues, opSL2Cvalues, spinorsSL2C_c, spSL2C_c_values, n_der
from . import PROJECTION_PATH, CONFIG_PATH, FORM_PATH, FORM_GENERAL_PATH, INPUT_PATH, LATEX_PATH, AUTOEFT_PATH
from .yProjection.read_write import get_terms, mathematica_to_form, writefile
from .yProjection.coefficient import Factor
from .sun_projection import equalize_field_indices
from .general import declaration_SL2C_sets
from .pyform.pyformfunction import pyForm

# logger_autoeft = logging.getLogger("autoeft")
# logger = logging.getLogger("autoeft.projection")
logger_autoeft = logging.getLogger("autoeft.projection")
logger = logger_autoeft.getChild(__name__)

# coupling constants: {lambdah,At,g1,g2,g3,mu,lambdaphi,kappa} and sometimes an "I"
# Ms = mass of the singlet
# Mu = \[Mu] = tachyonic mass of the Higgs
# eps = epsilon of DREG
# couplings that carry indices (only Yukawa couplings): {yu,yd,ye,conj[yu],conj[yd],conj[ye]}
# {H,e,u,b,l,q} & conj[H] and bar[psi], which are the Higgs,
#                                                 the right-handed leptons,
#                                                 the right-handed up-type quarks,
#                                                 the right-handed down-type quarks,
#                                                 the SU(2) lepton doublet and
#                                                 the SU(2) quark doublet.
# Name[{ind1},{ind2},{ind3},...,{indN}]
# {lor,spin,gauge,gaugeadj,colf,cola,flav} = {4d Lorentz index,
#                                             spinor index,
#                                             SU(2)-index in fundamental,
#                                             SU(2)-index in adjoint,
#                                             SU(3)-index in fundamental,
#                                             SU(3)-index in adjoint,
#                                             flavor index}
# Numbers attached to the index name are arbitrary, but all should be properly contracted
# field strength tensors of U(1),SU(2) and SU(3) named F, V and G or B, W and G
# covariant derivative cov[{lor},Field] of any field called Field
# completely antisymmetric symbol su2eps[{gauge1,gauge2}]
# group generators TT[{gaugeadj1},{gauge1, gauge2}]
# Gammamatrices, e.g.: gamma[{lor1}, {spin6497, spin6498}]
# Dictionary for names in form.


def findOpandCoeff(expression):
    """
    Give one long string "expression" and separate each term and for each term also the coefficient and coperators (i.e.
    fully contracted operators) of the term.

    Parameters
    ----------
    expression : str
        One long string from Benjamin Summ's output without whitespaces
        and curly braces at the beginning and end.

    Returns
    -------
    coefficient : [str, str, str, ...]
        Extracted coefficients.
    coperators : [str, str, str, ...]
        Extracted coperators.

    """
    # Idea: Use pattern0 to search for each summend and find afterwards fully contracted operators (=coperators) inside the terms.
    pattern0 = r"(/\(\d{0,3}\*?eps\*Ms(\^\d{1,2})?\)|/Ms(\^\d{1,2})?)"
    match = re.finditer(pattern0, expression)  # Find Terms by identifying always the quotient: /(6*eps*Ms^4).

    if match:
        matches = list(match)
        terms = []
        for i, v in enumerate(matches):
            if i > 0:
                terms.append(expression[matches[i - 1].end():v.end()])
            else:
                terms.append(expression[:v.end()])
                # print(v.group())

    n_terms = writefile("groups.txt", terms, write = False)  # number of terms found
    logger.debug(f"Number of terms = {n_terms:d}")
    pattern1 = r"(\*(H|F|V|G|TT|yu|yd|ye|conj\[yu\]|conj\[yd\]|conj\[ye\]|su2eps|e\[|u\[|b\[|l\[|q\[|conj\[H\]|bar\[e\]|bar\[u\]|bar\[b\]|bar\[l\]|bar\[q\]|cov\[|gamma\[)).*" + pattern0
    # terms = ["+(At*(-32*(1+eps)*kappa*Ms^4*mu+32*At^3*(2+eps)*(Ms^2+2*Mu^2)-32*At^2*mu*(2*Mu^2+eps*(Ms^2+2*Mu^2))+At*(3*(4+3*eps)*g2^2*Ms^4-8*(2+eps)*lambdaphi*Ms^4+16*(1+eps)*Ms^2*mu^2+64*(2+eps)*kappa*(Ms^4+Ms^2*Mu^2+Mu^4))-2*eps*(-8*kappa*Ms^4*mu-16*At^2*mu*Mu^2+16*At^3*(Ms^2+2*Mu^2)+At*(3*g2^2*Ms^4-4*lambdaphi*Ms^4+4*Ms^2*mu^2+32*kappa*(Ms^4+Ms^2*Mu^2+Mu^4)))*L[Ms^2,\[Mu]^2])*cov[{lor1},cov[{lor1},H[{gauge35226}]]]**H[{gauge35053}]**conj[H][{gauge35053}]**conj[H][{gauge35226}])/(32*eps*Ms^8)"]
    coperator = []
    coefficient = []
    for term in terms:
        match1 = re.search(pattern1, term)
        if match1:
            match2 = re.search(pattern0, match1.group())
            if match2:
                op = match1.group()[:match2.start()]
                coeff = match1.group()[match2.start():]
                if op[-1] == ")" and match1.group()[0] == "*":  # remove "*" in front of and ")" at the end of expression
                    coperator.append(op[1:-1])
                else:
                    coperator.append(op)
            else:
                coperator.append(match1.group())
                raise ImportError("Not all coperators and coefficient are successfully identified in the term")
            c = term[:match1.start()] + coeff
            if c[0] == "(":
                c = c[1:]
            elif c[0] == "-" and c[1] == "(":
                c = "-" + c[2:]
            elif c[0] == "+" and c[1] == "(":
                c = "+" + c[2:]
            coefficient.append(c)
        else:
            coperator.append("")
            raise ImportError("Not all terms are successfully identified in the expression")
    n_coperator = writefile("../coperators.txt", coperator, write = False)  # number of coperators successfully matched
    assert n_coperator == n_terms, "The program should find an fully contracted operator in EACH term."
    assert len(coefficient) == len(coperator), "The number of coefficients should match the number of coperators."
    logger.debug(f"Number of matched terms = {n_coperator:d}")

    return coefficient, coperator

def form(term, filename, cout=False, nDer=4):
    """
    Writes each term, i.e. summend in BS output in a FORM compatible way in a form file termI.frm which can be executed
    and further manipulated in FORM.
    Parameters
    ----------
    terms : Term
        Term object describing each term, i.e. summend in BS output.
    filename : str
        Specifies name of the directory in which FORM files are store and also the name of the main FORM file.

    Returns
    -------

    """
    TERM_PATH = FORM_PATH / filename
    TERM_PATH.mkdir(parents=True, exist_ok=True)  # Create directories if they don't exist.

    # Declarations:
    # To ensure that the indices for replacement of sigma2 and sigmabar2 by sigma and sigmabar are declared,
    # the method has to be called already now.
    fieldstrengthtensorDerivativeHandling1 = term.form_fieldstrengthtensorDerivativeHandling(nDer)
    if fieldstrengthtensorDerivativeHandling1:
        term.update_sortIndices()

    declarationsIndexGroups = "* Auxiliary sets for identity statements:\n"
    for indSet in ["lor", "spin", "spinA", "gauge", "gaugeadj", "colf", "cola", "flav"]:
        declarationsIndexGroups += "Set {0:s}: ".format(indSet) + ",".join(term.possible_indices[indSet]) + ";\n"
    simplifyFlavorMatrices1 = None  # term.form_simplifyFlavorMatrices()
    if simplifyFlavorMatrices1:
        simplifyFlavorMatrices, flMinitialize = simplifyFlavorMatrices1
        flMinitialize = "\n* FlavorMatrix Simplification\n" + flMinitialize
        declarationsIndexGroups += flMinitialize
    # Write declarations file:
    with open(TERM_PATH / "declarations.h", "w") as file:
        file.write(form_declarations() + declarationsIndexGroups)

    # To ensure that the indices for replacement of sigma2 and sigmabar2 by sigma and sigmabar are declared,
    # the method has to be called already now.
    simplifySigma2 = term.form_simplifySigma2()

    # SL2C Declarations
    declarationsSL2C = ""
    # Higgs part
    term.form_higgsDeclarations(nDer)
    # Declarations for derivative replacements
    term.form_spinorDeclarations(nDer)
    declarationsSL2C += "* Auxiliary variables for identity statements:\n"
    for ind in ["Usl", "Lsl", "Usldot", "Lsldot"]:
        declarationsSL2C += "Indices " + "=2,".join(term.possible_indices[ind]) + "=2;\n"
    declarationsSL2C += "\n"
    for indSet in ["Usl", "Lsl", "Usldot", "Lsldot"]:
        declarationsSL2C += "Set {0:s}: ".format(indSet) + ",".join(term.possible_indices[indSet]) + ";\n"
    alternatingList = {}
    alternatingList["ULsl"] = list(chain.from_iterable(zip(term.possible_indices["Usl"], term.possible_indices["Lsl"])))
    alternatingList["LUsl"] = list(chain.from_iterable(zip(term.possible_indices["Lsl"], term.possible_indices["Usl"])))
    alternatingList["ULsldot"] = list(chain.from_iterable(zip(term.possible_indices["Usldot"], term.possible_indices["Lsldot"])))
    alternatingList["LUsldot"] = list(chain.from_iterable(zip(term.possible_indices["Lsldot"], term.possible_indices["Usldot"])))
    for indSet in ["ULsl", "LUsl", "ULsldot", "LUsldot"]:
        declarationsSL2C += f"Set {indSet:s}: " + ",".join(alternatingList[indSet]) + ";\n"
    # Indices <Lsl1=2>,...,<Lsl100=2>;
    # Indices <Usl1=2>,...,<Usl100=2>;
    # Indices <Lsldot1}=2>,...,<Lsldot100=2>;
    # Indices <Usldot1=2>,...,<Usldot100=2>;
    # Autodeclare Indices Lsl=2, Usl=2, Lsldot=2, Usldot=2;\n
    # Set Usl: Usl1,...,Usl100;
    # Set Lsl: Lsl1, ..., Lsl100;
    # Set ULsl: <Usl1, Lsl1>, ..., <Usl100, Lsl100>;
    # Set LUsl: <Lsl1, Usl1>, ..., <Lsl100, Usl100>;
    # Set Usldot: Usldot1, ..., Usldot100;
    # Set Lsldot: Lsldot1, ..., Lsldot100;
    # Set ULsldot: <Usldot1, Lsldot1>, ..., <Usldot100, Lsldot100>;
    # Set LUsldot: <Lsldot1, Usldot1>, ..., <Lsldot100, Usldot100>;
    declarationsSL2C += "\n"
    with open(TERM_PATH / "declarationsSL2C.h", "w") as file:
        file.write(declarationsSL2C)

    formexpression = ""
    # General declarations of fields and so on:
    formexpression += "#-\n"
    formexpression += "#include declarations.h\n"
    formexpression += "#include declarationsSL2C.h\n"
    formexpression += "#+\n"

    # Expression of the Coperator:
    formexpression += term.get_form()

    # Rewrite derivatives of Higgsfields with placeholders:
    higgsDerivativetoCommutative1 = term.form_higgsDerivativetoCommutative(nDer)
    if higgsDerivativetoCommutative1:
        formexpression += "*\n* Derivative is rewritten with placeholders for the Higgsfield so that the Higgsfield replacements can be done with commuting fields.\n*\n"
        initialize, higgsDerivativetoCommutative = higgsDerivativetoCommutative1
        # Write file:
        with open(TERM_PATH / "higgsDerivativetoCommutative.prc", "w") as file:
            file.write(higgsDerivativetoCommutative)
        formexpression += initialize + "\n"
        formexpression += "#call higgsDerivativetoCommutative\n"
    # Replace derivatives by derivatives in SL2C notation:
    higgsDerivativetoSL2C = term.form_higgsDerivativetoSL2C()
    if higgsDerivativetoSL2C:
        with open(TERM_PATH / "higgsDerivativetoSL2C.prc", "w") as file:
            file.write(higgsDerivativetoSL2C)
        formexpression += "#call higgsDerivativetoSL2C\n"
        if cout: formexpression += "Print;\n"
        formexpression += ".sort\n\n"

    # Rewrite derivatives of spinors with placeholders:
    spinorDerivativetoCommutative1 = term.form_spinorDerivativetoCommutative(nDer)
    if spinorDerivativetoCommutative1:
        formexpression += "*\n* Derivative is rewritten with placeholders for the spinorfield so that the spinorfield replacements can be done with commuting fields.\n*\n"
        spinor_initialize, spinorDerivativetoCommutative = spinorDerivativetoCommutative1
        # Write file:
        with open(TERM_PATH / "spinorDerivativetoCommutative.prc", "w") as file:
            file.write(spinorDerivativetoCommutative)
        formexpression += spinor_initialize + "\n"
        formexpression += "#call spinorDerivativetoCommutative\n"
        formexpression += "\n"
    # Replace derivatives by derivatives in SL2C notation:
    spinorDerivativetoSL2C = term.form_spinorDerivativetoSL2C(nDer)
    if spinorDerivativetoSL2C:
        formexpression += "* Replace derivatives acting on a spinorfield by SL2C notated derivative and sigma matrices.\n"
        # Write file:
        with open(TERM_PATH / "spinorDerivativetoSL2C.prc", "w") as file:
            file.write(spinorDerivativetoSL2C)
        formexpression += "#call spinorDerivativetoSL2C\n"
        if cout: formexpression += "Print;\n"
        formexpression += ".sort\n\n"

    # Rewrite derivatives of fieldstrengthtensors with placeholders:
    if fieldstrengthtensorDerivativeHandling1:
        formexpression += "*\n* Derivative is rewritten with placeholders for the spinorfield so that the spinorfield replacements can be done with commuting fields.\n*\n"
        fieldstrengthtensor_initialize, fieldstrengthtensorDerivativetoCommutative, fieldstrengthtensorDerivativetoSL2C, fieldstrengthtensorCommutativetoDerivative = fieldstrengthtensorDerivativeHandling1
        # Write file:
        with open(TERM_PATH / "fieldstrengthtensorDerivativetoCommutative.prc", "w") as file:
            file.write(fieldstrengthtensorDerivativetoCommutative)
        with open(TERM_PATH / "fieldstrengthtensorDerivativetoSL2C.prc", "w") as file:
            file.write(fieldstrengthtensorDerivativetoSL2C)
        with open(TERM_PATH / "fieldstrengthtensorCommutativetoDerivative.prc", "w") as file:
            file.write(fieldstrengthtensorCommutativetoDerivative)
        formexpression += fieldstrengthtensor_initialize + "\n"
        formexpression += "#call fieldstrengthtensorDerivativetoCommutative\n"
        formexpression += "\n"
        formexpression += "#call fieldstrengthtensorDerivativetoSL2C\n"
        formexpression += "\n"

    # Replace two in the adjoint representation contracted SU2 Generators by Kronecker deltas in fundamental indices.
    replaceSU2Generators = term.form_replaceSU2Generators()
    if replaceSU2Generators:
        formexpression += replaceSU2Generators
        formexpression += "\n"
    # Replace two in the adjoint representation contracted SU3 Generators by Kronecker deltas in fundamental indices.
    replaceSU3Generators = term.form_replaceSU3Generators()
    if replaceSU3Generators:
        formexpression += replaceSU3Generators
        formexpression += "\n"
    # Substitute spinors and gamma matrices with the corresponding SL2C spinors:
    substituteSpinors = term.form_substituteSpinors()
    if substituteSpinors:
        with open(TERM_PATH / "substituteSpinors.prc", "w") as file:
            file.write(substituteSpinors)
        formexpression += "#call substituteSpinors\n"
        formexpression += "\n"

    higgsReplacebyEps = term.form_higgsReplacebyEps()
    if higgsReplacebyEps:
        formexpression += higgsReplacebyEps
        if cout: formexpression += "Print;\n"
        formexpression += ".sort\n\n"

    # Simplify contracted SU2 epsilon tensors:
    simplifyEpsSU2 = term.form_simplifyEpsSU2()
    if simplifyEpsSU2:
        with open(TERM_PATH / "simplifyEpsSU2.prc", "w") as file:
            file.write(simplifyEpsSU2)
        formexpression += "#call simplifyEpsSU2\n\n"
        if cout: formexpression += "Print;\n"

    # Replace sigma2 and sigmabar2 by sigma and sigmabar
    if simplifySigma2:
        formexpression += "* Replace sigma2 and sigmabar2 by sigma and sigmabar.\n"
        with open(TERM_PATH / "simplifySigma2.prc", "w") as file:
            file.write(simplifySigma2)
        formexpression += "#call simplifySigma2\n"

    # Replace contracted sigmas by SL2C epsilons
    formexpression += "* Replace contracted sigmas by epsilons.\n"
    with open(TERM_PATH / "replaceSigmabyEps.prc", "w") as file:
        file.write(term.form_replaceSigmabyEps())
    formexpression += "#call replaceSigmabyEps\n"

    # Replace contracted SL2C Epsilons by Kronecker deltas
    with open(TERM_PATH / "simplifyEps.prc", "w") as file:
        file.write(term.form_simplifyEps())
    formexpression += "#call simplifyEps\n"

    # Symmetric properties of Fieldstrengthtensors:
    # TODO: Only has to be done if there exists a fieldstrengthtensor -> write in function like fieldstrengthtensorDerivativeHandling
    formexpression += "*\n* Discard fieldstrengthtensor which are contracted with an SL2C epsilontensor, because the SL2C indices are symmetric\n*\n"
    formexpression += "id {0:s}CL?FieldS(?a, Lsl1?LUsl[k], Lsl2?LUsl[m])*[su2eps](Usl1?ULsl[k], Usl2?ULsl[m]) = 0;\n".format(opname["F"])
    formexpression += "id {0:s}CL?FieldS(?a, Lsl2?LUsl[m], Lsl1?LUsl[k])*[su2eps](Usl1?ULsl[k], Usl2?ULsl[m]) = 0;\n".format(opname["F"])
    formexpression += "id {0:s}CL?FieldS(?a, Lsldot1?LUsldot[k], Lsldot2?LUsldot[m])*[su2eps](Usldot1?ULsldot[k], Usldot2?ULsldot[m]) = 0;\n".format(opname["F"])
    formexpression += "id {0:s}CL?FieldS(?a, Lsldot2?LUsldot[m], Lsldot1?LUsldot[k])*[su2eps](Usldot1?ULsldot[k], Usldot2?ULsldot[m]) = 0;\n".format(opname["F"])
    formexpression += "\n"
    formexpression += "*\n* Discard fieldstrengthtensor of SU2 which are contracted with an epsilontensor, because the SU2 indices are symmetric\n*\n"
    formexpression += "id {0:s}CL?".format(opname["V"]) + "{" + "{0:s}CL,{0:s}CR".format(opname["V"]) + "}(op?, gauge1?gauge[k], gauge2?gauge[m], ?a)*[su2eps](gauge1?gauge[k],gauge2?gauge[m]) = 0;\n"
    formexpression += "id {0:s}CL?".format(opname["V"]) + "{" + "{0:s}CL,{0:s}CR".format(opname["V"]) + "}(op?, gauge2?gauge[m], gauge1?gauge[k], ?a)*[su2eps](gauge1?gauge[k],gauge2?gauge[m]) = 0;\n"
    formexpression += "\n"

    # Simplify Flavormatrices
    if simplifyFlavorMatrices1:
        formexpression += simplifyFlavorMatrices

    substituteSU3Spinors = term.form_substituteSU3Spinors()
    if substituteSU3Spinors:
        formexpression += "* Substitute hermitian conjugated fields wich have colour indices and transform in anti-fundamental representations by SU3 epsilon tensor and field with 2 color indices of the fundamental representation.\n"
        # Write file:
        with open(TERM_PATH / "substituteSU3Spinors.prc", "w") as file:
            file.write(substituteSU3Spinors)
        formexpression += "#call substituteSU3Spinors\n"

    # Simplify contracted SU3 epsilon tensors:
    simplifyEpsSU3 = term.form_simplifyEpsSU3()
    if simplifyEpsSU3:
        formexpression += "* Simplify SU3, i.e. 3 component epsilon tensor.\n"
        # Write file:
        with open(TERM_PATH / "simplifyEpsSU3.prc", "w") as file:
            file.write(simplifyEpsSU3)
        formexpression += "#call simplifyEpsSU3\n"

    formexpression += "\n"
    # Substitute placeholders by original fields inside of the derivatives:
    higgsCommutativetoDerivative = term.form_higgsCommutativetoDerivative(nDer)
    if higgsCommutativetoDerivative:
        # Write Higgs fields again inside of the derivative
        with open(TERM_PATH / "higgsCommutativetoDerivative.prc", "w") as file:
            file.write(higgsCommutativetoDerivative)
        formexpression += "#call higgsCommutativetoDerivative\n"
        formexpression += ".sort\n"
    spinorCommutativetoDerivative = term.form_spinorCommutativetoDerivative(nDer)
    if spinorCommutativetoDerivative:
        # Write spinor fields again inside of the derivative
        with open(TERM_PATH / "spinorCommutativetoDerivative.prc", "w") as file:
            file.write(spinorCommutativetoDerivative)
        formexpression += "#call spinorCommutativetoDerivative\n"
        formexpression += ".sort\n"
    formexpression += "\n"
    if fieldstrengthtensorDerivativeHandling1:
        formexpression += "#call fieldstrengthtensorCommutativetoDerivative\n"
        formexpression += "\n"

    formexpression += "* Write derivatives implicit with indices inside of fields.\n"
    with open(TERM_PATH / "derivativeasIndex.prc", "w") as file:
        file.write(Term.form_derivativeasIndex(nDer))
    formexpression += "#call derivativeasIndex\n"

    formexpression += "* Order fields by their helicity and then alpabetically.\n"
    # Write file:
    with open(TERM_PATH / "sortfields.prc", "w") as file:
        file.write(Term.form_sortfields(Term.extractOrder()))
    formexpression += "#call sortfields\n"

    formexpression += "* Write derivatives again outside of fields.\n"
    with open(TERM_PATH / "indexasDerivative.prc", "w") as file:
        file.write(Term.form_indexasDerivative(nDer))
    formexpression += "#call indexasDerivative\n"


    # Rearrange terms in a need way:
    commutingOps = list(opname_sorted["bosonfields"].values())  # opvalues[9:14]
    noncommutingOps = list(opname_sorted["fermionfields"].values())  # opvalues[14:]
    ops = ", ".join(commutingOps + noncommutingOps) + ", " + "{f:s}L, {f:s}R, {v:s}L, {v:s}R, {g:s}L, {g:s}R".format(f=opname["F"], v=opname["V"], g=opname["G"])
    formexpression += "Bracket " + ops + f", {opnameSL2C['[d_C]']}, {opnameSL2C['[e_C]']}, {opnameSL2C['L']}, {opnameSL2C['Q']}, {opnameSL2C['[u_C]']}" \
                                         f", {opnameSL2C['[d_C+]']}, {opnameSL2C['[e_C+]']}, {opnameSL2C['[L+]']}, {opnameSL2C['[Q+]']}, {opnameSL2C['[u_C+]']};\n"
    # formexpression += "Bracket sigma, sigmabar, [su2eps], [su2dK], [sl2Ceps], [sl2CdK];\n"
    formexpression += "Format 255;\n"  # Printing width is set to the maximum width to avoid linebreaks.
    formexpression += "Print +ss;\n"
    formexpression += ".end"

    # Write the main FORM file:
    with open(TERM_PATH / f"{filename}.frm", "w") as file:
        file.write(formexpression)

def form_declarations(nAuxIndices=100):
    """
    Contains all general declarations valid for any term, i.e. for example the declaration of all fields and indices.
    Returns
    -------
    form : str
        Content of the FORM file "declarations.h".
    """
    form = ""
    # write commuting and anticommuting operators of each term in separate list for initialization. Thus
    # Duplicated operators are removed.
    commutingOps = list(opname_sorted["bosonfields"].values())  # opvalues[9:14]
    noncommutingOps = list(opname_sorted["fermionfields"].values())  # opvalues[14:]
    paramOps = list(opname_sorted["coefficients"].values())  # opvalues[:9]
    commutingOps = "Function " + ", ".join(commutingOps) + ";\n"
    noncommutingOps = "Function " + ", ".join(noncommutingOps) + ";\n"
    paramOps = "CFunction " + ", ".join(paramOps) + ";\n"
    form += commutingOps  # Function H, [H+], F, V, G;
    form += "CFunction HC, [HC+];\n" # Auxiliary commuting Higgs
    form += noncommutingOps  # Function D, e, u, b, l, q, [ebar], [ubar], [bbar], [lbar], [qbar];
    form += "CFunction eC,uC,bC,lC,qC,[eCbar],[uCbar],[bCbar],[lCbar],[qCbar];\n"  # Auxiliary commuting spinors
    form += "\n"
    form += "Function " + ", ".join(opSL2Cvalues[0:10]) + ";\n"  # [e_C], [u_C], [d_C], L, Q, [e_C+], [u_C+], [d_C+], [L+], [Q+]
    form += "CFunction " + ", ".join(spSL2C_c_values) + ";\n"
    form += "Function {f:s}L, {f:s}R, {v:s}L, {v:s}R, {g:s}L, {g:s}R;\n".format(f=opname["F"], v=opname["V"], g=opname["G"])
    form += "CFunction {f:s}CL, {f:s}CR, {v:s}CL, {v:s}CR, {g:s}CL, {g:s}CR;\n".format(f=opname["F"], v=opname["V"], g=opname["G"])
    form += "Set FieldS: {f:s}CL, {f:s}CR, {v:s}CL, {v:s}CR, {g:s}CL, {g:s}CR;\n".format(f=opname["F"], v=opname["V"], g=opname["G"])
    form += paramOps  # CFunction yu, yd, ye, [yu+], [yd+], [ye+], [su2eps], T, gamma;
    form += "\n"
    form += "CFunction xi, [xi+], chi, [chi+];\n"
    form += "\n"
    form += "Set spinors: eC,uC,bC,lC,qC;\n"
    form += "Set spinorsAdj: [eCbar],[uCbar],[bCbar],[lCbar],[qCbar];\n"
    form += "Set spinorsAll: eC,uC,bC,lC,qC,[eCbar],[uCbar],[bCbar],[lCbar],[qCbar];\n"
    form += "\n"
    form += "Set spinorsN: " + ", ".join(spSL2C_c_values[0:5]) + ";\n"  # [e_C], [u_C], [d_C], L, Q
    form += "Set spinorsAdjN: " + ", ".join(spSL2C_c_values[5:10]) + ";\n"  # [e_C+], [u_C+], [d_C+], [L+], [Q+]
    form += "Set spinorsAllN: " + ", ".join(spSL2C_c_values[0:10]) + ";\n"  # [e_C], [u_C], [d_C], L, Q, [e_C+], [u_C+], [d_C+], [L+], [Q+]
    form += "\n"
    form += """AutoDeclare Indices lor      = 4; * 4d Lorentz index
AutoDeclare Indices lorA      = 4; * Auxiliary 4d Lorentz index
AutoDeclare Indices spin     = 4; * Index for Gamma matrices/ spinor index
AutoDeclare Indices spinA     = 2; * Auxiliary index for Gamma matrices/ spinor index
AutoDeclare Indices gauge    = 2; * SU(2)-index in fundamental
AutoDeclare Indices gaugeA    = 2; * Auxiliary SU(2)-index in fundamental
AutoDeclare Indices gaugeadj = 3; * SU(2)-index in adjoint
AutoDeclare Indices colf     = 3; * SU(3)-index in fundamental
AutoDeclare Indices colfA     = 3; * Auxiliary SU(3)-index in fundamental
AutoDeclare Indices cola     = 8; * SU(3)-index in adjoint
AutoDeclare Indices flav     = n; * flavor index"""
    form += "\n\n"
    form += "AutoDeclare Indices op; * auxiliary index for converting between commuting and noncommuting operators.\n\n"
    form += "* Declare some Symbols for pattern matching\n"
    form += "Symbols k,m;\n"
    form += "*\n* Indices and functions for derivatives in SL2C notation.\n*\n"
    form += "CFunction sigma, sigmabar, [sl2Ceps], [su3eps];\n"
    form += "CFunction sigma2, sigmabar2;\n"
    form += "* Auxiliary antisymmtric epsilons, used in combination with replace_.\n"
    form += "CFunction [sl2CepsA](antisymmetric), [su2epsA](antisymmetric), [su3epsA](antisymmetric);\n"
    form += "\n"
    form += "* Declare Kronecker Delta symbol for Sl2C Indices, because built in can not handle upper and lower (un-)dottet indices.\n"
    form += "* Since two indices are also symmetric when they are cyclic and vice versa and pattern matching is not allowed for symmetric function but for cyclic it is, [sl2CdK] is declared as cyclic.\n"
    form += "CFunction [sl2CdK](cyclic), [su2dK](cyclic), [su3dK](cyclic);\n" # [sl2CdK](cyclic)
    form += "\n\n"
    form += "Off Statistics;\n\n"

    return form

def run_form(filename):
    """
    Runs the FORM program given by the filename and returns the output, i.e. the translated expression.
    Parameters
    ----------
    filename : str
        Name of the FORM file and directory.
    Returns
    -------
        str: Converted expression.
    """
    TERM_PATH = FORM_PATH / filename
    try:
        formprocess = subprocess.run(
            ["form", "-p", TERM_PATH, TERM_PATH / f"{filename}.frm"],
            capture_output=True,  # If capture_output is true, stdout and stderr will be captured. When used,
            # the internal Popen object is automatically created with stdout=PIPE and stderr=PIPE. The
            # stdout and stderr arguments may not be supplied at the same time as capture_output. If
            # you wish to capture and combine both streams into one, use stdout=PIPE and stderr=STDOUT
            # instead of capture_output.
            text=True,  # output in stdout is now a string and not a byte sequence anymore
            check=True,  # If check is true, and the process exits with a non-zero exit code, a CalledProcessError
            # exception will be raised. Attributes of that exception hold the arguments, the exit code,
            # and stdout and stderr if they were captured.
        )
    except subprocess.CalledProcessError as exc:
        exc.cmd = list(map(str,list(exc.cmd)))
        logger.error(" ".join(exc.cmd) + "\n" + str(exc.stdout))
        logger.error(f"form returned non-zero exit status {exc.returncode}")
        sys.exit("STOP")
    else:
        # No Error occured
        logger.debug(" ".join(map(str,formprocess.args)))
        output = formprocess.stdout
        output = re.sub(r" *", "", output)  # Remove all whitespaces
        output = re.sub(r"\\", "",
                        output)  # Sometimes FORM splits indices in long expressions with an backslash "\" which is discarded.
        pattern = r"Print(\+s{1,2})?;\n{2}expr=\n{1,2}(?P<expression>(.|\n)*);"
        pattern_short = r"Print(\+s{1,2})?;\n{2}expr=(?P<expression>(.|\n)*);"  # Pattern for extremely short expressions, i.e. fitting in one line.
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

def print_form(filename, original, converted):
    p = "{0:s}: {1:s}\n".format(filename, original)
    p += "is converted into:\n"
    # Run FORM:
    p += converted
    logger.debug(p)
    # print(p)
    # print("======================")
    return p

def tensor_fields_indices(tensor_field, indices):
    """
    Returns latex expression of a tensor or field together with raised and lowered indices.
    Which indices are raised and lowered and in which order they are printed is specified here.
    Parameters
    ----------
    tensor_field: form_Operator
    indices: List[Index]

    Returns
    -------
        LaTex expression.
    """
    if tensor_field.expression == "0":
        return "0"
    if tensor_field.name in ["su2eps", "su3eps", "sl2Ceps", "su2dK", "su3dK", "sl2CdK"]:
        tex = f"{tensor_field.tex}"
    else:
        tex = f"({tensor_field.tex})"
    subscript_indices = {"Lsl": [], "gauge": [], "colf": []}
    superscript_indices = {"Usl": [], "gaugeadj": [], "cola": [], "flav": []}
    for index in indices:
        if index.typ in ["Lsldot","Lsl"]:
            subscript_indices["Lsl"].append(index)
        if index.typ in ["Usldot", "Usl"]:
            superscript_indices["Usl"].append(index)
    for index in indices:
        if index.typ == "gauge":
            subscript_indices["gauge"].append(index)
        if index.typ == "gaugeadj":
            superscript_indices["gaugeadj"].append(index)
    for index in indices:
        if index.typ == "colf":
            subscript_indices["colf"].append(index)
        if index.typ == "cola":
            superscript_indices["cola"].append(index)
    for index in indices:
        if index.typ == "flav":
            superscript_indices["flav"].append(index)
    subscript_indices_tex = ""
    for typ in subscript_indices.keys():
        if subscript_indices[typ] != []:
            for i in subscript_indices[typ]:
                subscript_indices_tex += i.tex_name + " "
            subscript_indices_tex = subscript_indices_tex[:-1] + ", "
    if subscript_indices_tex != "":
        tex += "_{" + f"{subscript_indices_tex[:-2]}" + "}"
    superscript_indices_tex = ""
    for typ in superscript_indices.keys():
        if superscript_indices[typ] != []:
            for i in superscript_indices[typ]:
                superscript_indices_tex += i.tex_name + " "
            superscript_indices_tex = superscript_indices_tex[:-1] + ", "
    if superscript_indices_tex != "":
        tex += "^{" + f"{superscript_indices_tex[:-2]}" + "}"

    return tex

def print_tex(Term, coeff=False):
    """
    Returns Latex output for one term which can consist of multiple summands.
    Parameters
    ----------
    Term : Term_form
        Object containing the term in SL2C notation and correctly assigned indices.
    Returns
    -------

    """
    tex = ""
    for summand in Term.terms:
        if summand.coeff == "0":
            continue
        tex += "\t"  # + r"\item $ "
        if coeff:
            # tex += "+\\" + "left(" + f"{summand.coeff.tex_expression:s}" + "\\" + "right)" + "*"
            tex += f"{summand.coeff.tex_expression:s}" + "*"
        for tensor in summand.tensors:
            tex += tensor_fields_indices(tensor, tensor.indices)
            tex += " "
        tex += r"\cdot "
        for field in summand.fields:
            tex += tensor_fields_indices(field, field.indices)
            tex += " "
            # tex += r"\cdot "
        tex = tex[:-1] + "\\"+ "\\" + "\n"  # " $\n"
    if tex == "":
        tex += "\t 0 \n"
        return tex
    # print(tex)

    return tex
    # else:
    #     return "\t\item $ 0 $\n"

def get_termobject(args):
    """
    Create Term object from given coefficient and contracted operator.
    Parameters
    ----------
    args
        args[0], args[1], args[2] = coeff, coperator, id
    Returns
    -------

    """
    coeff, coperator, id = args[0], args[1], args[2]
    return Term(coeff, coperator, id)

def convertviaform(term, n_der, max_dim):
    """
    Convert Term via FORM and return converted object.
    Parameters
    ----------
    term

    Returns
    -------

    """
    logger.info("Convert operators into SL2C notation.")
    # Write SL2C and set FORM-file:
    assert len(term) == 1
    summand = term[0]
    # FIXME: Correct index structure.
    sl2c_dirac_to_weyl, label = summand.form_convertDirac(max_dim)
    sl2c_derivative_in_SL2C = summand.form_convertDerivative(n_der)
    # TODO: Convert fieldstrength tensors.

    form_SL2C = declaration_SL2C_sets(term[0].possible_indices, n_der)
    send_to_form = [form_SL2C, f"{summand:c}", sl2c_dirac_to_weyl, sl2c_derivative_in_SL2C]
    expression = pyForm(FORM_GENERAL_PATH / "converttoSL2C.frm", send_to_form, input_dir=FORM_GENERAL_PATH, prompt= "READY")  # , debug=True, preprocessor_only=True
    expression = re.sub(r"(\s)*", "", expression)
    # Create FORM files:
    filename = term.name  # f"term{i:d}"
    # TODO: Update!
    form(term, filename)
    formoutput = run_form(filename)
    formoutput_formatted = print_form(filename,
                                      original=f"{term.cops_original:s}",
                                      converted=formoutput)
    ###
    TERM_PATH = FORM_PATH / filename
    with open(TERM_PATH / f"{filename}.h", "w") as file:
        formoutput=re.sub(r"(\s)*", "", formoutput)
        file.write(formoutput)

    new_term = get_terms(filepath=TERM_PATH / f"{filename}.h", as_one=True, name=term.name)
    for summand in new_term:
        summand.coeff *= Factor(term.coeff.expression)
    ###
    return formoutput_formatted, new_term  # Term_form(formoutput, term.coeff, term.name)

def converttoSL2C(terms, max_dim):
    """
    Output Terms of BSUOLEA are read in and formatted in SL2C Notation via FORM.
    Parameters
    ----------
    inputfile: str
        Filename of the inputfile in the input directory.
    max_dim
        Maximum mass dimension which occurs in the Lagrangian.
    pprint: bool
        Print formatted terms if True.

    Returns
    -------
    Array of formatted term_form objects.
    """


    # Write formfiles:
    logger.info("Run FORM")

    # with mp.Pool() as pool:
    #     terms_after_form = list(map(list, zip(*pool.map(convertviaform, terms))))
    # TODO:
    nDer = n_der(max_dim)
    terms_after_form = [list(convertviaform(term, nDer, max_dim)) for term in terms]  # list(map(list, zip(*map(convertviaform, terms))))
    ops = terms_after_form[0]
    form_terms = terms_after_form[1]

    # write a file containing all formatted operators written as a term to check the format:
    writefile(PROJECTION_PATH / "operators_formatted.txt", ops)
    del ops

    if pprint:
        for i, terms in enumerate(form_terms):
            print(f"Term {i:d}:")
            for term in terms.terms:
                print(f"{term:s}\n")
    return form_terms

# TODO: Implement progress bar: https://stackoverflow.com/questions/3160699/python-progress-bar
# import sys
#
# def progressbar(it, prefix="", size=60, file=sys.stdout):
#     count = len(it)
#     def show(j):
#         x = int(size*j/count)
#         file.write("%s[%s%s] %i/%i\r" % (prefix, "#"*x, "."*(size-x), j, count))
#         file.flush()
#     show(0)
#     for i, item in enumerate(it):
#         yield item
#         show(i+1)
#     file.write("\n")
#     file.flush()
# Usage:
#
# import time
#
# for i in progressbar(range(15), "Computing: ", 40):
#     time.sleep(0.1) # any code you need