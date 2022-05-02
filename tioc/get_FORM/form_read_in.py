import re
import sys
import logging.config

from pathlib import Path
from fractions import Fraction
from typing import Dict, List, Tuple
from copy import copy

from tioc import PROJECTION_PATH
from .class_index import Indices
from .class_form_operator import Operator
from tioc.get_BS.class_coefficient import Coefficient
from .. import opname, opnameSL2C

logger_autoeft = logging.getLogger("autoeft.projection")
logger = logger_autoeft.getChild(__name__)

class Term_s(Coefficient):
    """
    Single term consisting of an overall coefficient and a product of operators.
    """
    fields: List[str]
    tensors: List[str]
    coeff: Coefficient
    fieldcounter: Dict[str, int]

    def __init__(self, fields, tensors, coeff):
        fields_tmp = []
        for i, field in enumerate(fields):
            fields_tmp.append(Operator(field, type="field", numID=i+1))
        self.fields = tuple(fields_tmp)
        self.fieldcounter, self.n_D = self.get_fieldcounts()
        tensors_tmp = []
        for tensor in tensors:
            tensors_tmp.append(Operator(tensor, type="tensor"))
        self.tensors = tuple(tensors_tmp)
        # Gauge fields contain special projection index of the form idxF2I1, which is transmitted to the contracted tensors:
        self.gaugeIndicesforProjection_tensors()
        self.coeff = coeff
        #TODO: Get Name of SU2_W out of model file
        self.gaugeTensorsSUN = {"SU2_W": self.get_SUN_tensors(2), "SU3_C": self.get_SUN_tensors(3)}

    def __str__(self):
        """Specify the format for printing with str() or print() statement function. """
        tensor = "*".join(map(str,self.tensors))
        contractedOp = "*".join(map(str,self.fields))
        return f"{tensor:s}*{contractedOp:s}"  # {str(self.coeff)}*
        # return f"{self.cops:s}"  # {str(self.coeff)}*

    def __repr__(self):
        return self.__str__()

    def __format__(self, key):
        """Specify the format for "format" function in print statement: Here the same as the print statement itself."""
        return self.__str__()

    def get_projectionIndex_from_Field(self, indexID):
        """
        Returns the projection index of the contracted gauge field with index specified by indexID.
        Parameters
        ----------
        indexID

        Returns
        -------

        """
        for field in self.fields:
            for index in field.indices:
                if index.id == indexID:
                    return index.projection

    def gaugeIndicesforProjection_tensors(self):
        """
        Gauge fields contain special projection index of the form idxF2I1, which is transmitted to the contracted tensors.

        Returns
        -------

        """
        for tensor in self.tensors:
            for index in tensor.indices:
                for indtype in ["gauge", "colf"]:
                    if index.typ == indtype:
                        index.projection = self.get_projectionIndex_from_Field(index.id)

    def get_SUN_tensors(self, N):
        """
        Returns coefficient tensors (i.e. epsilon and Kronecker delta) of the group SU("N") where the epsilontensor is
        and the Kronecker delta are rewritten as the regarding tensors in FORM, i.e. e_(...) and d_(...,...).
        -------
        """
        epsilon = f"su{N:d}eps"
        delta = f"su{N:d}dK"
        suNeps = []
        suNdK = []
        for tensor in self.tensors:
            if tensor.name == epsilon:
                suNeps.append(f"e_({','.join([index.projection for index in tensor.indices])})")
            elif tensor.name == delta:
                suNdK.append(f"d_({','.join([index.projection for index in tensor.indices])})")

        form_tensors = "*".join(suNeps + suNdK)

        return form_tensors

    def get_fieldcounts(self):
        """
        Return dictionary which contains number of fields for each field type:
        Returns
        -------

        """
        #TODO: Translation before in FORM?
        translate = {f"{opname['F']}L": "BL",
                     f"{opname['G']}L": "GL",
                     f"{opname['V']}L": "WL",
                     opnameSL2C["[d_C]"]: "dC",
                     opnameSL2C["[e_C]"]: "eC",
                     opnameSL2C["L"]: "L",
                     opnameSL2C["Q"]: "Q",
                     opnameSL2C["[u_C]"]: "uC",
                     opname["H"]: "H",
                     opname["conj[H]"]: "H+",
                     opnameSL2C["[d_C+]"]: "dC+",
                     opnameSL2C["[e_C+]"]: "eC+",
                     opnameSL2C["[L+]"]: "L+",
                     opnameSL2C["[Q+]"]: "Q+",
                     opnameSL2C["[u_C+]"]: "uC+",
                     f"{opname['F']}R": "BL+",
                     f"{opname['G']}R": "GL+",
                     f"{opname['V']}R": "WL+"}
        fieldcount = {key: 0 for key in translate.values()}
        n_D = 0
        for field in self.fields:
            for sfield in translate.keys():
                if field.name == sfield:
                    fieldcount[translate[sfield]] += 1
            n_D += field.n_D # number of derivatives in a term
        return fieldcount, n_D

class Term_form(Term_s):
    """
    A term from the form output can consist of many summand. The term ist therefore split in Term_s objects
    consisting of only one summand.
    """
    cops: str
    coeff: Coefficient
    def __init__(self, cops, coeff, name):
        self.cops = cops
        self.coeff = coeff
        self.name = name
        self.terms = []
        for term in self.get_terms():
            self.terms.append(Term_s(term["cops"], term["tensors"], term["coeff"]))
        self.indices = self.get_indices()

    def __str__(self):
        """Specify the format for printing with str() or print() statement function. """
        return f"{self.terms}"  # {str(self.coeff)}*
    def __repr__(self):
        return self.__str__()
    def __format__(self, key):
        """Specify the format for "format" function in print statement: Here the same as the print statement itself."""
        return self.__str__()

    def get_terms(self):
        """
        The terms of the form output are extracted and written in individual Term_form objects. Since Terms may consist
        of multiple summand, each summand is stored as an Term_s object and the coefficients are separated for each
        Term_s object.
        For Example: A term like a*(c*A*B - d*C*D) is written in to different Term_s object:
        1. Term_s object: Coefficient = a*c, Fields = [A,B]
        2. Term_s object: Coefficient = - a*d, Fields = [C,D].
        Returns
        -------

        """
        # print(self.cops)
        zero = re.match(r"0", self.cops)
        match = re.finditer(r"\n", self.cops)
        terms = []
        term = {"cops": [""], "tensors": [""], "coeff": ""}
        if zero:
            operators = []
            return [{"cops": ["0"], "tensors": ["0"], "coeff": "0"}]
        elif match:
            matches = list(match)
            operators = []
            for i, m in enumerate(matches):
                if i > 0:
                    cop = self.cops[matches[i - 1].end():m.start()]
                else:
                    # Append first operator as the FIRST element in the list.
                    cop = self.cops[:m.start()]
                cop_match = re.match(r"\*", cop)
                if cop_match and len(cop) > 1:
                    cop = cop[1:]
                if cop: # Removes empty lines
                    operators.append(cop)
            # Append also last operator:
            cop = self.cops[matches[-1].end():]
            cop_match = re.match(r"\*", cop)
            if cop_match and len(cop) > 1:
                cop = cop[1:]
            operators.append(cop)
            # Extract coefficent and single terms
            coefficients = []  # Indices where the coefficients are
            start_coeff = None
            end_coeff = 0
            while end_coeff != len(operators)-1:
                for i in range(end_coeff, len(operators)):
                    match = re.search(r".+\*\(", operators[i])
                    if match:
                        operators[i] = operators[i][:-2]
                        start_coeff = i+1
                        break
                for i in range(start_coeff, len(operators)):
                    if operators[i] == ")":
                        end_coeff = i
                        break
                coefficients.append((start_coeff, end_coeff))

            for i,c in enumerate(coefficients):
                start_coeff, end_coeff = c
                if i == 0:
                    fields = operators[:start_coeff]
                else:
                    fields = operators[coefficients[i-1][1]+1:start_coeff]
                coeff_ops = operators[start_coeff:end_coeff]
                del start_coeff, end_coeff

                if fields[0] == "+":
                    overall_sign = +1
                    fields.pop(0)
                elif fields[0] == "-":
                    overall_sign = -1
                    fields.pop(0)
                # extract sumand in coefficients
                coeff_range = []
                start_coeff_range = None
                end_coeff_range = 0
                while end_coeff_range != len(coeff_ops) - 1:
                    for i in range(end_coeff_range, len(coeff_ops)):
                        if coeff_ops[i] == "+" or coeff_ops[i] == "-":
                            start_coeff_range = i
                            break
                    for i in range(start_coeff_range+1, len(coeff_ops)):
                        if coeff_ops[i] == "+" or coeff_ops[i] == "-":
                            end_coeff_range = i - 1
                            break
                        elif i == len(coeff_ops)-1:
                            end_coeff_range = len(coeff_ops)-1
                            break
                    coeff_range.append((start_coeff_range, end_coeff_range))

                coeff_range_ops = []
                for i, c in enumerate(coeff_range):
                    start_coeff, end_coeff = c
                    coeff_range_ops.append(coeff_ops[start_coeff:end_coeff+1])

                for tensors in coeff_range_ops:
                    for i, coeff in enumerate(tensors):
                        if coeff == "+" or coeff == "-":
                            if coeff == "+":
                                sign = overall_sign * (+1)
                            else:
                                sign = overall_sign * (-1)
                            match = re.match(r"\d{1,3}/\d{1,3}", tensors[i + 1])
                            if match:
                                # Sign and fraction as coefficient
                                factor = f"({str(Fraction(float(Fraction(tensors[i + 1])) * sign))})"
                                tensors.pop(0)
                                tensors.pop(0)
                                break
                            else:
                                # Only a sign coefficient
                                factor = f"({overall_sign})"
                                tensors.pop(0)
                                break
                    numerical_coefficient = copy(self.coeff) # Necessary to make a copy, to avoid multiple multiplications of factor
                    numerical_coefficient *= factor
                    numerical_coefficient.update_coefficient()
                    terms.append({"cops": fields, "tensors":tensors, "coeff": numerical_coefficient})  # numerical_coefficient

        return terms

    def get_indices(self):
        """
        Extract all indices occurring in a term, i.e. all indices in the coefficient tensors which may be contracted and
        all indices in the field tensors which may not be contracted. The indices are extracted as Index object and it is
        ensured that each index occurs at most ones, before they are store all together in an indices object which assigns
        them automatically unique LaTex indices. These indices are then reassigned to all the indices of the fields.
        Returns
        -------
        Returns a list of all occurring indices.
        """
        indices_check_tensor = []  # Different sumands should(?) contain the same indices.
        indices_check_field = []
        for sumand in self.terms:
            if sumand.coeff == "0" or sumand.fields[0].expression == "0":
                continue
            indices_field_tmp = ()
            indices_tensor_tmp = []
            for field in sumand.fields:
                # print(field)
                # print(field.indices)
                indices_field_tmp += field.indices
            for tensor in sumand.tensors:
                indices_tensor_tmp += tensor.indices
            tensor_check = [index.id for index in indices_tensor_tmp]
            field_check = [index.id for index in indices_field_tmp]
            if not set(field_check).issubset(set(tensor_check)):
                difference = set(field_check) - set(tensor_check)
                error = f"Indices of fields ({','.join(field_check)}) should be all also in the indices of the coefficients ({','.join(tensor_check)}) but ({','.join(list(difference))}) is/ are not in the tensor indices."
                logger.error(error)
                logger.debug(f"The corresponding expression is: {sumand}")
                raise AssertionError(error)
            indices_check_tensor.append(Indices(indices_tensor_tmp))  # Append the indices of the coefficient tensors, because there could be more tensor indices then field indices.
            indices_check_field.append(Indices(indices_field_tmp))

        if len(indices_check_tensor) == 0 or len(indices_check_field) == 0:
            return Indices([])

        if len(indices_check_tensor) > 1:
            pivot = indices_check_tensor[0]
            check = []
            for i in indices_check_tensor[1:]:
                check.append(pivot == i)
            if not any(check):
                # Since terms like BL(Lsl1,Lsl2)*(...+...) + BR(Usldot1,Usldot2)*(...+...) are possible, not all coefficient tensor have the same indices.
                # Thus, it isn't checked that all coefficient tensors should have the same indices, but at least one should have in this case.
                logger.error("There are not the same indices in the coefficient tensors, in all sumands of the term.")
                sys.exit("STOP")
        if len(indices_check_field) > 1:
            pivot = indices_check_field[0]
            check = []
            for i in indices_check_field[1:]:
                check.append(pivot == i)
            if not any(check):
                logger.error("There are not the same indices in the fields, in all sumands of the term.")
                sys.exit("STOP")

        indices_tensor = indices_check_tensor[0]
        for i in indices_check_tensor[1:]:
            indices_tensor += i
        indices_field = indices_check_field[0]
        for i in indices_check_field[1:]:
            indices_field += i
        indices_all = indices_field + indices_tensor
        for index in indices_all.indices:
            index.tex_name = indices_all.tex_indices[index.name]
        for sumand in self.terms:
            for field in sumand.fields:
                for ind in field.indices:
                    ind.tex_name = indices_all.indices[indices_all.indices.index(ind)].tex_name
            for tensor in sumand.tensors:
                for ind in tensor.indices:
                    ind.tex_name = indices_all.indices[indices_all.indices.index(ind)].tex_name

        # print("=========================")
        # for sumand in self.terms:
        #     for field in sumand.fields:
        #         for ind in field.indices:
        #             print(f"{ind:tex}")
        #     for tensor in sumand.tensors:
        #         for ind in tensor.indices:
        #             print(f"{ind:tex}")

        return indices_all
