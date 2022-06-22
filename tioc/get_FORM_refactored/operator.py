import re
import logging
import sys
from yaml import safe_load
from typing import Dict, List, Tuple
from abc import ABC, abstractmethod
from collections.abc import MutableMapping
from copy import copy

from tioc import CONFIG_PATH, op_config, escape_regex, model, index_config, op_pattern, index_pattern, dummy_index_pattern, op_name_pattern
from .index import Index, Dummy_Index
from .indices import Indices_Operator
from tioc import index_number_pattern as inp

logger_autoeft = logging.getLogger("autoeft.projection")
logger = logger_autoeft.getChild(__name__)

class Operator_Model(Index):
    """
    A class describing the operators in a term which is previously transformed via form.
    Possible Operators:
        {H, e, u, b, l, q} & {[H+], [e_C+], [u_C+], [d_C+], [L+], [Q+]}
        field strength tensors of U(1),SU(2) and SU(3) named B, W and G
        completely antisymmetric symbol su2eps(gauge1,gauge2), su3eps, ...
        completely symmetric symbol su2dK(gauge1,gauge2), su3dK, ...
        couplings that carry indices (only Yukawa couplings): {yu, yd, ye, [yu+], [yd+], [ye+]}
    """
    expr: str
    indices: Indices_Operator
    name: str
    isconj: bool
    non_conj_name: str
    tex: str
    description: str

    @abstractmethod
    def __init__(self, expr: str):
        self.expr = expr

    @abstractmethod
    def __repr__(self):
        """Specify the format the general string representation and for printing with repr()."""
        return f"{self.expr:s}"

    def __str__(self):
        """Specify the format for printing with str() or print() statement function: Here the same as the string representation repr() itself."""
        return self.__repr__()

    def __format__(self, key):
        """Specify the format for "format" function in print statement: Here the same as the string representation repr() itself."""
        if key == "tex":
            return self.tex
        else:
            return self.__repr__()

    @property
    @abstractmethod
    def tex(self):
        """Create tex expression of operator without derivatives."""
        return self.texed_op_wo_der(subscript_indices=("lor", "Lsldot", "Lsl", "gauge", "colf"), superscript_indices=("Usldot", "Usl", "gaugeadj", "cola", "flav"))

    def texed_op_wo_der(self, subscript_indices: Tuple[str], superscript_indices: Tuple[str]):
        """Create tex expression of operator without derivatives and with specified sub and superscript indices."""
        # Get texed name of operator
        tex_expr_op = self.tex_name

        if isinstance(self, Tensor):
            sub_indices = [[f"{index:tex}" for index in self.indices[ind_typ]] for ind_typ in subscript_indices]
            super_indices = [[f"{index:tex}" for index in self.indices[ind_typ]] for ind_typ in superscript_indices]
        elif isinstance(self, Field):
            sub_indices = [[f"{index:tex}" for index in self.indices[ind_typ] if not index.derIndex] for ind_typ in
                                 subscript_indices]
            super_indices = [[f"{index:tex}" for index in self.indices[ind_typ] if not index.derIndex] for ind_typ in
                             superscript_indices]
        else:
            logger.error("Class of operator unknown.")
            sys.exit("STOP")

        if r"_"  in tex_expr_op or r"^"  in tex_expr_op:
            tex_expr_op = f"({tex_expr_op})"
        tex_expr = tex_expr_op

        if any([len(sub_index_typ) for sub_index_typ in sub_indices]):
            tex_expr += f"_{{{', '.join([' '.join(sub_index_typ) for sub_index_typ in sub_indices if len(sub_index_typ) != 0])}}}"
        if any([len(super_index_typ) for super_index_typ in super_indices]):
            tex_expr += f"^{{{', '.join([' '.join(super_index_typ) for super_index_typ in super_indices if len(super_index_typ) != 0])}}}"
        return tex_expr

    @staticmethod
    def read_in_operator(expression:str) -> (Tuple[str, bool, str], Tuple[Index], int, Tuple[Index]):
        der_indices = []
        cov = escape_regex(op_config["fermionfields"]["D"]["mathematica"]["cov"])
        # ReadinMethod for cov(lor1234 derivative
        match_lor = re.finditer(r"(" + cov + r"\((?P<index>lor" + inp + r"),)", expression)
        matches_lor = list(match_lor)  # contain derivative indices
        # ReadinMethod for cov(Lsl, Usldot derivative
        ind_pat = f"(Lsl(?!dot)|Usl(?!dot)|Lsldot|Usldot)" + inp
        match_sl2C = re.finditer(r"(" + cov + r"\((?P<index1>" + ind_pat + r"),(?P<index2>" + ind_pat + r"))", expression)
        matches_sl2C = list(match_sl2C)
        if any(matches_lor):
            for i, v in enumerate(matches_lor):
                index = v.group("index")
                der_indices.append(Index(index, derIndex=i+1))
            nD = len(matches_lor)
            op = expression[matches_lor[-1].end(): (-1)*len(matches_lor)]
        elif any(matches_sl2C):
            for i, v in enumerate(matches_sl2C):
                index1 = v.group("index1")
                index2 = v.group("index2")
                der_indices.append(Index(index1, derIndex=i+1))
                der_indices.append(Index(index2, derIndex=i+1))
            nD = len(matches_sl2C)
            op = expression[matches_sl2C[-1].end()+1: (-1) * len(matches_sl2C)]
        else:
            nD = 0
            op = expression
        op_indices = []
        match = re.match(r"(?P<name>" + op_name_pattern + ")", op)  # \((" + index_pattern + r",?)+\)
        if match:
            name = match.group("name")
            # Test if operator is a conjugated one:
            for operator in {**op_config["tensors"], **op_config["bosonfields"], **op_config["fermionfields"]}.values():
                form_names = list(operator["mathematica"].values())
                if len(form_names) == 1 and name == form_names[0]:
                    non_conj_name = name
                    isconj = False
                    break
                elif len(form_names) == 2:
                    if name == form_names[1]:
                        non_conj_name =  form_names[0]
                        isconj = True
                        break
                    elif name == form_names[0]:
                        non_conj_name = name
                        isconj = False
                        break
            indices = op[match.end()+1:-1]
            indices = indices.split(",")
            for index in indices:
                # use match here in order to assure that one matches the start of the index.
                match = re.match(r"(?P<index>" + index_pattern + r")", index)
                if match:
                    index_matched = match.group("index")
                    op_indices.append(Index(index_matched))
                else:
                    match_dummy = re.match(dummy_index_pattern, index)
                    if match_dummy:
                        index_matched = Dummy_Index(int(match_dummy.group("number")))
                        op_indices.append(index_matched)
                    else:
                        logger.error("Index can not be identified.")
                        sys.exit("STOP")
        else:
            logger.error("No index found.")
            sys.exit("STOP")
        return (name, isconj, non_conj_name), tuple(op_indices), nD, tuple(der_indices)

    @property
    def expr(self):
        if type(self) == Tensor:
            expr = self.name
            expr += f"({','.join(map(str, self.indices))})"
        elif type(self) == Field:
            cov = op_config['fermionfields']['D']['mathematica']['cov']
            expr = ""
            if self.nD > 0:
                # Term contain derivatives
                for i in range(1, self.nD + 1):
                    expr += f"{cov}({','.join(map(str, [nonD_index for nonD_index in self.indices if nonD_index.derIndex == i]))},"
            n_brackets = self.nD * ")"
            expr += self.name
            non_Derivative_indices = [nonD_index for nonD_index in self.indices if not nonD_index.derIndex]
            expr += f"({','.join(map(str, non_Derivative_indices))})"

            expr += n_brackets

        # print(f"{self._expr} == {expr}")
        if self._expr != expr:
            logger.debug(f"Expression changed from {self._expr:s} to {expr:s}.")
        # assert self._expr == expr
        return expr
        # return self._expr

    @expr.setter
    def expr(self, fp_expr):
        """
        Set expression of the operator and extract the name and a string tuple of indices.
        Returns
        -------
        """
        names, op_indices, self.nD, der_indices = Operator_Model.read_in_operator(fp_expr)
        self.name, self.isconj, self.non_conj_name = names # names[0], names[1], names[2]
        # self.indices = Indices_Operator(der_indices + op_indices)
        # Indexstructure of the operator only without the derivative.
        ind_structure = [index.typ for index in op_indices]
        def assertion(der_indices, op_indices, fieldtype, non_conj_name):
            indices = Indices_Operator(der_indices + op_indices)
            # Indexstructure of the operator only without the derivative.
            ind_structure_op = [index.typ for index in op_indices]
            if fieldtype == "tensors":
                ind_structure = op_config[fieldtype][non_conj_name]["index_structure"]
            else:
                bosonsANDfermions = {**op_config["bosonfields"], **op_config["fermionfields"]}
                ind_structure = bosonsANDfermions[non_conj_name]["index_structure"]

            # Check that index structure matches:
            if ind_structure_op not in ind_structure:
                structure_match = [[1 if ind_structure_op[i] == should_index else 0 for i, should_index in enumerate(ind_struc)] for ind_struc in ind_structure]
                dummy_in_indices = [True if index_typ == "dummy" else False for index_typ in ind_structure_op]
                if all(dummy_in_indices):
                    # Since dummy indices should in general only occur in tensors and those are rewritten such that
                    # they only have fundamental indices, i.e. only indices of the same type, for each tensor which
                    # could occur here the 'ind_structure' should consist of only one unique possible structure.
                    # Thus, even for tensors with only dummy indices, the will_be_typ of the index can be determined.
                    # As a small check it is ensured that the type of the field is really a tensor.
                    assert fieldtype == "tensors"
                if any(dummy_in_indices):
                    best_matches = [sum(match) for match in structure_match]
                    if max(best_matches) > 0 or all(dummy_in_indices):
                        best_match = ind_structure[best_matches.index(max(best_matches))]
                        for i, should_be_typ in enumerate(best_match):
                            if op_indices[i].typ == "dummy":
                                op_indices[i].will_be_typ = should_be_typ

                        assertion(der_indices, op_indices, fieldtype, non_conj_name)
                else:
                    logger.error("Indexstructure doesn't match the required structure for this field.")
                    sys.exit("STOP")

            return indices

        if type(self) == Tensor:
            self.indices = assertion(der_indices, op_indices, "tensors", self.non_conj_name)
        elif type(self) == Field:
            self.indices = assertion(der_indices, op_indices, "field", self.non_conj_name)

        self._expr = fp_expr

    @property
    def tex_name(self):
        """Create tex expression of operator."""
        all_ops = {**op_config["tensors"], **op_config["bosonfields"], **op_config["fermionfields"]}
        if self.isconj:
            try:
                tex_expr = all_ops[self.non_conj_name]["tex_hc"]
            except KeyError:
                tex_expr = all_ops[self.non_conj_name]["tex"] + r"^{\dagger}"
        else:
            tex_expr = all_ops[self.non_conj_name]["tex"]

        return tex_expr


    @property
    def description(self):
        """Gives description to the operator."""
        all_ops = {**op_dict["tensors"], **op_dict["bosonfields"], **op_dict["fermionfields"]}
        description_expr = all_ops[self.non_conj_name]["description"]
        return description_expr

    @property
    def autoeft(self):
        """Returns autoeft name of the field/ tensor."""
        all_ops = {**op_dict["tensors"], **op_dict["bosonfields"], **op_dict["fermionfields"]}
        if self.isconj:
            form_field = list(all_ops[self.non_conj_name]["mathematica"].values())[1]
        else:
            form_field = list(all_ops[self.non_conj_name]["mathematica"].values())[0]
        autoeft_expr = all_ops[self.non_conj_name]["autoeft"][form_field]
        return autoeft_expr


class Tensor(Operator_Model):
    expr: str
    indices: Indices_Operator
    name: str
    tex: str
    description: str

    def __init__(self, expr: str):
        super().__init__(expr)

    def __repr__(self):
        """Specify the format the general string representation and for printing with repr()."""
        return super().__repr__()

    @property
    def tex(self):
        """Create tex expression of operator without derivatives."""
        return super().texed_op_wo_der(subscript_indices=("lor", "Lsldot", "Lsl", "gauge", "colf"),
                                    superscript_indices=("Usldot", "Usl", "gaugeadj", "cola", "flav"))

class Field(Operator_Model):
    expr: str
    field_pos: int  # Position of the Fiel starting at 1.
    indices: Indices_Operator
    name: str
    tex: str
    ac: bool
    description: str
    nD: int

    def __init__(self, expr: str, pos: int):
        super().__init__(expr)
        self.field_pos = pos
        self.gaugeIndicesforProjection()  # Write projection Indices for gauge indices of field.

    def __repr__(self):
        """Specify the format the general string representation and for printing with repr()."""
        return super().__repr__()

    @property
    def ac(self) -> bool:
        """Specififes whether Field commutes or anticommutes."""
        ac_expr = op_config["fermionfields"][self.non_conj_name]["ac"]
        return ac_expr

    @property
    def tex(self):
        """Create tex expression of operator with derivatives."""
        op_texed = super().texed_op_wo_der(subscript_indices=("lor", "Lsldot", "Lsl", "gauge", "colf"),
                                       superscript_indices=("Usldot", "Usl", "gaugeadj", "cola", "flav"))

        derIndices = Indices_Operator([index for index in self.indices if index.derIndex].copy())
        if not derIndices:
            # no derivatives in the operator
            return op_texed
        else:
            derIndices_sorted = {}
            for derivativeIndex in derIndices:
                try:
                    derIndices_sorted[derivativeIndex.derIndex].append(derivativeIndex)
                except KeyError:
                    derIndices_sorted[derivativeIndex.derIndex] = [derivativeIndex]
            tex_derivatives = ""
            cov_tex = op_config["fermionfields"]["D"]["tex"]
            for key, index in derIndices_sorted.items():
                if index[0].typ == "lor":
                    tex_derivatives += cov_tex + "_{" + f"{index[0]:tex}" + "}"
                elif index[0].typ in ["Lsldot", "Lsl"]:
                    tex_derivatives += cov_tex + "_{" + f"{index[0]:tex}" + "}" + "^{" + f"{index[1]:tex}" + "}"
                elif index[0].typ in ["Usldot", "Usl"]:
                    tex_derivatives += cov_tex + "_{" + f"{index[1]:tex}" + "}" + "^{" + f"{index[0]:tex}" + "}"


            return tex_derivatives + op_texed

    def gaugeIndicesforProjection(self):
        """
        The gauge indices of a field should be denoted by the pattern idxF2I1, where '2' denotes the second field
        (remember the unique ordering by helicity) and '1' the first index of this field, i.e. idxF2I1 denotes the first
        gauge index of the second field in this term.
        Note that the indices of the different gauge group are not distinguished and thus should always be part of the
        index object.

        Returns
        -------

        """
        pre = ""  # previously this was 'idx'
        counter = {"gauge": 1, "colf": 1}
        for indtype in ["gauge", "colf"]:
            for index in self.indices:
                if index.typ == indtype:
                    index.projection = f"{pre:s}F{self.field_pos:d}I{counter[indtype]:d}"
                    counter[indtype] += 1


