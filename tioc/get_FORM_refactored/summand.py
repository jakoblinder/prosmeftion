import logging.config
import re
import sys
from abc import ABC, abstractmethod
from copy import copy
from fractions import Fraction
from typing import Dict, List, Tuple

from tioc import model, op_config, index_config, get_SUN_name
from .operator import Tensor, Field
from .coefficient import Coefficient
from .indices import Indices_Summand, Indices_Operator, Possible_Indices
from .index import Dummy_Index
from .operators import Tensors, Fields

logger_autoeft = logging.getLogger("autoeft.projection")
logger = logger_autoeft.getChild(__name__)

class Summand_Model(ABC):  # Tensor, Field, Coefficient
    """Base class for single term consisting of an overall coefficient and a product of operators."""
    tensors: List[Tensor]
    fields: Tuple[Field]
    coeff: Coefficient
    fieldcounter: Dict[str, int]
    d: int  # massdimension of operator
    indices: Indices_Summand
    nD: int  # number of derivatives in the summand

    @abstractmethod
    def __init__(self, tensors: List[str], fields: Tuple[str], coeff: str, fp_name: str):
        """Load in fields, tensors and the coefficient."""
        self.name = fp_name
        self.tensors = tensors
        self.fields = fields
        self.coeff = coeff
        self.nD
        self.fieldcounter
        self.fieldstructure

        self.possible_indices = Possible_Indices(
            [index for index in list(self.indices).copy() if not isinstance(index, Dummy_Index)])

    @abstractmethod
    def __repr__(self):
        """Specify the format the general string representation and for printing with repr()."""
        tensor = "*".join(map(str, self.tensors))
        contractedOp = "*".join(map(str, self.fields))
        return f"{tensor:s}*{contractedOp:s}"

    def __str__(self):
        """Specify the format for printing with str() or print() statement function: Here the same as the string representation repr() itself."""
        return self.__repr__()

    def __format__(self, key):
        """Specify the format for "format" function in print statement: Here the same as the string representation repr() itself."""
        if key == "tex":
            return self.tex
        elif key == "complete" or key == "c":
            coeff = self.coeff
            tensor = "*".join(map(str, self.tensors))
            contractedOp = "*".join(map(str, self.fields))
            if tensor:
                return f"({coeff})*{tensor:s}*{contractedOp:s}"
            else:
                return f"({coeff})*{contractedOp:s}"
        else:
            return self.__repr__()

    @property
    def tex(self):
        """
        Tex expression of the entire summand. However, before this can be generated, the tex indices are generated
        using generators, with the range specified in the index configuration file.
        Returns
        -------
            LaTex expression.
        """
        # Instantiate an index generator for each type of index
        generator = {}
        for index in self.indices:
            indrange = Indices_Summand.get_tex_range(index_config[index.typ]["tex_indices"])
            if indrange:
                generator[index.typ] = Indices_Summand.infinite_Indices(indrange)
            else:
                generator[index.typ] = Indices_Summand.infinite_numIndices(index_config[index.typ]["tex_indices"])

        tex_indices = {}  # Dictionary for unique tex names of indices.
        for index in self.indices:
            try:
                if index.typ == "sbasis":
                    ind = index_config
                    tex_indices[index.indname] = f"{index_config['sbasis']['tex_indices']}_{{{index.id}}}"
                else:
                    tex_indices[index.indname] = next(generator[index.typ])
            except StopIteration:
                logger.error("Specified index range is to small to map all occurring indices.")
                sys.exit("STOP")

        #write latex index for all indices:
        for operator in self.tensors:
            for index in operator.indices:
                index.tex = tex_indices[index.indname]
        for operator in self.fields:
            for index in operator.indices:
                index.tex = tex_indices[index.indname]

        # Build tex expression:
        mult_sign = ''  # '*'
        tex_expr = ""
        tensors = mult_sign.join([f"{op:tex}" for op in self.tensors])
        fields = mult_sign.join([f"{op:tex}" for op in self.fields])
        tex_expr += rf"{self.coeff:tex}\\" + mult_sign
        if tensors:
            tex_expr += f"{tensors}{mult_sign}{fields}"
        else:
            tex_expr += f"{fields}"

        return tex_expr

    @property
    def tensors(self):
        return self._tensors

    @tensors.setter
    def tensors(self, value: List):
        self._tensors = Tensors(tuple(map(Tensor, value)))

    @property
    def fields(self):
        return self._fields

    @fields.setter
    def fields(self, fp_fields):
        self._fields = Fields(tuple(map(Field, fp_fields, range(1, len(fp_fields)+1))))

    @property
    def coeff(self):
        return self._coeff

    @coeff.setter
    def coeff(self, value):
        # TODO: Abstraction into operators
        self._coeff = Coefficient(value, self.name)

    @property
    def d(self):
        """
        Mass dimension of operator.
        """
        return self.coeff.mdim

    @property
    def fieldcounter(self):
        """
        Returns dictionary, e.g.
        {'BL': 0, 'WL': 0, 'L': 0, 'Q': 0, 'dC': 0, 'eC': 0, 'uC': 0, 'H': 1,
        'H+': 1, 'L+': 0, 'Q+': 0, 'dC+': 0, 'eC+': 0, 'uC+': 0, 'BL+': 0, 'WL+': 0}
        which tells that there are one ordinary and one conjugated Higgs field in the Summand.
        Returns
        -------
        """
        fieldcount = {key.name: 0 for key in model.fields.values()}
        autoeft_to_projection = {autoeft: projection for field_name, field in {**op_config["bosonfields"], **op_config["fermionfields"]}.items() if field_name != "D" for projection, autoeft in field["autoeft"].items()}
        for field in self.fields:
            try:
                for counted_field in fieldcount.keys():
                    if autoeft_to_projection[counted_field] == field.name:
                        fieldcount[counted_field] += 1
                        break
            except KeyError:
                logger.error(f"Field {field.name:s} doesn't have an autoeft translation.")
                sys.exit("STOP")

        return fieldcount

    @property
    def fieldcounter_stripped(self):
        """
        Remove all 0 entries, e.g.:
        {'BL': 0, 'WL': 0, 'L': 0, 'Q': 0, 'dC': 0, 'eC': 0, 'uC': 0, 'H': 1,
        'H+': 1, 'L+': 0, 'Q+': 0, 'dC+': 0, 'eC+': 0, 'uC+': 0, 'BL+': 0, 'WL+': 0}
        -> {'H': 1, 'H+': 1}
        """
        fieldcount_s = {field: count for field, count in self.fieldcounter.items() if count != 0}

        return fieldcount_s

    @property
    def fieldstructure(self):
        """
        Get the exact field and derivative structure, e.g. term D(uC)*H*H*d+ has structure:
        ((uC,1),(H,0),(H,0),(d+,0))
        I.e. a tuple of as many tuple as there are fields in the Summand is return, where each tuple give first,
        the kind of field at this position and second, the number of derivatives acting on it.
        """
        projection_to_autoeft = {projection: autoeft for field_name, field in
                              {**op_config["bosonfields"], **op_config["fermionfields"]}.items() if field_name != "D"
                              for projection, autoeft in field["autoeft"].items()}
        fieldstructure = []
        for field in self.fields:
            try:
                structure = (projection_to_autoeft[field.name] , field.nD)
                fieldstructure.append(structure)
            except KeyError:
                logger.error(f"Field {field.name:s} doesn't have an autoeft translation.")
                sys.exit("STOP")

        return tuple(fieldstructure)

    @property
    def nD(self):
        """ Number of derivatives in a summand."""
        n_D = 0
        for field in self.fields:
            n_D += field.nD
        return n_D

    @property
    def indices(self):
        """All indices of a summand."""
        indices_expr = Indices_Operator([])
        for tensor in self.tensors:
            indices_expr += tensor.indices
        for field in self.fields:
            indices_expr += field.indices
        if any(isinstance(index, Dummy_Index) for index in indices_expr):
            # Allow uncontracted indices when there is a Dummy_index in the indices.
            return Indices_Summand(indices_expr.indices, allow_uncontracted=True)
        else:
            return Indices_Summand(indices_expr.indices)


    def replace_SUN_indices_by_projection_indices(self):
        """
        Replace all SUN-indices by the corresponding projection indices, which should be saved at this time in the
        projection attribute of each index.
        E.g.: [su2eps](gauge1,gauge2)*H(gauge1)*[H+](gauge2) -> [su2eps](gaugeF1I1,gaugeF2I1)*H(gaugeF1I1)*[H+](gaugeF2I1)
        Returns
        -------
        """
        for sun_group in ["gauge", "colf"]:
            for tensor in self.tensors:
                for index in tensor.indices[sun_group]:
                    index.expr = index.projection
            for field in self.fields:
                for index in field.indices[sun_group]:
                    index.expr = index.projection


class Summand(Summand_Model):
    """Single term consisting of an overall coefficient and a product of operators."""
    fieldcounter: Dict[str, int]

    def __init__(self, tensors: List[str], fields: List[str], coeff: str, fp_name: str):
        super().__init__(tensors, fields, coeff, fp_name)

        # Gauge fields contain special projection index of the form idxF2I1, which is transmitted to the contracted tensors:
        self.gaugeIndicesforProjection_tensors()

        su2 = get_SUN_name(2)
        su3 = get_SUN_name(3)
        self.gaugeTensorsSUN = {su2: self.tensors[su2], su3: self.tensors[su3]}

    def __repr__(self):
        """Specify the format the general string representation and for printing with repr()."""
        return super().__repr__()

    def get_projectionIndex_from_Field(self, indexID):
        """
        Returns the projection index of the contracted gauge field with index specified by indexID.
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
