import logging.config
import re
import sys
from abc import ABC, abstractmethod
from copy import copy
from fractions import Fraction
from typing import Dict, List, Tuple

from tioc import model, op_config
from .operator import Tensor, Field
from .coefficient import Coefficient
from .indices import Indices_Summand, Indices_Operator

logger_autoeft = logging.getLogger("autoeft.projection")
logger = logger_autoeft.getChild(__name__)

class Summand_Model(Tensor, Field, Coefficient):
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
            mult_sign = '*'
            tex_expr = ""
            tensors = mult_sign.join([f"{op:tex}" for op in self.tensors])
            fields = mult_sign.join([f"{op:tex}" for op in self.fields])
            tex_expr += rf"{self.coeff:tex}\\" + "*"
            if tensors:
                tex_expr += f"{tensors}{mult_sign}{fields}"
            else:
                tex_expr += f"{fields}"
            return tex_expr
        else:
            return self.__repr__()

    @property
    def tensors(self):
        return self._tensors

    @tensors.setter
    def tensors(self, value: List):
        self._tensors = list(map(Tensor, value))

    @property
    def fields(self):
        return self._fields

    @fields.setter
    def fields(self, fp_fields):
        self._fields = tuple(map(Field, fp_fields, range(1, len(fp_fields)+1)))

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
        fieldcount = {key.name: 0 for key in model.fields.values()}
        # op = {**op_config["bosonfields"], **op_config["fermionfields"]}
        autoeft_projection = {autoeft: projection for field_name, field in {**op_config["bosonfields"], **op_config["fermionfields"]}.items() if field_name != "D" for projection, autoeft in field["autoeft"].items()}
        for field in self.fields:
            try:
                for counted_field in fieldcount.keys():
                    if autoeft_projection[counted_field] == field.name:
                        fieldcount[counted_field] += 1
                        break
            except KeyError:
                logger.error(f"Field {field.name:s} doesn't have an autoeft translation.")
                sys.exit("STOP")

        return fieldcount

    @property
    def fieldcounter_stripped(self):
        # Remove all 0 entries:
        fieldcount_s = {field: count for field, count in self.fieldcounter.items() if count != 0}

        return fieldcount_s


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
        return Indices_Summand(indices_expr.indices)

class Summand(Summand_Model):
    """Single term consisting of an overall coefficient and a product of operators."""
    fieldcounter: Dict[str, int]

    def __init__(self, tensors: List[str], fields: List[str], coeff: str, fp_name: str):
        super().__init__(tensors, fields, coeff, fp_name)

        # Gauge fields contain special projection index of the form idxF2I1, which is transmitted to the contracted tensors:
        self.gaugeIndicesforProjection_tensors()

        def get_SUN_name(N):
            """Get Name of SU2_W out of model file."""
            for group_name, group_properties in model.sun_groups.items():
                if group_properties.N == N:
                    return group_name

        self.gaugeTensorsSUN = {get_SUN_name(2): self.get_SUN_tensors(2), get_SUN_name(3): self.get_SUN_tensors(3)}

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

    def get_SUN_tensors(self, N):
        """
        Returns coefficient tensors (i.e. epsilon and Kronecker delta) of the group SU("N") where the epsilontensor
        and the Kronecker delta are rewritten as the regarding tensors in FORM, i.e. e_(...) and d_(...,...).
        -------
        """
        epsilon = f"[su{N:d}eps]"
        delta = f"[su{N:d}dK]"
        suNeps = []
        suNdK = []
        for tensor in self.tensors:
            if tensor.name == epsilon:
                suNeps.append(f"e_({','.join([index.projection for index in tensor.indices])})")
            elif tensor.name == delta:
                suNdK.append(f"d_({','.join([index.projection for index in tensor.indices])})")

        form_tensors = "*".join(suNeps + suNdK)

        return form_tensors
