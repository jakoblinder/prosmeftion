import logging.config
import re
import sys
from abc import ABC, abstractmethod
from copy import copy
from typing import Dict, List, Tuple
from collections.abc import MutableMapping

from tioc import model, op_config, index_config, get_SUN_name
from .operator import Tensor, Field, Operator_Model
from .operators import Tensors
from .coefficient import Coefficient
from .indices import Indices_Summand, Indices_Operator
from .index import Dummy_Index

logger_autoeft = logging.getLogger("autoeft.projection")
logger = logger_autoeft.getChild(__name__)

class MonBasisTensor():
    def __init__(self, group: str, tensors: Tensors, coeff: Coefficient):
        self.tensors = tensors
        self.coeff = coeff
        self.group = group

    @property
    def fieldcontent(self):
        pass

    def __repr__(self):
        """Specify the format the general string representation and for printing with repr()."""
        repr_operators = "*".join(map(repr, self.tensors))
        return f"{repr_operators:s}"

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
    def tex(self):
        # TODO
        tex_expr = ""
        return tex_expr

class SymBasisTensor():
    def __init__(self, monoms: List[MonBasisTensor]):
        self.monoms = monoms

    def __repr__(self):
        """Specify the format the general string representation and for printing with repr()."""
        repr_monoms = "+".join(map(repr, self.monoms))
        return f"{repr_monoms:s}"

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
    def tex(self):
        tex_expr = "+".join([f"{monom:tex}" for monom in self.monoms])
        return tex_expr

class Model_BasisTensors(MutableMapping):
    def __init__(self, basisTensors: Tuple):
        self.basisTensors = basisTensors
    # TODO

class SymBasisTensors(Model_BasisTensors):
    def __init__(self, basisTensors: Tuple[SymBasisTensor]):
        super().__init__(basisTensors)

class MonBasisTensors(Model_BasisTensors):
    def __init__(self, basisTensors: Tuple[MonBasisTensor]):
        super().__init__(basisTensors)