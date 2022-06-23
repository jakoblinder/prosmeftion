import logging.config
import re
import sys
from abc import ABC, abstractmethod
from copy import copy
from fractions import Fraction
from typing import Dict, List, Tuple
from collections.abc import MutableMapping

from tioc import model, op_config, index_config, get_SUN_name
from .operator import Tensor, Field
from .coefficient import Coefficient
from .indices import Indices_Summand, Indices_Operator
from .index import Dummy_Index

logger_autoeft = logging.getLogger("autoeft.projection")
logger = logger_autoeft.getChild(__name__)

class Operators_Model(Tensor, Field, MutableMapping):
    """Bas class for the list of tensors and fields which is saved in a Summand object"""
    @abstractmethod
    def __init__(self, operators: List[Tensor, Field]):
        """Load in fields, tensors and the coefficient."""
        self.operators = operators
        # self.nD
        # self.fieldcounter
        # self.fieldstructure

    @abstractmethod
    def __repr__(self):
        """Specify the format the general string representation and for printing with repr()."""
        repr_operators = "*".join(map(str, self.operators))
        return f"{repr_operators:s}"

    def __str__(self):
        """Specify the format for printing with str() or print() statement function: Here the same as the string representation repr() itself."""
        return self.__repr__()

    def __format__(self, key):
        """Specify the format for "format" function in print statement: Here the same as the string representation repr() itself."""
        if key == "tex":
            return self.tex
        elif key == "complete" or key == "c":
            texed_operators = "*".join(map(str, self.operators))
            return f"{texed_operators:s}"
        else:
            return self.__repr__()

     @property
     def tex(self):
        mult_sign = ''  # '*'
        tex_expr = mult_sign.join([f"{op:tex}" for op in self.operators])
        return tex_expr

    def __len__(self):
        """List length"""
        return len(self.operators)

    def __getitem__(self, key):
        """
        If 'key' is of type slice or int we just get a list item. If 'key' is of type str, all indices or the type indicated
        by 'key' are returned in their occurring order.
        Parameters
        ----------
        key: int, slice, str

        Returns
        -------
        """
        if isinstance(key, int):
            return self.operators[key]
        elif isinstance(key, slice):
            return type(self)(self.operators[key])
        elif isinstance(key, str):
            op_types = [get_SUN_name(i) for i in [2,3]] ##[indextyp for indextyp in index_config.keys()]
            if key not in (op_types):
                logger.error(f"Key {key:s} is not a possible operator typ.")
                sys.exit("STOP")
            else:
                operator_of_typ = []
                for operator in self.operators:
                    # TODO
                    if operator.typ == key:
                        operator_of_typ.append(index)
            return type(self)(tuple(operator_of_typ))
        else:
            logger.error(f"Key is of type {type(key)}, but it should be of type int, str or slice.")
            sys.exit("STOP")

    def __setitem__(self, key, value):
        if not isinstance(key, int) or (isinstance(key, int) and key < 0):
            logger.error("Key has to be of typ int and >= 0.")
            sys.exit("STOP")
        if not isinstance(value, Index):
            logger.error("The value which will be set has to be of type Index.")
            sys.exit("STOP")
        logger.debug(f"The index {self.indices[key]:s} will be rewritten with {value:s}.")
        indices = list(self.indices)
        indices[key] = value
        self.indices = tuple(indices)

    def __delitem__(self, key):
        if not isinstance(key, int) or (isinstance(key, int) and key < 0):
            logger.error("Key has to be of typ int and >= 0.")
            sys.exit("STOP")
        logger.warning(f"The index {self.indices[key]:s} will be deleted.")
        indices = list(self.indices)
        del indices[key]
        self.indices = tuple(indices)

    def __iter__(self):
        return iter(self.indices)

    def insert(self, ii, val):
        if not isinstance(value, Index):
            logger.error("The value which will be inserted has to be of type Index.")
            sys.exit("STOP")
        logger.warning(f"The index {self.indices[key]:s} will be inserted.")
        indices = list(self.indices)
        indices.insert(ii, val)
        self.indices = tuple(indices)

    def append(self, val):
        if not isinstance(value, Index):
            logger.error("The value which will be appended has to be of type Index.")
            sys.exit("STOP")
        logger.warning(f"The index {self.indices[key]:s} will be appended.")
        self.insert(len(self.indices), val)

    def clear(self):
        return self.indices.clear()

    def copy(self):
        return self.indices.copy()

class Tensors(Operators_Model):
    @abstractmethod
    def __init__(self, operators: List[Tensor]):
        self.operators = operators

class Fields(Operators_Model):
    @abstractmethod
    def __init__(self, operators: List[Field]):
        self.operators = operators