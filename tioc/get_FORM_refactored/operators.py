import logging.config
import re
import sys
from abc import ABC, abstractmethod
from copy import copy
from fractions import Fraction
from typing import Dict, List, Tuple
from collections.abc import MutableMapping

from tioc import model, op_config, index_config, get_SUN_name
from .operator import Tensor, Field, Operator_Model
from .coefficient import Coefficient
from .indices import Indices_Summand, Indices_Operator
from .index import Dummy_Index

logger_autoeft = logging.getLogger("autoeft.projection")
logger = logger_autoeft.getChild(__name__)

class Operators_Model(MutableMapping):
    """Bas class for the list of tensors and fields which is saved in a Summand object"""
    @abstractmethod
    def __init__(self, operators: Tuple[Operator_Model]):
        """Load in fields, tensors and the coefficient."""
        self.operators = tuple(operators)
        # self.nD
        # self.fieldcounter
        # self.fieldstructure

    @abstractmethod
    def __repr__(self):
        """Specify the format the general string representation and for printing with repr()."""
        repr_operators = "*".join(map(repr, self.operators))
        return f"{repr_operators:s}"

    def __str__(self):
        """Specify the format for printing with str() or print() statement function: Here the same as the string representation repr() itself."""
        return self.__repr__()

    def __format__(self, key):
        """Specify the format for "format" function in print statement: Here the same as the string representation repr() itself."""
        if key == "tex":
            return self.tex
        elif key == "complete" or key == "c":
            operators = "*".join(map(repr, self.operators))
            return f"{operators:s}"
        elif key == "autoeft" or key == "a":
            return "*".join([f"{op:a}" for op in self])
        elif key == "projection" or key == "p":
            return "*".join([f"{op:p}" for op in self])
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
            op_SUN_typ = {get_SUN_name(group): index_name for group, index_name in zip([2,3], ["gauge", "colf"])} ##[indextyp for indextyp in index_config.keys()]
            op_SUN_typ_check = list(op_SUN_typ.keys())
            op_other_typ = ["sl2C", "yukawa"]
            if key not in (op_SUN_typ_check + op_other_typ):
                logger.error(f"Key {key:s} is not a possible operator typ.")
                sys.exit("STOP")
            elif key in op_SUN_typ_check:
                # look for SUN typ tensors.
                operator_of_typ = []
                for operator in self.operators:
                    if operator.indices[op_SUN_typ[key]]:
                        # Indices of the requested typ exist in the operator.
                        operator_of_typ.append(operator)
            elif key == "yukawa":
                # Look for Yukawa matrices.
                operator_of_typ = []
                for operator in self.operators:
                    if operator.indices["flav"]:
                        # Indices of the typ "flav" exist in the operator, which only occur (for tensors) in the Yukawa matrices.
                        operator_of_typ.append(operator)
            elif key == "sl2C":
                # Look for sl2C-tensors.
                operator_of_typ = []
                for operator in self.operators:
                    if type(operator) == Tensor:
                        if operator.indices["sl"] or operator.indices["sldot"]:
                            # Sl2C Indices exist in the operator.
                            operator_of_typ.append(operator)
                    elif type(operator) == Field:
                        sl2C_indices = operator.indices["sl"] + operator.indices["sldot"]
                        if [index for index in sl2C_indices if not index.derIndex]:
                            # Sl2C Indices exist in the operator, which don't come from a derivative.
                            operator_of_typ.append(operator)
            return type(self)(tuple(operator_of_typ))
        else:
            logger.error(f"Key is of type {type(key)}, but it should be of type int, str or slice.")
            sys.exit("STOP")

    def __setitem__(self, key, value):
        type_check1 = [isinstance(key, int), key > 0 if isinstance(key, int) else False, isinstance(key, slice)]
        if not any(type_check1):
            logger.error("Key has to be of typ int and >= 0 or of type slice.")
            sys.exit("STOP")
        type_check2 = [isinstance(value, Tensor), isinstance(value, Field)]
        if not any(type_check2):
            logger.error("The value which will be set has to be of type Tensor or Field.")
            sys.exit("STOP")
        logger.debug(f"The operator {self.operators[key]:s} will be rewritten with {value:s}.")
        operators = list(self.operators)
        operators[key] = value
        self.operators = tuple(operators)

    def __delitem__(self, key):
        type_check = [isinstance(key, int), key > 0 if isinstance(key, int) else False, isinstance(key, slice),  isinstance(key, str)]
        if not any(type_check):
            logger.error("Key has to be of typ int and >= 0 or of type slice or of type string.")
            sys.exit("STOP")
        logger.debug(f"The operator(s) {self[key]:s} will be deleted.")
        operators = list(self.operators)
        if isinstance(key, str) or isinstance(key, slice):
            for op in self[key]:
                del operators[type(self).index(operators, op)]
        else:
            del operators[key]
        self.operators = tuple(operators)

    def __iter__(self):
        return iter(self.operators)

    @staticmethod
    def index(operators, op):
        """Returns index of searched operator op."""
        for i, ref_op in enumerate(operators):
            if op == ref_op:
                return i

    def insert(self, ii, value):
        type_check = [isinstance(value, Tensor), isinstance(value, Field)]
        if not any(type_check):
            logger.error("The value which will be inserted has to be of type Tensor or Field.")
            sys.exit("STOP")
        logger.debug(f"The operator {value:s} will be inserted at the position {ii:d}.")
        operators = list(self.operators)
        operators.insert(ii, value)
        self.operators = tuple(operators)

    def append(self, value):
        type_check = [isinstance(value, Tensor), isinstance(value, Field), isinstance(value, self)]
        if not any(type_check):
            logger.error("The value which will be appended has to be of type Index.")
            sys.exit("STOP")
        logger.debug(f"The operator(s) {self.operators[key]:s} will be appended.")
        self.insert(len(self.operators), value)

    def clear(self):
        return self.operators.clear()

    def copy(self):
        return type(self)(list(self.operators).copy())

class Tensors(Operators_Model):
    def __init__(self, operators: List[Tensor]):
        super().__init__(operators)

    def __repr__(self):
        """Specify the format the general string representation and for printing with repr()."""
        return super().__repr__()

class Fields(Operators_Model):
    def __init__(self, operators: List[Field]):
        super().__init__(operators)

    def __repr__(self):
        """Specify the format the general string representation and for printing with repr()."""
        return super().__repr__()