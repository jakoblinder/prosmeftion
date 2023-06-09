import logging.config
import re
import sys
from abc import ABC, abstractmethod
from copy import copy
from typing import Dict, List, Tuple
from collections.abc import MutableMapping

from prosmeftion import op_config, index_config
from .operator import Tensor, Field, Operator_Model
from .operators import Tensors
from .coefficient import Factor
from .indices import Indices_Summand, Indices_Operator
from .index import Index

logger_autoeft = logging.getLogger("autoeft.projection")
logger = logger_autoeft.getChild(__name__)

class MonBasisTensor():
    def __init__(self, group, tableau, coeff: Factor):
        self.group = group.name
        self.tableau = tableau

        if group.N == 2:
            ind_prefix = "gauge"
        elif group.N == 3:
            ind_prefix = "colf"
        else:
            logger.error("No index prefix defined for this group.")
            sys.exit("STOP")

        tensors = [Tensor(f"[su{group.N:d}eps](" + ",".join(f"{ind_prefix}F{index[0]}I{index[1]}" for index in column) + ")")
                     for column in self.tableau.transposed()
                  ]
        for tensor in tensors:
            for index in tensor.indices:
                index.projection = index.expr
        self.tensors = Tensors(tensors)

        self.coeff = Factor(coeff)

    def __repr__(self):
        """Specify the format the general string representation and for printing with repr()."""
        repr_operators = f"{self.coeff:s}" + "*" + "*".join(map(repr, self.tensors))
        return f"{repr_operators:s}"

    def __str__(self):
        """Specify the format for printing with str() or print() statement function: Here the same as the string representation repr() itself."""
        return self.__repr__()

    def __format__(self, key):
        """Specify the format for "format" function in print statement: Here the same as the string representation repr() itself."""
        if key == "tex":
            return self.tex
        elif key == "autoeft" or key == "a":
            return f"({self.coeff:s})*" + f"{self.tensors:a}"
        elif key == "projection" or key == "p":
            return f"({self.coeff:s})*" + f"{self.tensors:p}"
        else:
            return self.__repr__()

    @property
    def tex(self):
        mult_sign = ""
        self.coeff.name = "factor"
        if self.coeff.expr == "1":
            tex_expr = f"{mult_sign} ".join([f"{tensor:tex}" for tensor in self.tensors])
        else:
            tex_expr = f"{self.coeff:tex} {mult_sign} " + f"{mult_sign} ".join([f"{tensor:tex}" for tensor in self.tensors])
        return tex_expr


class SymBasisTensor():
    monoms: MonBasisTensor  # Monomial basis tensors in the correct permutation of indices with symmetrisation factor.
    tensorIndex: Index  # Index which determines the number of the basis tensor.
    def __init__(self, monoms: List[MonBasisTensor], tensorIndex: int):
        self.monoms = monoms
        self.tensorIndex = Index(f"sbasis{tensorIndex:d}")

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
        elif key == "autoeft" or key == "a":
            return f"{'+'.join([f'{monom:a}' for monom in self.monoms])}"
        elif key == "projection" or key == "p":
            return f"{'+'.join([f'{monom:p}' for monom in self.monoms])}"
        elif key == "abbreviation" or key == "abb":
            if "2" in self.group:
                N = 2
            elif "3" in self.group:
                N = 3
            else:
                logger.error("Group not defined.")
                sys.exit("STOP")
            return f"TSU{N:d}({','.join([repr(index) for index in self.indices])})"
        else:
            return self.__repr__()

    @property
    def tex(self):
        if "2" in self.group:
            N = 2
        elif "3" in self.group:
            N = 3
        else:
            logger.error("Group not defined.")
            sys.exit("STOP")

        tex_expr = f"({op_config['tensors'][f'TSU{N:d}']['tex']})"
        sub_indices = [f"{index:tex}" for ind_typ in ["gauge", "colf"] for index in self.indices[ind_typ]]
        # TODO: Tex index for sbasis index and equalize sbasis index in indices with tensorIndex.
        tex_sbasis = index_config["sbasis"]["tex_indices"]
        # just write only the number
        super_indices = [f"{index.id}" for ind_typ in ["sbasis"] for index in self.indices[ind_typ]]
        assert super_indices and sub_indices
        tex_expr += f"^{{{', '.join(super_indices)}}}"
        tex_expr += f"_{{{' '.join(sub_indices)}}}"

        tex_expr += " &= "

        tex_expr += " ".join([f"{monom:tex}" for monom in self.monoms])
        return tex_expr

    @property
    def group(self):
        return self.monoms[0].group

    @property
    def fieldcontent(self):
        return self.monoms[0].fieldcontent

    @property
    def indices(self):
        basis_tensor_index = [self.tensorIndex]
        def sort_indices(indices: List[Index]):
            """
            Sort indices with id F1I2 in given list of indices by their order in fields, i.e. for example:
            gaugeF2I1, gaugeF1I1 -> gaugeF1I1, gaugeF2I1
            Parameters
            ----------
            indices

            Returns
            -------
            """
            indices_aux = []
            for index in indices:
                match = re.match(r"(gauge|colf)F(?P<field>\d{1,2})I(?P<index>\d{1,2})", index.projection)
                f_number = int(match.group("field"))
                i_number = int(match.group("index"))
                indices_aux.append((f_number + 0.1*i_number, index))
            indices_aux = sorted(indices_aux, key=lambda tup: tup[0])
            return [index[1] for index in indices_aux]

        tensor_indices = sort_indices([index for tensor in self.monoms[0].tensors for index in tensor.indices])

        return Indices_Operator(basis_tensor_index + tensor_indices)


class Model_BasisTensors(ABC, MutableMapping):
    """
    Contain symmetrized and unsymmetrized basis tensors of SUN.
    """
    basisTensors: Tuple[SymBasisTensor]
    @abstractmethod
    def __init__(self, basisTensors: Tuple):
        self.basisTensors = tuple(basisTensors)

    @abstractmethod
    def __repr__(self):
        """Specify the format the general string representation and for printing with repr()."""
        return f"{', '.join([repr(tensor) for tensor in self])}"

    def __str__(self):
        """Specify the format for printing with str() or print() statement function: Here the same as the string representation repr() itself."""
        return self.__repr__()

    def __format__(self, key):
        """Specify the format for "format" function in print statement: Here the same as the string representation repr() itself."""
        if key == "tex":
            return "\\\\ \n".join([f"{basisTensor:tex}" for basisTensor in self]) + "\n"
        elif key == "abbreviation" or key == "abb":
            return "\n".join([f"{basisTensor:abb}" for basisTensor in self])
        else:
            return self.__repr__()

    def __len__(self):
        """List length"""
        return len(self.basisTensors)

    def __getitem__(self, key):
        """
        If 'key' is of type slice or int we just get a list item.
        Parameters
        ----------
        key: int, slice

        Returns
        -------
        """
        if isinstance(key, int):
            return self.basisTensors[key]
        elif isinstance(key, slice):
            return type(self)(self.basisTensors[key])
        else:
            logger.error(f"Key is of type {type(key)}, but it should be of type int or slice.")
            sys.exit("STOP")

    def __setitem__(self, key, value):
        type_check1 = [isinstance(key, int), key > 0 if isinstance(key, int) else False, isinstance(key, slice)]
        if not any(type_check1):
            logger.error("Key has to be of typ int and >= 0 or of type slice.")
            sys.exit("STOP")
        type_check2 = [isinstance(value, MonBasisTensor), isinstance(value, SymBasisTensor)]
        if not any(type_check2):
            logger.error("The value which will be set has to be of type MonBasisTensor or SymBasisTensor.")
            sys.exit("STOP")
        logger.debug(f"The basis Tensor {self.basisTensors[key]:s} will be rewritten with {value:s}.")
        basisTensors = list(self.basisTensors)
        basisTensors[key] = value
        self.basisTensors = tuple(basisTensors)

    def __delitem__(self, key):
        type_check = [isinstance(key, int), key > 0 if isinstance(key, int) else False, isinstance(key, slice)]
        if not any(type_check):
            logger.error("Key has to be of typ int and >= 0 or of type slice.")
            sys.exit("STOP")
        logger.debug(f"The basis tensor(s) {self[key]:s} will be deleted.")
        basisTensors = list(self.basisTensors)
        if isinstance(key, slice):
            for i in range(key.start,key.stop,key.step):
                del basisTensors[i]
        else:
            del basisTensors[key]
        self.basisTensors = tuple(basisTensors)

    def __iter__(self):
        return iter(self.basisTensors)

    def clear(self):
        return self.basisTensors.clear()

    def copy(self):
        return type(self)(list(self.basisTensors).copy())

    @property
    def group(self):
        return self.basisTensors[0].group

    @property
    def fieldcontent(self):
        return self.basisTensors[0].fieldcontent


class SymBasisTensors(Model_BasisTensors):
    basisTensors: Tuple[SymBasisTensor]

    def __init__(self, basisTensors: Tuple[SymBasisTensor]):
        super().__init__(basisTensors)

    def __repr__(self):
        """Specify the format the general string representation and for printing with repr()."""
        return super().__repr__()


class MonBasisTensors(Model_BasisTensors):
    basisTensors: Tuple[MonBasisTensor]
    def __init__(self, basisTensors: Tuple[MonBasisTensor]):
        super().__init__(basisTensors)

    def __repr__(self):
        """Specify the format the general string representation and for printing with repr()."""
        return super().__repr__()