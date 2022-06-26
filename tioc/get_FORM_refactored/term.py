import re
import logging
import sys
import sage.all
import sage.matrix as mx
from sage.rings.rational_field import QQ
from yaml import safe_load
from typing import Dict, List, Tuple
from abc import ABC, abstractmethod

from tioc import CONFIG_PATH, opname, escape_regex, model, index_config
from tioc.get_FORM_refactored.coefficient import Coefficient
from tioc.get_FORM_refactored.basisTensors import MonBasisTensors, SymBasisTensors
from .summand import Summand
from .indices import Indices_Term, Indices_Summand
from .index import Index

logger_autoeft = logging.getLogger("autoeft.projection")
logger = logger_autoeft.getChild(__name__)

class Term_Model():  # Summand
    """Bas class for single term consisting of an overall coefficient and a product of operators."""
    terms: List[Summand]
    name: str
    # indices: Indices
    @abstractmethod
    def __init__(self, terms: List[Dict[str, List[str]]], name: str):
        """Load in fields, tensors and the coefficient."""
        self.terms = [Summand(term["tensors"], term["fields"], term["coefficient"], name) for term in terms]
        self.name = name
        # self.indices

    @abstractmethod
    def __repr__(self):
        """Specify the format the general string representation and for printing with repr()."""
        if len(self.terms) == 0:
            return "0"
        else:
            return f"{'+'.join(map(repr, self.terms))}"

    def __str__(self):
        """Specify the format for printing with str() or print() statement function: Here the same as the string representation repr() itself."""
        return self.__repr__()

    def __format__(self, key):
        """Specify the format for "format" function in print statement: Here the same as the string representation repr() itself."""
        if key == "tex":
            if len(self.terms) == 0:
                return "0"
            else:
                return self.tex
        else:
            return self.__repr__()

    def __len__(self):
        """List length"""
        return len(self.terms)

    def __getitem__(self, ii):
        """Get a list item"""
        if isinstance(ii, slice):
            return type(self)(self.terms[ii])
        else:
            return self.terms[ii]

    def __setitem__(self, key, value):
        if not isinstance(key, int) or (isinstance(key, int) and key < 0):
            logger.error("Key has to be of typ int and >= 0.")
            sys.exit("STOP")
        if not isinstance(value, Summand):
            logger.error("The value which will be set has to be of type Summand.")
            sys.exit("STOP")
        logger.error("A Summand cannot be overwritten afterwards.")
        sys.exit("STOP")
        # logger.warning(f"The Summand {self.terms[key]:s} will be rewritten with {value:s}.")
        # self.indices[key] = value

    def __delitem__(self, key):
        if not isinstance(key, int) or (isinstance(key, int) and key < 0):
            logger.error("Key has to be of typ int and >= 0.")
            sys.exit("STOP")
        logger.error("A Summand cannot be deleted.")
        sys.exit("STOP")
        # logger.warning(f"The Summand {self.terms[key]:s} will be deleted.")
        # del self.indices[key]

    def __iter__(self):
        return iter(self.terms)


    def insert(self, ii, val):
        assert type(val) == Summand, f"Inserted expression '{val}' doesn't has the type Summand."
        self.terms.insert(ii, val)

    def append(self, val):
        # insert value at the end
        self.insert(len(self.terms), val)

    @property
    def indices(self):
        """
        Extract all indices occurring in a term, i.e. all indices in the coefficient tensors which may be contracted and
        all indices in the field tensors which may not be contracted. The indices are extracted as Index object and it is
        ensured that each index occurs at most ones, before they are stored all together in an indices object which assigns
        them automatically unique LaTex indices. These indices are then reassigned to all the indices of the fields.
        Returns
        -------
        Returns a list of all occurring indices.
        """
        indices_expr = Indices_Summand([])
        for term in self.terms:
            indices_expr += term.indices

        return Indices_Term(indices_expr.indices)

    @property
    def tex(self):
        """
        Returns
        -------
            Tex expression of whole term.
        """
        tex_expr = ""
        for summand in self:  # .terms
            # tex_expr += "+"  # Sign is already include in the coefficient and thus not necessary here.
            tex_expr += f"{summand:tex}"

        return tex_expr

class Term(Term_Model):
    """
    A term from the form output can consist of many summand. The term ist therefore split in Term_s objects
    consisting of only one summand.
    """
    name: str
    terms: List[Summand]
    # indices: Indices
    def __init__(self,  terms: List[Dict[str, List[str]]], name: str):
        if all(isinstance(term, dict) for term in terms):
            super().__init__(terms, name)
        elif all(isinstance(term, Summand) for term in terms):
            self.terms = terms
            self.name = name
        self.indices

    def __repr__(self):
        """Specify the format the general string representation and for printing with repr()."""
        return super().__repr__()

class TermType(Term_Model):
    """
    Summand objects sorted by their type.
    """
    terms: List[Summand]
    name: str
    nD: int  # Number of derivatives in the term.
    field_content: Dict
    d: int  # mass dimension
    sun_projection_tensors: Dict[str,Dict[str,SymBasisTensors]]  # Sun_projection tensors
    sun_projection_matrix: Dict[str, sage.matrix.matrix_rational_dense.Matrix_rational_dense]
    def __init__(self, summand : Summand, field_content: Dict):
        self.terms = [summand]  # [Summand([tensor.expr for tensor in summand.tensors], [field.expr for field in summand.fields], summand.coeff.expr, name) for term in terms]
        self.name = str({"nD": self.nD, **field_content})
        self.field_content = field_content

        # SUN-projection:
        self.sun_projection_tensors = None
        self.sun_projection_matrix = {key: mx.constructor.matrix(QQ, 0, 0, []) for key in model.sun_groups.keys()}

    def __repr__(self):
        """Specify the format the general string representation and for printing with repr()."""
        return super().__repr__()

    @property
    def nD(self):
        # Number of derivatives in the Term
        return self.terms[0].nD

    @property
    def d(self):
        # Mass dimension of the terms
        return self.terms[0].d

