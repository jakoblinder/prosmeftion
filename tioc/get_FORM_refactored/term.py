import re
import logging
import sys
from yaml import safe_load
from typing import Dict, List, Tuple
from abc import ABC, abstractmethod
# from functools import cache  # to cache properties
from tioc import cached_property

from tioc import CONFIG_PATH, opname, escape_regex, model, index_config
from tioc.get_FORM_refactored.coefficient import Coefficient
from .summand import Summand
from .indices import Indices_Term
from .index import Index

logger_autoeft = logging.getLogger("autoeft.projection")
logger = logger_autoeft.getChild(__name__)

class Term_Model(Summand):
    """Bas class for single term consisting of an overall coefficient and a product of operators."""
    terms: List[Summand]
    name: str
    # indices: Indices
    @abstractmethod
    def __init__(self, terms: List[Dict[str, str]], name: str):
        """Load in fields, tensors and the coefficient."""
        self.terms = [Summand(term["tensors"], term["fields"], term["coefficient"], name) for term in terms]
        self.name = name
        # # The following has to stay in order to generate the tex expression for the indices on operator level
        # self.indices

    @abstractmethod
    def __repr__(self):
        """Specify the format the general string representation and for printing with repr()."""
        return f"{'+'.join(map(repr, self.terms))}"

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
        indices_expr = Indices_Term([])
        for term in self.terms:
            indices_expr += term.indices

        return indices_expr

    @property
    def tex(self):
        """
        Tex expression of whole term, but before
        tex indices are generated with created generators from the range specified in config files.
        Returns
        -------
            Dictionary with name of index and corresponding tex expression.
        """
        # Instantiate an index generator for each type of index
        generator = {}
        for index in self.indices:
            indrange = Indices_Term.get_tex_range(index_config[index.typ]["tex_indices"])
            if indrange:
                generator[index.typ] = Indices_Term.infinite_Indices(indrange)
            else:
                generator[index.typ] = Indices_Term.infinite_numIndices(index_config[index.typ]["tex_indices"])

        tex_indices = {}  # Dictionary for unique tex names of indices.
        for index in self.indices:
            try:
                tex_indices[index.indname] = next(generator[index.typ])
            except StopIteration:
                logger.error("Specified index range is to small to map all occurring indices.")
                sys.exit("STOP")

        for summand in self.terms:
            for operator in summand.tensors:
                for index in operator.indices:
                    index.tex = tex_indices[index.indname]
            for operator in summand.fields:
                for index in operator.indices:
                    index.tex = tex_indices[index.indname]

        tex_expr = ""
        for summand in self.terms:
            # tex_expr += "+"
            tex_expr += f"{summand:tex}"

        return tex_expr

    # TODO: Tex expression for índices.

class Term(Term_Model):
    """
    A term from the form output can consist of many summand. The term ist therefore split in Term_s objects
    consisting of only one summand.
    """
    name: str
    terms: List[Summand]
    # indices: Indices
    def __init__(self,  terms: List[Dict[str, str]], name: str):
        super().__init__(terms, name)
        self.indices

    def __repr__(self):
        """Specify the format the general string representation and for printing with repr()."""
        return super().__repr__()

