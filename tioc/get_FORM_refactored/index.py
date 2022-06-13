import re
from yaml import safe_load
import logging
import sys
from typing import Dict, List, Tuple
from abc import ABC, abstractmethod

from tioc import CONFIG_PATH, index_pattern

logger_autoeft = logging.getLogger("autoeft.projection")
logger = logger_autoeft.getChild(__name__)

#TODO: Index range could be read from the modelfile
with open(CONFIG_PATH / "index.yml", "r") as file:
    index_config = safe_load(file)

class Index_Model(ABC):
    """
    E.g.: lor1234, gauge345, ...
    """
    expr: str
    typ: str  # gauge, colf, lorentz, etc.
    indname: str  # Ensures clear contractions
    id: str  # cryptic id, e.g1234A in combination with the type
    tex: str  # readable name which should be printed in LaTex for the index
    d: int  # dimension of the index
    description: str
    derIndex: bool or int #  False if it is not an index of a derivative. If True, number specifies the number of the derivative.

    @abstractmethod
    def __init__(self, expr, derIndex=False):
        self.expr = expr
        self.derIndex = derIndex

    @abstractmethod
    def __repr__(self):
        return f"{self.expr:s}"

    def __str__(self):
        """Specify the format for printing with str() or print() statement function: Here the same as the string representation repr() itself."""
        return self.__repr__()

    def __format__(self, key):
        """Specify the format for "format" function in print statement: Here the same as the string representation repr() itself."""
        if key == "tex":
            return f"{self.tex:s}"
        elif key == "s":
            return f"{self.expr:s}"
        else:
            return self.__repr__()

    def __eq__(self, other):
        """
        Note: Only the non unique "name" is compared, which is for Sl2C Indices Lsl, Usl, Lsldot, Usldot written like:
        sl-1234 or sldot-1234.
        Since in order to check that indices are contracted an upper and lower SL2C index has to be treated the same.
        """
        return self.expr == other.expr
        # return self.indname == other.indname

    @property
    def typ(self) -> str:
        """
        E.g.: lor1234, gauge345, ... have type lor, gauge, ...
        """
        match = re.match(index_pattern, self.expr)
        if match:
            typ_expr = match.group("typ")
        else:
            logger.error(
                f"The index {self.expr} is not written in a valid index format. An index should look like e.g.: flav67Av89.")
            sys.exit("STOP")
        return typ_expr

    @property
    def id(self) -> str:
        """
        E.g.: lor1234, gauge345, ... have id 1234, 345, ...
        """
        match = re.match(index_pattern, self.expr)
        if match:
            id_expr = match.group("id")
        else:
            logger.error(
                f"The index {self.expr} is not written in a valid index format. An index should look like e.g.: "
                f"flav67Av89.")
            sys.exit("STOP")
        return id_expr

    @property
    def indname(self) -> str:
        """
        E.g.: lor1234, gauge345, Lsldot123, Usldot123 ... have indname: lor-1234, gauge-345, sldot-123, sldot-123 ...
        """
        match = re.match(index_pattern, self.expr)
        if match:
            typ_expr = match.group("typ")
            num = match.group("id")
            if typ_expr == "Lsl" or typ_expr == "Usl":
                indname_expr = f"sl-{num}"
            elif typ_expr == "Lsldot" or typ_expr == "Usldot":
                indname_expr = f"sldot-{num}"
            else:
                indname_expr = f"{typ_expr}-{num}"
        else:
            logger.error(
                f"The index {self.expr} is not written in a valid index format. An index should look like e.g.: flav67Av89.")
            sys.exit("STOP")
        return indname_expr

    @property
    def d(self) -> int:
        """Dimension of the index."""
        d_expr = index_config[self.typ]["dimension"]
        return d_expr

    @property
    def description(self) -> str:
        """Some information about the index."""
        description_expr = index_config[self.typ]["description"]
        return description_expr

    def is_in(self,listofIndices):
        """
        Checks if index(self) exists in the given list of Indices.
        Parameters
        ----------
        listofIndices: List[Index]
            given list of Indices.
        Returns
        -------
            #index occurs in list - otherwise False.
        """
        exist = 0
        for i in listofIndices:
            if self == i:
                # Equal index exists in list
                exist += 1
        if exist:
            return exist
        else:
            return False

class Index(Index_Model):
    name: str
    typ: str  # gauge, colf, lorentz, etc.
    indname: str  # Ensures clear contractions
    id: str  # cryptic id, e.g1234A in combination with the type
    tex: str  # readable name which should be printed in LaTex for the index
    d: int
    description: str
    projection: str  # Index for projection of gauge indices.

    def __init__(self, expr, derIndex=False):
        super().__init__(expr, derIndex)
        self.indname
        self.typ
        self.d
        self.description

    def __repr__(self):
        return super().__repr__()

    def __format__(self, key):
        """Specify the format for "format" function in print statement: Here the same as the print statement itself."""
        if key == "tex":
            return f"{self.tex}"
        if key == "debug":
            return f"{self.tex}({self.name})"
        else:
            return self.__repr__()