import re
from yaml import safe_load
import logging
import sys
from typing import Dict, List, Tuple
from abc import ABC, abstractmethod

from tioc import CONFIG_PATH, index_pattern

logger_autoeft = logging.getLogger("autoeft.projection")
logger = logger_autoeft.getChild(__name__)

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

        # return self.expr == other.dual_index.expr

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
    def description(self) -> str:
        """Some information about the index."""
        description_expr = index_config[self.typ]["description"]
        return description_expr

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
    def dual_index(self):
        """
        Returns dual index, i.e. for a lowered index, the upper one and for an upper the lower index. If they are the same
        it just returns the same index.
        """
        if self.typ == "Lsl":
            expr = f"Usl{self.id}"
        elif self.typ == "Usl":
            expr = f"Lsl{self.id}"
        elif self.typ == "Lsldot":
            expr = f"Usldot{self.id}"
        elif self.typ == "Usldot":
            expr = f"Lsldot{self.id}"
        else:
            expr = self.expr
        dual_ind = type(self)(expr, derIndex=self.derIndex)
        # Here, other attributes could possibly be set equal to the attributes of the non-dual index.
        try:
            dual_ind.lp = self.lp
        except AttributeError:
            pass
        return dual_ind

    def copy(self):
        new_ind = type(self)(self.expr, derIndex=self.derIndex)
        # Here, other attributes could possibly be set equal to the attributes of the non-dual index.
        try:
            new_ind.projection = self.projection
        except AttributeError:
            pass
        try:
            new_ind.lp = self.lp.copy()
        except AttributeError:
            pass
        return new_ind

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
        elif key == "projection" or key == "p":
            return f"{self.projection:s}"
        else:
            return self.__repr__()

    @property
    def d(self) -> int:
        """Dimension of the index."""
        d_expr = index_config[self.typ]["dimension"]
        return d_expr

class Dummy_Index(Index_Model):
    expr: str
    typ: str  # gauge, colf, lorentz, etc.
    id: str  # Number of the dummy index.
    indname: str  # Ensures clear contractions
    description: str
    derIndex: bool or int

    def __init__(self, number: int, derIndex=False):
        super().__init__(f"N{number:d}_?", derIndex)
        self.number = number
        self.will_be_typ = None
        self.will_be_expr = None

    def __repr__(self):
        return super().__repr__()

    @property
    def typ(self) -> str:
        """
        E.g.: lor1234, gauge345, ... have type lor, gauge, ...
        """
        if self.will_be_typ:
            return self.will_be_typ
        else:
            return "dummy"

    @property
    def id(self) -> str:
        """
        E.g.: lor1234, gauge345, ... have id 1234, 345, ...
        """
        return self.number

    @property
    def indname(self) -> str:
        """
        E.g.: lor1234, gauge345, Lsldot123, Usldot123 ... have indname: lor-1234, gauge-345, sldot-123, sldot-123 ...
        """
        return f"N-{self.id:d}"

    @property
    def description(self) -> str:
        """Some information about the index."""
        return "Dummy index of FORM."

class LP_Index():
    """
    Index for lorentz projection. This index is the possibly ambiguous attribute lp of a 'normal' index and specifies
    first of all the position of the field.
    """
    fp: int  # field position
    derIndex: int  # Number of the derivative. This can be for some time also an empty object().
    typ: str  # gauge, colf, lorentz, etc.
    indname: str  # Ensures clear contractions
    id: int  # cryptic id, e.g1234A in combination with the type
    # tex: str  # readable name which should be printed in LaTex for the index
    description: str
    # derIndex: bool or int  # False if it is not an index of a derivative. If True, number specifies the number of the derivative.

    def __init__(self, field_pos: int, derIndex: int):
        self.fp = int(field_pos)
        self.derIndex = derIndex

    def __repr__(self):
        return f"{self.fp:d}-{str(self.derIndex):s}"

    def __str__(self):
        """Specify the format for printing with str() or print() statement function: Here the same as the string representation repr() itself."""
        return self.__repr__()

    def __format__(self, key):
        """Specify the format for "format" function in print statement: Here the same as the string representation repr() itself."""
        return self.__repr__()

    def __eq__(self, other):
        if self.fp is other.fp and self.derIndex is other.derIndex:
            return True
        else:
            return False

    def copy(self):
        return LP_Index(self.fp, self.derIndex)

    @property
    def typ(self) -> str:
        return "lp"

    @property
    def id(self) -> str:
        return f"{self.fp:d}-{str(self.derIndex):s}"

    @property
    def indname(self) -> str:
        """
        E.g.: lor1234, gauge345, Lsldot123, Usldot123 ... have indname: lor-1234, gauge-345, sldot-123, sldot-123 ...
        """
        return f"lp-{self.id:s}"

    @property
    def description(self) -> str:
        """Some information about the index."""
        return "Index for lorentz projection, specifying the position of the field."