import logging
import sys
from fractions import Fraction
from pathlib import Path
from typing import List, Dict, Union

from .tableau import Young_Tableau
from .index import Index, LP_Index

logger_autoeft = logging.getLogger("autoeft.projection")
logger = logger_autoeft.getChild(__name__)

class LR_Tableaux():

    l_tab: Young_Tableau
    r_tab: Young_Tableau
    factor: Fraction
    def __init__(self, l_tab: Young_Tableau, r_tab: Young_Tableau, factor: Fraction):
        self.l_tab = l_tab
        self.r_tab = r_tab
        self.factor = factor if isinstance(factor, Fraction) else Fraction(factor).limit_denominator()

    def __repr__(self):
        return f"factor: {repr(self.factor)}; l_tab: {repr(self.l_tab)}; r_tab:\n{repr(self.r_tab)}"

    def __format__(self, key):
        if key == "lp" or key == "fp":
            if not self.l_tab:
                # empty tableau
                l_tab_output = self.l_tab.tableau.__format__("nice")
            elif isinstance(self.l_tab[0, 0], Index):
                # create copy which contains lp entries instead of indices
                if key == "lp":
                    l_tab = type(self.l_tab)([[ele.lp for ele in row] for row in self.l_tab.tableau])
                else:
                    l_tab = type(self.l_tab)([[ele.lp.fp for ele in row] for row in self.l_tab.tableau])
                l_tab_output = l_tab.tableau.__format__("nice")

            if not self.r_tab:
                # empty tableau
                r_tab_output = self.r_tab.tableau.__format__("nice")
            elif isinstance(self.r_tab[0, 0], Index):
                # create copy which contains lp entries instead of indices
                if key == "lp":
                    r_tab = type(self.r_tab)([[ele.lp for ele in row] for row in self.r_tab.tableau])
                else:
                    r_tab = type(self.r_tab)([[ele.lp.fp for ele in row] for row in self.r_tab.tableau])
                r_tab_output = r_tab.tableau.__format__("nice")

            return f"factor:\n{repr(self.factor)}\nl_tab:\n{l_tab_output:s}\nr_tab:\n{r_tab_output:s}"
        else:
            return f"factor:\n{repr(self.factor)}\nl_tab:\n{self.l_tab:{key}}\nr_tab:\n{self.r_tab:{key}}"

    def __str__(self):
        return self.__str__()

    def __mul__(self, other):
        if isinstance(other, LR_Tableaux):
            if self.l_tab == other.l_tab and self.r_tab == other.r_tab:
                factor = self.factor * other.factor
                return LR_Tableaux(self.l_tab, self.r_tab, factor)
            else:
                logger.error("l_tab and r_tab are not the same and thus cannot be multiplied.")
                sys.exit("STOP")
        elif isinstance(other, int) or isinstance(other, float) or isinstance(other, Fraction):
            factor = self.factor * Fraction(other).limit_denominator()
            return LR_Tableaux(self.l_tab, self.r_tab, factor)
        else:
            logger.error("Factor is of the wrong type.")
            sys.exit("STOP")

    def __add__(self, other):
        if isinstance(other, LR_Tableaux):
            if self.l_tab == other.l_tab and self.r_tab == other.r_tab:
                factor = self.factor + other.factor
                return LR_Tableaux(self.l_tab, self.r_tab, factor)
            else:
                logger.error("l_tab and r_tab are not the same and thus cannot be multiplied.")
                sys.exit("STOP")
        elif isinstance(other, int) or isinstance(other, float) or isinstance(other, Fraction):
            factor = self.factor + Fraction(other).limit_denominator()
            return LR_Tableaux(self.l_tab, self.r_tab, factor)
        else:
            logger.error("Factor is of the wrong type.")
            sys.exit("STOP")

    def __sub__(self, other):
        if isinstance(other, LR_Tableaux):
            other.factor *= Fraction(-1)
            return (self + other)
        elif isinstance(other, int) or isinstance(other, float) or isinstance(other, Fraction):
            summand = (-1)*other
            return (self + summand)

    def __truediv__(self, other):
        if isinstance(other, LR_Tableaux):
            other.factor = 1 / other.factor
            return (self * other)
        elif isinstance(other, int) or isinstance(other, float) or isinstance(other, Fraction):
            factor = 1 / other
            return (self * factor)




