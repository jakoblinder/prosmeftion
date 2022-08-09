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

    l_tab: Young_Tableau  # l_tab and r_tab are normal order, i.e. the field position indices are
    r_tab: Young_Tableau  # increasing in each column from top to bottom.
    factor: Fraction
    def __init__(self, l_tab: Young_Tableau, r_tab: Young_Tableau, factor: Fraction):
        self.l_tab, l_sign = l_tab.normal_order_tableau_col()
        self.r_tab, r_sign = r_tab.normal_order_tableau_col()
        # order also the row in increasing order
        self.l_tab = self.l_tab.normal_order_tableau_row()
        self.r_tab = self.r_tab.normal_order_tableau_row()

        self.factor = factor if isinstance(factor, Fraction) else Fraction(factor).limit_denominator()
        self.factor *= l_sign * r_sign

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

    def copy(self):
        l_tab = self.l_tab.copy()
        r_tab = self.r_tab.copy()
        return LR_Tableaux(l_tab, r_tab, self.factor)

    def ibp(self, derivative: LP_Index, N: int):
        """

        Parameters
        ----------
        lr = LR_tableaux(l_tab, r_tab, factor)
        derivative: LP_Index
            This index specifies the derivative, which should be integrated by parts:
        N
            Number of fields in the operator.

        Returns
        -------

        """
        lr = self.copy()  # lr has to be deeply copied in order to not change the given lr.
        assert lr.l_tab.ncols() == 1
        assert lr.r_tab.ncols() == 1

        found = 0

        l_tab_index = lr.l_tab.index_lp(derivative)
        if l_tab_index:
            l_row_index_der, colum = l_tab_index
            assert colum == 0
            del l_tab_index
            found += 1

        r_tab_index = lr.r_tab.index_lp(derivative)
        if r_tab_index:
            r_row_index_der, colum = r_tab_index
            assert colum == 0
            del r_tab_index
            found += 1
        # Note: l_row_index_der and r_row_index_der denote the row in LH and RH-tableau where the derivative Index sits.

        if found != 2:
            logger.warning("No index was found for the given LP_Index. Therefore, no ibp was made.")
            return lr
        else:
            # Indices of the derivative Indices where found and the integration by parts can be done.
            tabs = []
            spec_derivative = object()  # specify uniquely the derivative
            for i in range(1, N + 1):  # iterate over all possible field positions
                if i == derivative.fp:  # except the one which should be integrated
                    continue
                else:
                    # Change the field position and specify the derivative 'position' uniquely
                    lr.l_tab[l_row_index_der, 0].lp = LP_Index(field_pos=i, derIndex=spec_derivative)
                    lr.l_tab[l_row_index_der, 0].derIndex = spec_derivative

                    lr.r_tab[r_row_index_der, 0].lp = LP_Index(field_pos=i, derIndex=spec_derivative)
                    lr.r_tab[r_row_index_der, 0].derIndex = spec_derivative
                    # Sign changes due to ibp
                    lr *= (-1)

                    tabs.append(lr)
            return tabs





