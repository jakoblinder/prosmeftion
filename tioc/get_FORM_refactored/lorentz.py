import logging
import sys
from fractions import Fraction
from pathlib import Path
from typing import List, Dict, Union, Tuple
from itertools import zip_longest

from .tableau import Young_Tableau, Lorentz_Tableau, OpClass
from .index import Index, LP_Index

logger_autoeft = logging.getLogger("autoeft.projection")
logger = logger_autoeft.getChild(__name__)

default = object()  # Default value for functions
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

    @classmethod
    def permutation_sign(cls, l: list):
        """
        Determines the permutations of sign of the given list of digits:
        For each entry, count how many smaller numbers are to the right of it. This indicates how many adjacent
        transpositions are needed to put the permutation in order.
        This is counted by 'cnt' and thus the permutation sign is given by (-1)**(cnt % 2)
        Parameters
        ----------
        l

        Returns
        -------

        """
        if len(l) <= 1:
            return 1
        assert all([isinstance(i, int) for i in l]), f"Entries have to be integers not {type(l[0])}"
        # check that all entries are different
        n_entry= [l.count(i) for i in l]
        if any([i > 1 for i in n_entry]):
            logger.warning(f"Some numbers appear at least twice in {l}.")
            return
        n = len(l)
        cnt = 0
        for i in range(n):
            for j in range(i + 1, n):
                if (l[i] > l[j]):
                    cnt += 1
        return (-1) ** (cnt % 2)

    def from_lr_tableaux(self, op_class: OpClass):
        """
        Construct lorentz tableau from left- and (conjugated) right-handed tableaux for given class.
        Note: For undotted epsilon tensors the normal order, i.e. the order of the rows in the corresponding tableaux
        is ascending w.r.t. to the first entry (and the field index in this entry) in each row. Further also from
        top to bottom in each row the field indices are ascending.
        The same is done for the dotted epsilon tensors and the corresponding right-handed tableau. If one construct
        from such an ordered right-handed tableau like for example
        -------
        | 1| 2|
        -------
        | 3| 4|
        -------
        an Lorentz tableau, this would look like
        -------
        | 2| 1|
        -------
        | 4| 3|
        -------.
        I.e. the Lorentz tableau is not SSYT only due to the normal ordering. If one changes the to columns in the
        Lorentz tableau it becomes a SSYT.
        Therefore, in the following construction of a Lorentz tableau, the ordering of the columns in the right-handed
        tableau is reversed!
        Parameters
        ----------
        op_class

        Returns
        -------
            Lorentz tableau and the overall sign.
        """
        N = op_class.N
        lorentz_tab, sign = Lorentz_Tableau([], op_class), 1
        # print(f"{self.r_tab:lp}")
        # print(f"{self.l_tab:lp}")
        for r_col in self.r_tab.iter_col():
            r_col_fp = [ele[0].lp.fp for ele in r_col.iter_row()]
            col = [j for j in range(1, N + 1) if j not in r_col_fp]
            lt = Lorentz_Tableau([[i] for i in col], op_class)
            lorentz_tab.append_col(lt)
            new_sign = LR_Tableaux.permutation_sign(col + r_col_fp)  # Note: This is not the same as r_col_fp + col! (see dim 8 Paper p. 20)
            if new_sign:
                sign *= new_sign
            else:
                # There appeared twice the same number and thus the tableau cannot be SSYT:
                return False

        lorentz_tab = lorentz_tab.reverse_cols()

        for l_col in self.l_tab.iter_col():
            l_col_fp = [ele[0].lp.fp for ele in l_col.iter_row()]
            lt = Lorentz_Tableau([[i] for i in l_col_fp], op_class)
            lorentz_tab.append_col(lt)

        # print(f"{lorentz_tab:nice}")

        return lorentz_tab, sign

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
            l_row_index_der, column = l_tab_index
            assert column == 0
            del l_tab_index
            found += 1

        r_tab_index = lr.r_tab.index_lp(derivative)
        if r_tab_index:
            r_row_index_der, column = r_tab_index
            assert column == 0
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
                    lr_cp = lr.copy()  # !
                    # Change the field position and specify the derivative 'position' uniquely
                    lr_cp.l_tab[l_row_index_der, 0].lp = LP_Index(field_pos=i, derIndex=spec_derivative)
                    lr_cp.l_tab[l_row_index_der, 0].derIndex = spec_derivative

                    lr_cp.r_tab[r_row_index_der, 0].lp = LP_Index(field_pos=i, derIndex=spec_derivative)
                    lr_cp.r_tab[r_row_index_der, 0].derIndex = spec_derivative
                    # Sign changes due to ibp
                    lr_cp *= (-1)

                    tabs.append(lr_cp)
            return tabs

    def der_23(self, der2: LP_Index, der3: LP_Index, N: int):
        """
        Replace a derivative on field 2 and a derivative on field 3 which are completely contracted, by
        a sum of derivatives on other fields.
        Parameters
        ----------
        lr = LR_tableaux(l_tab, r_tab, factor)
        der2: LP_Index
            This index specifies the derivative on field 2, which should replaced.
        der3: LP_Index
            This index specifies the derivative on field 3, which should replaced.
        N
            Number of fields in the operator.

        Returns
        -------

        """
        lr = self.copy()  # lr has to be deeply copied in order to not change the given lr.
        assert lr.l_tab.ncols() == 1
        assert lr.r_tab.ncols() == 1

        print(f"{lr:lp}")

        found = 0

        l_tab_index = lr.l_tab.index_lp(der2)
        if l_tab_index:
            l_row_index_der2, column = l_tab_index
            assert column == 0
            del l_tab_index
            del column
            found += 1

        r_tab_index = lr.r_tab.index_lp(der2)
        if r_tab_index:
            r_row_index_der2, column = r_tab_index
            assert column == 0
            del r_tab_index
            del column
            found += 1

        l_tab_index = lr.l_tab.index_lp(der3)
        if l_tab_index:
            l_row_index_der3, column = l_tab_index
            assert column == 0
            del l_tab_index
            del column
            found += 1

        r_tab_index = lr.r_tab.index_lp(der3)
        if r_tab_index:
            r_row_index_der3, column = r_tab_index
            assert column == 0
            del r_tab_index
            del column
            found += 1

        # Note: l_row_index_der2 and r_row_index_der2 denote the row in LH and RH-tableau where the derivative Index sits.

        if found != 4:
            logger.warning("No index was found for the given LP_Index. Therefore, no ibp was made.")
            return lr
        else:
            # Indices of the derivative Indices where found and the integration by parts can be done.
            tabs = []
            spec_derivative2 = object()  # specify uniquely the derivative
            spec_derivative3 = object()
            # specify an additional sign if the initial indices where in the wrong order:
            # Wished order would be:
            # -----   -----
            # |2-1|   |2-1|
            # ----- X -----
            # |3-1|   |3-1|
            # -----   -----
            # i.e.: l_row_index_der2 = 0, r_row_index_der2 = 0 & l_row_index_der3 = 1, r_row_index_der3 = 1
            sign = 1
            if l_row_index_der2 > l_row_index_der3:
                sign *= (-1)
            if r_row_index_der2 > r_row_index_der3:
                sign *= (-1)

            # Contracted derivatives on first field (-> see formula)

            lr_cp = lr.copy()  # !
            # Change the field position and specify the derivative 'position' uniquely
            lr_cp.l_tab[l_row_index_der2, 0].lp = LP_Index(field_pos=1, derIndex=spec_derivative2)
            lr_cp.l_tab[l_row_index_der2, 0].derIndex = spec_derivative2

            lr_cp.r_tab[r_row_index_der2, 0].lp = LP_Index(field_pos=1, derIndex=spec_derivative2)
            lr_cp.r_tab[r_row_index_der2, 0].derIndex = spec_derivative2

            lr_cp.l_tab[l_row_index_der3, 0].lp = LP_Index(field_pos=1, derIndex=spec_derivative3)
            lr_cp.l_tab[l_row_index_der3, 0].derIndex = spec_derivative3

            lr_cp.r_tab[r_row_index_der3, 0].lp = LP_Index(field_pos=1, derIndex=spec_derivative3)
            lr_cp.r_tab[r_row_index_der3, 0].derIndex = spec_derivative3
            # Factor 1/2 (-> see formula) and no ibp sign!
            lr_cp *= (1 / 2)
            # Possible additional sign due to wrong order in the initial indices:
            lr_cp *= sign

            tabs.append(lr_cp)

            for i in range(2, N + 1):    # iterate over all possible field positions
                for j in range(2, N + 1):
                    if (i, j) == (2,3) or (i, j) == (3,2):  # except the one which should be integrated
                        continue
                    else:
                        lr_cp = lr.copy()  # !
                        # Change the field position and specify the derivative 'position' uniquely
                        lr_cp.l_tab[l_row_index_der2, 0].lp = LP_Index(field_pos=i, derIndex=spec_derivative2)
                        lr_cp.l_tab[l_row_index_der2, 0].derIndex = spec_derivative2

                        lr_cp.r_tab[r_row_index_der2, 0].lp = LP_Index(field_pos=i, derIndex=spec_derivative2)
                        lr_cp.r_tab[r_row_index_der2, 0].derIndex = spec_derivative2

                        lr_cp.l_tab[l_row_index_der3, 0].lp = LP_Index(field_pos=j, derIndex=spec_derivative3)
                        lr_cp.l_tab[l_row_index_der3, 0].derIndex = spec_derivative3

                        lr_cp.r_tab[r_row_index_der3, 0].lp = LP_Index(field_pos=j, derIndex=spec_derivative3)
                        lr_cp.r_tab[r_row_index_der3, 0].derIndex = spec_derivative3
                        # Factor 1/2 (-> see formula)
                        lr_cp *= (1/2)
                        # Sign changes due to ibp
                        lr_cp *= (-1)
                        # Possible additional sign due to wrong order in the initial indices:
                        lr_cp *= sign

                        tabs.append(lr_cp)
            return tabs

    def is_in_fp(self, num_l: Union[int, Tuple], num_r: Union[int, Tuple]) -> Union[bool, Tuple[int]]:
        """
        Check whether the field indices num_l and num_r are in the l_tab and in the r_tab.
        Parameters
        ----------
        num_l
        num_r

        Returns
        -------

        """
        check = (self.l_tab.is_in(num_l), self.r_tab.is_in(num_r))

        if all([type(i) is int for i in check]):  # Note: all(check) doesn't work, because (0,0) would yield a False, but only (False, False) or (False, True), ... should.
            return check
        else:
            return False

    def is_in(self, i: Union[LP_Index, int]=default, j: Union[LP_Index, int]=default, k: Union[LP_Index, int]=default, l: Union[LP_Index, int]=default) -> Union[bool, Tuple[int]]:
        """
        TODO: Check
        Search for specific rows in l_tab and r_tab numbered by
        l_tab = [[i], [j]]
        r_tab = [[k], [l]]
        if specific entry is not given it is assumed that this entry is arbitrary.
        i, j, k, l are either given in terms of LP_Indices to specify that they are a derivative, or they are just an
        integer specifying the field position only.
        If they are LP_Indices the 'derIndex' can also be used to specify which of the i, j, k or l belong to which
        derivative, by giving the corresponding i, j, k or l the same 'derIndex' attribute.

        Parameters
        ----------
        i = l_tab[0,0]
        j = l_tab[1,0]
        k = r_tab[0,0]
        l = r_tab[1,0]
        =>  -----   -----
            |i  |   |k  |
            ----- X -----
            |j  |   |l  |
            -----   -----
        Returns
        -------
            False if nothing found
            or the numbers of the two rows and the LP_index of the derivative which should be ibp.
        """
        found = False
        # At least one index in each row has to be given:
        if i is default:
            assert j is not default
        elif j is default:
            assert i is not default

        if k is default:
            assert l is not default
        elif l is default:
            assert k is not default
        # If i or j, i.e. an index in the l_tab is and LP_Index, there also has to be an LP_Index in k and l
        # and vice versa. In order to have a "full" derivative.
        if type(i) is LP_Index or type(j) is LP_Index:
            assert type(k) is LP_Index or type(l) is LP_Index
        if type(k) is LP_Index or type(k) is LP_Index:
            assert type(i) is LP_Index or type(j) is LP_Index

        if all([type(m) is LP_Index for m in (i, j, k, l)]) and (i != k and j != l) and (i != l and j != k):
            logger.error(f"Four different derivative indices are meaningless in a sense.")
            sys.exit("STOP")
        if [type(m) is LP_Index for m in (i,j)].count(True) > 0:
            assert [type(m) is LP_Index for m in (k,l)].count(True) > 0


        # Find candidate rows in left tableau:
        der_i = []  # if
        der_j = []
        left = False

        if type(i) is int and type(j) is int:
            # ordinary, i.e. easy case:
            left = self.l_tab.is_in((i, j))
        elif i is default and type(j) is int:
            # ordinary, i.e. easy case:
            left = self.l_tab.is_in(j)
        elif type(i) is int and j is default:
            # ordinary, i.e. easy case:
            left = self.l_tab.is_in(i)
        elif type(i) is LP_Index and j is default:
            # find all rows where the LP_index occurs
            for m, col in enumerate(self.l_tab.iter_col()):
                if col[0,0].lp.fp == i.fp and col[0,0].lp.derIndex:
                    der_i.append((m, col[0,0].lp))
                elif col[1,0].lp.fp == i.fp and col[1,0].lp.derIndex:
                    der_i.append((m, col[1, 0].lp))
        elif type(i) is LP_Index and type(j) is int:
            for m, col in enumerate(self.l_tab.iter_col()):
                if col[0, 0].lp.fp == i.fp and col[0,0].lp.derIndex and col[1, 0].lp.fp == j:
                    der_i.append((m, col[0, 0].lp))
                elif col[1, 0].lp.fp == i.fp and col[1,0].lp.derIndex and col[0, 0].lp.fp == j:
                    der_i.append((m, col[1, 0].lp))
        elif type(i) is int and type(j) is LP_Index:
            for m, col in enumerate(self.l_tab.iter_col()):
                if col[0, 0].lp.fp == i and col[1, 0].lp.fp == j.fp and col[1,0].lp.derIndex:
                    der_j.append((m, col[1, 0].lp))
                elif col[1, 0].lp.fp == i and col[0, 0].lp.fp == j.fp and col[0,0].lp.derIndex:
                    der_j.append((m, col[0, 0].lp))
        elif i is default and type(j) is LP_Index:
            for m, col in enumerate(self.l_tab.iter_col()):
                if col[1, 0].lp.fp == j.fp and col[1,0].lp.derIndex:
                    der_j.append((m, col[1, 0].lp))
                elif col[0, 0].lp.fp == j.fp and col[0,0].lp.derIndex:
                    der_j.append((m, col[0, 0].lp))
        elif type(i) is LP_Index and type(j) is LP_Index:
            for m, col in enumerate(self.l_tab.iter_col()):
                if col[0, 0].lp.fp == i.fp and col[0,0].lp.derIndex and col[1, 0].lp.fp == j.fp and col[1,0].lp.derIndex:
                    der_i.append((m, col[0, 0].lp))
                    der_j.append((m, col[1, 0].lp))
                elif col[1, 0].lp.fp == i.fp and col[1,0].lp.derIndex and col[0, 0].lp.fp == j.fp and col[0,0].lp.derIndex:
                    der_i.append((m, col[1, 0].lp))
                    der_j.append((m, col[0, 0].lp))

        # If der_i has entries, i is an LP_Index and rows with corresponding entries where found.
        # I.e. already at this place it we know that if left == False and both der_i and der_j are empty, nothing was found
        if not left and not der_i and not der_j:
            return False

        # Find candidate rows in right tableau:
        der_k = []  # if
        der_l = []
        right = False

        if type(k) is int and type(l) is int:
            # ordinary, i.e. easy case:
            right = self.r_tab.is_in((k, l))
        elif k is default and type(l) is int:
            # ordinary, i.e. easy case:
            right = self.r_tab.is_in(l)
        elif type(k) is int and l is default:
            # ordinary, i.e. easy case:
            right = self.r_tab.is_in(k)
        elif type(k) is LP_Index and l is default:
            # find all rows where the LP_index occurs
            for m, col in enumerate(self.r_tab.iter_col()):
                if col[0, 0].lp.fp == k.fp and col[0,0].lp.derIndex:
                    der_k.append((m, col[0, 0].lp))
                elif col[1, 0].lp.fp == k.fp and col[1,0].lp.derIndex:
                    der_k.append((m, col[1, 0].lp))
        elif type(k) is LP_Index and type(l) is int:
            for m, col in enumerate(self.r_tab.iter_col()):
                if col[0, 0].lp.fp == k.fp and col[0,0].lp.derIndex and col[1, 0].lp.fp == l:
                    der_k.append((m, col[0, 0].lp))
                elif col[1, 0].lp.fp == k.fp and col[1,0].lp.derIndex and col[0, 0].lp.fp == l:
                    der_k.append((m, col[1, 0].lp))
        elif type(k) is int and type(l) is LP_Index:
            for m, col in enumerate(self.r_tab.iter_col()):
                if col[0, 0].lp.fp == k and col[1, 0].lp.fp == l.fp and col[1,0].lp.derIndex:
                    der_l.append((m, col[1, 0].lp))
                elif col[1, 0].lp.fp == k and col[0, 0].lp.fp == l.fp and col[0,0].lp.derIndex:
                    der_l.append((m, col[0, 0].lp))
        elif k is default and type(l) is LP_Index:
            for m, col in enumerate(self.r_tab.iter_col()):
                if col[1, 0].lp.fp == l.fp and col[1,0].lp.derIndex:
                    der_l.append((m, col[1, 0].lp))
                elif col[0, 0].lp.fp == l.fp and col[0,0].lp.derIndex:
                    der_l.append((m, col[0, 0].lp))
        elif type(k) is LP_Index and type(l) is LP_Index:
            for m, col in enumerate(self.r_tab.iter_col()):
                if col[0, 0].lp.fp == k.fp and col[0,0].lp.derIndex and col[1, 0].lp.fp == l.fp and col[1,0].lp.derIndex:
                    der_k.append((m, col[0, 0].lp))
                    der_l.append((m, col[1, 0].lp))
                elif col[1, 0].lp.fp == k.fp and col[1,0].lp.derIndex and col[0, 0].lp.fp == l.fp and col[0,0].lp.derIndex:
                    der_k.append((m, col[1, 0].lp))
                    der_l.append((m, col[0, 0].lp))

        if not right and not der_k and not der_l:
            return False

        # At this point candidates for both rows are found.
        # Further, the following cases are finished up to the return statement
        # i, j = default or int
        # k, l = default or int
        if left and right:
            # all non default given indices where integers
            assert not der_i and not der_j and not der_k and not der_l
            return left, right

        if der_i and der_j:
            # Both i and j are derivatives and for both candidates where found.
            # Find a possible match in der_k and der_l.
            if der_k and der_l:
                if i == k and j != l:
                    # i and k should be indices of the same derivative, but j and l not, or it is at least not demanded for them.
                    for cand_i in der_i:
                        for cand_k in der_k:
                            if cand_i[1] == cand_k[1]:
                                return cand_i[0], cand_k[0], cand_i[1]
                elif i == k and j == l:
                    # i and k should be indices of the same derivative and j and l should be indices of the same derivative
                    for cand_i, cand_j in zip(der_i, der_j):
                        for cand_k, cand_l in zip(der_k, der_l):
                            if cand_i[1] == cand_k[1] and cand_j[1] == cand_l[1]:
                                return cand_i[0], cand_k[0], cand_i[1], cand_j[1]
                elif i != k and j == l:
                    # j and l should be indices of the same derivative, but i and k not, or it is at least not demanded for them.
                    for cand_j in der_j:
                        for cand_l in der_l:
                            if cand_j[1] == cand_l[1]:
                                return cand_j[0], cand_l[0], cand_j[1]

                elif i == l and j != k:
                    # i and l should be indices of the same derivative, but j and k not, or it is at least not demanded for them.
                    for cand_i in der_i:
                        for cand_l in der_l:
                            if cand_i[1] == cand_l[1]:
                                return cand_i[0], cand_l[0], cand_i[1]
                elif i == l and j == k:
                    # i and l should be indices of the same derivative and j and k should be indices of the same derivative
                    for cand_i, cand_j in zip(der_i, der_j):
                        for cand_k, cand_l in zip(der_k, der_l):
                            if cand_i[1] == cand_l[1] and cand_j[1] == cand_k[1]:
                                return cand_i[0], cand_l[0], cand_i[1], cand_j[1]
                elif i != l and j == k:
                    # j and k should be indices of the same derivative, but i and l not, or it is at least not demanded for them.
                    for cand_j in der_j:
                        for cand_k in der_k:
                            if cand_j[1] == cand_k[1]:
                                return cand_j[0], cand_k[0], cand_j[1]
            elif der_k:
                if i == k:
                    # i and k should be indices of the same derivative, but j and l not, or it is at least not demanded for them.
                    for cand_i in der_i:
                        for cand_k in der_k:
                            if cand_i[1] == cand_k[1]:
                                return cand_i[0], cand_k[0], cand_i[1]
                elif j == k:
                    # j and k should be indices of the same derivative, but i and l not, or it is at least not demanded for them.
                    for cand_j in der_j:
                        for cand_k in der_k:
                            if cand_j[1] == cand_k[1]:
                                return cand_j[0], cand_k[0], cand_j[1]
            elif der_l:
                if i == l:
                    # i and l should be indices of the same derivative, but j and k not, or it is at least not demanded for them.
                    for cand_i in der_i:
                        for cand_l in der_l:
                            if cand_i[1] == cand_l[1]:
                                return cand_i[0], cand_l[0], cand_i[1]
                elif j == l:
                    # j and l should be indices of the same derivative, but i and k not, or it is at least not demanded for them.
                    for cand_j in der_j:
                        for cand_l in der_l:
                            if cand_j[1] == cand_l[1]:
                                return cand_j[0], cand_l[0], cand_j[1]
        elif der_i:
            # Remember the given candidate also respects that j was default or int.
            if der_k and der_l:
                if i == k:
                    # i and k should be indices of the same derivative, but j and l not, or it is at least not demanded for them.
                    for cand_i in der_i:
                        for cand_k in der_k:
                            if cand_i[1] == cand_k[1]:
                                return cand_i[0], cand_k[0], cand_i[1]
                elif i == l:
                    # i and l should be indices of the same derivative, but j and k not, or it is at least not demanded for them.
                    for cand_i in der_i:
                        for cand_l in der_l:
                            if cand_i[1] == cand_l[1]:
                                return cand_i[0], cand_l[0], cand_i[1]
            elif der_k:
                if i == k:
                    # i and k should be indices of the same derivative, but j and l not, or it is at least not demanded for them.
                    for cand_i in der_i:
                        for cand_k in der_k:
                            if cand_i[1] == cand_k[1]:
                                return cand_i[0], cand_k[0], cand_i[1]
            elif der_l:
                if i == l:
                    # i and l should be indices of the same derivative, but j and k not, or it is at least not demanded for them.
                    for cand_i in der_i:
                        for cand_l in der_l:
                            if cand_i[1] == cand_l[1]:
                                return cand_i[0], cand_l[0], cand_i[1]
        elif der_j:
            if der_k and der_l:
                if j == k:
                    # j and k should be indices of the same derivative, but i and l not, or it is at least not demanded for them.
                    for cand_j in der_j:
                        for cand_k in der_k:
                            if cand_j[1] == cand_k[1]:
                                return cand_j[0], cand_k[0], cand_j[1]
                elif j == l:
                    # j and l should be indices of the same derivative, but i and k not, or it is at least not demanded for them.
                    for cand_j in der_j:
                        for cand_l in der_l:
                            if cand_j[1] == cand_l[1]:
                                return cand_j[0], cand_l[0], cand_j[1]
            elif der_k:
                if j == k:
                    # j and k should be indices of the same derivative, but i and l not, or it is at least not demanded for them.
                    for cand_j in der_j:
                        for cand_k in der_k:
                            if cand_j[1] == cand_k[1]:
                                return cand_j[0], cand_k[0], cand_j[1]
            elif der_l:
                if j == l:
                    # j and l should be indices of the same derivative, but i and k not, or it is at least not demanded for them.
                    for cand_j in der_j:
                        for cand_l in der_l:
                            if cand_j[1] == cand_l[1]:
                                return cand_j[0], cand_l[0], cand_j[1]

        return False

    def eom(self, N: int):
        """
        Determines whether there could be an equation of motion in the term.
        This can be relatively easily done, because eoms only may occur when the SL2C indices are antisymmetric,
        i.e. if 2 indices in one column of the l_tab and r_tab have the same field index, this could yield an EOM
        and this property would return True.
        E.g.:
        -----   -----
        |1-1|   |1-1|
        ----- X -----
        |1  |   |l  |
        -----   -----

        Parameters
        ----------
        N
            Number of field in the corresponding Summand.

        Returns
        -------

        """
        l_tab = self.l_tab.copy()
        r_tab = self.r_tab.copy()

        for i in range(1, N+1):
            left  = l_tab.is_in((i, i))
            right = r_tab.is_in((i, i))
            if type(left) is not bool or type(right) is not bool:
                return True

        return False


    def remove_ibp(self, N:int):
        """

        Returns
        -------

        """
        lr = self.copy()  # lr has to be deeply copied in order to not change the given lr.
        tabs = []  # List[LR_Tableaux, bool]  -> The bool specifies whether EOMs might occur.
        def apply_ibp(lr_tab: LR_Tableaux, N:int):
            """
            Recursive method applying all ibp relations. FIXME: Should not be recursive, since after each ibp the may occur EOMs in the first field.
            Parameters
            ----------
            lr

            Returns
            -------

            """
            print(f"{lr_tab:lp}")
            der1 = object()
            der2 = object()
            # if lr := lr_tab.is_in(i=LP_Index(2, der1), j=LP_Index(3, der2), k=LP_Index(2, der1), l=LP_Index(3, der2)):
            #     l = lr[0]
            #     r = lr[1]
            #     der_LP = lr[2]
            #     extracted_lr_tab = LR_Tableaux(lr_tab.l_tab[:, l], lr_tab.r_tab[:, r], lr_tab.factor)
            #     print(f"{extracted_lr_tab:lp}")
            #     print("----------------")

            if lr := lr_tab.is_in(i=LP_Index(1,der1), k=LP_Index(1,der1)):
                # -----   -----
                # |1-1|   |1-1|
                # ----- X -----
                # |j  |   |l  |
                # -----   -----
                l = lr[0]
                r = lr[1]
                der_LP = lr[2]
                # Get the rows corresponding rows from lr_tab and do the ibp with them
                extracted_lr_tab = LR_Tableaux(lr_tab.l_tab[:,l], lr_tab.r_tab[:,r], lr_tab.factor)
                # Integrate the selected rows by parts
                tabs_tmp = extracted_lr_tab.ibp(der_LP, N)
                # Reinsert the new rows again
                for new_lr_tab in tabs_tmp:
                    lr_tab_cp = lr_tab.copy()
                    lr_tab_cp.l_tab[:, l] = new_lr_tab.l_tab[:, 0]
                    lr_tab_cp.r_tab[:, r] = new_lr_tab.r_tab[:, 0]
                    lr_tab_cp.factor = new_lr_tab.factor
                    # Check whether there could occur EOMs
                    if lr_tab_cp.eom(N):
                        # There might be an EOM
                        # Do not apply to each lr_tableaux again all possible ibp relations, because the 'ibp' method also
                        # write derivatives again on the first field. These could be EOMs, but the algorithm here
                        # doesn't know that. Thus, the EOMs need to be removed first and afterwards the algorithm
                        # can be applied again.
                        tabs.append((lr_tab_cp, True))
                    else:
                        # No EOM occur in this step. The next iteration can already be done
                        apply_ibp(lr_tab_cp, N)

            elif lr := lr_tab.is_in(i=LP_Index(2,der1), k=LP_Index(2,der1), l=1):
                # -----   -----
                # |2-1|   |2-1|
                # ----- X -----
                # |j  |   |1  |
                # -----   -----
                l = lr[0]
                r = lr[1]
                der_LP = lr[2]
                # Get the rows corresponding rows from lr_tab and do the ibp with them
                extracted_lr_tab = LR_Tableaux(lr_tab.l_tab[:, l], lr_tab.r_tab[:, r], lr_tab.factor)
                # Integrate the selected rows by parts
                tabs_tmp = extracted_lr_tab.ibp(der_LP, N)
                # Reinsert the new rows again
                for new_lr_tab in tabs_tmp:
                    lr_tab_cp = lr_tab.copy()
                    lr_tab_cp.l_tab[:, l] = new_lr_tab.l_tab[:, 0]
                    lr_tab_cp.r_tab[:, r] = new_lr_tab.r_tab[:, 0]
                    lr_tab_cp.factor = new_lr_tab.factor
                    # Check whether there could occur EOMs
                    if lr_tab_cp.eom(N):
                        # Due to possible EOMs, do not apply the ibp relations directly again -> remove EOMs first.
                        tabs.append((lr_tab_cp, True))
                    else:
                        # No EOM occur in this step. The next iteration can already be done
                        apply_ibp(lr_tab_cp, N)

            elif lr := lr_tab.is_in(i=LP_Index(2,der1), j=1, k=LP_Index(2,der1)):
                # -----   -----
                # |2-1|   |2-1|
                # ----- X -----
                # |1  |   |l  |
                # -----   -----
                l = lr[0]
                r = lr[1]
                der_LP = lr[2]
                # Get the rows corresponding rows from lr_tab and do the ibp with them
                extracted_lr_tab = LR_Tableaux(lr_tab.l_tab[:, l], lr_tab.r_tab[:, r], lr_tab.factor)
                # Integrate the selected rows by parts
                tabs_tmp = extracted_lr_tab.ibp(der_LP, N)
                # Reinsert the new rows again
                for new_lr_tab in tabs_tmp:
                    lr_tab_cp = lr_tab.copy()
                    lr_tab_cp.l_tab[:, l] = new_lr_tab.l_tab[:, 0]
                    lr_tab_cp.r_tab[:, r] = new_lr_tab.r_tab[:, 0]
                    lr_tab_cp.factor = new_lr_tab.factor
                    # Check whether there could occur EOMs
                    if lr_tab_cp.eom(N):
                        # Due to possible EOMs, do not apply the ibp relations directly again -> remove EOMs first.
                        tabs.append((lr_tab_cp, True))
                    else:
                        # No EOM occur in this step. The next iteration can already be done
                        apply_ibp(lr_tab_cp, N)

            elif lr := lr_tab.is_in(i=LP_Index(3,der1), j=2, k=LP_Index(3,der1), l=1):
                # -----   -----
                # |3-1|   |3-1|
                # ----- X -----
                # |2  |   |1  |
                # -----   -----
                l = lr[0]
                r = lr[1]
                der_LP = lr[2]
                # Get the rows corresponding rows from lr_tab and do the ibp with them
                extracted_lr_tab = LR_Tableaux(lr_tab.l_tab[:, l], lr_tab.r_tab[:, r], lr_tab.factor)
                # Integrate the selected rows by parts
                tabs_tmp = extracted_lr_tab.ibp(der_LP, N)
                # Reinsert the new rows again
                for new_lr_tab in tabs_tmp:
                    lr_tab_cp = lr_tab.copy()
                    lr_tab_cp.l_tab[:, l] = new_lr_tab.l_tab[:, 0]
                    lr_tab_cp.r_tab[:, r] = new_lr_tab.r_tab[:, 0]
                    lr_tab_cp.factor = new_lr_tab.factor
                    # Check whether there could occur EOMs
                    if lr_tab_cp.eom(N):
                        # Due to possible EOMs, do not apply the ibp relations directly again -> remove EOMs first.
                        tabs.append((lr_tab_cp, True))
                    else:
                        # No EOM occur in this step. The next iteration can already be done
                        apply_ibp(lr_tab_cp, N)

            elif lr := lr_tab.is_in(i=LP_Index(3,der1), j=1, k=LP_Index(3,der1), l=2):
                # -----   -----
                # |3-1|   |3-1|
                # ----- X -----
                # |1  |   |2  |
                # -----   -----
                l = lr[0]
                r = lr[1]
                der_LP = lr[2]
                # Get the rows corresponding rows from lr_tab and do the ibp with them
                extracted_lr_tab = LR_Tableaux(lr_tab.l_tab[:, l], lr_tab.r_tab[:, r], lr_tab.factor)
                # Integrate the selected rows by parts
                tabs_tmp = extracted_lr_tab.ibp(der_LP, N)
                # Reinsert the new rows again
                for new_lr_tab in tabs_tmp:
                    lr_tab_cp = lr_tab.copy()
                    lr_tab_cp.l_tab[:, l] = new_lr_tab.l_tab[:, 0]
                    lr_tab_cp.r_tab[:, r] = new_lr_tab.r_tab[:, 0]
                    lr_tab_cp.factor = new_lr_tab.factor
                    # Check whether there could occur EOMs
                    if lr_tab_cp.eom(N):
                        # Due to possible EOMs, do not apply the ibp relations directly again -> remove EOMs first.
                        tabs.append((lr_tab_cp, True))
                    else:
                        # No EOM occur in this step. The next iteration can already be done
                        apply_ibp(lr_tab_cp, N)

            elif lr := lr_tab.is_in(i=LP_Index(2,der1), j=LP_Index(3,der2), k=LP_Index(2,der1), l=LP_Index(3,der2)):
                # -----   -----
                # |2-1|   |2-1|
                # ----- X -----
                # |3-1|   |3-1|
                # -----   -----
                l = lr[0]
                r = lr[1]
                der_LP2 = lr[2]  # i
                der_LP3 = lr[3]  # j
                # Get the rows corresponding rows from lr_tab and do the ibp with them
                extracted_lr_tab = LR_Tableaux(lr_tab.l_tab[:, l], lr_tab.r_tab[:, r], lr_tab.factor)
                # Integrate the selected rows by parts
                tabs_tmp = extracted_lr_tab.der_23(der_LP2, der_LP3, N)
                # Reinsert the new rows again
                for new_lr_tab in tabs_tmp:
                    lr_tab_cp = lr_tab.copy()
                    lr_tab_cp.l_tab[:, l] = new_lr_tab.l_tab[:, 0]
                    lr_tab_cp.r_tab[:, r] = new_lr_tab.r_tab[:, 0]
                    lr_tab_cp.factor = new_lr_tab.factor
                    # Do not apply to each lr_tableaux again all possible ibp relation, because the 'der_23' method also
                    # write derivatives again on the first field. These are of course EOMs, but the algorithm here
                    # doesn't know that. Thus, the EOMs need to be removed first and afterwards the algorithm
                    # can be applied again.
                    # Check whether there could occur EOMs
                    if lr_tab_cp.eom(N):
                        # Due to possible EOMs, do not apply the ibp relations directly again -> remove EOMs first.
                        tabs.append((lr_tab_cp, True))
                    else:
                        # No EOM occur in this step. The next iteration can already be done
                        apply_ibp(lr_tab_cp, N)
            else:
                tabs.append((lr_tab, False))

        apply_ibp(lr, N)

        return tabs

    def remove_schouten(self):
        """
        Remove Schouten identities, i.e. Substitute in l_tab and r_tab 2 epsilon tensors with indices of the type:
        ---------     ---------   ---------
        |1  |2  |     |1  |3  |   |1  |2  |
        --------- = - --------- + ---------
        |4  |3  |     |2  |4  |   |3  |4  |
        ---------     ---------   ---------


        Returns
        -------

        """
        def schouten(tab: Young_Tableau):
            """Schouten identities for one epsilon tensor."""
            if tab.ncols() < 2:
                return [(tab.copy(), 1)]
            else:
                tabs = []
                max_col = tab.ncols()  # number of columns
                for m in range(max_col):
                    for n in range(m + 1, max_col): # n > m
                        i = tab[0, m]
                        l = tab[1, m]
                        j = tab[0, n]
                        k = tab[1, n]
                        if i.lp.fp < j.lp.fp and j.lp.fp < k.lp.fp and k.lp.fp < l.lp.fp:
                            # i < j < k < l -> Schouten identity needs to be applied
                            new_tab1 = tab.copy()
                            new_tab1[0, m] = i
                            new_tab1[1, m] = j
                            new_tab1[0, n] = k
                            new_tab1[1, n] = l
                            tabs.append((new_tab1, -1))

                            new_tab2 = tab.copy()
                            new_tab2[0, m] = i
                            new_tab2[1, m] = k
                            new_tab2[0, n] = j
                            new_tab2[1, n] = l
                            tabs.append((new_tab2, +1))
                if not tabs:
                    # No Schouten id -> just return the original element
                    return [(tab.copy(), 1)]
                else:
                    return tabs

        new_l_tabs = schouten(self.l_tab)
        new_r_tabs = schouten(self.r_tab)

        lr_tableaux = []

        for new_l_tab in new_l_tabs:
            for new_r_tab in new_r_tabs:
                coeff = self.factor * new_l_tab[1] * new_r_tab[1]
                lr_tableaux.append(LR_Tableaux(new_l_tab[0],new_r_tab[0], coeff))

        return lr_tableaux


