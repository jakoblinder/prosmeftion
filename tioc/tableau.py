import sys
from collections.abc import MutableSequence
from typing import Dict, List, Tuple
from .sun_projection import get_basis
from autoeft.combinat import Tableau
from autoeft.invariants import OpClass, LorentzTableau

class Young_List(MutableSequence):
    def __init__(self, list_input):
        self.list = list_input

    def __repr__(self):
        return repr(self.list)

    def __format__(self, key):
        if key == "nice":
            rows = []
            for row in self.list:
                rows.append(f"|{'|'.join([f'{repr(entry):>2s}' for entry in row])}|")
            output = "-"*len(rows[0]) + "\n"
            for row in rows:
                output += row + "\n"
                output += "-"*len(row) + "\n"
            return output
        else:
            return repr(self)

    @property
    def list(self):
        return self._list

    @list.setter
    def list(self, value: List[List]):
        try:
            self.diagram = [len(row) for row in value]
        except TypeError:
            print("The input list has to be of typ List[List].")
            sys.exit("STOP")
        if not sorted(self.diagram)[::-1] == self.diagram:
            print("The number of columns has to not to be increase when going from the top row to the bottom row.")
            sys.exit("STOP")
        self._list = [list(row) for row in value]

    @property
    def nrows(self):
        """Number of rows."""
        return len(self.list)

    def __getitem__(self, ii):
        max_dim = 2
        if isinstance(ii, int):
            return self.list[ii]
        elif isinstance(ii, slice):
            return type(self)(self.list[ii])
        elif isinstance(ii, tuple) and all([isinstance(i, int) for i in ii]):
            assert len(ii) == max_dim, "To many keys given."
            return self.list[ii[0]][ii[1]]
        elif isinstance(ii, tuple):
            # If any of the keys is an integer, rewrite them as a slice returning exactly the one item the number would have returned.
            ii = list(ii)
            for j in range(max_dim):
                if isinstance(ii[j], int):
                    ii[j] = slice(ii[j], ii[j] + 1, 1)
                elif isinstance(ii[j], slice):
                    if ii[j].start == None:
                        ii[j] = slice(0, ii[j].stop, ii[j].step)
                    if ii[j].step == None:
                        ii[j] = slice(ii[j].start, ii[j].stop, 1)
                    if j == 0 and ii[0].stop == None:
                        ii[0] = slice(ii[0].start, len(self.list), ii[0].step)
                    # ii[1].stop is possibly different for each list
            ii = tuple(ii)
            new_list = []
            # row:
            for i in range(ii[0].start, ii[0].stop, ii[0].step):
                start = ii[1].start
                if ii[1].stop == None:
                    stop = len(self.list[i])
                else:
                    stop = ii[1].stop
                step = ii[1].step

                row = self.list[i][start:stop:step]
                # Note that for slices there never is an index error, but it is rather returned just an empty list.
                if row:
                    # row not empty
                    new_list.append(row)

            return type(self)(new_list)

    def __delitem__(self, ii):
        del self.list[ii]

    def __setitem__(self, ii, value):
        if isinstance(ii, int):
            # substitution of one row
            assert isinstance(value, list) or isinstance(value, tuple)
            assert len(self.list[ii]) == value, "Rows should have the same lengths."
            value = list(value)
            self.list[ii] = value
        elif isinstance(ii, slice):
            # substitution of multiple rows -> inserted diagram also has to be a diagram.
            assert isinstance(value, type(self))
            assert value.diagram == self.diagram[ii], "Diagram should have the same shape afterwards."
            self.list[ii] = value
            type(self)(self.list)
        elif isinstance(ii, tuple):
            assert isinstance(value, type(self))
            max_dim = 2
            assert len(ii) == max_dim, "To many keys given."
            # If any of the keys is an integer, rewrite them as a slice returning exactly the one item the number would have returned.
            ii = list(ii)
            for j in range(max_dim):
                if isinstance(ii[j], int):
                    ii[j] = slice(ii[j], ii[j] + 1, 1)
                elif isinstance(ii[j], slice):
                    if ii[j].start == None:
                        ii[j] = slice(0, ii[j].stop, ii[j].step)
                    if ii[j].step == None:
                        ii[j] = slice(ii[j].start, ii[j].stop, 1)
                    if j == 0 and ii[0].stop == None:
                        ii[0] = slice(ii[0].start, len(self.list), ii[0].step)
                    # ii[1].stop is possibly different for each list
            ii = tuple(ii)
            # Makes sure that shapes match:
            if not self[ii].diagram == value.diagram:
                print(f"Chosen range {str(self[ii].diagram)} doesn't match th substitute change {str(value.diagram)}.")
                sys.exit("STOP")

            new_list = []
            # row:
            for i_val, i_old in enumerate(range(ii[0].start, ii[0].stop, ii[0].step)):
                if i_val >= value.nrows:
                    break
                start = ii[1].start
                if ii[1].stop == None:
                    stop = len(self.list[i_old])
                else:
                    stop = ii[1].stop
                step = ii[1].step

                self.list[i_old][start:stop:step] = value[i_val]

    def __len__(self):
        return len(self.list)

    def insert(self, ii, value) -> None:
        self.list.insert(ii, value)

# Examples:
# b = [[1,2,4,5],[3,5,4],[2,3],[3,4]]
# t = Young_List(b)
# print(f"{t:nice}")
# # print(f"First row:\n{t[0,:]:nice}")
# # print(f"First column:\n{t[:,0]:nice}")
# t[3,:] = Young_List([[1,3]])
# print(f"First column:\n{t:nice}")
# print(f"{t[0,0]}")

class Young_Tableau(Tableau):
    def __init__(self, tableau):
        self.tableau = Young_List(tableau)

    def __repr__(self):
        return self.tableau.__repr__()

    def __format__(self, key):
        return self.tableau.__format__(key)

    def __str__(self):
        return super().__str__()

# Examples:
# b = [[1,2,4,5],[3,1,4],[2,3],[3,4]]
# yt = Young_Tableau(b)
# print(repr(yt))
# print(str(yt))
# print(f"{yt:nice}")

# print(yt.contains(4))
# print(yt.contains(6))

# print(yt.get_number(1,1))
# print(yt.get_number(1,2))
# for i in yt.replace_entry(1,2):
# FIXME: Wofür wird diese Methode gebraucht? Warum werden die indices nur iterative getauscht?
#     print(type(i))
#     print(f"{i:nice}")
# print(f"{yt:nice}")

class Lorentz_Tableau(LorentzTableau):
    """
    Represent a Lorentz tensor as tableau.

    The tableau is filled with integer numbers from 1 to N,
    where N is the total number of fields in an operator.
    The entry i refers to the SL(2,C) indices of the i-th field.
    The position of the index on the field is not relevant,
    since all building blocks are assumed to be totally symmetric
    in their SL(2,C) indices.
    Note l_tabs are the undotted epsilon tensors which are contracted with the left-handed fields, whereas
         r_tabs denote the dotted epsilon tensors, contracted with right-handed fields.
    """

    tableau: Tuple[Tuple[int, ...], ...]
    op_class: OpClass
    def __init__(self, tableau: Tuple[Tuple[int, ...], ...], op_class: OpClass):
        self.tableau = Young_List(tableau)
        self.op_class = op_class

    def __repr__(self):
        return self.tableau.__repr__()

    def __format__(self, key):
        return self.tableau.__format__(key)

    def __str__(self):
        return super().__str__()

    @property
    def lr_tableaux(self) -> Tuple[Tableau, Tableau, int]:
        """Return the left- and the (conjugated) right-handed tableaux as well as the overall sign."""

        l_tab, r_tab, sign = LorentzTableau.lr_tableaux.fget(self)

        return Young_Tableau(l_tab), Young_Tableau(r_tab), sign


def test(basis, max_dim):
    # Examples:
    c = [[1,2], [3,4]]
    # Properties of the operator
    mass_dim = 6
    field_content = {"H":2, "H+":2}
    derivatives = 2
    # get all the invariants associated with the operator
    op_class, op_subclass, op_type, operator = basis[mass_dim].get_operator(field_content, derivatives)

    lt = Lorentz_Tableau(c, op_class)
    print(f"{lt:nice}")

    l_tab, r_tab, sign = lt.lr_tableaux
    print(f"LH:\n{l_tab:nice}")
    print(f"RH:\n{r_tab:nice}")
    print(f"sign: {sign:d}")
    lorentzt, sign2 = Lorentz_Tableau.from_lr_tableaux(l_tab, r_tab, op_class)
    print(f"Lorentz tableau:\n{lorentzt:nice}")
    print(f"sign: {sign2:d}")

    # FIXME: Remove FrozenInstanceError:
    #  University/Physik/9Physik/ws21/Masterthesis/autoeft/projection/tioc/tableau.py", line 208, in __init__
    #     self.tableau = Young_List(tableau)
    #   File "<string>", line 4, in __setattr__
    # dataclasses.FrozenInstanceError: cannot assign to field 'tableau'

    print("TEST")