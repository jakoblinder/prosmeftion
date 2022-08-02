import sys
from collections.abc import MutableSequence
from dataclasses import dataclass
from itertools import zip_longest
from typing import Iterator, List, Tuple, Dict
from sage.combinat.permutation import Permutation
from sage.combinat.skew_tableau import SkewTableau

from tioc import n_der, op_config, model

# from autoeft.combinat import Tableau
# from autoeft.invariants import LorentzTableau
from autoeft.combinat import Partition
# from autoeft.invariants import OpClass


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
        wrong_type = False
        if not (isinstance(value,list) or isinstance(value,tuple)):
            wrong_type = True
        if not all([(isinstance(row,list) or isinstance(row,tuple)) for row in value]):
            wrong_type = True
        if wrong_type:
            print("The input list has to be of typ List[List].")
            sys.exit("STOP")

        self._list = [list(row) for row in value]
        # check that the tableau is a Young tableau:
        self.assert_yt()

    @property
    def diagram(self):
        """Shape of the tableau."""
        return [len(row) for row in self.list]

    def nrows(self, col: int=0):
        """
        Number of rows in general or number of row in the i-th colum.
        Note: The number of rows in general is equal to the number of rows in the first column
        """
        if col == 0:
            return len(self.list)
        else:
            n = 0
            for ncols in self.diagram:
                if ncols >= col + 1:
                    n += 1
            return n
    def assert_yt(self):
        """
        Assert that the tableau saved in list is a Young tableau, i.e. the number of columns is not increasing
        when going from the top row to the bottom row.
        Returns
        -------

        """
        if not sorted(self.diagram)[::-1] == self.diagram:
            print("The number of columns, i.e. boxes in one row does not need to be increased when going from the top to the bottom row in order to be a Young tableau.")
            sys.exit("STOP")
        return True
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
            # row:
            for i in range(ii[0].start, ii[0].stop, ii[0].step):
                start = ii[1].start
                if ii[1].stop == None:
                    stop = len(self.list[i])
                else:
                    stop = ii[1].stop
                step = ii[1].step

                del self.list[i][start:stop:step]

            # Delete possibly empty lists, i.e. rows.
            for i, row in enumerate(self.list):
                if not row:
                    del self.list[i]

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
            assert isinstance(value, type(self)), f"The inserted value should be of type {type(self)} and not {type(value)}."
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
                if i_val >= value.nrows():
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
        # Since effectively a row is inserted, the 'value' has to be a row.
        assert len(value.diagram) == 1
        self.list.insert(ii, value[0])
        # check that the tableau is still a Young tableau:
        self.assert_yt()
    def insert_row(self, ii, value) -> None:
        """Insert a row."""
        self.insert(ii, value)

    def append_row(self, value) -> None:
        """Append a row."""
        self.list.insert_row(len(self.list), value)

    def insert_col(self, ii, value) -> None:
        """Insert a column."""
        # Can only insert one column in one step:
        assert isinstance(ii, int), "Can only insert one column in one step."
        # If the index is smaller than 0 it will be inserted before the first column of the old diagram
        if ii < 0:
            before_diagram = True
        else:
            before_diagram = False
        # If the index is larger than the diagram has columns, the new column will be inserted after the old diagram:
        if ii >= max(self.diagram):
            after_diagram = True
        else:
            after_diagram = False

        # Ensure the correct shape -> self.diagram = [1,1,1,...]
        assert value.nrows() == sum(value.diagram)
        # Ensure that also after inserting the column, the YT has the correct shape:
        nrows_inserted = value.nrows()
        if before_diagram:
            assert nrows_inserted >= self.nrows(0), "The number of boxes in the column which is inserted before the other columns has to be larger (or equal) then the number of boxes in the old first column."
            insert_at = 0
        elif after_diagram:
            assert nrows_inserted <= self.nrows(self.diagram[0] - 1), "The number of boxes in the column which is inserted after the other columns has to be smaller (or equal) then the number of boxes in the old last column."
            insert_at = self.diagram[0]
        else:
            assert nrows_inserted >= self.nrows(ii + 1), "The number of boxes in the column which is inserted has to be larger (or equal) then the number of boxes in the following column."
            assert nrows_inserted <= self.nrows(ii - 1), "The number of boxes in the column which is inserted has to be smaller (or equal) then the number of boxes in the preceding column."
            insert_at = ii

        # Insert one empty column at the correct place.
        for i, row in enumerate(self.list):
            if i >= nrows_inserted: break
            row.insert(insert_at, None)

        self[:nrows_inserted,insert_at] = value

        # check that the tableau is still a Young tableau:
        self.assert_yt()

    def append_col(self, value) -> None:
        self.insert_col(self.diagram[0], value)

# Examples:
# b = [[1,2,4,5],[3,5,4],[2,3],[3,4]]
# t = Young_List(b)
# print(f"{t:nice}")
# # print(f"First row:\n{t[0,:]:nice}")
# # print(f"First column:\n{t[:,0]:nice}")
# t[3,:] = Young_List([[1,3]])
# print(f"First column:\n{t:nice}")
# print(f"{t[0,0]}")
# Number of rows in each column
# print(f"{t:nice}")
# print(t.nrows)
# print(t.n_rows(0))
# print(t.n_rows(2))
# print(t.n_rows(3))
# print(t.n_rows(4))

# Copy Tableau class to unfreeze it
@dataclass(order=True, frozen=False)
class Tableau:
    """Represent a tableau."""

    tableau: Tuple[Tuple[int, ...], ...]

    def __str__(self):
        return (
            "T["
            + ",".join(
                "[" + ",".join(str(num) for num in row) + "]" for row in self.tableau
            )
            + "]"
        )

    def __len__(self):
        return len(self.tableau)

    def __iter__(self):
        return iter(self.tableau)

    @property
    def shape(self) -> Partition:
        """Return partition corresponding to the tableau."""
        return Partition(tuple(len(row) for row in self.tableau))

    def contains(self, num: int) -> bool:
        """Return whether the tableau contains the given number."""
        return any(num in row for row in self.tableau)

    def get_number(self, row_num: int, col_num: int) -> int:
        """Return number in given row and column, starting from 0."""
        return self.tableau[row_num][col_num]

    def replace_entry(self, num: int, replacement: int) -> Iterator["Tableau"]:
        """Replace, one by one, each occurrence of the given number by the replacement number."""
        for i, row in enumerate(self.tableau):
            for j, entry in enumerate(row):
                if entry == num:
                    replaced_tab = [list(replaced_row) for replaced_row in self.tableau]
                    replaced_tab[i][j] = replacement
                    yield self.from_list(replaced_tab)

    @classmethod
    def from_list(cls, tableau: List[List[int]]) -> "Tableau":
        """Construct from list."""
        return cls(tuple(tuple(row) for row in tableau))

    @classmethod
    def from_sage_tableau(cls, sage_tableau: List[list]) -> "Tableau":
        """Construct from Sage tableau."""
        return cls(
            tuple(tuple(int(num) if num else 0 for num in row) for row in sage_tableau)
        )

class Young_Tableau(Tableau):
    def __init__(self, tableau):
        self.tableau = Young_List(tableau)

    def __repr__(self):
        return self.tableau.__repr__()

    def __format__(self, key):
        return self.tableau.__format__(key)

    def __str__(self):
        return super().__str__()

    def __getitem__(self, ii):
        self.tableau[ii]

    def __delitem__(self, ii):
        del self.tableau[ii]

    def __setitem__(self, ii, value):
        self.tableau[ii] = value.tableau

    def __len__(self):
        return len(self.tableau)

    def insert_row(self, ii, value) -> None:
        self.tableau.insert_row(ii, value.tableau)

    def append_row(self, value) -> None:
        self.insert_row(len(self.tableau), value)

    def insert_col(self, ii, value) -> None:
        self.tableau.insert_col(ii, value.tableau)

    def append_col(self, value) -> None:
        self.insert_col(self.tableau.diagram[0], value)

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

# Copy OpClass class to unfreeze it AND rearrange some stuff.
@dataclass(order=True, frozen=False)
class OpClass:
    """Represent an operator by the total number of fields N, the number of SU(2) left indices 2nl, and the number of SU(2) right indices 2nr."""

    N: int
    nl: int
    nr: int

    def __str__(self):
        N, nl, nr = self.N, self.nl, self.nr
        return f"(N={N},nl={nl},nr={nr})"

    @property
    def dir_name(self) -> str:
        return str(self.N)

    @property
    def primary_partition(self) -> Partition:
        """
        Return the primary partition.

        The primary partition is defined by the Young diagram with
        nr columns of length N - 2 and nl columns of length 2.
        """
        return Partition.from_list([self.nl + self.nr] * 2 + [self.nr] * (self.N - 4))

# Copy LorentzTableau class to unfreeze it
@dataclass(order=True, frozen=False)
class LorentzTableau:
    """
    Represent a Lorentz tensor as tableau.

    The tableau is filled with integer numbers from 1 to N,
    where N is the total number of fields in an operator.
    The entry i refers to the SL(2,C) indices of the i-th field.
    The position of the index on the field is not relevant,
    since all building blocks are assumed to be totally symmetric
    in their SL(2,C) indices.
    """

    tableau: Tuple[Tuple[int, ...], ...]
    op_class: OpClass

    def __str__(self):
        return (
            "LT["
            + ",".join(
                "[" + ",".join(str(i) for i in row) + "]" for row in self.tableau
            )
            + "]"
        )

    def __iter__(self):
        return iter(self.tableau)

    @property
    def shape(self) -> Partition:
        """Return partition of the tableau."""
        return Partition(tuple(len(row) for row in self.tableau))

    @property
    def is_ssyt(self) -> bool:
        """Return whether this tableau is semi-standard."""
        return SkewTableau(self.tableau).is_semistandard()

    @property
    def lr_tableaux(self) -> Tuple[Tableau, Tableau, int]:
        """Return the left- and the (conjugated) right-handed tableaux as well as the overall sign."""
        N, nl, nr = self.op_class.N, self.op_class.nl, self.op_class.nr
        l_tab, r_tab, sign = [], [], 1
        for i, col in enumerate(zip_longest(*self.tableau), start=1):
            if i <= nr:

                """v3.8
                r_tab.append(r_col := [j for j in range(1, N + 1) if j not in col])
                """
                r_col = [j for j in range(1, N + 1) if j not in col]
                r_tab.append(r_col)

                sign *= Permutation(list(col) + r_col).sign()
            else:
                l_tab.append([j for j in col if j])
        return Tableau(tuple(zip(*l_tab))), Tableau(tuple(zip(*reversed(r_tab)))), sign  # FIXME: Why is the r_tab reversed? (It changes nothing, but shouldn't it be unnecessary?)

    @classmethod
    def from_lr_tableaux(cls, l_tab, r_tab, op_class: OpClass):
        """Construct from left- and (conjugated) right-handed tableaux for given class and return the overall sign."""
        N = op_class.N
        tab, sign = [], 1
        for r_col in reversed(list(zip(*r_tab))):  # FIXME: Reversing necessary due to reverse statement of r_tab above.
            """v3.8
            tab.append(col := [j for j in range(1, N + 1) if j not in r_col])
            """
            col = [j for j in range(1, N + 1) if j not in r_col]
            tab.append(col)

            sign *= Permutation(col + list(r_col)).sign()
        for l_col in zip(*l_tab):
            tab.append(l_col)
        lorentz_tab = cls(
            tuple(tuple(num for num in row if num) for row in zip_longest(*tab)),
            op_class,
        )
        return lorentz_tab, sign

    @classmethod
    def from_sage_tableau(
        cls, sage_tableau: List[list], op_class: OpClass
    ) -> "LorentzTableau":
        """Construct from Sage tableau."""
        lorentz_tab = cls(
            tuple(tuple(int(num) for num in row) for row in sage_tableau), op_class
        )
        return lorentz_tab

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
        # self.lt = LorentzTableau(Young_List(tableau), op_class)
        self.tableau = Young_List(tableau)
        self.op_class = op_class
        assert len(self) == 2  # each epsilon tensor only has 2 indices

    def __repr__(self):
        return self.tableau.__repr__()

    def __format__(self, key):
        return self.tableau.__format__(key)

    def __str__(self):
        return super().__str__()

    def __getitem__(self, ii):
        self.tableau[ii]

    def __delitem__(self, ii):
        del self.tableau[ii]

    def __setitem__(self, ii, value):
        self.tableau[ii] = value.tableau

    def __len__(self):
        return len(self.tableau)

    def insert_row(self, ii, value) -> None:
        self.tableau.insert_row(ii, value.tableau)

    def append_row(self, value) -> None:
        self.insert_row(len(self.tableau), value)

    def insert_col(self, ii, value) -> None:
        self.tableau.insert_col(ii, value.tableau)

    def append_col(self, value) -> None:
        self.insert_col(self.tableau.diagram[0], value)

    @property
    def lr_tableaux(self) -> Tuple[Tableau, Tableau, int]:
        """Return the left- and the (conjugated) right-handed tableaux as well as the overall sign."""

        l_tab, r_tab, sign = LorentzTableau.lr_tableaux.fget(self)

        return Young_Tableau(l_tab.tableau), Young_Tableau(r_tab.tableau), sign

def get_op_class(field_content: Dict[str,int], derivatives: int, mass_dim: int):
    """
    Returns the operator class object OpClass for given:
    ----------
    field_content
        e.g.: {"H":2, "H+":2}
    derivatives
        e.g.: 4
    mass_dim
        e.g.: 8

    Returns
    -------
        operator class: OpClass(N, nl, nr)
    """
    N = sum(field_content.values())

    nl = derivatives / 2
    nr = derivatives / 2
    for field, multiplicity in field_content.items():
        helicity = model.fields[field].helicity
        nl += multiplicity * (abs(helicity) - helicity) / 2
        nr += multiplicity * (abs(helicity) + helicity) / 2

    nl = int(nl)
    nr = int(nr)

    assert N + nl + nr == mass_dim

    return OpClass(N, nl, nr)


# Examples:
def test():
    tab = [[2,1,1,2], [4,3,3,4]]
    # Properties of the operator
    mass_dim = 8
    field_content = {"H":2, "H+":2}
    derivatives = 4

    op_class = get_op_class(field_content, derivatives, mass_dim)

    # The lorentz tableau is a tableau from which the dotted and undotted epsilon tensors can be derived.
    lt = Lorentz_Tableau(tab, op_class)
    print(f"{lt:nice}")
    # Change one column of the lorentz tableau.
    lt2 = Lorentz_Tableau([[3], [1]], op_class)
    lt[:, 1] = lt2
    print(f"{lt:nice}")

    # Derive the tableaus for the dotted and undotted epsilon tensors.
    l_tab, r_tab, sign = lt.lr_tableaux
    print(f"LH:\n{l_tab:nice}")
    print(f"RH:\n{r_tab:nice}")
    print(f"sign: {sign:d}")

    # Change one column of the "right-handed tableau", i.e. the tableau for the dotted epsilon tensor.
    yt_r = Young_Tableau([[3], [1]])
    r_tab[:,1] = yt_r
    print(f"RH:\n{r_tab:nice}")

    # Derive from the left- and right-handed tableaus the lorentz tableau. This is in a first step the most relevant
    # part for the projection.
    lorentz, sign2 = Lorentz_Tableau.from_lr_tableaux(l_tab, r_tab, op_class)
    print(f"Lorentz tableau:\n{lorentz:nice}")
    print(f"sign: {sign2:d}")

    # Can Young_tableau also handle strings?
    yt = Young_Tableau([["a"], ["b"]])
    print(f"\n{yt:nice}")
    # Append a row
    yt_arow = Young_Tableau([["c"]])
    print(f"\n{yt_arow:nice}")
    yt.append_row(yt_arow)
    print(f"\n{yt:nice}")
    # Append a column
    yt_acol = Young_Tableau([["a"], ["b"]])
    print(f"\n{yt_acol:nice}")
    yt.append_col(yt_acol)
    print(f"\n{yt:nice}")
    del yt[:,0]
    print(f"\n{yt:nice}")


    print("TEST")