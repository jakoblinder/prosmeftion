import re
from yaml import safe_load
import logging
import sys
from typing import Dict, List, Tuple, Union
from collections.abc import MutableMapping
from abc import ABC, abstractmethod
from copy import copy

from tioc import CONFIG_PATH, index_pattern, index_config
from .index import Index, Dummy_Index

logger_autoeft = logging.getLogger("autoeft.projection")
logger = logger_autoeft.getChild(__name__)

class Indices_Model(Index, MutableMapping):
    indices: Tuple[Index]  # Tuple of indices in one term.

    abc = "abcdefghijklmnopqrstuvwxyz"
    ABC = abc.upper()
    alpha_beta_gamma = [r"\alpha", r"\beta", r"\gamma", r"\delta", r"\varepsilon", r"\zeta", r"\eta", r"\theta",
                        r"\iota", r"\kappa", r"\lambda", r"\mu", r"\nu", r"\omicron", r"\pi", r"\rho", r"\sigma",
                        r"\tau", r"\upsilon", r"\phi", r"\chi", r"\psi", r"\omega",
                        ]
    alpha_beta_gamma_dot = list(map(lambda ind : rf"\dot{{{ind}}}", alpha_beta_gamma))

    @abstractmethod
    def __init__(self, indices: Tuple[Index], allow_uncontracted=False):
        # assert type(indices) == tuple
        if allow_uncontracted:
            # Remove double occuring indices
            reduced_indices = list(indices)

            # Remove more than one time occurring indices from self.indices.
            for index in reduced_indices.copy():
                while index.is_in(reduced_indices) > 1:
                    reduced_indices.remove(index)

            self._indices = tuple(reduced_indices)
            # self._indices = indices
        else:
            self.indices = indices

    @abstractmethod
    def __repr__(self):
        return f"{', '.join(map(repr, self.indices))}"

    def __str__(self):
        """Specify the format for printing with str() or print() statement function. """
        return self.__repr__()

    def __format__(self, key):
        """Specify the format for "format" function in print statement: Here the same as the print statement itself."""
        if key == "projection" or key == "p":
            return f"{','.join([f'{index:p}'for index in self]):s}"
        else:
            return self.__repr__()

    def __eq__(self, other):
        if len(self.indices) != len(other.indices):
            # logger.debug("Different number of indices in each object.")
            return False
        equal = []
        for i in self.indices:
            if i.is_in(other.indices):
                equal.append(True)
            else:
                # logger.debug(f"Index {i} doesn't occur in {other.indices}.")
                return False
        if len(equal) == len(self.indices) and all(equal):
            return True

    def __add__(self, other):
        """Combine two sets of indices and remove all indices occurring more than one time."""
        indices = list(self.indices)
        # Remove more than one time occurring indices from self.indices.
        for i in indices.copy():
            while i.is_in(indices) > 1:
                indices.remove(i)
        other_indices = list(other.indices)
        for i in other_indices:
            if not i.is_in(indices):
                indices.append(i)
        return type(self)(indices, allow_uncontracted=True)

    def __len__(self):
        """List length"""
        return len(self.indices)

    def __getitem__(self, key):
        """
        If 'key' is of type slice or int we just get a list item. If 'key' is of type str, all indices of the type indicated
        by 'key' are returned in their occurring order, except for Summand_indices where the indices are sorted by their id.
        Parameters
        ----------
        key: int, slice, str

        Returns
        -------
        """
        if isinstance(key, int):
            return self.indices[key]
        elif isinstance(key, slice):
            return type(self)(self.indices[key])
        elif isinstance(key, str):
            index_types = [indextyp for indextyp in index_config.keys()]
            index_special_types = ["sl", "sldot", "sl2C", "dummy"]
            # index_types = ['lor', 'Lsl', 'Usl', 'Lsldot', 'Usldot', 'spin', 'gauge', 'gaugeadj', 'colf', 'cola', 'flav']
            if key not in (index_types + index_special_types):
                logger.error(f"Key {key:s} is not a possible index typ.")
                sys.exit("STOP")
            elif key in index_special_types:
                index_of_typ = []
                if key == "sl":
                    for index in self.indices:
                        if index.typ == "Lsl" or index.typ == "Usl":
                            index_of_typ.append(index)
                elif key == "sldot":
                    for index in self.indices:
                        if index.typ == "Lsldot" or index.typ == "Usldot":
                            index_of_typ.append(index)
                elif key == "sl2C":
                    for index in self.indices:
                        if index.typ in ("Lsl", "Usl", "Lsldot", "Usldot"):
                            index_of_typ.append(index)
                elif key == "dummy":
                    for index in self.indices:
                        if isinstance(index, Dummy_Index):
                            index_of_typ.append(index)
            else:
                index_of_typ = []
                for index in self.indices:
                    if index.typ == key:
                        index_of_typ.append(index)
            # if index_of_typ:
            return type(self)(tuple(index_of_typ), allow_uncontracted=True)
            # else:
            #     return []
        else:
            logger.error(f"Key is of type {type(key)}, but it should be of type int, str or slice.")
            sys.exit("STOP")

    def __setitem__(self, key, value):
        if not isinstance(key, int) or (isinstance(key, int) and key < 0):
            logger.error("Key has to be of typ int and >= 0.")
            sys.exit("STOP")
        if not isinstance(value, Index):
            logger.error("The value which will be set has to be of type Index.")
            sys.exit("STOP")
        logger.debug(f"The index {self.indices[key]:s} will be rewritten with {value:s}.")
        indices = list(self.indices)
        indices[key] = value
        self.indices = tuple(indices)

    def __delitem__(self, key):
        if not isinstance(key, int) or (isinstance(key, int) and key < 0):
            logger.error("Key has to be of typ int and >= 0.")
            sys.exit("STOP")
        logger.debug(f"The index {self.indices[key]:s} will be deleted.")
        indices = list(self.indices)
        del indices[key]
        self.indices = tuple(indices)

    def __iter__(self):
        return iter(self.indices)

    def insert(self, key, val):
        if not isinstance(val, Index):
            logger.error("The value which will be inserted has to be of type Index.")
            sys.exit("STOP")
        if not isinstance(key, int):
            logger.error("The key where the value is inserted has to be of type integer.")
            sys.exit("STOP")
        if key == len(self):
            logger.debug(f"The index {val:s} will be appended.")
        else:
            logger.debug(f"The index {val:s} will be inserted.")
        indices = list(self.indices)
        indices.insert(key, val)
        self.indices = tuple(indices)

    def append(self, val):
        if not isinstance(val, Index):
            logger.error("The value which will be appended has to be of type Index.")
            sys.exit("STOP")
        self.insert(len(self.indices), val)

    def clear(self):
        return type(self)(tuple(self.indices.clear()), allow_uncontracted=True)

    def copy(self):
        return type(self)(tuple(list(self.indices).copy()), allow_uncontracted=True)

    @staticmethod
    def infinite_Indices(finite_list, max=5):
        """
        Generator from a given finite tuple of indices, an infinite generator which yields an index in the order
        of the indices occurring in the tuple and attaches primes if all indices are used until the maximum number
        of primes is exhausted.
        ==> Maximum number of indices is therefore len(finite_list)*max
        Parameters
        ----------
        finite_list : tuple
            Contains str with indices.
        max
            Maximum number of primes which should be added to an index.
        Returns
        -------
        """
        for j in range(max):
            prime = r"^{" + r"\prime" * j + r"}"
            for i in finite_list:
                if j:  # j > 0
                    yield i + prime
                else:  # j == 0
                    yield i

    @staticmethod
    def infinite_numIndices(single_index, max=50):
        """
        Gives generator which returns from a given single index, e.g. \alpha indices \alpha_{1} to \alpha_{max}.
        Parameters
        ----------
        single_index
        max

        Returns
        -------
        """
        for i in range(1, max + 1):
            yield f"{single_index}_{{{i:d}}}"

    @staticmethod
    def get_tex_range(tex_indices):
        r"""
        Create finite list of indices out of specified range or by numbering a single given index.
        Parameters
        ----------
        tex_indices (list or str)
            String like "\mu...\sigma", "u...z", "\alpha" or "\dot{\alpha}" or list of strings which will
            be the tex indices.
        Returns
        -------

        """
        if type(tex_indices) == list or type(tex_indices) == tuple:
            return tuple(tex_indices)
        match_range = re.match(r"(?P<start>[a-zA-Z\\\{\}]+)\.\.\.(?P<end>[a-zA-Z\\\{\}]+)", tex_indices)
        match_single = re.match(r"(?P<ind>[a-zA-Z\\\{\}]+)$", tex_indices)
        if match_range:
            if match_range.group("start") == match_range.group("end"):
                logger.warning(f"The defined range of indices {match_range.group('start')}...{match_range.group('end')} has the same beginning and end.")
            if match_range.group("start") in Indices_Model.abc and match_range.group("end") in Indices_Model.abc:
                start = Indices_Model.abc.index(match_range.group("start"))
                end = Indices_Model.abc.index(match_range.group("end"))
                if start > end:
                    indices = Indices_Model.abc[end: start+1]
                    indices = indices[::-1]
                else:
                    indices = Indices_Model.abc[start: end+1]
            elif match_range.group("start") in Indices_Model.ABC and match_range.group("end") in Indices_Model.ABC:
                start = Indices_Model.ABC.index(match_range.group("start"))
                end = Indices_Model.ABC.index(match_range.group("end"))
                if start > end:
                    indices = Indices_Model.ABC[end: start+1]
                    indices = indices[::-1]
                else:
                    indices = Indices_Model.ABC[start: end + 1]
            elif (match_range.group("start") in Indices_Model.alpha_beta_gamma and match_range.group("end") in Indices_Model.alpha_beta_gamma) or (match_range.group("start") in Indices_Model.alpha_beta_gamma and match_range.group("end") == r"\epsilon") or (match_range.group("start") == r"\epsilon") and match_range.group("end") in Indices_Model.alpha_beta_gamma:
                if match_range.group("start") == r"\epsilon":
                    start = Indices_Model.alpha_beta_gamma.index(r"\varepsilon")
                else:
                    start = Indices_Model.alpha_beta_gamma.index(match_range.group("start"))
                if match_range.group("end") == r"\epsilon":
                    end = Indices_Model.alpha_beta_gamma.index(r"\varepsilon")
                else:
                    end = Indices_Model.alpha_beta_gamma.index(match_range.group("end"))
                if start > end:
                    indices = Indices_Model.alpha_beta_gamma[end: start+1]
                    indices = indices[::-1]
                else:
                    indices = Indices_Model.alpha_beta_gamma[start: end + 1]
            elif (match_range.group("start") in Indices_Model.alpha_beta_gamma_dot and match_range.group("end") in Indices_Model.alpha_beta_gamma_dot) or (match_range.group("start") in Indices_Model.alpha_beta_gamma_dot and match_range.group("end") == r"\epsilon") or (match_range.group("start") == r"\epsilon") and match_range.group("end") in Indices_Model.alpha_beta_gamma_dot:
                if match_range.group("start") == r"\epsilon":
                    start = Indices_Model.alpha_beta_gamma_dot.index(r"\varepsilon")
                else:
                    start = Indices_Model.alpha_beta_gamma_dot.index(match_range.group("start"))
                if match_range.group("end") == r"\epsilon":
                    end = Indices_Model.alpha_beta_gamma_dot.index(r"\varepsilon")
                else:
                    end = Indices_Model.alpha_beta_gamma_dot.index(match_range.group("end"))
                if start > end:
                    indices = Indices_Model.alpha_beta_gamma_dot[end: start+1]
                    indices = indices[::-1]
                else:
                    indices = Indices_Model.alpha_beta_gamma_dot[start: end + 1]
            else:
                logger.error("Start and end index of the index range are not part of the same naming range or at least one of them is not in any naming range.")
                sys.exit("STOP")
        elif match_single:
            # single Indices matched, which gets numbers as subscript
            return False
        else:
            logger.error("Indexrange does not fit the required pattern.")
            sys.exit("STOP")
        return tuple(indices)

class Indices_Operator(Indices_Model):
    indices: Tuple[Index]  # Tuple of indices in one term.

    def __init__(self, indices: Tuple[Index], allow_uncontracted=False):
        super().__init__(indices, allow_uncontracted)

    @property
    def indices(self):
        return self._indices

    @indices.setter
    def indices(self, fp_indices):
        self._indices = fp_indices

    def __repr__(self):
        return super().__repr__()

    def __add__(self, other):
        """Combine two sets of indices and KEEP all indices which occur more than one time."""
        indices_expr = tuple(list(self.indices) + list(other.indices))
        return type(self)(indices_expr)

class Possible_Indices(Indices_Model):
    indices: Tuple[Index]  # Tuple of all possible indices in one Summand/ Term or whatever.

    def __init__(self, indices: Tuple[Index], allow_uncontracted=True):
        super().__init__(tuple(indices), allow_uncontracted)

    def __repr__(self):
        return super().__repr__()

    def __format__(self, key):
        """Specify the format for "format" function in print statement: Here the same as the print statement itself."""
        return super().__format__(key)

    @property
    def indices(self):
        return self._indices

    @indices.setter
    def indices(self, fp_indices):
        reduced_indices = list(fp_indices)

        # Remove more than one time occurring indices from self.indices.
        for index in reduced_indices.copy():
            while index.is_in(reduced_indices) > 1:
                reduced_indices.remove(index)

        self._indices = tuple(reduced_indices)

    def __len__(self):
        """List length"""
        return len(self.indices)

    def __getitem__(self, key: Union[int, slice, str]):
        """
        If 'key' is of type slice or int we just get a list item. If 'key' is of type str, all indices of the type indicated
        by 'key' are sorted by their 'id' (i.e. for example flav2, flav1 -> flav1, flav2 and so on) and returned.

        The reason behind this sorting is that for further manipulation of the expression there are Sets of SL2C-indices
        in FORM necessary, which specify the contraction of those indices. Especially these are Sets like
        Set ULsl: Usl1, Lsl1, Usl2, Lsl2, ... ;
        Set LUsl: Lsl1, Usl1, Lsl2, Usl2, ... ;
        which can then be used in id-statements like for example:
        id [sl2Ceps](Usl1?Usl, Usl2?ULsl[k]) * [sl2Ceps](Lsl1?LUsl[k], Lsl3?Lsl) = + [sl2CdK](Usl1,Lsl3);
        where Usl2 has to be a superscript SL2C-index and Lsl1 has to be exactly the same subscript SL2C-index.

        Note: The index order in the self.indices attribute is completely unimportant at this point and is only relevant
              in the indices of an operator.

        Parameters
        ----------
        key: int, slice, str

        Returns
        -------
        """
        indices = super().__getitem__(key)
        if isinstance(indices, Possible_Indices) and isinstance(key, str):
            # Note: sorted of a string of letters capital and non-capital and of numbers is sorted like:
            # 0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz
            sort_ind = sorted(indices, key=lambda x: str(x.id))
            return type(self)(sort_ind, allow_uncontracted=True)
        else:
            return indices

    def __setitem__(self, key, value):
        super().__setitem__(key, value)

    def __delitem__(self, key):
        super().__delitem__(key)

    def __iter__(self):
        return super().__iter__()

    def get_sl2C_sets(self):
        """
        Return multiple lists of SL2C-indices, which specify the Usl, Lsl, Usldot, Lsldot, ULsl, LUsl, ULsldot, LUsldot
        lists of indices, necessary for contraction of SL2C-indices.

        Ensure that for each index in the Lsl list there is the exact dual index in the same spot in the Usl list and vice versa.

        Returns
        -------

        """
        lsl    = self["Lsl"]
        usl    = self["Usl"]
        lsldot = self["Lsldot"]
        usldot = self["Usldot"]
        def complete_indices(indices, dual_indices):
            """

            Parameters
            ----------
            indices
                Indices which are checked, i.e. it is checked that for each index there is the correct
            dual_indices
                dual index given by dual_indices in the correct position.
            Returns
            -------

            """
            for index in indices:
                # check that the dual usl index is in the usl list and if not add him
                if not index.dual_index.is_in(dual_indices):
                    dual_indices.append(index.dual_index)
                    # sort them again
                    dual_indices = type(self)(sorted(dual_indices, key=lambda x: str(x.id)))
            # check, that position are correct
            for i, index in enumerate(indices):
                if not index.dual_index == dual_indices[i]:
                    logger.error(f"The index {index} doesn't have the correct dual index at the right position, but rather {dual_indices[i]}.")
                    sys.exit("STOP")

            return indices, dual_indices

        lsl,    usl    = complete_indices(lsl,    usl)
        usl,    lsl    = complete_indices(usl,    lsl)
        lsldot, usldot = complete_indices(lsldot, usldot)
        usldot, lsldot = complete_indices(usldot, lsldot)

        return lsl, usl, lsldot, usldot

    def generate_index(self, typ: str, exclude_indices: List[Index] = [], derIndex: Union[bool, int, str]=False, fp_min: int=1, fp_max: int=None):
        """
        Call with
            gen_index = summand.possible_indices.generate_index("flav")
            new_index = next(gen_index)
        for non SL2C-indices and with
            gen_index = summand.possible_indices.generate_index("Lsl")
            new_indices = next(gen_index)
        for the generation of two new SL2C indices, e.g. Lsl1 and Usl1. The indices have to be generated directly both,
        due to the specific declaration and Set construction of the FORM indices.

        Parameters
        ----------
        typ
            typ of index
        derIndex
            Is this an index of a derivative or not and if yes on which derivative (int).
            Write "new" to assign unique "object()" to SL2C-indices.
        fp_min: int
            Minimum number on an index.
        fp_max: int
            Possible maximum number on an index.
            =>E.g.: For fp_min = 2 and fp_max=12345: flav2, flav3, ..., flav12345 would be returned.
        exclude_indices
            Optional parameter, which allows to specify additional indices which should not be taken in the generation of new indices.

        Returns
        -------
        Returns a generator object, which generates unused indices of the specified typ.
        """
        assert type(derIndex) in [bool, int, str], 'derIndex is of the wrong type.'
        index_types = [indextyp for indextyp in index_config.keys()]
        index_types += ["sl2C"]  # generate all 4 types of SL2C-indices at ones -> necessary for derivative index generation in ibp relations.
        if typ not in index_types:
            logger.error(f"The type {typ} is not one of the possible types {', '.join(index_types)}")
            sys.exit("STOP")

        if derIndex == "new": assert  typ == "sl2C", "For generation of new indices it is necessary to generate all 4 indices at ones with the 'type' sl2C."

        sentinel = object()
        def count(min, max=None):
            """
            count(10) --> 10 11 12 13 14 ...
            count(2.5, 7.5) -> 2.5 3.5 4.5 5.5 6.5 -> max itself is not printed
            """
            n = min
            while True:
                if max:
                    if n >= max:
                        yield sentinel
                yield n
                n += 1

        continue_object = object()
        def gen_new_index(typ: str, i: int, derIndex: Union[bool, int, object]):
            index = Index(f"{typ}{i}", derIndex)
            if exclude_indices:
                index_list = self[typ] + exclude_indices[typ]
            else:
                index_list = self[typ]
            if index.is_in(index_list.indices):
                return continue_object
                # continue
            elif i is sentinel:
                logger.error("Not possible to generate a new index, since generator is out of range.")
                sys.exit("STOP")
            else:
                self.append(index)
                if index.dual_index != index:
                    self.append(index.dual_index)
                    # yield index, index.dual_index
                    return index, index.dual_index
                else:
                    # yield index
                    return index

        for i in count(fp_min, fp_max):
            spec_derivative = object()  # specify derivative uniquely
            if derIndex == "new" or type(derIndex) == object: derIndex = spec_derivative
            if typ == "sl2C":
                new_index_sl = gen_new_index("Lsl", i, derIndex)
                new_index_sldot = gen_new_index("Lsldot", i, derIndex)
                if new_index_sl is continue_object or new_index_sldot is continue_object:
                    continue
                else:
                    lsl, usl = new_index_sl
                    lsldot, usldot = new_index_sldot
                    yield lsl, lsldot, usl, usldot
            else:
                new_index = gen_new_index(typ, i, derIndex)
                if new_index is continue_object:
                    continue
                else:
                    yield new_index

class Indices_Summand(Indices_Model):
    indices: Tuple[Index]  # Tuple of indices in one term.
    tex_indices: Dict  # Dictionary for each index containing the unique tex name.

    def __init__(self, indices: Tuple[Index], allow_uncontracted=False):
        super().__init__(indices, allow_uncontracted)
        # self.possible_indices = Possible_Indices([index for index in list(self.indices).copy() if not isinstance(index, Dummy_Index) ])

    @property
    def indices(self):
        return self._indices

    @indices.setter
    def indices(self, fp_indices: Indices_Operator):
        """Check that indices are contracted."""
        contracted = False
        contract = []
        uncontractedInd = []
        for i in fp_indices:
            if fp_indices.count(i) == 2:
                contract.append(True)
            elif fp_indices.count(i) > 2:
                logger.error(f"At least the index {i} occurs more then two times.")
                sys.exit("STOP")
            else:
                uncontractedInd.append(i)
        # Copy list for iteration in such away that only elements from "uncontractedInd" are removed in an
        # "uncontractedInd.remove("something")" order.
        uncontractedInd_tmp = uncontractedInd.copy()
        # uncontractedInd_tmp = list(map(str, uncontractedInd_tmp))
        # Consider now the possible uncontracted indices which can only be SL2C-indices:
        # TODO: Could be simplified a little bit by the use of the dual_index of an index.
        for index in uncontractedInd_tmp:
            if (index.typ == "Usl" and Index(f"Lsl{index.id}") in uncontractedInd) or (index.typ == "Lsl" and Index(f"Usl{index.id}") in uncontractedInd):
                contract += [True, True]
                uncontractedInd.remove(Index(f"Usl{index.id}"))
                uncontractedInd.remove(Index(f"Lsl{index.id}"))
            if (index.typ == "Usldot" and Index(f"Lsldot{index.id}") in uncontractedInd) or (index.typ == "Lsldot" and Index(f"Usldot{index.id}") in uncontractedInd):
                contract += [True, True]
                uncontractedInd.remove(Index(f"Usldot{index.id}"))
                uncontractedInd.remove(Index(f"Lsldot{index.id}"))
        if len(uncontractedInd) == 0 and all(contract):
            contracted = True
        elif all([index.typ == "sbasis" for index in uncontractedInd]):
            contracted = True

        if contracted:
            # 'Remove double occuring indices'
            reduced_indices = list(fp_indices)

            # Remove more than one time occurring indices from self.indices.
            for index in reduced_indices.copy():
                while index.is_in(reduced_indices) > 1:
                    reduced_indices.remove(index)

            self._indices = tuple(reduced_indices)
        else:
            logger.error(f"Not all indices are contracted. The indices {', '.join(map(str,uncontractedInd_tmp))} are not contracted.")
            sys.exit("STOP")

    def __repr__(self):
        return super().__repr__()

    def __len__(self):
        """List length"""
        return len(self.indices)

    def __getitem__(self, key: Union[int, slice, str]):
        """
        If 'key' is of type slice or int we just get a list item. If 'key' is of type str, all indices of the type indicated
        by 'key' are sorted by their 'id' (i.e. for example flav2, flav1 -> flav1, flav2 and so on) and returned.

        The reason behind this sorting is that for further manipulation of the expression there are Sets of SL2C-indices
        in FORM necessary, which specify the contraction of those indices. Especially these are Sets like
        Set ULsl: Usl1, Lsl1, Usl2, Lsl2, ... ;
        Set LUsl: Lsl1, Usl1, Lsl2, Usl2, ... ;
        which can then be used in id-statements like for example:
        id [sl2Ceps](Usl1?Usl, Usl2?ULsl[k]) * [sl2Ceps](Lsl1?LUsl[k], Lsl3?Lsl) = + [sl2CdK](Usl1,Lsl3);
        where Usl2 has to be a superscript SL2C-index and Lsl1 has to be exactly the same subscript SL2C-index.

        Note: The index order in the self.indices attribute is completely unimportant at this point and is only relevant
              in the indices of an operator.

        Parameters
        ----------
        key: int, slice, str

        Returns
        -------
        """
        indices = super().__getitem__(key)
        if isinstance(indices, Indices_Summand) and isinstance(key, str):
            # Note: sorted of a string of letters capital and non-capital and of numbers is sorted like:
            # 0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz
            sort_ind = sorted(indices, key=lambda x: str(x.id))
            return type(self)(sort_ind, allow_uncontracted=True)
        else:
            return indices

    def __setitem__(self, key, value):
        super().__setitem__(key, value)

    def __delitem__(self, key):
        super().__delitem__(key)

    def __iter__(self):
        return super().__iter__()

class Indices_Term(Indices_Model):
    indices: Tuple[Index]  # Tuple of indices in one term.

    def __init__(self, indices: Tuple[Index], allow_uncontracted=False):
        super().__init__(indices, allow_uncontracted)

    def __repr__(self):
        return super().__repr__()

    def __format__(self, key):
        """Specify the format for "format" function in print statement: Here the same as the print statement itself."""
        if key == "debug":
            return ", ".join(f"{self.tex_indices[index.name]}({index.name})" for index in self.indices)
        else:
            return self.__repr__()

    @property
    def indices(self):
        return self._indices

    @indices.setter
    def indices(self, fp_indices):
        reduced_indices = list(fp_indices)

        # Remove more than one time occurring indices from self.indices.
        for index in reduced_indices.copy():
            while index.is_in(reduced_indices) > 1:
                reduced_indices.remove(index)

        self._indices = tuple(reduced_indices)


