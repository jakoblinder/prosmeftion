import re
from yaml import safe_load
import logging
import sys
from typing import Dict, List, Tuple
from abc import ABC, abstractmethod
from copy import copy

from tioc import CONFIG_PATH, index_pattern
from .index import Index

logger_autoeft = logging.getLogger("autoeft.projection")
logger = logger_autoeft.getChild(__name__)

class Indices_Model(Index):
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
            self._indices = indices
        else:
            self.indices = indices

    @abstractmethod
    def __repr__(self):
        return f"{','.join(map(str, self.indices))}"

    def __str__(self):
        """Specify the format for printing with str() or print() statement function. """
        return self.__repr__()

    def __format__(self, key):
        """Specify the format for "format" function in print statement: Here the same as the print statement itself."""
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

    def __getitem__(self, ii):
        """Get a list item"""
        if isinstance(ii, slice):
            return type(self)(self.indices[ii])
        else:
            return self.indices[ii]

    def __setitem__(self, ii, val):
        indices = list(self.indices)
        assert type(val) == Index
        indices[ii] = val
        self.indices = tuple(indices)

    # def __delitem__(self, ii):
    #     """Delete an item"""
    #     del self.indices[ii]
    #

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

    def __init__(self, indices: Tuple[Index]):
        super().__init__(indices)

    def __repr__(self):
        return super().__repr__()

    def __add__(self, other):
        """Combine two sets of indices and KEEP all indices which occur more than one time."""
        indices_expr = tuple(list(self.indices) + list(other.indices))
        return type(self)(indices_expr)

class Indices_Summand(Indices_Model):
    indices: Tuple[Index]  # Tuple of indices in one term.

    def __init__(self, indices: Tuple[Index], allow_uncontracted=False):
        super().__init__(indices, allow_uncontracted)

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
            else:
                uncontractedInd.append(i)
        # Copy list for iteration in such away that only elements from "uncontractedInd" are removed in an
        # "uncontractedInd.remove("something")" order.
        uncontractedInd_tmp = uncontractedInd.copy()
        # uncontractedInd_tmp = list(map(str, uncontractedInd_tmp))
        # Consider now the possible uncontracted indices which can only be SL2C-indices:
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

        if contracted:
            'Remove double occuring indices'
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

class Indices_Term(Indices_Model):
    indices: Tuple[Index]  # Tuple of indices in one term.
    tex_indices: Dict  # Dictionary for each index containing the unique tex name.
    def __init__(self, indices: Tuple[Index]):
        super().__init__(indices)

    def __repr__(self):
        return ", ".join(f"{self.tex_indices[index.name]}({index.name})" for index in self.indices)

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


