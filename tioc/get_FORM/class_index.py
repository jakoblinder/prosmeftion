import re
from yaml import safe_load
import logging
import sys
from typing import Dict, List, Tuple

from tioc import CONFIG_PATH

logger_autoeft = logging.getLogger("autoeft.projection")
logger = logger_autoeft.getChild("index")

#TODO: Index range could be read from the modelfile
with open(CONFIG_PATH / "index.yml", "r") as file:
    index_config = safe_load(file)

# print(index_config["Usldot"]["tex_indices"])
class Index():
    name: str
    typ: str  # gauge, colf, lorentz, etc.
    id: str  # cryptic id, e.g1234A in combination with the type it ensures clear contractions
    tex_name: str  # readable name which should be printed in LaTex for the index
    dimension: int
    description: str

    def __init__(self, name, **data):
        self.name = name
        if "tex_name" in data.keys():
            self.tex_name = data["tex_name"]
        else:
            self.tex_name = None
        self.typ, self.number, self.id = Index.get_typ_id(self.name)
        self.dimension = index_config[self.typ]["dimension"]
        self.description = index_config[self.typ]["description"]

    def __repr__(self):
        if self.tex_name:
            return f"{self.tex_name:s}"
        else:
            return f"{self.name:s}"
    def __str__(self):
        """Specify the format for printing with str() or print() statement function. """
        return self.__repr__()
    def __format__(self, key):
        """Specify the format for "format" function in print statement: Here the same as the print statement itself."""
        if key == "tex":
            return f"{self.tex_name}({self.name})"
        else:
            return self.__repr__()

    def __eq__(self, other):
        return self.name == other.name

    @staticmethod
    def get_typ_id(index_expr):
        index = index_expr  # e.g.: "flav67Av89"
        match = re.match(r"(?P<typ>lor|Lsl(?!dot)|Usl(?!dot)|Lsldot|Usldot|spin|gauge(?!adj)|gaugeadj|colf|cola|flav)(?P<number>[a-zA-Z\d]+)",
                         index)
        if match:
            typ = match.group("typ")
            num = match.group("number")
            if typ == "Lsl" or typ == "Usl":
                ID = f"sl-{num}"
            elif typ == "Lsldot" or typ == "Usldot":
                ID = f"sldot-{num}"
            else:
                ID = f"{typ}-{num}"
        else:
            logger.error(f"The index {index} is not written in a valid index format. An index should look like e.g.: flav67Av89.")
            sys.exit("STOP")
        return typ, num, ID

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

class Indices(Index):
    indices: Tuple[Index]  # Tuple of indices in one term.
    # tex_indices: List  # List of allowed indices for naming
    tex_indices: Dict  # Dictionary for each index containing the unique tex name.

    abc = "abcdefghijklmnopqrstuvwxyz"
    ABC = abc.upper()
    alpha_beta_gamma = [r"\alpha", r"\beta", r"\gamma", r"\delta", r"\varepsilon", r"\zeta", r"\eta", r"\theta",
                        r"\iota", r"\kappa", r"\lambda", r"\mu", r"\nu", r"\omicron", r"\pi", r"\rho", r"\sigma",
                        r"\tau", r"\upsilon", r"\phi", r"\chi", r"\psi", r"\omega",
                        ]
    alpha_beta_gamma_dot = list(map(lambda ind : f"\dot{{{ind}}}", alpha_beta_gamma))
    def __init__(self, indices):
        # assert type(indices) == tuple
        self.indices = indices
        self.tex_indices = self.get_tex_indices()

    def __repr__(self):
        output = []
        for index in self.indices:
            output.append(f"{self.tex_indices[index.name]}({index.name})")
        return ", ".join(output)
        # return f"{','.join(map(str, self.indices))}"
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
        indices = list(self.indices)
        # Remove frequently occuring indices from self.indices.
        for i in indices.copy():
            while i.is_in(indices) > 1:
                indices.remove(i)
        other_indices = list(other.indices)
        for i in other_indices:
            if not i.is_in(indices):
                indices.append(i)
        return Indices(indices)

    @staticmethod
    def infinte_Indices(finite_list, max=5):
        """
        Generator from a given finite tuple of indices, which yields and index and in the order of the tuple index and
        attaches primes if all indices are used until the maximum number of primes is exhausted.
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
            prime = r"^{" + "\prime" * j + r"}"
            for i in finite_list:
                if j:
                    yield i + prime
                else:
                    yield i

    @staticmethod
    def infinte_numIndices(single_index, max=50):
        for i in range(1, max + 1):
            yield f"{single_index}_{{{i:d}}}"

    @staticmethod
    def get_tex_range(tex_indices):
        match_range = re.match(r"(?P<start>[a-zA-Z\\\{\}]+)\.\.\.(?P<end>[a-zA-Z\\\{\}]+)", tex_indices)
        match_single = re.match(r"(?P<ind>[a-zA-Z\\\{\}]+)$", tex_indices)
        if match_range:
            if match_range.group("start") == match_range.group("end"):
                logger.warning(f"The defined range of indices {match_range.group('start')}...{match_range.group('end')} has the same beginning and end.")
            if match_range.group("start") in Indices.abc and match_range.group("end") in Indices.abc:
                start = Indices.abc.index(match_range.group("start"))
                end = Indices.abc.index(match_range.group("end"))
                if start > end:
                    indices = Indices.abc[end: start+1]
                    indices = indices[::-1]
                else:
                    indices = Indices.abc[start: end+1]
            elif match_range.group("start") in Indices.ABC and match_range.group("end") in Indices.ABC:
                start = Indices.ABC.index(match_range.group("start"))
                end = Indices.ABC.index(match_range.group("end"))
                if start > end:
                    indices = Indices.ABC[end: start+1]
                    indices = indices[::-1]
                else:
                    indices = Indices.ABC[start: end + 1]
            elif (match_range.group("start") in Indices.alpha_beta_gamma and match_range.group("end") in Indices.alpha_beta_gamma) or (match_range.group("start") in Indices.alpha_beta_gamma and match_range.group("end") == r"\epsilon") or (match_range.group("start") == r"\epsilon") and match_range.group("end") in Indices.alpha_beta_gamma:
                if match_range.group("start") == r"\epsilon":
                    start = Indices.alpha_beta_gamma.index(r"\varepsilon")
                else:
                    start = Indices.alpha_beta_gamma.index(match_range.group("start"))
                if match_range.group("end") == r"\epsilon":
                    end = Indices.alpha_beta_gamma.index(r"\varepsilon")
                else:
                    end = Indices.alpha_beta_gamma.index(match_range.group("end"))
                if start > end:
                    indices = Indices.alpha_beta_gamma[end: start+1]
                    indices = indices[::-1]
                else:
                    indices = Indices.alpha_beta_gamma[start: end + 1]
            elif (match_range.group("start") in Indices.alpha_beta_gamma_dot and match_range.group("end") in Indices.alpha_beta_gamma_dot) or (match_range.group("start") in Indices.alpha_beta_gamma_dot and match_range.group("end") == r"\epsilon") or (match_range.group("start") == r"\epsilon") and match_range.group("end") in Indices.alpha_beta_gamma_dot:
                if match_range.group("start") == r"\epsilon":
                    start = Indices.alpha_beta_gamma_dot.index(r"\varepsilon")
                else:
                    start = Indices.alpha_beta_gamma_dot.index(match_range.group("start"))
                if match_range.group("end") == r"\epsilon":
                    end = Indices.alpha_beta_gamma_dot.index(r"\varepsilon")
                else:
                    end = Indices.alpha_beta_gamma_dot.index(match_range.group("end"))
                if start > end:
                    indices = Indices.alpha_beta_gamma_dot[end: start+1]
                    indices = indices[::-1]
                else:
                    indices = Indices.alpha_beta_gamma_dot[start: end + 1]
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

    def get_tex_indices(self):
        # Instantiate an index generator for each type of index
        generator = {}
        for index in self.indices:
            indrange = Indices.get_tex_range(index_config[index.typ]["tex_indices"])
            if indrange:
                generator[index.typ] = Indices.infinte_Indices(indrange)
            else:
                generator[index.typ] = Indices.infinte_numIndices(index_config[index.typ]["tex_indices"])

        tex_indices = {}  # Dictionary for unique tex names of indices.
        for index in self.indices:
            try:
                tex_indices[index.name] = next(generator[index.typ])
            except StopIteration:
                logger.error("Specified index range is to small to map all occurring indices.")
                sys.exit("STOP")

        return tex_indices


# ind = Indices.infinte_Indices(["a","b"])
# ind2 = Indices.infinte_Indices(["a","b"])
# indnum = Indices.infinte_numIndices(r"\alpha")
# dot = Indices.alpha_beta_gamma_dot
# print(next(ind))
# print(next(ind))
# print(next(ind))
# print("Pause")
# print(next(ind2))
# indices = Indices.get_tex_indices(r"\dot{\omega}...\dot{\omicron}")
# print(indices)

# alpha = Index("Lsldot12345", tex_name = r"\alpha")
# alpha = Index("Lsldot123")
# beta = Index("Lsldot124")
# gamma = Index("Lsldot125")
# gauge1 = Index("gauge3")
# gauge2 = Index("gauge1")
# flav1 = Index("flav2345")
# flav2 = Index("flav987")
#
# indices = Indices((alpha, gauge2, beta, flav2, gamma, flav1, gauge1))
#
# print(indices)


