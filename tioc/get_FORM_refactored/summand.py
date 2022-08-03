import logging.config
import sys
from abc import ABC, abstractmethod
from typing import Dict, List, Tuple
from pathlib import Path

from tioc import model, op_config, index_config, get_SUN_name
from tioc.general import create_procedure

from .coefficient import Coefficient
from .index import Dummy_Index
from .indices import Indices_Summand, Indices_Operator, Possible_Indices
from .operator import Tensor, Field
from .operators import Tensors, Fields
from tioc.tableau import get_op_class

logger_autoeft = logging.getLogger("autoeft.projection")
logger = logger_autoeft.getChild(__name__)

class Summand_Model(ABC):  # Tensor, Field, Coefficient
    """Base class for single term consisting of an overall coefficient and a product of operators."""
    tensors: List[Tensor]
    fields: Tuple[Field]
    coeff: Coefficient
    fieldcounter: Dict[str, int]
    d: int  # massdimension of operator
    indices: Indices_Summand
    nD: int  # number of derivatives in the summand

    @abstractmethod
    def __init__(self, tensors: List[str], fields: Tuple[str], coeff: str, fp_name: str):
        """Load in fields, tensors and the coefficient."""
        self.name = fp_name
        self.tensors = tensors
        self.fields = fields
        self.coeff = coeff
        self.nD
        self.fieldcounter
        self.fieldstructure

        self.possible_indices = Possible_Indices(
            [index for index in list(self.indices).copy() if not isinstance(index, Dummy_Index)])

    @abstractmethod
    def __repr__(self):
        """Specify the format the general string representation and for printing with repr()."""
        tensor = "*".join(map(str, self.tensors))
        contractedOp = "*".join(map(str, self.fields))
        return f"{tensor:s}*{contractedOp:s}"

    def __str__(self):
        """Specify the format for printing with str() or print() statement function: Here the same as the string representation repr() itself."""
        return self.__repr__()

    def __format__(self, key):
        """Specify the format for "format" function in print statement: Here the same as the string representation repr() itself."""
        if key == "tex":
            return self.tex
        elif key == "complete" or key == "c":
            coeff = self.coeff
            tensor = "*".join(map(str, self.tensors))
            contractedOp = "*".join(map(str, self.fields))
            if tensor:
                return f"({coeff})*{tensor:s}*{contractedOp:s}"
            else:
                return f"({coeff})*{contractedOp:s}"
        else:
            return self.__repr__()

    @property
    def tex(self):
        """
        Tex expression of the entire summand. However, before this can be generated, the tex indices are generated
        using generators, with the range specified in the index configuration file.
        Returns
        -------
            LaTex expression.
        """
        # Instantiate an index generator for each type of index
        generator = {}
        for index in self.indices:
            indrange = Indices_Summand.get_tex_range(index_config[index.typ]["tex_indices"])
            if indrange:
                generator[index.typ] = Indices_Summand.infinite_Indices(indrange)
            else:
                generator[index.typ] = Indices_Summand.infinite_numIndices(index_config[index.typ]["tex_indices"])

        tex_indices = {}  # Dictionary for unique tex names of indices.
        for index in self.indices:
            try:
                if index.typ == "sbasis":
                    ind = index_config
                    tex_indices[index.indname] = f"{index_config['sbasis']['tex_indices']}_{{{index.id}}}"
                else:
                    tex_indices[index.indname] = next(generator[index.typ])
            except StopIteration:
                logger.error("Specified index range is to small to map all occurring indices.")
                sys.exit("STOP")

        #write latex index for all indices:
        for operator in self.tensors:
            for index in operator.indices:
                index.tex = tex_indices[index.indname]
        for operator in self.fields:
            for index in operator.indices:
                index.tex = tex_indices[index.indname]

        # Build tex expression:
        mult_sign = ''  # '*'
        tex_expr = ""
        tensors = mult_sign.join([f"{op:tex}" for op in self.tensors])
        fields = mult_sign.join([f"{op:tex}" for op in self.fields])
        tex_expr += rf"{self.coeff:tex}\\" + mult_sign
        if tensors:
            tex_expr += f"{tensors}{mult_sign}{fields}"
        else:
            tex_expr += f"{fields}"

        return tex_expr

    @property
    def tensors(self):
        return self._tensors

    @tensors.setter
    def tensors(self, value: List):
        self._tensors = Tensors(tuple(map(Tensor, value)))

    @property
    def fields(self):
        return self._fields

    @fields.setter
    def fields(self, fp_fields):
        self._fields = Fields(tuple(map(Field, fp_fields, range(1, len(fp_fields)+1))))

    @property
    def coeff(self):
        return self._coeff

    @coeff.setter
    def coeff(self, value):
        # TODO: Abstraction into operators
        self._coeff = Coefficient(value, self.name)

    @property
    def d(self):
        """
        Mass dimension of operator.
        """
        return self.coeff.mdim

    @property
    def fieldcounter(self):
        """
        Returns dictionary, e.g.
        {'BL': 0, 'WL': 0, 'L': 0, 'Q': 0, 'dC': 0, 'eC': 0, 'uC': 0, 'H': 1,
        'H+': 1, 'L+': 0, 'Q+': 0, 'dC+': 0, 'eC+': 0, 'uC+': 0, 'BL+': 0, 'WL+': 0}
        which tells that there are one ordinary and one conjugated Higgs field in the Summand.
        Returns
        -------
        """
        fieldcount = {key.name: 0 for key in model.fields.values()}
        autoeft_to_projection = {autoeft: projection for field_name, field in {**op_config["bosonfields"], **op_config["fermionfields"]}.items() if field_name != "D" for projection, autoeft in field["autoeft"].items()}
        for field in self.fields:
            try:
                for counted_field in fieldcount.keys():
                    if autoeft_to_projection[counted_field] == field.name:
                        fieldcount[counted_field] += 1
                        break
            except KeyError:
                logger.error(f"Field {field.name:s} doesn't have an autoeft translation.")
                sys.exit("STOP")

        return fieldcount

    @property
    def fieldcounter_stripped(self):
        """
        Remove all 0 entries, e.g.:
        {'BL': 0, 'WL': 0, 'L': 0, 'Q': 0, 'dC': 0, 'eC': 0, 'uC': 0, 'H': 1,
        'H+': 1, 'L+': 0, 'Q+': 0, 'dC+': 0, 'eC+': 0, 'uC+': 0, 'BL+': 0, 'WL+': 0}
        -> {'H': 1, 'H+': 1}
        """
        fieldcount_s = {field: count for field, count in self.fieldcounter.items() if count != 0}

        return fieldcount_s

    @property
    def fieldstructure(self):
        """
        Get the exact field and derivative structure, e.g. term D(uC)*H*H*d+ has structure:
        ((uC,1),(H,0),(H,0),(d+,0))
        I.e. a tuple of as many tuple as there are fields in the Summand is return, where each tuple give first,
        the kind of field at this position and second, the number of derivatives acting on it.
        """
        projection_to_autoeft = {projection: autoeft for field_name, field in
                              {**op_config["bosonfields"], **op_config["fermionfields"]}.items() if field_name != "D"
                              for projection, autoeft in field["autoeft"].items()}
        fieldstructure = []
        for field in self.fields:
            try:
                structure = (projection_to_autoeft[field.name] , field.nD)
                fieldstructure.append(structure)
            except KeyError:
                logger.error(f"Field {field.name:s} doesn't have an autoeft translation.")
                sys.exit("STOP")

        return tuple(fieldstructure)

    @property
    def nD(self):
        """ Number of derivatives in a summand."""
        n_D = 0
        for field in self.fields:
            n_D += field.nD
        return n_D

    @property
    def indices(self):
        """All indices of a summand."""
        indices_expr = Indices_Operator([])
        for tensor in self.tensors:
            indices_expr += tensor.indices
        for field in self.fields:
            indices_expr += field.indices
        if any(isinstance(index, Dummy_Index) for index in indices_expr):
            # Allow uncontracted indices when there is a Dummy_index in the indices.
            return Indices_Summand(indices_expr.indices, allow_uncontracted=True)
        else:
            return Indices_Summand(indices_expr.indices)


    def replace_SUN_indices_by_projection_indices(self):
        """
        Replace all SUN-indices by the corresponding projection indices, which should be saved at this time in the
        projection attribute of each index.
        E.g.: [su2eps](gauge1,gauge2)*H(gauge1)*[H+](gauge2) -> [su2eps](gaugeF1I1,gaugeF2I1)*H(gaugeF1I1)*[H+](gaugeF2I1)
        Returns
        -------
        """
        for sun_group in ["gauge", "colf"]:
            for tensor in self.tensors:
                for index in tensor.indices[sun_group]:
                    index.expr = index.projection
            for field in self.fields:
                for index in field.indices[sun_group]:
                    index.expr = index.projection


class Summand(Summand_Model):
    """Single term consisting of an overall coefficient and a product of operators."""
    fieldcounter: Dict[str, int]

    def __init__(self, tensors: List[str], fields: List[str], coeff: str, fp_name: str):
        super().__init__(tensors, fields, coeff, fp_name)

        # Gauge fields contain special projection index of the form idxF2I1, which is transmitted to the contracted tensors:
        self.gaugeIndicesforProjection_tensors()

        su2 = get_SUN_name(2)
        su3 = get_SUN_name(3)
        self.gaugeTensorsSUN = {su2: self.tensors[su2], su3: self.tensors[su3]}

    def __repr__(self):
        """Specify the format the general string representation and for printing with repr()."""
        return super().__repr__()

    def get_projectionIndex_from_Field(self, ref_index):  # indexID
        """
        Returns the projection index of the contracted gauge field with index specified by indexID.
        """
        for field in self.fields:
            for index in field.indices:
                # if index.id == indexID:
                if index == ref_index.dual_index:
                    return index.projection

    def gaugeIndicesforProjection_tensors(self):
        """
        Gauge fields contain special projection index of the form idxF2I1, which is transmitted to the contracted tensors.

        Returns
        -------

        """
        for tensor in self.tensors:
            for index in tensor.indices:
                for indtype in ["gauge", "colf"]:
                    if index.typ == indtype:
                        # index.projection = self.get_projectionIndex_from_Field(index.id)
                        index.projection = self.get_projectionIndex_from_Field(index)

    @property
    def op_class(self):
        """Operator class of the operator."""
        return get_op_class(self.fieldcounter_stripped, self.nD, self.d)

    #########################
    ### EOM Substitutions ###
    #########################

    def form_higgsEOM(self, form_path: Path):
        """
        Replace D2(H(gauge1)) and D2([H+](gauge1)) by the SM equation of motion.
        Parameters
        -------
        form_path
            Path where the form file is stored
        Returns
        -------

        """
        sl2Ceps = op_config["tensors"]["[sl2Ceps]"]["mathematica"]["sl2Ceps"]  # '[sl2Ceps]'
        su2eps = op_config["tensors"]["[su2eps]"]["mathematica"]["su2eps"]
        su3eps = op_config["tensors"]["[su3eps]"]["mathematica"]["su3eps"]
        # Higgs field
        higgs = op_config["bosonfields"]["H"]["mathematica"]["H"]  # 'H'
        higgs_dagger = op_config["bosonfields"]["H"]["mathematica"]["conj[H]"]  # 'H'
        # Fermions
        l = op_config["fermionfields"]["L"]["mathematica"]["L"]               # L
        l_dagger = op_config["fermionfields"]["L"]["mathematica"]["conj[L]"]  # [L+]
        q = op_config["fermionfields"]["Q"]["mathematica"]["Q"]
        q_dagger = op_config["fermionfields"]["Q"]["mathematica"]["conj[Q]"]
        d = op_config["fermionfields"]["[d_C]"]["mathematica"]["dC"]               # [d_C]
        d_dagger = op_config["fermionfields"]["[d_C]"]["mathematica"]["conj[dC]"]  # [d_C+]
        u = op_config["fermionfields"]["[u_C]"]["mathematica"]["uC"]
        u_dagger = op_config["fermionfields"]["[u_C]"]["mathematica"]["conj[uC]"]
        e = op_config["fermionfields"]["[e_C]"]["mathematica"]["eC"]
        e_dagger = op_config["fermionfields"]["[e_C]"]["mathematica"]["conj[eC]"]
        # Coefficients
        Mu = op_config["coefficients"]["Mu"]["mathematica"]["Mu"]
        lambdah= op_config["coefficients"]["lambdah"]["mathematica"]["lambdah"]
        # Yukawa matrices
        yu = op_config["tensors"]["yu"]["mathematica"]["yu"]
        yu_dagger = op_config["tensors"]["yu"]["mathematica"]["conj[yu]"]
        yd = op_config["tensors"]["yd"]["mathematica"]["yd"]
        yd_dagger = op_config["tensors"]["yd"]["mathematica"]["conj[yd]"]
        ye = op_config["tensors"]["ye"]["mathematica"]["ye"]
        ye_dagger = op_config["tensors"]["ye"]["mathematica"]["conj[ye]"]

        @create_procedure(form_path, "higgsEOM")
        def create_form():
            form = ""
            id_statements = []
            gen_index_gauge = self.possible_indices.generate_index("gauge")  # get a new, i.e. unused SU2 index with 'next(gen_index_gauge)'
            gen_index_colf = self.possible_indices.generate_index("colf")
            gen_index_flav = self.possible_indices.generate_index("flav")
            gen_index_sl = self.possible_indices.generate_index("Lsl")  # next(gen_index_sl) will generate new Lsl and new Usl index
            gen_index_sldot = self.possible_indices.generate_index("Lsldot")
            gauge1 = next(gen_index_gauge)
            gauge2 = next(gen_index_gauge)
            gauge3 = next(gen_index_gauge)
            colf1 = next(gen_index_colf)
            colf2 = next(gen_index_colf)
            colf3 = next(gen_index_colf)
            flav1 = next(gen_index_flav)
            flav2 = next(gen_index_flav)
            lsl1, usl1 = next(gen_index_sl)
            lsl2, usl2 = next(gen_index_sl)
            lsldot1, usldot1 = next(gen_index_sldot)
            lsldot2, usldot2 = next(gen_index_sldot)
            #  EOM of the ordinary Higgs field
            id_statement1 = f"D2({higgs}({gauge1}?gauge)) ="
            id_statement1 += f" + ({Mu}^2) * {higgs}({gauge1})"

            id_statement1 += f" - {lambdah} * {su2eps}({gauge2},{gauge3})*{higgs_dagger}({gauge3}) * {higgs}({gauge2}) * {higgs}({gauge1})"

            id_statement1 += f" + (1/2) * {yu}({flav1},{flav2}) * {su3eps}({colf1},{colf2},{colf3}) * {sl2Ceps}({lsldot1},{lsldot2})"
            id_statement1 += f" * {q_dagger}({usldot2},{gauge1},{colf2},{colf3},{flav1}) * {u_dagger}({usldot1},{colf1},{flav2})"

            id_statement1 += f" - (1/2) * {yd_dagger}({flav2},{flav1}) * {su3eps}({colf1},{colf2},{colf3}) * {sl2Ceps}({usl1},{usl2})"
            id_statement1 += f" * {d}({lsl2},{colf2},{colf3},{flav2}) * {q}({lsl1},{gauge1},{colf1},{flav1})"

            id_statement1 += f" - {ye_dagger}({flav2},{flav1}) * {sl2Ceps}({usl1},{usl2})"
            id_statement1 += f" * {e}({lsl2},{flav2}) * {l}({lsl1},{gauge1},{flav1})"

            id_statements.append(f"id once ifmatch -> 3 {id_statement1:s};\n")

            # EOM of the conjugated Higgs field
            id_statement2 = f"D2({higgs_dagger}({gauge1}?gauge)) ="
            id_statement2 += f" + ({Mu}^2) * {higgs_dagger}({gauge1})"

            id_statement2 += f" - {lambdah} * {su2eps}({gauge2},{gauge3})*{higgs_dagger}({gauge3}) * {higgs}({gauge2}) * {higgs_dagger}({gauge1})"

            id_statement2 += f" + (1/2) * {yu_dagger}({flav2},{flav1}) * {su3eps}({colf1},{colf2},{colf3}) * {sl2Ceps}({usl1},{usl2})" # {sl2Ceps}({usl1},{usl2})
            id_statement2 += f" * {u}({lsl2},{colf2},{colf3},{flav2}) * {q}({lsl1},{gauge1},{colf1},{flav1})"

            id_statement2 += f" - (1/2) * {yd}({flav1},{flav2}) * {su3eps}({colf1},{colf2},{colf3}) * {sl2Ceps}({lsldot1},{lsldot2})" # {sl2Ceps}({lsldot1},{lsldot2})
            id_statement2 += f" * {q_dagger}({usldot2},{gauge1},{colf2},{colf3},{flav1}) * {d_dagger}({usldot1},{colf1},{flav2})"

            id_statement2 += f" - {ye}({flav1},{flav2}) * {sl2Ceps}({lsldot1},{lsldot2})"
            id_statement2 += f" * {l_dagger}({usldot2},{gauge1},{flav1}) * {e_dagger}({usldot1},{flav2})"

            id_statements.append(f"id once ifmatch -> 3 {id_statement2:s};\n")

            for id in id_statements:
                form += f"{id:s}"

            return form

        create_form()

    def form_fieldstrengthtensor(self, form_path: Path):
        """
        # FL(Lsl1, Lsl2, H(?a))
        # FR(Usldot1, Usldot2, H(?a))
        Replace
            FL(Lsl1, Lsl2, H(gauge1)), FL(Lsl1, Lsl2, [H+](gauge1)),
            FR(Usldot1, Usldot2, H(gauge1)), FR(Usldot1, Usldot2, [H+](gauge1))
        by a sum of the U(1) and SU(2) field strength tensors B and W.

        Parameters
        -------
        form_path
            Path where the form file is stored
        Returns
        -------

        """
        su2eps = op_config["tensors"]["[su2eps]"]["mathematica"]["su2eps"]
        # Higgs field
        higgs = op_config["bosonfields"]["H"]["mathematica"]["H"]  # 'H'
        higgs_dagger = op_config["bosonfields"]["H"]["mathematica"]["conj[H]"]  # 'H'
        # field strength tensors
        bl = op_config["bosonfields"]["BL"]["mathematica"]["BL"]
        br = op_config["bosonfields"]["BL"]["mathematica"]["conj[BL]"]
        wl = op_config["bosonfields"]["WL"]["mathematica"]["WL"]
        wr = op_config["bosonfields"]["WL"]["mathematica"]["conj[WL]"]
        # Coefficients
        g1 = op_config["coefficients"]["g1"]["mathematica"]["g1"]
        g2 = op_config["coefficients"]["g2"]["mathematica"]["g2"]

        @create_procedure(form_path, "fieldstrengthtensor")
        def create_form():
            form = "" # #do field={" + f"{higgs}, {higgs_dagger}" + "}\n
            id_statements = []
            gen_index_gauge = self.possible_indices.generate_index("gauge")  # get a new, i.e. unused SU2 index with 'next(gen_index_gauge)'
            gen_index_sl = self.possible_indices.generate_index("Lsl")  # next(gen_index_sl) will generate new Lsl and new Usl index
            gen_index_sldot = self.possible_indices.generate_index("Lsldot")
            gauge1 = next(gen_index_gauge)
            gauge2 = next(gen_index_gauge)
            gauge3 = next(gen_index_gauge)
            gauge4 = next(gen_index_gauge)
            lsl1, usl1 = next(gen_index_sl)
            lsl2, usl2 = next(gen_index_sl)
            lsldot1, usldot1 = next(gen_index_sldot)
            lsldot2, usldot2 = next(gen_index_sldot)

            # Note: The Hypercharges of H(gauge1) is (1/2) and that of [H+](gauge1) is (-1/2).
            # Note: Epsilon in W(gauge, gauge) H(gauge) are there, because two fundamental indices are contracted.

            #  Left handed field strength tensor: FL(Lsl1, Lsl2, H(gauge1))
            id_statement1 = f"FL({lsl1}, {lsl2}, {higgs}({gauge1}?gauge)) ="
            id_statement1 += f" + ({g1}/2) * {bl}({lsl1},{lsl2}) * {higgs}({gauge1})"
            id_statement1 += f" + ({g2}/2) * {su2eps}({gauge1},{gauge2}) * {wl}({lsl1},{lsl2},{gauge2},{gauge3}) * {su2eps}({gauge3},{gauge4}) * {higgs}({gauge4})"

            id_statements.append(f"id once ifmatch -> 4 {id_statement1:s};\n")

            #  Left handed field strength tensor: FL(Lsl1, Lsl2, [H+](gauge1))
            id_statement2 = f"FL({lsl1}, {lsl2}, {higgs_dagger}({gauge1}?gauge)) ="
            id_statement2 += f" + (- {g1}/2) * {bl}({lsl1},{lsl2}) * {higgs_dagger}({gauge1})"
            id_statement2 += f" + ({g2}/2) * {su2eps}({gauge1},{gauge2}) * {wl}({lsl1},{lsl2},{gauge2},{gauge3}) * {su2eps}({gauge3},{gauge4}) * {higgs_dagger}({gauge4})"

            id_statements.append(f"id once ifmatch -> 4 {id_statement2:s};\n")

            #  Right handed field strength tensor: FR(Usldot1, Usldot2, H(gauge1))
            id_statement1 = f"FR({usldot1}, {usldot2}, {higgs}({gauge1}?gauge)) ="
            id_statement1 += f" + ({g1}/2) * {br}({usldot1},{usldot2}) * {higgs}({gauge1})"
            id_statement1 += f" + ({g2}/2) * {su2eps}({gauge1},{gauge2}) * {wr}({usldot1},{usldot2},{gauge2},{gauge3}) * {su2eps}({gauge3},{gauge4}) * {higgs}({gauge4})"

            id_statements.append(f"id once ifmatch -> 4 {id_statement1:s};\n")

            #  Right handed field strength tensor: FR(Usldot1, Usldot2, [H +](gauge1))
            id_statement2 = f"FR({usldot1}, {usldot2}, {higgs_dagger}({gauge1}?gauge)) ="
            id_statement2 += f" + (- {g1}/2) * {br}({usldot1},{usldot2}) * {higgs_dagger}({gauge1})"
            id_statement2 += f" + ({g2}/2) * {su2eps}({gauge1},{gauge2}) * {wr}({usldot1},{usldot2},{gauge2},{gauge3}) * {su2eps}({gauge3},{gauge4}) * {higgs_dagger}({gauge4})"

            id_statements.append(f"id once ifmatch -> 4 {id_statement2:s};\n")


            for id in id_statements:
                form += f"{id:s}"

            # form += "#enddo"
            return form

        create_form()

    def form_fermionEOM(self, form_path: Path):
        """
        Replace EOM(`field', Lsl1, ?a) and EOM(`field', Usldot1, ?a) by the SM equation of motion for fermions.
        Parameters
        -------
        form_path
            Path where the form file is stored
        Returns
        -------

        """
        su2eps = op_config["tensors"]["[su2eps]"]["mathematica"]["su2eps"]
        # Higgs field
        higgs = op_config["bosonfields"]["H"]["mathematica"]["H"]  # 'H'
        higgs_dagger = op_config["bosonfields"]["H"]["mathematica"]["conj[H]"]  # 'H'
        # Fermions
        l = op_config["fermionfields"]["L"]["mathematica"]["L"]               # L
        l_dagger = op_config["fermionfields"]["L"]["mathematica"]["conj[L]"]  # [L+]
        q = op_config["fermionfields"]["Q"]["mathematica"]["Q"]
        q_dagger = op_config["fermionfields"]["Q"]["mathematica"]["conj[Q]"]
        d = op_config["fermionfields"]["[d_C]"]["mathematica"]["dC"]               # [d_C]
        d_dagger = op_config["fermionfields"]["[d_C]"]["mathematica"]["conj[dC]"]  # [d_C+]
        u = op_config["fermionfields"]["[u_C]"]["mathematica"]["uC"]
        u_dagger = op_config["fermionfields"]["[u_C]"]["mathematica"]["conj[uC]"]
        e = op_config["fermionfields"]["[e_C]"]["mathematica"]["eC"]
        e_dagger = op_config["fermionfields"]["[e_C]"]["mathematica"]["conj[eC]"]
        # Yukawa matrices
        yu = op_config["tensors"]["yu"]["mathematica"]["yu"]
        yu_dagger = op_config["tensors"]["yu"]["mathematica"]["conj[yu]"]
        yd = op_config["tensors"]["yd"]["mathematica"]["yd"]
        yd_dagger = op_config["tensors"]["yd"]["mathematica"]["conj[yd]"]
        ye = op_config["tensors"]["ye"]["mathematica"]["ye"]
        ye_dagger = op_config["tensors"]["ye"]["mathematica"]["conj[ye]"]

        @create_procedure(form_path, "fermionEOM")
        def create_form():
            form = ""
            id_statements = []
            gen_index_gauge = self.possible_indices.generate_index("gauge")  # get a new, i.e. unused SU2 index with 'next(gen_index_gauge)'
            gen_index_colf = self.possible_indices.generate_index("colf")
            gen_index_flav = self.possible_indices.generate_index("flav")
            gen_index_sl = self.possible_indices.generate_index("Lsl")  # next(gen_index_sl) will generate new Lsl and new Usl index
            gen_index_sldot = self.possible_indices.generate_index("Lsldot")
            gauge1 = next(gen_index_gauge)
            gauge2 = next(gen_index_gauge)
            colf1 = next(gen_index_colf)
            colf2 = next(gen_index_colf)
            colf3 = next(gen_index_colf)
            flav1 = next(gen_index_flav)
            flav2 = next(gen_index_flav)
            lsl1, usl1 = next(gen_index_sl)
            lsldot1, usldot1 = next(gen_index_sldot)
            #  EOM of Q
            id_statement1 = f"EOM({q},{usldot1}?Usldot,{gauge1}?gauge,{colf1}?colf,{flav1}?flav) ="
            id_statement1 += f" + i_ * {yu}({flav1},{flav2}) * {higgs_dagger}({gauge1}) * {u_dagger}({usldot1},{colf1},{flav2})"
            id_statement1 += f" + i_ * {yd}({flav1},{flav2}) * {higgs}({gauge1}) * {d_dagger}({usldot1},{colf1},{flav2})"
            id_statements.append(f"id once ifmatch -> 4 {id_statement1:s};\n")

            #  EOM of Q+
            id_statement2 = f"EOM({q_dagger},{lsl1}?Lsl,{gauge1}?gauge,{colf2}?colf,{colf3}?colf,{flav1}?flav) ="
            id_statement2 += f" + i_ * {yu_dagger}({flav2},{flav1}) * {u}({lsl1},{colf2},{colf3},{flav2}) * {higgs}({gauge1})"
            id_statement2 += f" + i_ * {yd_dagger}({flav2},{flav1}) * {d}({lsl1},{colf2},{colf3},{flav2}) * {higgs_dagger}({gauge1})"
            id_statements.append(f"id once ifmatch -> 4 {id_statement2:s};\n")
            id_statements.append("\n")

            #  EOM of L
            id_statement1 = f"EOM({l},{usldot1}?Usldot,{gauge1}?gauge,{flav1}?flav) ="
            id_statement1 += f" + i_ * {ye}({flav1},{flav2}) * {higgs}({gauge1}) * {e_dagger}({usldot1},{flav2})"
            id_statements.append(f"id once ifmatch -> 4 {id_statement1:s};\n")

            #  EOM of L+
            id_statement2 = f"EOM({l_dagger},{lsl1}?Lsl,{gauge1}?gauge,{flav1}?flav) ="
            id_statement2 += f" + i_ * {ye_dagger}({flav2},{flav1}) * {e}({lsl1},{flav2}) * {higgs_dagger}({gauge1})"
            id_statements.append(f"id once ifmatch -> 4 {id_statement2:s};\n")
            id_statements.append("\n")

            #  EOM of [d_C+]
            id_statement1 = f"EOM({d_dagger},{lsl1}?Lsl,{colf1}?colf,{flav1}?flav) ="
            id_statement1 += f" + i_ * {yd_dagger}({flav1},{flav2}) * {su2eps}({gauge1},{gauge2}) * {higgs_dagger}({gauge2}) * {q}({lsl1},{gauge1},{colf1},{flav2})"
            id_statements.append(f"id once ifmatch -> 4 {id_statement1:s};\n")

            #  EOM of [d_C]
            id_statement2 = f"EOM({d},{usldot1}?Usldot,{colf2}?colf,{colf3}?colf,{flav1}?flav) ="
            id_statement2 += f" + i_ * {yd}({flav2},{flav1}) * {su2eps}({gauge1},{gauge2}) * {q_dagger}({usldot1},{gauge2},{colf2},{colf3},{flav2}) * {higgs}({gauge1})"
            id_statements.append(f"id once ifmatch -> 4 {id_statement2:s};\n")
            id_statements.append("\n")

            #  EOM of [u_C+]
            id_statement1 = f"EOM({u_dagger},{lsl1}?Lsl,{colf1}?colf,{flav1}?flav) ="
            id_statement1 += f" + (- i_) * {yu_dagger}({flav1},{flav2}) * {su2eps}({gauge1},{gauge2}) * {higgs}({gauge1}) * {q}({lsl1},{gauge2},{colf1},{flav2})"
            id_statements.append(f"id once ifmatch -> 4 {id_statement1:s};\n")

            #  EOM of [u_C]
            id_statement2 = f"EOM({u},{usldot1}?Usldot,{colf2}?colf,{colf3}?colf,{flav1}?flav) ="
            id_statement2 += f" + i_ * {yu}({flav2},{flav1}) * {su2eps}({gauge1},{gauge2}) * {q_dagger}({usldot1},{gauge2},{colf2},{colf3},{flav2}) * {higgs_dagger}({gauge1})"
            id_statements.append(f"id once ifmatch -> 4 {id_statement2:s};\n")
            id_statements.append("\n")

            #  EOM of [e_C+]
            id_statement1 = f"EOM({e_dagger},{lsl1}?Lsl,{flav1}?flav) ="
            id_statement1 += f" + i_ * {ye_dagger}({flav1},{flav2}) * {su2eps}({gauge1},{gauge2}) * {higgs_dagger}({gauge2}) * {l}({lsl1},{gauge1},{flav2})"
            id_statements.append(f"id once ifmatch -> 4 {id_statement1:s};\n")

            #  EOM of [e_C]
            id_statement2 = f"EOM({e},{usldot1}?Usldot,{flav1}?flav) ="
            id_statement2 += f" + i_ * {ye}({flav2},{flav1}) * {su2eps}({gauge1},{gauge2}) * {l_dagger}({usldot1},{gauge2},{flav2}) * {higgs}({gauge1})"
            id_statements.append(f"id once ifmatch -> 4 {id_statement2:s};\n")

            for id in id_statements:
                form += f"{id:s}"

            return form

        create_form()

    def form_fieldstrengthtensorEOM(self, form_path: Path):
        """
        Replace EOM(BL, Lsl1, Usldot1), EOM(WL, Lsl1, Usldot1, gauge1, gauge2)
        and
        EOM(BR, Lsl1, Usldot1), EOM(WR, Lsl1, Usldot1, gauge1, gauge2)
        by the SM equation of motion.
        Parameters
        -------
        form_path
            Path where the form file is stored
        Returns
        -------

        """
        cov = op_config["fermionfields"]["D"]["mathematica"]["cov"]  # 'D'
        su2eps = op_config["tensors"]["[su2eps]"]["mathematica"]["su2eps"]
        su3eps = op_config["tensors"]["[su3eps]"]["mathematica"]["su3eps"]
        # Higgs field
        higgs = op_config["bosonfields"]["H"]["mathematica"]["H"]  # 'H'
        higgs_dagger = op_config["bosonfields"]["H"]["mathematica"]["conj[H]"]  # 'H'
        # Fermions
        l = op_config["fermionfields"]["L"]["mathematica"]["L"]               # L
        l_dagger = op_config["fermionfields"]["L"]["mathematica"]["conj[L]"]  # [L+]
        q = op_config["fermionfields"]["Q"]["mathematica"]["Q"]
        q_dagger = op_config["fermionfields"]["Q"]["mathematica"]["conj[Q]"]
        d = op_config["fermionfields"]["[d_C]"]["mathematica"]["dC"]               # [d_C]
        d_dagger = op_config["fermionfields"]["[d_C]"]["mathematica"]["conj[dC]"]  # [d_C+]
        u = op_config["fermionfields"]["[u_C]"]["mathematica"]["uC"]
        u_dagger = op_config["fermionfields"]["[u_C]"]["mathematica"]["conj[uC]"]
        e = op_config["fermionfields"]["[e_C]"]["mathematica"]["eC"]
        e_dagger = op_config["fermionfields"]["[e_C]"]["mathematica"]["conj[eC]"]
        # field strength tensors
        bl = op_config["bosonfields"]["BL"]["mathematica"]["BL"]
        br = op_config["bosonfields"]["BL"]["mathematica"]["conj[BL]"]
        wl = op_config["bosonfields"]["WL"]["mathematica"]["WL"]
        wr = op_config["bosonfields"]["WL"]["mathematica"]["conj[WL]"]
        # Coefficients
        g1 = op_config["coefficients"]["g1"]["mathematica"]["g1"]
        g2 = op_config["coefficients"]["g2"]["mathematica"]["g2"]

        @create_procedure(form_path, "fieldstrengthtensorEOM")
        def create_form():
            form = ""
            id_statements = []
            gen_index_gauge = self.possible_indices.generate_index("gauge")  # get a new, i.e. unused SU2 index with 'next(gen_index_gauge)'
            gen_index_colf = self.possible_indices.generate_index("colf")
            gen_index_flav = self.possible_indices.generate_index("flav")
            gen_index_sl = self.possible_indices.generate_index("Lsl")  # next(gen_index_sl) will generate new Lsl and new Usl index
            gen_index_sldot = self.possible_indices.generate_index("Lsldot")
            gauge1 = next(gen_index_gauge)
            gauge2 = next(gen_index_gauge)
            gauge3 = next(gen_index_gauge)
            gauge4 = next(gen_index_gauge)
            colf1 = next(gen_index_colf)
            colf2 = next(gen_index_colf)
            colf3 = next(gen_index_colf)
            flav1 = next(gen_index_flav)
            lsl1, usl1 = next(gen_index_sl)
            lsldot1, usldot1 = next(gen_index_sldot)
            #  EOM of the left-handed field strength tensor of U(1), BL
            id_statement1 = f"EOM({bl},{lsl1}?Lsl,{usldot1}?Usldot) ="
            rHS1 = f" {g1}*("

            rHS1 += f" + (i_/2) * {su2eps}({gauge1},{gauge2}) * ( {cov}({lsl1},{usldot1},{higgs}({gauge1})) * {higgs_dagger}({gauge2})"
            rHS1 += f" - {higgs}({gauge1}) * {cov}({lsl1},{usldot1},{higgs_dagger}({gauge2})) )"

            rHS1 += f" - (1/6) * {su3eps}({colf1},{colf2},{colf3}) * {su2eps}({gauge1},{gauge2}) * {q}({lsl1},{gauge1},{colf1},{flav1}) * {q_dagger}({usldot1},{gauge2},{colf2},{colf3},{flav1})"

            rHS1 += f" + (2/3) * {su3eps}({colf1},{colf2},{colf3}) * {u}({lsl1},{colf2},{colf3},{flav1}) * {u_dagger}({usldot1},{colf1},{flav1})"

            rHS1 += f" - (1/3) * {su3eps}({colf1},{colf2},{colf3}) * {d}({lsl1},{colf2},{colf3},{flav1}) * {d_dagger}({usldot1},{colf1},{flav1})"

            rHS1 += f" + {su2eps}({gauge1},{gauge2}) * {l}({lsl1},{gauge1},{flav1}) * {l_dagger}({usldot1},{gauge2},{flav1})"

            rHS1 += f" + 2 * {e}({lsl1},{flav1}) * {e_dagger}({usldot1},{flav1})"

            rHS1 += ")"

            id_statement1 += rHS1
            id_statements.append(f"id once ifmatch -> 4 {id_statement1:s};\n")

            #  EOM of the right-handed field strength tensor of U(1), BR
            id_statement2 = f"EOM({br},{lsl1}?Lsl,{usldot1}?Usldot) ="

            id_statement2 += rHS1
            id_statements.append(f"id once ifmatch -> 4 {id_statement2:s};\n")


            #  EOM of the left-handed field strength tensor of SU(2), WL
            # Note: The SU2-indices are explicitly symmetrized.
            id_statement3 = f"EOM({wl},{lsl1}?Lsl,{usldot1}?Usldot,{gauge1}?gauge,{gauge2}?gauge) ="
            rHS2 = f" - ({g2}/2) * {su2eps}({gauge1},{gauge3}) * {su2eps}({gauge2},{gauge4}) * ("

            def gauge_indices(gaugek, gaugel):
                """For the symmetrization of the SU2-indices."""
                output = f" + i_ * ( {cov}({lsl1},{usldot1},{higgs}({gaugel})) * {higgs_dagger}({gaugek})"
                output += f" - {higgs}({gaugel}) * {cov}({lsl1},{usldot1},{higgs_dagger}({gaugek})) )"
                output += f" - (1/3) * {su3eps}({colf1},{colf2},{colf3}) * {q}({lsl1},{gaugel},{colf1},{flav1}) * {q_dagger}({usldot1},{gaugek},{colf2},{colf3},{flav1})"
                output += f" + {l}({lsl1},{gaugel},{flav1}) * {l_dagger}({usldot1},{gaugek},{flav1})"

                return output

            # 'Normal' order of SU2 indices:
            rHS2 += gauge_indices(gauge3, gauge4)
            # SU2 indices swapped:
            rHS2 += gauge_indices(gauge4, gauge3)

            rHS2 += ")"

            id_statement3 += rHS2
            id_statements.append(f"id once ifmatch -> 4 {id_statement3:s};\n")

            #  EOM of the right-handed field strength tensor of SU(2), WR
            id_statement4 = f"EOM({wr},{lsl1}?Lsl,{usldot1}?Usldot,{gauge1}?gauge,{gauge2}?gauge) ="
            id_statement4 += rHS2
            id_statements.append(f"id once ifmatch -> 4 {id_statement4:s};\n")

            for id in id_statements:
                form += f"{id:s}"

            return form

        create_form()


