import logging.config
import sys
from abc import ABC, abstractmethod
from typing import Dict, List, Tuple, Union
from pathlib import Path
from copy import copy

from prosmeftion import model, op_config, index_config, get_SUN_name, get_commuting_op
from prosmeftion.general import create_procedure

from .coefficient import Coefficient, Factor
from .index import Index, Dummy_Index, LP_Index
from .indices import Indices_Summand, Indices_Operator, Possible_Indices
from .operator import Tensor, Field
from .operators import Tensors, Fields
from .tableau import get_op_class, Young_Tableau
from .lorentz import LR_Tableaux

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
        tensors = "*".join(map(str, self.tensors))
        contractedOp = "*".join(map(str, self.fields))

        return f"{tensors:s}*{contractedOp:s}" if tensors else f"{contractedOp:s}"

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
        tensors = f"\n{mult_sign}".join([f"{op:tex}" for op in self.tensors])
        fields = f"\n{mult_sign}".join([f"{op:tex}" for op in self.fields])
        tex_expr += rf"{self.coeff:tex}" + f"\n{mult_sign}"
        if tensors:
            tex_expr += f"{tensors}\n{mult_sign}{fields}"
        else:
            tex_expr += f"{fields}"

        return tex_expr

    @property
    def tensors(self):
        return self._tensors

    @tensors.setter
    def tensors(self, fp_tensors: List):
        if isinstance(fp_tensors, Tensors):
            self._tensors = fp_tensors
        else:
            self._tensors = Tensors(tuple(map(Tensor, fp_tensors)))

    @property
    def fields(self):
        return self._fields

    @fields.setter
    def fields(self, fp_fields):
        if isinstance(fp_fields, Fields):
            self._fields = fp_fields
        else:
            self._fields = Fields(tuple(map(Field, fp_fields, range(1, len(fp_fields)+1))))

    @property
    def coeff(self):
        return self._coeff

    @coeff.setter
    def coeff(self, fp_coeff):
        if isinstance(fp_coeff, Coefficient):
            self._coeff = fp_coeff
        else:
            self._coeff = Coefficient(fp_coeff, self.name)

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
                logger.warning(f"Field {field.name:s} doesn't have an autoeft translation.")
                return
                # sys.exit("STOP")

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
                structure = (projection_to_autoeft[field.name], field.nD)
                fieldstructure.append(structure)
            except KeyError:
                logger.warning(f"Field {field.name:s} doesn't have an autoeft translation.")
                return
                # sys.exit("STOP")

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

    @property
    def possible_indices(self):
        """
        All possibly in addition occurring indices which have to be declared for FORM, before FORM is run are saved in
        the attribute _possible_indices, except the ones which occur for sure in the expression, since those are the
        indices of the operators.
        However, since the attribute possible_indices should return also the "for sure" indices, those are added
        manually. The philosophy is that the declaration of to many indices is better than to less.
        Returns
        -------

        """
        indices_of_summand = Possible_Indices([index for index in list(self.indices).copy() if not isinstance(index, Dummy_Index)])
        self._possible_indices += indices_of_summand
        return self._possible_indices

    @possible_indices.setter
    def possible_indices(self, fp_possible_indices):
        self._possible_indices = fp_possible_indices

    def copy(self):
        """
        Returns a copy of the object.
        -------

        """
        return type(self)(self.tensors.copy(), self.fields.copy(), self.coeff.copy(), self.name)

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

    def __init__(self, tensors: List[Union[str,Tensor]], fields: List[Union[str,Field]], coeff: Union[str,Coefficient], fp_name: str):
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

    def copy(self):
        summand = super().copy()
        su2 = get_SUN_name(2)
        su3 = get_SUN_name(3)
        summand.gaugeTensorsSUN = {su2: summand.tensors[su2], su3: summand.tensors[su3]}

        return summand

    #######################
    ### SL2C Conversion ###
    #######################

    def form_convertDirac(self):
        """
        Convert Dirac spinors into Weyl spinors
        #call convertDirac(leftspinor, rightspinor)
        Do replacements in general way, but replace xi, [xi+], chi, [chi+] by the stuff written in the diracspinor op_config entry,
         like [0, "[e_C+]"] -> xi = 0 and [chi+] = [e_C+] for the electron.
        convertDirac then has to be called for every combination of spinors (not two charge conjugated ones).
        in convertDirac is written.

        Returns
        -------

        """
        # Spinors
        all_spinors = {name: field for name, field in op_config["fermionfields"].items() if
                       "diracspinor" in field.keys()}
        all_spinors = {name: weyl for field in all_spinors.values() for name, weyl in field["diracspinor"].items()}
        all_spinors = {name: [get_commuting_op(entry) if entry else 0 for entry in weyl] for name, weyl in all_spinors.items()}
        # CFunction xi1, [xi1+], chi1, [chi1+], xi2, [xi2+], chi2, [chi2+]
        # Tensors
        sl2Ceps = op_config["tensors"]["[sl2Ceps]"]["mathematica"]["sl2Ceps"]
        gamma = op_config["tensors"]["gamma"]["mathematica"]["gamma"]
        sigma2lor = op_config["tensors"]["sigma2lor"]["mathematica"]["sigma2lor"]
        sigma = op_config["tensors"]["sigma"]["mathematica"]["sigma"]
        sigmabar = op_config["tensors"]["sigmabar"]["mathematica"]["sigmabar"]
        sigma2 = op_config["tensors"]["sigma2"]["mathematica"]["sigma2"]
        sigmabar2 = op_config["tensors"]["sigmabar2"]["mathematica"]["sigmabar2"]
        # Indices
        gen_index_sl = self.possible_indices.generate_index("Lsl")  # next(gen_index_sl) will generate new Lsl and new Usl index
        gen_index_sldot = self.possible_indices.generate_index("Lsldot")

        def for_unique_spinorpair(lsbar, rs, indices, repetition):
            lsl1, usl2, lsl2, usl2, lsldot1, usldot1, lsldot2, usldot2 = indices
            id_statements = []

            if lsbar[-3:] == "bar" and rs[-1:] != "C":
                ######################
                ### psibar*...*psi ###
                ######################
                chi1, xi1_dagger = all_spinors[lsbar]
                xi2, chi2_dagger = all_spinors[rs]
                lsbar = get_commuting_op(lsbar)
                rs = get_commuting_op(rs)
                # psibar*psi
                id = ""
                id += f"{lsbar}(op1?op,spin1?spin,?a)*{rs}(op2?op,spin1?spin,?b) = "
                if chi1 == 0 or xi2 == 0:
                    id += "0"
                else:
                    id += f"{sl2Ceps}({usl1},{usl2})*{chi1}(op1,{lsl2},?a)*{xi2}(op2,{lsl1},?b)"
                id += " + "
                if xi1_dagger == 0 or chi2_dagger == 0:
                    id += "0"
                else:
                    id += f"{sl2Ceps}({lsldot1},{lsldot2})*{xi1_dagger}(op1,{usldot2},?a)*{chi2_dagger}(op2,{usldot1},?b)"
                id_statements.append(id)

                # psibar*gamma*psi
                id = ""
                id += f"{lsbar}(op1?op,spin1?spin,?a)*{gamma}(lor1?lor,spin1?spin,spin2?spin)*{rs}(op2?op,spin2?spin,?b) = "
                if chi1 == 0 or chi2_dagger == 0:
                    id += "0"
                else:
                    id += f"- {chi1}(op1,{lsl1},?a)*{sigma}(lor1,{usl1},{lsldot1})*{chi2_dagger}(op2,{usldot1},?b)"
                id += " - "
                if xi1_dagger == 0 or xi2 == 0:
                    id += "0"
                else:
                    id += f"{xi1_dagger}(op1,{usldot1},?a)*{sigmabar}(lor1,{lsldot1},{usl1})*{xi2}(op2,{lsl1},?b)"
                id_statements.append(id)

                # psibar*simga2lor*psi
                id = ""
                id += f"{lsbar}(op1?op,spin1?spin,?a)*{sigma2lor}(lor1?lor,lor2?lor,spin1?spin,spin2?spin)*{rs}(op2?op,spin2?spin,?b) = "
                if chi1 == 0 or xi2 == 0:
                    id += "0"
                else:
                    id += f"- {chi1}(op1,{lsl1},?a)*{sigma2}(lor1,lor2,{usl1},{usl2})*{xi2}(op2,{lsl2},?b)"
                id += " - "
                if xi1_dagger == 0 or chi2_dagger == 0:
                    id += "0"
                else:
                    id += f"{xi1_dagger}(op1,{usldot1},?a)*{sigmabar2}(lor1,lor2,{lsldot1},{lsldot2})*{chi2_dagger}(op2,{usldot2},?b)"
                id_statements.append(id)

            elif lsbar[-4:] == "barC" and rs[-1:] != "C":
                ######################
                ### psibarC*...*psi ###
                ######################
                xi1, chi1_dagger = all_spinors[lsbar]
                xi2, chi2_dagger = all_spinors[rs]
                lsbar = get_commuting_op(lsbar)
                rs = get_commuting_op(rs)

                # psibar*psi
                id = ""
                id += f"{lsbar}(op1?op,spin1?spin,?a)*{rs}(op2?op,spin1?spin,?b) = "
                if xi1 == 0 or xi2 == 0:
                    id += "0"
                else:
                    id += f"{sl2Ceps}({usl1},{usl2})*{xi1}(op1,{lsl2},?a)*{xi2}(op2,{lsl1},?b)"
                id += " + "
                if chi1_dagger == 0 or chi2_dagger == 0:
                    id += "0"
                else:
                    id += f"{sl2Ceps}({lsldot1},{lsldot2})*{chi1_dagger}(op1,{usldot2},?a)*{chi2_dagger}(op2,{usldot1},?b)"
                id_statements.append(id)

                # psibar*gamma*psi
                id = ""
                id += f"{lsbar}(op1?op,spin1?spin,?a)*{gamma}(lor1?lor,spin1?spin,spin2?spin)*{rs}(op2?op,spin2?spin,?b) = "
                if xi1 == 0 or chi2_dagger == 0:
                    id += "0"
                else:
                    id += f"- {xi1}(op1,{lsl1},?a)*{sigma}(lor1,{usl1},{lsldot1})*{chi2_dagger}(op2,{usldot1},?b)"
                id += " - "
                if chi1_dagger == 0 or xi2 == 0:
                    id += "0"
                else:
                    id += f"{chi1_dagger}(op1,{usldot1},?a)*{sigmabar}(lor1,{lsldot1},{usl1})*{xi2}(op2,{lsl1},?b)"
                id_statements.append(id)

                # psibar*simga2lor*psi
                id = ""
                id += f"{lsbar}(op1?op,spin1?spin,?a)*{sigma2lor}(lor1?lor,lor2?lor,spin1?spin,spin2?spin)*{rs}(op2?op,spin2?spin,?b) = "
                if xi1 == 0 or xi2 == 0:
                    id += "0"
                else:
                    id += f"- {xi1}(op1,{lsl1},?a)*{sigma2}(lor1,lor2,{usl1},{usl2})*{xi2}(op2,{lsl2},?b)"
                id += " - "
                if chi1_dagger == 0 or chi2_dagger == 0:
                    id += "0"
                else:
                    id += f"{chi1_dagger}(op1,{usldot1},?a)*{sigmabar2}(lor1,lor2,{lsldot1},{lsldot2})*{chi2_dagger}(op2,{usldot2},?b)"
                id_statements.append(id)

            elif lsbar[-3:] == "bar" and rs[-1:] == "C":
                ######################
                ### psibar*...*psiC ###
                ######################
                lsbar = leftspinor
                rs = rightspinor + "C"
                chi1, xi1_dagger = all_spinors[lsbar]
                chi2, xi2_dagger = all_spinors[rs]
                lsbar = get_commuting_op(lsbar)
                rs = get_commuting_op(rs)

                # psibar*psi
                id = ""
                id += f"{lsbar}(op1?op,spin1?spin,?a)*{rs}(op2?op,spin1?spin,?b) = "
                if chi1 == 0 or chi2 == 0:
                    id += "0"
                else:
                    id += f"{sl2Ceps}({usl1},{usl2})*{chi1}(op1,{lsl2},?a)*{chi2}(op2,{lsl1},?b)"
                id += " + "
                if xi1_dagger == 0 or xi2_dagger == 0:
                    id += "0"
                else:
                    id += f"{sl2Ceps}({lsldot1},{lsldot2})*{xi1_dagger}(op1,{usldot2},?a)*{xi2_dagger}(op2,{usldot1},?b)"
                id_statements.append(id)

                # psibar*gamma*psi
                id = ""
                id += f"{lsbar}(op1?op,spin1?spin,?a)*{gamma}(lor1?lor,spin1?spin,spin2?spin)*{rs}(op2?op,spin2?spin,?b) = "
                if chi1 == 0 or xi2_dagger == 0:
                    id += "0"
                else:
                    id += f"- {chi1}(op1,{lsl1},?a)*{sigma}(lor1,{usl1},{lsldot1})*{xi2_dagger}(op2,{usldot1},?b)"
                id += " - "
                if xi1_dagger == 0 or chi2 == 0:
                    id += "0"
                else:
                    id += f"{xi1_dagger}(op1,{usldot1},?a)*{sigmabar}(lor1,{lsldot1},{usl1})*{chi2}(op2,{lsl1},?b)"
                id_statements.append(id)

                # psibar*simga2lor*psi
                id = ""
                id += f"{lsbar}(op1?op,spin1?spin,?a)*{sigma2lor}(lor1?lor,lor2?lor,spin1?spin,spin2?spin)*{rs}(op2?op,spin2?spin,?b) = "
                if chi1 == 0 or chi2 == 0:
                    id += "0"
                else:
                    id += f"- {chi1}(op1,{lsl1},?a)*{sigma2}(lor1,lor2,{usl1},{usl2})*{chi2}(op2,{lsl2},?b)"
                id += " - "
                if xi1_dagger == 0 or xi2_dagger == 0:
                    id += "0"
                else:
                    id += f"{xi1_dagger}(op1,{usldot1},?a)*{sigmabar2}(lor1,lor2,{lsldot1},{lsldot2})*{xi2_dagger}(op2,{usldot2},?b)"
                id_statements.append(id)

            form = ""
            for id in id_statements:
                form += f"id once ifmatch->{repetition} {id};\n"

            return form

        barfields = []
        fields = []
        for field in self.fields:
            if field.non_conj_name in op_config["fermionfields"].keys():
                if "bar" in field.name:
                    barfields.append(field.name)
                else:
                    fields.append(field.name)


        n_bilinears = (len(fields) + len(barfields)) // 2

        # Group the occurring fermions into pair in which they could have occurred in the summand.
        possible_bilinears = []
        for barfield in barfields:
            for field in fields:
                if barfield[-4:] == "barC" and field[-1:] == "C":
                    continue
                possible_bilinears.append((barfield, field))

        possible_bilinears = list(set(possible_bilinears))

        form = []
        for i in range(1, n_bilinears + 1):
            lsl1, usl1 = next(gen_index_sl)
            lsl2, usl2 = next(gen_index_sl)
            lsldot1, usldot1 = next(gen_index_sldot)
            lsldot2, usldot2 = next(gen_index_sldot)
            indices = [lsl1, usl2, lsl2, usl2, lsldot1, usldot1, lsldot2, usldot2]
            for leftspinor, rightspinor in possible_bilinears:
                    form.append(for_unique_spinorpair(leftspinor, rightspinor, indices, i))
            form.append(f"\nlabel {i};\n\n")

        return "".join(form), n_bilinears, n_bilinears

    def form_convertDerivative(self, n_der:int=4):
        sigma = op_config["tensors"]["sigma"]["mathematica"]["sigma"]
        D = op_config["fermionfields"]["D"]["mathematica"]["cov"]
        # Indices
        gen_index_sl = self.possible_indices.generate_index("Lsl")  # next(gen_index_sl) will generate new Lsl and new Usl index
        gen_index_sldot = self.possible_indices.generate_index("Lsldot")

        def derivatives(nD:int, name:str):
            left = ""
            right_derivatives = ""
            right_sigmas = f"((-1/2)^{nD})*"

            for i in range(1, nD + 1):
                left += f"{D}(lor{i}?lor,"
                lsl, usl = next(gen_index_sl)
                lsldot, usldot = next(gen_index_sldot)
                right_derivatives += f"{D}({lsl},{usldot},"
                right_sigmas += f"{sigma}(lor{i},{usl},{lsldot})*"

            left += f"{name}(op1?op)"
            left += nD*")"

            right_derivatives += f"{name}(op1)"
            right_derivatives += nD*")"
            right = right_sigmas + right_derivatives

            return f"{left} = {right}"


        id_statements = []
        for field in self.fields:
            if field.nD:
                id_statements.append(derivatives(field.nD, field.name))

        form = ""
        for id in id_statements:
            form += f"id once {id};\n"

        return form

    def form_convertFieldstrengthTensor(self):
        fieldstrength_tensors = {name: field for name, field in op_config["bosonfields"].items() if "helicity" not in field.keys()}
        fieldstrength_tensors_names = [list(fsT["mathematica"].values())[0] for fsT in fieldstrength_tensors.values()]
        handed_fsTs = {list(fsT["mathematica"].values())[0]: op_config["bosonfields"][list(fsT["autoeft"].keys())[0]]["mathematica"] for fsT in fieldstrength_tensors.values()}
        fLs = {name: list(mathematica.values())[0] for name, mathematica in handed_fsTs.items()}
        fRs = {name: list(mathematica.values())[1] for name, mathematica in handed_fsTs.items()}

        # Tensors
        sigma2 = op_config["tensors"]["sigma2"]["mathematica"]["sigma2"]
        sigmabar2 = op_config["tensors"]["sigmabar2"]["mathematica"]["sigmabar2"]
        # imaginary unit
        im = op_config["coefficients"]["I"]["mathematica"]["I"]

        fsTs = []
        for field in self.fields:
            if field.name in fieldstrength_tensors_names:
                fsTs.append(field.name)

        # Indices
        if fsTs:
            gen_index_sl = self.possible_indices.generate_index("Lsl")  # next(gen_index_sl) will generate new Lsl and new Usl index
            gen_index_sldot = self.possible_indices.generate_index("Lsldot")

        form = ""
        for fsT in fsTs:
            fL = get_commuting_op(fLs[fsT])
            fR = get_commuting_op(fRs[fsT])
            fsT = get_commuting_op(fsT)

            lsl1, usl1 = next(gen_index_sl)
            lsl2, usl2 = next(gen_index_sl)
            lsldot1, usldot1 = next(gen_index_sldot)
            lsldot2, usldot2 = next(gen_index_sldot)

            # ordinary fieldstrength tensor
            id = f"{fsT}(op1?op,?a,lor1?lor,lor2?lor,?b)"
            id += " = "
            id += f"(+ {im}/4)*({fR}(op1,?a,?b,{usldot1},{usldot2})*{sigmabar2}(lor1,lor2,{lsldot1},{lsldot2}) - {fL}(op1,?a,?b,{lsl1},{lsl2})*{sigma2}(lor1,lor2,{usl1},{usl2}))"
            form += f"id once {id};\n"

            # For dual fieldstrength tensors, i.e. eps(lor1,lor2,lor3,lor4)*F(lor3,lor4), substitute F(lor1,lor2) and
            # simplify the sigma2 matrices with '#call lorepsandSigma2'.

        return form, len(fsTs)

    def form_simplifySigma2(self, n_sigma2):
        """

        Parameters
        ----------
        n_sigma2
            Number of possible sigma2 and sigma2bar matrices: Since they (can) only occur due to fieldstrength tensors
            and Dirac bilinears, the number of the maximum possible ones is given by the number of bilinears plus
            the number of fieldstrenght tensors.

        Returns
        -------

        """
        # Tensors
        sl2Ceps = op_config["tensors"]["[sl2Ceps]"]["mathematica"]["sl2Ceps"]
        sigma = op_config["tensors"]["sigma"]["mathematica"]["sigma"]
        sigmabar = op_config["tensors"]["sigmabar"]["mathematica"]["sigmabar"]
        sigma2 = op_config["tensors"]["sigma2"]["mathematica"]["sigma2"]
        sigmabar2 = op_config["tensors"]["sigmabar2"]["mathematica"]["sigmabar2"]
        # imaginary unit
        im = op_config["coefficients"]["I"]["mathematica"]["I"]

        # Indices
        if n_sigma2:
            gen_index_sl = self.possible_indices.generate_index("Lsl")  # next(gen_index_sl) will generate new Lsl and new Usl index
            gen_index_sldot = self.possible_indices.generate_index("Lsldot")

        id_statements = []
        for i in range(n_sigma2):
            lsl3, usl3 = next(gen_index_sl)
            lsldot3, usldot3 = next(gen_index_sldot)
            id = f"{sigma2}(lor1?lor,lor2?lor,UslA1?Usl,UslA2?Usl)"
            id += " = "
            id += f"({im}/2)*{sl2Ceps}(UslA1,{usl3})*("
            id += f"{sigma}(lor1,{lsl3},{lsldot3})*{sigmabar}(lor2,{usldot3},UslA2)"
            id += " - "
            id += f"{sigma}(lor2,{lsl3},{lsldot3})*{sigmabar}(lor1,{usldot3},UslA2)"
            id += ")"
            id_statements.append(id)

            lsl4, usl4 = next(gen_index_sl)
            lsldot4, usldot4 = next(gen_index_sldot)
            id = f"{sigmabar2}(lor1?lor,lor2?lor,LsldotA1?Lsldot,LsldotA2?Lsldot)"
            id += " = "
            id += f"({im}/2)*{sl2Ceps}(LsldotA1,{lsldot4})*("
            id += f"{sigmabar}(lor1,{usldot4},{usl4})*{sigma}(lor2,{lsl4},LsldotA2)"
            id += " - "
            id += f"{sigmabar}(lor2,{usldot4},{usl4})*{sigma}(lor1,{lsl4},LsldotA2)"
            id += ")"
            id_statements.append(id)


        form = ""
        for id in id_statements:
            form += f"id once {id};\n"

        return form

    ########################
    ### Gauge Conversion ###
    ########################

    def form_antifundamentalIndices(self, used_label):
        # Tensors
        su2eps = op_config["tensors"]["[su2eps]"]["mathematica"]["su2eps"]
        su3eps = op_config["tensors"]["[su3eps]"]["mathematica"]["su3eps"]
        gen = op_config["tensors"]["T"]["mathematica"]["TT"]
        # Fields
        # Fixme: Decide somehow automatically which fields transform originally in the antifundamental representation
        su2anti = [list(op_config["bosonfields"]["H"]["mathematica"].values())[1],
                   list(op_config["fermionfields"]["L"]["mathematica"].values())[1],
                   list(op_config["fermionfields"]["Q"]["mathematica"].values())[1]]
        su3anti = [list(op_config["fermionfields"]["Q"]["mathematica"].values())[1],
                   list(op_config["fermionfields"]["[u_C]"]["mathematica"].values())[0],
                   list(op_config["fermionfields"]["[d_C]"]["mathematica"].values())[0]]

        # Translation
        spinors_dirac_weyl = {name: field for name, field in op_config["fermionfields"].items() if "diracspinor" in field.keys()}
        spinors_dirac_weyl = {dirac: weyl[0] if weyl[0] else weyl[1]  for field in spinors_dirac_weyl.values() for dirac, weyl in field["diracspinor"].items()}
        # Indices
        gen_index_gauge = self.possible_indices.generate_index("gauge")
        gen_index_colf = self.possible_indices.generate_index("colf")

        form = ".sort\n"
        form += "CFunction f1;\n"

        id_statements = {}
        for field in self.fields:
            if field.name == list(op_config["bosonfields"]["H"]["mathematica"].values())[1]:
                # Higgs
                name = field.name
            else:
                if field.name not in spinors_dirac_weyl.keys(): continue
                name = spinors_dirac_weyl[field.name]
            if name in su2anti:
                cname = get_commuting_op(name)
                gauge2 = next(gen_index_gauge)

                used_label += 1
                label1 = used_label
                id_statements[label1] = []

                id = f"f1?AllcconFields(op1?op,?a,gauge1A?gauge,?b)*"
                id += f"{cname}(op2?op,?c,gauge1A?gauge,?d)"
                id += " = "
                id += f"f1(op1,?a,gauge1A,?b)*"
                id += f"{cname}(op2,?c,{gauge2},?d)*{su2eps}({gauge2},gaugeA1)"
                id_statements[label1].append(id)

                id = f"{gen}(gaugeadjA1?gaugeadj,gaugeA1?gauge,gaugeA2?gauge)*"
                id += f"{cname}(op2?op,?c,gaugeA1?gauge,?d)"
                id += " = "
                id += f"{gen}(gaugeadjA1,gaugeA1,gaugeA2)*"
                id += f"{cname}(op2,?c,{gauge2},?d)*{su2eps}({gauge2},gaugeA1)"
                id_statements[label1].append(id)

            if name in su3anti:
                cname = get_commuting_op(name)
                colf2 = next(gen_index_colf)
                colf3 = next(gen_index_colf)

                used_label += 1
                label2 = used_label
                id_statements[label2] = []

                id = f"f1?AllcconFields(op1?op,?a,colfA1?colf,?b)*"
                id += f"{cname}(op2?op,?c,colfA1?colf,?d)"
                id += " = "
                id += f"f1(op1,?a,colfA1,?b)*"
                id += f"(1/2)*{su3eps}(colfA1,{colf2},{colf3})*{cname}(op2,?c,{colf2},{colf3},?d)"
                id_statements[label2].append(id)

                id = f"{gen}(colaA1?cola,colfA1?colf,colfA2?colf)*"
                id += f"{cname}(op2?op,?c,colfA1?colf,?d)"
                id += " = "
                id += f"{gen}(colaA1,colfA1,colfA2)*"
                id += f"(1/2)*{su3eps}(colfA1,{colf2},{colf3})*{cname}(op2,?c,{colf2},{colf3},?d)"
                id_statements[label2].append(id)

        for label in id_statements.keys():
            for id in id_statements[label]:
                form += f"id once ifmatch->{label} {id};\n"
            form += f"\nlabel {label};\n\n"

        return form, used_label

    def form_adjointIndices(self):
        # Tensors
        su2eps = op_config["tensors"]["[su2eps]"]["mathematica"]["su2eps"]
        su3eps = op_config["tensors"]["[su3eps]"]["mathematica"]["su3eps"]
        gen = op_config["tensors"]["T"]["mathematica"]["TT"]
        # Fields
        # Fixme: Decide somehow automatically which fields transform originally in the adjoint representation
        su2adj = {field: op_config["bosonfields"]["WL"]["non_singlet"] for field in op_config["bosonfields"]["WL"]["mathematica"].values()}
        su3adj = {field: op_config["bosonfields"]["GL"]["non_singlet"] for field in op_config["bosonfields"]["GL"]["mathematica"].values()}

        # Indices
        gen_index_gauge = self.possible_indices.generate_index("gauge")
        gen_index_colf = self.possible_indices.generate_index("colf")

        form = ".sort\n"
        form += "CFunction f1, f2;\n"

        id_statements = []
        for field in self.fields:
            if field.name in su2adj.values():
                gauge1 = next(gen_index_gauge)
                gauge2 = next(gen_index_gauge)
                gauge3 = next(gen_index_gauge)
                for cname in [get_commuting_op(field) for field in su2adj.keys()]:
                    id = f"{cname}(op1?op,?a,gaugeadj1?gaugeadj,?b)"
                    id += " = "
                    id += f"{gen}(gaugeadj1,{gauge3},{gauge1})*{cname}(op1,?a,{gauge1},{gauge2},?b)"
                    id += f"*{su2eps}({gauge2},{gauge3})"
                    id_statements.append(id)
            elif field.name in su3adj.values():
                for cname in [get_commuting_op(field) for field in su3adj.keys()]:
                    colf1 = next(gen_index_colf)
                    colf2 = next(gen_index_colf)
                    colf3 = next(gen_index_colf)
                    colf4 = next(gen_index_colf)
                    id = f"{cname}(op1?op,?a,cola1?cola,?b)"
                    id += " = "
                    id += f"(1/2)*{gen}(cola1,{colf4},{colf2})*{su3eps}({colf1},{colf3},{colf4})"
                    id += f"*{cname}(op1,?a,{colf1},{colf2},{colf3},?b)"
                    id_statements.append(id)
            else:
                continue

        for id in id_statements:
            form += f"id once {id};\n"

        return form


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

            id_statement1 += f" - {lambdah} * {su2eps}({gauge3},{gauge2})*{higgs_dagger}({gauge3}) * {higgs}({gauge2}) * {higgs}({gauge1})"

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

            id_statement2 += f" - {lambdah} * {su2eps}({gauge3},{gauge2})*{higgs_dagger}({gauge3}) * {higgs}({gauge2}) * {higgs_dagger}({gauge1})"

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
        # TODO: Do the same for all fields not just the Higgs.

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

            #  Left-handed field strength tensor: FL(Lsl1, Lsl2, H(gauge1))
            id_statement1 = f"FL({lsl1}?Lsl, {lsl2}?Lsl, {higgs}({gauge1}?gauge)) ="
            id_statement1 += f" + ({g1}/2) * {bl}({lsl1},{lsl2}) * {higgs}({gauge1})"
            id_statement1 += f" + ({g2}/2) * {wl}({gauge1},{gauge2},{lsl1},{lsl2}) * {su2eps}({gauge2},{gauge3}) * {higgs}({gauge3})"

            id_statements.append(f"id once ifmatch -> 4 {id_statement1:s};\n")

            #  Left-handed field strength tensor: FL(Lsl1, Lsl2, [H+](gauge1))
            id_statement2 = f"FL({lsl1}?Lsl, {lsl2}?Lsl, {higgs_dagger}({gauge1}?gauge)) ="
            id_statement2 += f" + (- {g1}/2) * {bl}({lsl1},{lsl2}) * {higgs_dagger}({gauge1})"
            id_statement2 += f" + ({g2}/2) * {wl}({gauge1},{gauge2},{lsl1},{lsl2}) * {su2eps}({gauge2},{gauge3}) * {higgs_dagger}({gauge3})"

            id_statements.append(f"id once ifmatch -> 4 {id_statement2:s};\n")

            #  Right-handed field strength tensor: FR(Usldot1, Usldot2, H(gauge1))
            id_statement1 = f"FR({usldot1}?Usldot, {usldot2}?Usldot, {higgs}({gauge1}?gauge)) ="
            id_statement1 += f" + ({g1}/2) * {br}({usldot1},{usldot2}) * {higgs}({gauge1})"
            id_statement1 += f" + ({g2}/2) * {wr}({gauge1},{gauge2},{usldot1},{usldot2}) * {su2eps}({gauge2},{gauge3}) * {higgs}({gauge3})"

            id_statements.append(f"id once ifmatch -> 4 {id_statement1:s};\n")

            #  Right-handed field strength tensor: FR(Usldot1, Usldot2, [H +](gauge1))
            id_statement2 = f"FR({usldot1}?Usldot, {usldot2}?Usldot, {higgs_dagger}({gauge1}?gauge)) ="
            id_statement2 += f" + (- {g1}/2) * {br}({usldot1},{usldot2}) * {higgs_dagger}({gauge1})"
            id_statement2 += f" + ({g2}/2) * {wr}({gauge1},{gauge2},{usldot1},{usldot2}) * {su2eps}({gauge2},{gauge3}) * {higgs_dagger}({gauge3})"

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
            id_statement2 += f" - i_ * {yu_dagger}({flav2},{flav1}) * {u}({lsl1},{colf2},{colf3},{flav2}) * {higgs}({gauge1})"
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
            id_statement1 += f" + i_ * {yd_dagger}({flav1},{flav2}) * {su2eps}({gauge2},{gauge1}) * {higgs_dagger}({gauge2}) * {q}({lsl1},{gauge1},{colf1},{flav2})"
            id_statements.append(f"id once ifmatch -> 4 {id_statement1:s};\n")

            #  EOM of [d_C]
            id_statement2 = f"EOM({d},{usldot1}?Usldot,{colf2}?colf,{colf3}?colf,{flav1}?flav) ="
            id_statement2 += f" + i_ * {yd}({flav2},{flav1}) * {su2eps}({gauge2},{gauge1}) * {q_dagger}({usldot1},{gauge2},{colf2},{colf3},{flav2}) * {higgs}({gauge1})"
            id_statements.append(f"id once ifmatch -> 4 {id_statement2:s};\n")
            id_statements.append("\n")

            #  EOM of [u_C+]
            id_statement1 = f"EOM({u_dagger},{lsl1}?Lsl,{colf1}?colf,{flav1}?flav) ="
            id_statement1 += f" + (- i_) * {yu_dagger}({flav1},{flav2}) * {su2eps}({gauge1},{gauge2}) * {higgs}({gauge1}) * {q}({lsl1},{gauge2},{colf1},{flav2})"
            id_statements.append(f"id once ifmatch -> 4 {id_statement1:s};\n")

            #  EOM of [u_C]
            id_statement2 = f"EOM({u},{usldot1}?Usldot,{colf2}?colf,{colf3}?colf,{flav1}?flav) ="
            id_statement2 += f" + i_ * {yu}({flav2},{flav1}) * {su2eps}({gauge2},{gauge1}) * {q_dagger}({usldot1},{gauge2},{colf2},{colf3},{flav2}) * {higgs_dagger}({gauge1})"
            id_statements.append(f"id once ifmatch -> 4 {id_statement2:s};\n")
            id_statements.append("\n")

            #  EOM of [e_C+]
            id_statement1 = f"EOM({e_dagger},{lsl1}?Lsl,{flav1}?flav) ="
            id_statement1 += f" + i_ * {ye_dagger}({flav1},{flav2}) * {su2eps}({gauge2},{gauge1}) * {higgs_dagger}({gauge2}) * {l}({lsl1},{gauge1},{flav2})"
            id_statements.append(f"id once ifmatch -> 4 {id_statement1:s};\n")

            #  EOM of [e_C]
            id_statement2 = f"EOM({e},{usldot1}?Usldot,{flav1}?flav) ="
            id_statement2 += f" + i_ * {ye}({flav2},{flav1}) * {su2eps}({gauge2},{gauge1}) * {l_dagger}({usldot1},{gauge2},{flav2}) * {higgs}({gauge1})"
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
        flavdK = op_config["tensors"]["[flavdK]"]["mathematica"]["flavdK"]
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
            flav1 = next(gen_index_flav)  # * {flavdK}({flav1},{flav2})
            flav2 = next(gen_index_flav)
            lsl1, usl1 = next(gen_index_sl)
            lsldot1, usldot1 = next(gen_index_sldot)
            #  EOM of the left-handed field strength tensor of U(1), BL
            id_statement1 = f"EOM({bl},{lsl1}?Lsl,{usldot1}?Usldot) ="
            rHS1 = f" {g1}*("

            rHS1 += f" + (i_/2) * {su2eps}({gauge2},{gauge1}) * ( {cov}({lsl1},{usldot1},{higgs}({gauge1})) * {higgs_dagger}({gauge2})"
            rHS1 += f" - {higgs}({gauge1}) * {cov}({lsl1},{usldot1},{higgs_dagger}({gauge2})) )"

            rHS1 += f" - (1/6) * {su3eps}({colf1},{colf2},{colf3}) * {su2eps}({gauge2},{gauge1}) * {flavdK}({flav1},{flav2}) * {q}({lsl1},{gauge1},{colf1},{flav1}) * {q_dagger}({usldot1},{gauge2},{colf2},{colf3},{flav2})"

            rHS1 += f" + (2/3) * {su3eps}({colf1},{colf2},{colf3}) * {flavdK}({flav1},{flav2}) * {u}({lsl1},{colf2},{colf3},{flav1}) * {u_dagger}({usldot1},{colf1},{flav2})"

            rHS1 += f" - (1/3) * {su3eps}({colf1},{colf2},{colf3}) * {flavdK}({flav1},{flav2}) * {d}({lsl1},{colf2},{colf3},{flav1}) * {d_dagger}({usldot1},{colf1},{flav2})"

            rHS1 += f" + {su2eps}({gauge2},{gauge1}) * {flavdK}({flav1},{flav2}) * {l}({lsl1},{gauge1},{flav1}) * {l_dagger}({usldot1},{gauge2},{flav2})"

            rHS1 += f" + 2 * {flavdK}({flav1},{flav2}) * {e}({lsl1},{flav1}) * {e_dagger}({usldot1},{flav2})"

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
            rHS2 = f" - ({g2}/2) * ("  # Factor 1/2 necessary for explicit symmetrization of SU2 indices.
            def gauge_indices(gaugei, gaugek):
                """For the symmetrization of the SU2-indices."""
                output = f" + i_ * ( {cov}({lsl1},{usldot1},{higgs}({gaugei})) * {higgs_dagger}({gaugek})"
                output += f" - {higgs}({gaugei}) * {cov}({lsl1},{usldot1},{higgs_dagger}({gaugek})) )"
                output += f" - (1/3) * {su3eps}({colf1},{colf2},{colf3}) * {flavdK}({flav1},{flav2}) * {q}({lsl1},{gaugei},{colf1},{flav1}) * {q_dagger}({usldot1},{gaugek},{colf2},{colf3},{flav2})"
                output += f" + {flavdK}({flav1},{flav2}) * {l}({lsl1},{gaugei},{flav1}) * {l_dagger}({usldot1},{gaugek},{flav2})"

                # output += f" - (1/3) * {su3eps}({colf1},{colf2},{colf3}) * {q}({lsl1},{gaugek},{colf1},{flav1}) * {q_dagger}({usldot1},{gaugei},{colf2},{colf3},{flav1})"
                # output += f" + {l}({lsl1},{gaugek},{flav1}) * {l_dagger}({usldot1},{gaugei},{flav1})"

                return output

            # 'Normal' order of SU2 indices:
            rHS2 += gauge_indices(gauge1, gauge2)
            # SU2 indices swapped:
            rHS2 += gauge_indices(gauge2, gauge1)

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

    ############################
    ### SL2C epsilon tensors ###
    ############################
    @property
    def lr(self):
        def set_lp_in_tensors(tensors: Tensors, ref_index: Index, lp: LP_Index, derIndex: Union[bool, int]):
            """
            Set in the index in tensors which is conjugated with the index 'ref_index', the 'lp' attribute to lp
            and the derIndex attribute to 'derIndex'.
            Parameters
            ----------
            tensors
            ref_index
            lp
            derIndex

            Returns
            -------
                True if index successfully replace - otherwise False.
            """
            found = False
            for tensor in tensors:
                for sl_index in tensor.indices["sl"]:
                    if sl_index == ref_index.dual_index:
                        found = True
                        sl_index.lp = lp
                        # derivative index has to be set in order to be able to infer later on the derivative structure
                        # only from the tensors.
                        sl_index.derIndex = derIndex
                        break
                for sldot_index in tensor.indices["sldot"]:
                    if sldot_index == ref_index.dual_index:
                        found = True
                        sldot_index.lp = lp
                        sldot_index.derIndex = derIndex
                        break
                if found: break

            if found:
                return True
            else:
                return False

        for field in self.fields:
            fp = field.field_pos
            for sl_index in field.indices["sl"]:
                lp = LP_Index(fp, sl_index.derIndex)
                sl_index.lp = lp
                # equalize all properties in the tensor indices
                assert set_lp_in_tensors(self.tensors, sl_index, lp, sl_index.derIndex)
            for sldot_index in field.indices["sldot"]:
                lp = LP_Index(fp, sldot_index.derIndex)
                sldot_index.lp = lp
                assert set_lp_in_tensors(self.tensors, sldot_index, lp, sldot_index.derIndex)

        # Infer from field and derivative structure the l_tab and r_tab, specifying the epsilon tensors:
        # First, get SL2C-eps tensors for undotted and dotted indices:
        eps_lh = self.tensors["sl"]
        eps_rh = self.tensors["sldot"]
        if eps_lh:
            eps0 = eps_lh[0]  # first epsilon tensor
            l_tab = Young_Tableau([[eps0.indices[0].copy()], [eps0.indices[1].copy()]])
            # Note: The copy is necessary in order get new indices, which are not linked to the old Summand indices and can thus be changed safely.
            for eps in eps_lh[1:]:
                # Append a column for the next epsilon tensor
                l_tab.append_col(Young_Tableau([[eps.indices[0].copy()], [eps.indices[1].copy()]]))

        else:
            l_tab = Young_Tableau([])
        if eps_rh:
            eps0 = eps_rh[0]  # first epsilon tensor
            r_tab = Young_Tableau([[eps0.indices[0].copy()], [eps0.indices[1].copy()]])
            for eps in eps_rh[1:]:
                # Append a column for the next epsilon tensor
                r_tab.append_col(Young_Tableau([[eps.indices[0].copy()], [eps.indices[1].copy()]]))
        else:
            r_tab = Young_Tableau([])

        # print(f"{eps_lh}\n->\n{l_tab:nice}\n{l_tab:lp}")
        # print(f"{eps_rh}\n->\n{r_tab:nice}\n{r_tab:lp}")
        # print("=========================================================================")
        return LR_Tableaux(l_tab, r_tab, 1)

    def get_term_from_lr_tabs(self, lr: LR_Tableaux, ignore_op_class: bool = False):
        """
        Infer from the l_tab and r_tab the derivative structure for a given field structure.
        For this purpose, the field content and the other indices are necessary, which is why the term with a predominantly
        incorrect derivative structure must also be specified.
        -> In a first step, all SL2C-indices and thus also all derivatives of the fields are removed and also the
           SL2C-epsilon tensors are removed.
        -> Since each row in the l_tab and r_tab specifies uniquely an epsilon tensor and those specify uniquely the indices
           of themselves and also those of the field.
        Parameters
        ----------
        term
        lr = LR_tableaux(l_tab, r_tab, factor)
            The lr tableaux have to be given explicitely, because a new term with the SAME fields whould be constructed
            for the GIVEN tableaux.
        ignore_op_class
            If True, the operator class of the Summand will not be checked, i.e. it is possible to have more derivatives then allowed by the operator class.

        Returns
        -------

        """
        if not ignore_op_class:
            assert lr.l_tab.ncols() == self.op_class.nl, "The number of columns in l_tab doesn't match the the one necessary for the operator class."
            assert lr.r_tab.ncols() == self.op_class.nr, "The number of columns in r_tab doesn't match the the one necessary for the operator class."

        def get_eps_from_tab(tab: Young_Tableau):
            sl2Ceps = op_config["tensors"]["[sl2Ceps]"]["mathematica"]["sl2Ceps"]  # '[sl2Ceps]'
            epss = []
            for i in range(tab.ncols()):
                col = tab[:, i]
                tensor = Tensor(f"{sl2Ceps}({col[0, 0]},{col[1, 0]})")
                for j, index in enumerate(tensor.indices):
                    index.lp = col[j, 0].lp
                    index.derIndex = col[j, 0].derIndex

                epss.append(tensor)

            dual_indices = []
            for eps in epss:
                for index in eps.indices:
                    dual_indices.append(index.dual_index)

            return Tensors(epss), Indices_Operator(dual_indices)

        eps_undotted, lsl_indices = get_eps_from_tab(lr.l_tab)
        eps_dotted, usldot_indices = get_eps_from_tab(lr.r_tab)
        # sort indices by fields:
        field_indices = {i: [] for i in range(1, self.op_class.N + 1)}

        for index in lsl_indices:
            field_indices[index.lp.fp].append(index)
        for index in usldot_indices:
            field_indices[index.lp.fp].append(index)

        #New (identical) Summand object:
        summand = self.copy()

        # Delete old SL2C-tensors
        del summand.tensors["sl2C"]
        # Append new SL2C-tensors
        for eps in eps_undotted:
            summand.tensors.append(eps)
        for eps in eps_dotted:
            summand.tensors.append(eps)

        def delete_sl2c_index(indices: Indices_Operator):
            """
            Recursively delete all SL2C-indices in 'indices'.
            Parameters
            ----------
            indices

            Returns
            -------

            """
            sl2c = [index.typ in ("Lsl", "Usl", "Lsldot", "Usldot") for index in indices]
            if any(sl2c):
                for i, index in enumerate(indices):
                    if index.typ in ("Lsl", "Usl", "Lsldot", "Usldot"):
                        del field.indices[i]
                        break
                return delete_sl2c_index(indices)
            else:
                return indices

        # Delete old SL2C-indices
        for field in summand.fields:
            delete_sl2c_index(field.indices)
        # Append new SL2C-indices
        for field in summand.fields:
            for new_index in field_indices[field.field_pos]:
                field.indices.append(new_index)
            # Recalculate the number of derivatives 'derIndex', which stand on every index.
            # field.reset_derIndex()

        # Adjust the coefficient
        summand.coeff *= Factor(lr.factor)

        return summand



