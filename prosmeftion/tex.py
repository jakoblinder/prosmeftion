# -*- coding: utf-8 -*-
import os
import subprocess
import sys
import re
import logging.config
import multiprocessing as mp

from itertools import chain
from pathlib import Path
from typing import List

from prosmeftion.sl2c.class_term import Term
from . import coeffvalues, opname_sorted, opname, opnameSL2C, opvalues, opSL2Cvalues, spinorsSL2C_c, spSL2C_c_values
from . import PROJECTION_PATH, CONFIG_PATH, FORM_PATH, INPUT_PATH, LATEX_PATH, AUTOEFT_PATH
from .yProjection.read_write import get_terms
from .yProjection.coefficient import Factor
from .sun_projection import equalize_field_indices

logger_autoeft = logging.getLogger("autoeft.projection")
logger = logger_autoeft.getChild(__name__)

def write_texfile_and_create_pdf(pdfname="terms"):
    """
    Define Decorator to generat tex files of given latex expression and name of pdf.
    Parameters
    ----------
    pdfname: str
        name of the pdf without postfix.
    Returns
    -------
    """
    def decorator(func):
        def wrapper_with_func_args(*args, **kwargs):
            pdf = f"{args[0]:s}_{pdfname:s}"
            logger.info(f"Create texed PDF {pdf:s}.pdf of all terms.")
            logger.debug("Create Tex file of all terms.")

            # call wrapped function:
            latex = func(*args, **kwargs)
            #  args[0]: filename
            #  args[1]: list of terms

            with open(LATEX_PATH / "terms_all.tex", "w") as file:
                file.write(latex)
            try:
                logger.debug("Construct PDF.")
                # print(subprocess.list2cmdline(["pdflatex", f"-output-directory={LATEX_PATH}", LATEX_PATH / "terms.tex"]))
                subprocess.run(["pdflatex", f"-output-directory={LATEX_PATH}", LATEX_PATH / "terms.tex"],
                               capture_output=True, text=True, check=True)
            except subprocess.CalledProcessError as exc:
                exc.cmd = list(map(str, list(exc.cmd)))
                logger.error(" ".join(exc.cmd) + "\n" + str(exc.stdout))
                logger.error(f"Error {exc.returncode}")
                sys.exit("STOP")
            else:
                # Move created pdf to the main projection folder.
                file_path = LATEX_PATH / "terms.pdf"
                file_path.rename(PROJECTION_PATH / f"{pdf:s}.pdf")
        return wrapper_with_func_args

    return decorator

@write_texfile_and_create_pdf("irreps")
def tex_unsorted_terms(inputfilename, terms):
    latex = ""

    for term in terms:
        latex += r"\paragraph{" + f"{term.name:s}" + "}\n"
        for summand in term.terms:
            latex += r"\begin{align}" + "\n"
            latex += "\t" + r"\begin{autobreak}" + "\n"
            latex += f"{summand:tex}" + "\n"
            latex += "\t" + r"\end{autobreak}" + "\n"
            latex += r"\end{align}" + "\n"

    return latex

@write_texfile_and_create_pdf("sorted")
def tex_sorted_terms(inputfilename, single_terms):
    single_terms = equalize_field_indices(single_terms)
    latex = ""
    for term_type in single_terms.values():
        for term_type_nD in term_type.values():
            latex += r"\section*{" + re.sub(r"'", "", term_type_nD.name) + "}\n"
            for term in term_type_nD:
                latex += r"\paragraph{" + f"{term.name:s}" + "}\n"
                latex += r"\begin{align}" + "\n"
                latex += "\t" + r"\begin{autobreak}" + "\n"
                latex += f"{term:tex}\n"
                latex += "\t" + r"\end{autobreak}" + "\n"
                latex += r"\end{align}" + "\n"

    return latex

@write_texfile_and_create_pdf("sorted_wo_doubles")
def tex_sorted_terms_wo_doubles(inputfilename, single_terms):
    latex = ""
    for term_type in single_terms.values():
        for term_type_nD in term_type.values():
            latex += r"\section*{" + re.sub(r"'", "", term_type_nD.name) + "}\n"
            for term in term_type_nD:
                latex += r"\paragraph{" + f"{term.name:s}" + "}\n"
                latex += r"\begin{align}" + "\n"
                latex += "\t" + r"\begin{autobreak}" + "\n"
                latex += f"{term:tex}\n"
                latex += "\t" + r"\end{autobreak}" + "\n"
                latex += r"\end{align}" + "\n"

    return latex

@write_texfile_and_create_pdf("before_sun_projection")
def tex_sorted_terms_before_sun(inputfilename, single_terms):
    latex = ""
    for term_type in single_terms.values():
        for term_type_nD in term_type.values():
            latex += r"\section*{" + re.sub(r"'", "", term_type_nD.name) + "}\n"
            for term in term_type_nD:
                latex += r"\paragraph{" + f"{term.name:s}" + "}\n"
                latex += r"\begin{align}" + "\n"
                latex += "\t" + r"\begin{autobreak}" + "\n"
                latex += f"{term:tex}\n"
                latex += "\t" + r"\end{autobreak}" + "\n"
                latex += r"\end{align}" + "\n"

    return latex

@write_texfile_and_create_pdf("sorted_wo_eoms")
def tex_sorted_terms_wo_eoms(inputfilename, single_terms):
    latex = ""
    for term_type in single_terms.values():
        for term_type_nD in term_type.values():
            latex += r"\section*{" + re.sub(r"'", "", term_type_nD.name) + "}\n"
            for term in term_type_nD:
                latex += r"\paragraph{" + f"{term.name:s}" + "}\n"
                latex += r"\begin{align}" + "\n"
                latex += "\t" + r"\begin{autobreak}" + "\n"
                latex += f"{term:tex}\n"
                latex += "\t" + r"\end{autobreak}" + "\n"
                latex += r"\end{align}" + "\n"

    return latex

@write_texfile_and_create_pdf("yBasis")
def tex_terms_sorted_sun_projection(inputfilename, single_terms, max_dim:int=6):
    latex = ""
    for mass_dim in range(2, max_dim + 1):
        if any([True if term_mass_dim.d == mass_dim else False for term_type in single_terms.values() for term_mass_dim in term_type.values()]):
            latex += r"\section*{" + f"Mass dimension: {mass_dim:d}" + "}\n"
        for term_type in single_terms.values():
            for term_mass_dim in term_type.values():
                if term_mass_dim.d is not mass_dim:
                    continue
                latex += r"\subsection*{" + re.sub(r"'", "", term_mass_dim.name) + "}\n"
                # get terms with specific field structure
                term_with_specific_field_structure = {}
                for term in term_mass_dim:
                    try:
                        term_with_specific_field_structure[term.fieldstructure].append(term)
                    except KeyError:
                        term_with_specific_field_structure[term.fieldstructure] = [term]

                latex_terms = ""
                for name_of_term, terms_specific in term_with_specific_field_structure.items():
                    name_form = "".join([f"{name}{nD}" for name, nD in name_of_term])
                    latex_terms += r"\paragraph{" + f"{name_form:s}" + "}\n"
                    for term in terms_specific:
                        latex_terms += r"\begin{align}" + "\n"
                        latex_terms += "\t" + r"\begin{autobreak}" + "\n"
                        latex_terms += f"{term:tex}\n"
                        latex_terms += "\t" + r"\end{autobreak}" + "\n"
                        latex_terms += r"\end{align}" + "\n"

                latex_sun_tensors = ""
                if term_mass_dim.sun_projection_tensors:
                    for sun_group, sun_tensors in term_mass_dim.sun_projection_tensors.items():
                        if not sun_tensors: continue
                        # get tex indices for the tensor indices
                        # Note 1: This has to be done after the generation of the tex expression of the term even if
                        # the latter is texed after the SUN-Tensors, since in this generation the tex indices for each
                        # operator are assigned.
                        # Note 2: Since the indices on all fields are always named in ascending order, exactly as it is
                        # the case for the SUN-Tensors. The SUN-Tensor indices should have all assigned the same LaTex
                        # indices and can thus be written one time, for all terms.
                        if "2" in sun_group:
                            index_name = "gauge"
                        elif "3" in sun_group:
                            index_name = "colf"
                        ref_indices = {index.expr: index.tex for index in term_mass_dim[0].indices[index_name]}
                        for term in term_mass_dim[1:]:
                            for index in term.indices[index_name]:
                                assert ref_indices[index.expr] == index.tex, f"Tex indices for indices of type {index_name} are not the same in all terms."

                        for basis_tensor in sun_tensors:
                            for monom in basis_tensor.monoms:
                                for tensor in monom.tensors:
                                    for index in tensor.indices:
                                        index.tex = ref_indices[index.expr]
                        latex_sun_tensors += r"\paragraph{" + "".join((map(lambda x: rf"\{x}" if x == "_" else x, f"{sun_group}"))) + "-Basis Tensors" + "}\n"  # escape '_' in paragraph.
                        latex_sun_tensors += r"\begin{align}" + "\n"
                        latex_sun_tensors += f"{sun_tensors:tex}"
                        latex_sun_tensors += r"\end{align}" + "\n"

                if latex_sun_tensors: latex += latex_sun_tensors
                latex += latex_terms
    return latex
