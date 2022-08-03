import logging.config
import pickle  # For development purposes only, to save time in debugging.
import timeit
from argparse import ArgumentParser
from datetime import datetime
from pathlib import Path

from yaml import safe_load

from tioc import CONFIG_PATH, FORM_GENERAL_PATH, n_der, get_basis
from tioc.general import form_declarations
from tioc.lor_projection import rearrange_derivatives, replace_eoms, ibp_and_schouten_ids
from tioc.projection import converttoSL2C, tex_unsorted_terms, tex_sorted_terms, tex_sorted_terms_wo_doubles, \
    tex_terms_sorted_sun_projection  # , main
from tioc.sun_projection import get_type, remove_doubles, sun_projection, replace_sun_tensors_by_projected_ones
from tioc.tableau import test

# configure logger
timestamp = datetime.now()
logger_autoeft = logging.getLogger("autoeft")
logger = logging.getLogger("autoeft.projection")
with open(CONFIG_PATH / "logger.yml", "r") as file:
    logconfig = safe_load(file)

parser = ArgumentParser(
    description="Project BSUOLEA output onto AUTOEFT-basis.",
)
parser.add_argument("matched", type=Path, help="path to file with matched lagrangian")
parser.add_argument("-s", "--skip", type=int, default=0, help="number of lines to skip in the matched file")
parser.add_argument(
    "-t", "--tex", action="store_true", help="save supplementary tex files"
)
parser.add_argument("-d", "--dimension", type=int, default=6, help="maximum mass dimension for the projection")
parser.add_argument(
    "-l",
    "--log",
    dest="logLevel",
    choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
    default="INFO",
    help="set the logging level (default: %(default)s)",
)
args = parser.parse_args()

logconfig['handlers']['console']['level'] = logging.getLevelName(args.logLevel)
logging.config.dictConfig(logconfig)

logger.info(timestamp.replace(microsecond=0).isoformat())

def number_terms(single_terms):
    nterms = 0
    for terms in single_terms.values():
        for terms_n in terms.values():
            nterms += len(terms_n)
    return nterms

def main(max_dim = 6):
    # get basis up to max_dim mass dimension
    basis = get_basis(max_dim)
    # bs_file = args.matched.resolve()  # "exampleOutputBS.m"
    #
    # with open(FORM_GENERAL_PATH / "declarations_general.h", "w") as file:
    #     file.write(form_declarations(n_der))
    #
    # terms = converttoSL2C(bs_file, header=args.skip, pprint=False)
    # del bs_file
    # # pickle.dump(terms, open(CONFIG_PATH / "terms.p", "wb"))
    # # terms = pickle.load(open(CONFIG_PATH / "terms.p", "rb"))
    #
    # if args.tex:
    #     tex_unsorted_terms(terms)
    #
    # # Sort terms by type for the projection:
    # single_terms = get_type(terms)
    # del terms
    #
    # if args.tex:
    #     tex_sorted_terms(single_terms)
    #
    # nterms_before = number_terms(single_terms)
    # # Remove double terms:
    # single_terms = remove_doubles(single_terms)
    #
    # nterms_after = number_terms(single_terms)
    #
    # logger.info(f"#Terms with doubles: {nterms_before:d} <-> #Terms without doubles: {nterms_after:d}")
    #
    # # print terms
    # if args.tex:
    #     tex_sorted_terms_wo_doubles(single_terms)
    #
    # pickle.dump(single_terms, open(CONFIG_PATH / "single_terms_with_eoms.p", "wb"))
    # # single_terms = pickle.load(open(CONFIG_PATH / "single_terms_with_eoms.p", "rb"))
    #
    # nterms_before = number_terms(single_terms)
    #
    # single_terms = replace_eoms(single_terms)
    # single_terms = replace_eoms(single_terms)
    #
    #
    # nterms_after = number_terms(single_terms)
    #
    # logger.info(f"#Terms with EOMs: {nterms_before:d} <-> #Terms with less EOMs: {nterms_after:d}")
    #
    # pickle.dump(single_terms, open(CONFIG_PATH / "single_terms_with_less_eoms.p", "wb"))
    # # single_terms = pickle.load(open(CONFIG_PATH / "single_terms_with_less_eoms.p", "rb"))
    #
    # # rearrange derivative on term of type: nD: 4 & {"H": 1, "H+": 1}
    # single_terms = rearrange_derivatives(single_terms)
    #
    # single_terms = replace_eoms(single_terms)
    #
    # # pickle.dump(single_terms, open(CONFIG_PATH / "single_terms.p", "wb"))
    # # single_terms = pickle.load(open(CONFIG_PATH / "single_terms.p", "rb"))
    #
    # nterms_before = number_terms(single_terms)
    # # Remove double terms:
    # single_terms = remove_doubles(single_terms)
    #
    # nterms_after = number_terms(single_terms)
    #
    # logger.info(f"#Terms with doubles: {nterms_before:d} <-> #Terms without doubles: {nterms_after:d}")

    # pickle.dump(single_terms, open(CONFIG_PATH / "single_terms.p", "wb"))
    single_terms = pickle.load(open(CONFIG_PATH / "single_terms.p", "rb"))

    # TODO: Apply the integration by parts and Schouten identities to the lorentz structure.
    # test()
    single_terms = ibp_and_schouten_ids(single_terms)

    # SUN_Projection of terms without doubles:
    single_terms = sun_projection(single_terms, basis, max_dim)

    single_terms = replace_sun_tensors_by_projected_ones(single_terms)

    nterms_after_sun_projection = 0
    for terms in single_terms.values():
        for terms_n in terms.values():
            nterms_after_sun_projection += len(terms_n)

    logger.info(f"Terms after SUN projection: {nterms_after_sun_projection:d}")

    if args.tex:
        tex_terms_sorted_sun_projection(single_terms)

    # pickle.dump(single_terms, open(CONFIG_PATH / "single_terms.p", "wb"))
    # single_terms = pickle.load(open(CONFIG_PATH / "single_terms.p", "rb"))


    return single_terms


start_time = timeit.default_timer()

stuff = main(args.dimension)

stop_time = timeit.default_timer()
logger.info(f"Done in {stop_time - start_time:.2f} sec.")

