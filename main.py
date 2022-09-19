import logging.config
import pickle  # For development purposes only, to save time in debugging.
import timeit
from argparse import ArgumentParser
from datetime import datetime
from pathlib import Path

from yaml import safe_load

from tioc import CONFIG_PATH, FORM_GENERAL_PATH, PROJECTION_PATH, n_der, get_basis, run_form
from tioc.general import form_declarations

from tioc.lor_projection import ibp_and_schouten_ids, replace_eoms, rearrange_derivatives
from tioc.converttoSL2C import converttoSL2C
from tioc.tex import tex_unsorted_terms, tex_sorted_terms, tex_sorted_terms_wo_doubles, tex_terms_sorted_sun_projection, tex_sorted_terms_wo_eoms, tex_sorted_terms_before_sun
from tioc.sun_projection import get_type, remove_doubles, sun_projection, replace_sun_tensors_by_projected_ones
from tioc.get_FORM_refactored.tableau import test
from tioc.get_FORM_refactored.indices import Possible_Indices
from tioc.get_FORM_refactored.rfr import rfr
from tioc.general import declaration_SL2C_sets

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

def main(max_dim = 6, debug=None):
    # get basis up to max_dim mass dimension
    basis = get_basis(max_dim)
    bs_file = args.matched.resolve()  # "exampleOutputBS.m"

    with open(FORM_GENERAL_PATH / "declarations_general.h", "w") as file:
        file.write(form_declarations(n_der))

    terms = converttoSL2C(bs_file, header=args.skip, pprint=False)
    del bs_file
    # pickle.dump(terms, open(CONFIG_PATH / "terms.p", "wb"))
    # terms = pickle.load(open(CONFIG_PATH / "terms.p", "rb"))

    if args.tex:
        tex_unsorted_terms(terms)

    # Sort terms by type for the projection:
    single_terms = get_type(terms)
    del terms

    if args.tex:
        tex_sorted_terms(single_terms)

    nterms_before = number_terms(single_terms)
    # Remove double terms:
    single_terms = remove_doubles(single_terms)

    nterms_after = number_terms(single_terms)

    logger.info(f"#Terms with doubles: {nterms_before:d} <-> #Terms without doubles: {nterms_after:d}")

    # print terms
    if args.tex:
        tex_sorted_terms_wo_doubles(single_terms)

    pickle.dump(single_terms, open(CONFIG_PATH / "single_terms_with_eoms.p", "wb"))
    # single_terms = pickle.load(open(CONFIG_PATH / "single_terms_with_eoms.p", "rb"))

    nterms_before = number_terms(single_terms)

    single_terms = replace_eoms(single_terms)
    single_terms = replace_eoms(single_terms)

    nterms_after = number_terms(single_terms)

    logger.info(f"#Terms with EOMs: {nterms_before:d} <-> #Terms with less EOMs: {nterms_after:d}")

    pickle.dump(single_terms, open(CONFIG_PATH / "single_terms_with_less_eoms.p", "wb"))
    # single_terms = pickle.load(open(CONFIG_PATH / "single_terms_with_less_eoms.p", "rb"))

    # rearrange derivative on term of type: nD: 4 & {"H": 1, "H+": 1}
    single_terms = rearrange_derivatives(single_terms)

    single_terms = replace_eoms(single_terms)  # (D^2H) * (D^2H+) -> ... + ~ H+ * H * H * (D^2H+) + ... => need to replace eom again.
    single_terms = replace_eoms(single_terms)

    # pickle.dump(single_terms, open(CONFIG_PATH / "single_terms.p", "wb"))
    # single_terms = pickle.load(open(CONFIG_PATH / "single_terms.p", "rb"))

    nterms_before = number_terms(single_terms)
    # Remove double terms:
    single_terms = remove_doubles(single_terms)

    nterms_after = number_terms(single_terms)

    logger.info(f"#Terms with doubles: {nterms_before:d} <-> #Terms without doubles: {nterms_after:d}")

    if args.tex:
        tex_sorted_terms_wo_eoms(single_terms)

    pickle.dump(single_terms, open(CONFIG_PATH / "single_terms.p", "wb"))
    # single_terms = pickle.load(open(CONFIG_PATH / "single_terms.p", "rb"))

    logger.info("Symmetries terms.")
    single_terms = rfr(single_terms)

    nterms_before = number_terms(single_terms)
    single_terms = remove_doubles(single_terms)
    nterms_after = number_terms(single_terms)
    logger.info(f"#Terms with doubles: {nterms_before:d} <-> #Terms without doubles: {nterms_after:d}")
    pickle.dump(single_terms, open(CONFIG_PATH / "single_terms_without_rfr.p", "wb"))
    # single_terms = pickle.load(open(CONFIG_PATH / "single_terms_without_rfr.p", "rb"))

    # test()
    # Apply ibp and schouten ids up to the point where all tableaus, are SSYT:
    # Note due to the replacement of contracted derivative on the second and third field,
    # this has to be done iteratively while removing always the EOMs before the next iteration
    if debug:
        single_terms = {debug: single_terms[debug]}

    number_iterations = 0

    while True:
        single_terms, ssyt = ibp_and_schouten_ids(single_terms, max_dim)
        number_iterations += 1
        single_terms = replace_eoms(single_terms)
        # pickle.dump(single_terms, open(CONFIG_PATH / "single_terms_test.p", "wb"))
        # single_terms = pickle.load(open(CONFIG_PATH / "single_terms_test.p", "rb"))
        # Remove double terms:
        single_terms = remove_doubles(single_terms)
        nterms_after = number_terms(single_terms)
        logger.info(f"#Terms after {number_iterations}. ibp iteration {nterms_after:d}.")
        if ssyt:
            break

    logger.info(f"After {number_iterations:d} iterations of the ibp algorithm all ibp relations and Schouten identities are applied.")

    if args.tex:
        tex_sorted_terms_before_sun(single_terms)
    # SUN_Projection of terms without doubles:
    single_terms = sun_projection(single_terms, basis, max_dim)

    single_terms = replace_sun_tensors_by_projected_ones(single_terms)

    nterms_after_sun_projection = number_terms(single_terms)

    logger.info(f"Terms after SUN projection: {nterms_after_sun_projection:d}")

    if args.tex:
        tex_terms_sorted_sun_projection(single_terms,max_dim)

    # pickle.dump(single_terms, open(CONFIG_PATH / "single_terms.p", "wb"))
    # single_terms = pickle.load(open(CONFIG_PATH / "single_terms.p", "rb"))


    return single_terms


start_time = timeit.default_timer()

terms = main(args.dimension)
# fieldstructure = tuple(map(int,list("0100000110000000")))
# terms = main(args.dimension, debug=fieldstructure)

output = ""
for term_type in terms.values():
    for term_mass_dim in term_type.values():
        for term in term_mass_dim:
            output += f"+{term:c}\n"
n_terms = number_terms(terms)

with open(PROJECTION_PATH / f"all_terms{n_terms}.h", "w") as file:
    file.write(output)

stop_time = timeit.default_timer()
logger.info(f"Done in {stop_time - start_time:.2f} sec.")

