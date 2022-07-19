import logging.config
import logging.config
import timeit
import pickle  # For development purposes only, to save time in debugging.
from argparse import ArgumentParser
from datetime import datetime
from pathlib import Path

from yaml import safe_load

from tioc.projection import converttoSL2C, tex_unsorted_terms, tex_sorted_terms, tex_sorted_terms_wo_doubles, tex_terms_sorted_sun_projection  # , main
from tioc.sun_projection import get_type, remove_doubles, sun_projection, replace_sun_tensors_by_projected_ones
from tioc import form_declarations, declaration_SL2C_sets, CONFIG_PATH, FORM_GENERAL_PATH, op_config

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


def main():
    # bs_file = args.matched.resolve()  # "exampleOutputBS.m"
    #
    # with open(FORM_GENERAL_PATH / "declarations_general.h", "w") as file:
    #     file.write(form_declarations())
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
    # # single_terms_debug = single_terms
    # # del terms
    #
    # if args.tex:
    #     tex_sorted_terms(single_terms)
    #
    # nterms_before = 0
    # for terms in single_terms.values():
    #     for terms_n in terms.values():
    #         nterms_before += len(terms_n)
    #
    # # Remove double terms:
    # single_terms = remove_doubles(single_terms)
    #
    # nterms_after = 0
    # for terms in single_terms.values():
    #     for terms_n in terms.values():
    #         nterms_after += len(terms_n)
    #
    # logger.info(f"#Terms with doubles: {nterms_before:d} <-> #Terms without doubles: {nterms_after:d}")
    #
    # # print terms
    # if args.tex:
    #     tex_sorted_terms_wo_doubles(single_terms)
    #
    # # SUN_Projection of terms without doubles:
    # single_terms = sun_projection(single_terms, max_dim=6)
    # single_terms = replace_sun_tensors_by_projected_ones(single_terms)
    #
    # nterms_after_sun_projection = 0
    # for terms in single_terms.values():
    #     for terms_n in terms.values():
    #         nterms_after_sun_projection += len(terms_n)
    #
    # logger.info(f"Terms after SUN projection: {nterms_after_sun_projection:d}")
    #
    # if args.tex:
    #     tex_terms_sorted_sun_projection(single_terms)
    #
    # pickle.dump(single_terms, open(CONFIG_PATH / "single_terms.p", "wb"))
    single_terms = pickle.load(open(CONFIG_PATH / "single_terms.p", "rb"))

    for type in single_terms.values():
        for term_mass_dim in type.values():
            for summand in term_mass_dim:
                form = declaration_SL2C_sets(summand.possible_indices)
                print(f"{summand:c}")
                print("---------------")

    return single_terms


start_time = timeit.default_timer()

stuff = main()

stop_time = timeit.default_timer()
logger.info(f"Done in {stop_time - start_time:.2f} sec.")

