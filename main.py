import logging.config
import logging.config
import timeit
from argparse import ArgumentParser
from datetime import datetime
from pathlib import Path

from yaml import safe_load

from tioc.projection import CONFIG_PATH, converttoSL2C, tex_unsorted_terms, tex_sorted_terms  # , main
from tioc.sun_projection import get_type, sun_projection
from tioc import form_declarations, FORM_GENERAL_PATH, op_config

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
    op = op_config
    bs_file = args.matched.resolve()  # "exampleOutputBS.m"

    with open(FORM_GENERAL_PATH / "declarations_general.h", "w") as file:
        file.write(form_declarations())

    terms = converttoSL2C(bs_file, header=args.skip, pprint=False)

    if args.tex:
        tex_unsorted_terms(terms)
    # Sort terms by type and do projection again.
    # for term in terms:
    #     for summand in term.terms:
    #         print(f"Fieldcounter: {summand.fieldcounter}")
    #         print(f"Summand: {summand}")

    single_terms = get_type(terms)
    del terms
    if args.tex:
        tex_sorted_terms(single_terms)

    # SUN_Projection:
    single_terms = sun_projection(single_terms, max_dim=6)

    return single_terms


start_time = timeit.default_timer()

stuff = main()

stop_time = timeit.default_timer()
logger.info(f"Done in {stop_time - start_time:.2f} sec.")

