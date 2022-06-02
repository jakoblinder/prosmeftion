import logging.config
import logging.config
import timeit
from argparse import ArgumentParser
from datetime import datetime
from pathlib import Path

from yaml import safe_load

from tioc.projection import CONFIG_PATH, converttoSL2C, write_texfile  # , main
from tioc.sun_projection import get_type, sun_projection
from tioc import form_declarations, FORM_GENERAL_PATH

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
    bs_file = args.matched.resolve()  # "exampleOutputBS.m"

    with open(FORM_GENERAL_PATH / "declarations_general.h", "w") as file:
        file.write(form_declarations())

    terms = converttoSL2C(bs_file, header=args.skip, pprint=False)

    # TODO: Write tex:
    if args.tex:
        write_texfile(terms)

    single_terms = get_type(terms)
    del terms
    # SUN_Projection:
    single_terms = sun_projection(single_terms, max_dim = 6)

    # print(type(single_terms))

    return single_terms

# for i in Term.extractOrder():
#     print(f"Name: {i.name:2s} - Form: {i.form_name:6s} - Helicity: {i.helicity:.1f} - Fermion: {i.fermion}")
# print(Term.form_sortfields(Term.extractOrder()))

start_time = timeit.default_timer()

main()

stop_time = timeit.default_timer()
logger.info(f"Done in {stop_time - start_time:.2f} sec.")

