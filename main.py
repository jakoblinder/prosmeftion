import logging
import logging.config
import timeit

from argparse import ArgumentParser
from pathlib import Path
from yaml import safe_load
from datetime import datetime

import tioc.get_FORM.class_index

from tioc.projection import CONFIG_PATH, converttoSL2C, write_texfile, get_type  # , main
from tioc.get_BS.class_term import Term

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

    terms = converttoSL2C(bs_file, header=args.skip, pprint=False)

    if args.tex:
        write_texfile(terms)

    single_terms = get_type(terms)
    print(len(list(single_terms.keys())[0]))
    return single_terms

# for i in Term.extractOrder():
#     print(f"Name: {i.name:2s} - Form: {i.form_name:6s} - Helicity: {i.helicity:.1f} - Fermion: {i.fermion}")
# print(Term.form_sortfields(Term.extractOrder()))

start_time = timeit.default_timer()

main()

stop_time = timeit.default_timer()
logger.info(f"Done in {stop_time - start_time:.2f} sec.")

