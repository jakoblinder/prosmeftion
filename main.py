import logging.config
import pickle  # For development purposes only, to save time in debugging.
import timeit
from argparse import ArgumentParser
from datetime import datetime
from pathlib import Path
from yaml import safe_load
import autoeft.io.basis as io_basis

from prosmeftion import CONFIG_PATH, FORM_GENERAL_PATH, PROJECTION_PATH, AUTOEFT_PATH, n_der, get_basis, run_form, op_config
from prosmeftion.general import form_declarations, declaration_SL2C_sets
from prosmeftion.config import sort_model_fields
from prosmeftion.lor_projection import ibp_and_schouten_ids, replace_eoms, rearrange_derivatives
from prosmeftion.yProjection.read_write import remove_header_and_spaces, mathematica_to_form, get_terms, get_type
from prosmeftion.yProjection.utils import equalize_field_indices, remove_doubles
from prosmeftion.converttoSL2C import converttoSL2C
from prosmeftion.tex import tex_unsorted_terms, tex_sorted_terms, tex_sorted_terms_wo_doubles, tex_terms_sorted_sun_projection, tex_sorted_terms_wo_eoms, tex_sorted_terms_before_sun
# FIXME: from prosmeftion.sun_projection import sun_projection, replace_sun_tensors_by_projected_ones
# from prosmeftion.yProjection.tableau import test
from prosmeftion.yProjection.indices import Possible_Indices
from prosmeftion.yProjection.rfr import rfr, rfr_sun

# configure logger
timestamp = datetime.now()
logger_autoeft = logging.getLogger("autoeft")
logger = logger_autoeft.getChild("projection")
# logger = logging.getLogger("autoeft.projection")
with open(CONFIG_PATH / "logger.yml", "r") as file:
    logconfig = safe_load(file)

parser = ArgumentParser(
    description="Project BSUOLEA output onto AUTOEFT-basis.",
)
parser.add_argument("lagrangian", type=Path, help="path to file with matched lagrangian")
parser.add_argument("basis", type=Path, help="path to file a tar containin the corresponding autoeft basis")

parser.add_argument("-s", "--skip", type=int, default=0, help="number of lines to skip in the matched file")
parser.add_argument(
    "-t", "--tex", action="store_true", help="save supplementary tex files"
)
parser.add_argument("-d", "--dimension", type=int, default=6, help="maximum mass dimension for the projection")
parser.add_argument(
    "-f",
    "--format",
    choices=["form", "mathematica"],
    default="form",
    help="Specify in which format the input is given.",
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

input_file = args.lagrangian.resolve()  # "exampleOutputBS.m"
basis_path = args.basis.resolve()
#  AUTOEFT_PATH / Path(f"{model.path}/"

logconfig['handlers']['console']['level'] = logging.getLevelName(args.logLevel)
logging.config.dictConfig(logconfig)

logger.info(timestamp.replace(microsecond=0).isoformat())

def number_terms(single_terms):
    nterms = 0
    for terms in single_terms.values():
        for terms_n in terms.values():
            nterms += len(terms_n)
    return nterms

# basis_path = AUTOEFT_PATH / Path("efts", "ssm-eft")  # "6", "basis"


def main(input_file, basis_path, max_dim = 6, debug=None):
    # Maximum number of derivatives
    nDer = n_der(max_dim)

    # get basis up to max_dim mass dimension
    model, basis = get_basis(basis_path, max_dim)

    fields_sorted_ac, model = sort_model_fields(model)

    inputfilename = input_file.stem

    with open(FORM_GENERAL_PATH / "declarations_general.h", "w") as file:
        file.write(form_declarations(nDer, fields_sorted_ac))

    expression = remove_header_and_spaces(input_file, header=args.skip)

    if args.format == "mathematica":
        expression = mathematica_to_form(expression)

    # Extract coefficient and Operator from the output and pack them into the desired object structure:
    terms = get_terms(expression, model)
    del expression

    # Conversion in y-Basis notation:
    terms = converttoSL2C(terms, nDer, model)

    if args.tex:
        tex_unsorted_terms(inputfilename, terms)

    # Sort terms by type for the projection:
    single_terms = get_type(terms, model)
    del terms

    if args.tex:
        tex_sorted_terms(inputfilename, single_terms)

    nterms_before = number_terms(single_terms)
    # Remove double terms:
    single_terms = remove_doubles(single_terms, model)

    nterms_after = number_terms(single_terms)

    logger.info(f"#Terms with doubles: {nterms_before:d} <-> #Terms without doubles: {nterms_after:d}")

    # print terms
    if args.tex:
        tex_sorted_terms_wo_doubles(inputfilename, single_terms)

    pickle.dump(single_terms, open(CONFIG_PATH / "single_terms_with_eoms.p", "wb"))
    # single_terms = pickle.load(open(CONFIG_PATH / "single_terms_with_eoms.p", "rb"))

    nterms_before = number_terms(single_terms)

    single_terms = replace_eoms(single_terms, nDer, model)
    single_terms = replace_eoms(single_terms, nDer, model)

    nterms_after = number_terms(single_terms)

    logger.info(f"#Terms with EOMs: {nterms_before:d} <-> #Terms with less EOMs: {nterms_after:d}")

    pickle.dump(single_terms, open(CONFIG_PATH / "single_terms_with_less_eoms.p", "wb"))
    # single_terms = pickle.load(open(CONFIG_PATH / "single_terms_with_less_eoms.p", "rb"))

    # rearrange derivative on term of type: nD: 4 & {"H": 1, "H+": 1}
    single_terms = rearrange_derivatives(single_terms, nDer, model)

    single_terms = replace_eoms(single_terms, nDer, model)  # (D^2H) * (D^2H+) -> ... + ~ H+ * H * H * (D^2H+) + ... => need to replace eom again.
    single_terms = replace_eoms(single_terms, nDer, model)

    # pickle.dump(single_terms, open(CONFIG_PATH / "single_terms.p", "wb"))
    # single_terms = pickle.load(open(CONFIG_PATH / "single_terms.p", "rb"))

    nterms_before = number_terms(single_terms)
    # Remove double terms:
    single_terms = remove_doubles(single_terms, model)

    nterms_after = number_terms(single_terms)

    logger.info(f"#Terms with doubles: {nterms_before:d} <-> #Terms without doubles: {nterms_after:d}")

    if args.tex:
        tex_sorted_terms_wo_eoms(inputfilename, single_terms)

    pickle.dump(single_terms, open(CONFIG_PATH / "single_terms.p", "wb"))
    # single_terms = pickle.load(open(CONFIG_PATH / "single_terms.p", "rb"))

    logger.info("Symmetries terms.")
    single_terms = rfr(single_terms, model, fields_sorted_ac)

    nterms_before = number_terms(single_terms)
    single_terms = remove_doubles(single_terms, model)
    nterms_after = number_terms(single_terms)
    logger.info(f"#Terms with doubles: {nterms_before:d} <-> #Terms without doubles: {nterms_after:d}")
    pickle.dump(single_terms, open(CONFIG_PATH / "single_terms_without_rfr.p", "wb"))
    # single_terms = pickle.load(open(CONFIG_PATH / "single_terms_without_rfr.p", "rb"))

    # test()
    # Apply ibp and schouten ids up to the point where all tableaux, are SSYT:
    # Note due to the replacement of contracted derivative on the second and third field,
    # this has to be done iteratively while removing always the EOMs before the next iteration
    if debug:
        single_terms = {debug: single_terms[debug]}

    number_iterations = 0

    while True:
        single_terms, ssyt = ibp_and_schouten_ids(single_terms, max_dim, model)
        number_iterations += 1
        single_terms = replace_eoms(single_terms, nDer, model)
        # pickle.dump(single_terms, open(CONFIG_PATH / "single_terms_test.p", "wb"))
        # single_terms = pickle.load(open(CONFIG_PATH / "single_terms_test.p", "rb"))
        # Remove double terms:
        single_terms = remove_doubles(single_terms, model)
        nterms_after = number_terms(single_terms)
        logger.info(f"#Terms after {number_iterations}. ibp iteration {nterms_after:d}.")
        if ssyt:
            break

    logger.info(f"After {number_iterations:d} iterations of the ibp algorithm all ibp relations and Schouten identities are applied.")

    if args.tex:
        tex_sorted_terms_before_sun(inputfilename, single_terms)

    pickle.dump(single_terms, open(CONFIG_PATH / "single_terms_before_sun.p", "wb"))
    # single_terms = pickle.load(open(CONFIG_PATH / "single_terms_before_sun.p", "rb"))

    single_terms = rfr_sun(single_terms, model)
    # single_terms = remove_doubles(single_terms, model)

    # # SUN_Projection of terms without doubles:
    # single_terms = sun_projection(single_terms, basis, max_dim, model)
    #
    # single_terms = replace_sun_tensors_by_projected_ones(single_terms, model)

    nterms_after_sun_projection = number_terms(single_terms)

    logger.info(f"Terms after SUN projection: {nterms_after_sun_projection:d}")

    if args.tex:
        tex_terms_sorted_sun_projection(inputfilename, single_terms, max_dim)

    # pickle.dump(single_terms, open(CONFIG_PATH / "single_terms.p", "wb"))
    # single_terms = pickle.load(open(CONFIG_PATH / "single_terms.p", "rb"))


    return single_terms


start_time = timeit.default_timer()

terms = main(input_file, basis_path, args.dimension)
# fieldstructure = tuple(map(int,list("0100000110000000")))
# terms = main(args.dimension, debug=fieldstructure)

output = ""
for term_type in terms.values():
    for term_mass_dim in term_type.values():
        for term in term_mass_dim:
            output += f"+{term:c}\n"
n_terms = number_terms(terms)

with open(PROJECTION_PATH / f"{input_file.stem:s}_yBasis_{n_terms}.h", "w") as file:
    file.write(output)

stop_time = timeit.default_timer()
logger.info(f"Done in {stop_time - start_time:.2f} sec.")

