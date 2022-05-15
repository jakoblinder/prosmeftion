import logging
import logging.config
import timeit

import sage.all
import sage.matrix as mx
from sage.rings.rational_field import QQ
from autoeft.sun_projection import tensor_projection

from argparse import ArgumentParser
from pathlib import Path
from yaml import safe_load
from datetime import datetime

import tioc.get_FORM.class_index

from tioc.projection import CONFIG_PATH, converttoSL2C, write_texfile, model  # , main
from tioc.sun_projection import get_type, get_basis, get_basis_tensors
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
    basis = get_basis(6)
    for type in single_terms.values():
        for term_mass_dim in type.values():
            field_content = term_mass_dim.field_content
            derivatives = term_mass_dim.n_D
            mass_dim = term_mass_dim.mass_dim
            try:
                term_mass_dim.sun_projection_tensors = get_basis_tensors(basis, field_content, derivatives, mass_dim)
                # There exists a term in the basis which matches the type and the projection can be done:
                # sun_tensors = list of all su2 tensors contained in all terms of the same type.
                for sun_group, tensors in term_mass_dim.sun_projection_tensors.items():
                    if tensors:
                        N = model.sun_groups[sun_group].N
                        sun_monom  = tensors["sun_monom"]
                        sun_tensor = tensors["sun_tensor"]
                    else:
                        # Since the basis doesn't contain basis tensors for this group, also the tensor itself shouldn't
                        # and thus the next possible SUN group is tried.
                        continue
                    # collect all sun tensors of the group and the tensor in one ordered list for the projection.
                    sun_field_tensors = []
                    for term in term_mass_dim.terms:
                        sun_field_tensors.append(term.gaugeTensorsSUN[sun_group])

                    # project all tensors simultaneously
                    P_map, G_map = tensor_projection(N, sun_monom, sun_tensor, sun_field_tensors)
                    dim = len(sun_monom)
                    n_projection_op = len(sun_field_tensors)
                    G = mx.constructor.matrix(QQ, dim, dim, G_map)
                    P = mx.constructor.matrix(QQ, n_projection_op, dim, P_map)

                    print(f"P: (\n{P})")
                    print("----")
                    print(f"G: (\n{G})")
                    print("----")
                    # P*G^(-1) gives the projection matrix
                    print(f"P*G^-1: (\n{P * G.inverse()})")
                    # print(P * G.inverse())
            except KeyError:
                # Term contains still redundancies and therefore their exist no basis tensor.
                term_mass_dim.sun_projection_tensors = False
            except ValueError:
                pass

    return single_terms

# for i in Term.extractOrder():
#     print(f"Name: {i.name:2s} - Form: {i.form_name:6s} - Helicity: {i.helicity:.1f} - Fermion: {i.fermion}")
# print(Term.form_sortfields(Term.extractOrder()))

start_time = timeit.default_timer()

main()

stop_time = timeit.default_timer()
logger.info(f"Done in {stop_time - start_time:.2f} sec.")

