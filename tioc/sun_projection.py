from fractions import Fraction
from pathlib import Path

from . import AUTOEFT_PATH, model
from autoeft.invariants import SUNTableau, field_projection_operator, symmetrize_tensors
from autoeft.io import load_basis
from autoeft.sun_projection import tensor_projection

def sun_internal_projector(op_type, group):
    """Construct a projector to impose the 'internal symmetries'."""
    field_id = 1
    projection_operator = []
    for field, multiplicity in op_type.fields:
        for i in range(multiplicity):
            field_id += i
            field_tableau = SUNTableau.fill_tableau(
                field_id, field.sun_reprs[group].partition
            )
            projection_operator.append(field_projection_operator(field_tableau))
    return projection_operator


def form_basis(op_type, tensors, group):
    """Return monomial and symmetrized tensors in FORM-readable output."""
    internal_projector = sun_internal_projector(op_type, group)

    monomial_basis = []
    tensor_basis = []
    for tensor in tensors:
        monomial_basis.append(tensor.to_form())
        sym_tensors = symmetrize_tensors([(tensor, Fraction(1))], internal_projector)
        tensor_basis.append(SUNTableau.tensors_to_form(sym_tensors))
    return monomial_basis, tensor_basis

# load basis
eft_path = AUTOEFT_PATH / Path("eft/ssmeft/")  # path to eft output

def sun_projection(max_dim: int):
    basis = {}  # dictionary with basis for each mass dimension from 4 to 6.
    for dim in range(4, max_dim + 1):
        # load_basis also returns some counters and the Hilbert series, which we don't need here...
        basis[dim], _, _ = load_basis(eft_path, dim)
    assert basis[max_dim].model.name == model.name, f"The model used for construction of the basis is {basis[max_dim].model.name}, but it should be {model.name}."

    return basis


# example for L L H H, SU(2) tensor
su2_example = [
    "e_(idxF1I1,idxF2I1)*e_(idxF3I1,idxF4I1)/2 + e_(idxF1I1,idxF3I1)*e_(idxF2I1,idxF4I1)/2"
]

su2 = model.sun_groups["SU2_W"]  # group name must match definition in model file
field_content = {"L": 2, "H": 2}  # field names must match definition in model file
derivatives = 0

# get all the invariants associated with the operator
op_class, op_subclass, op_type, operator = basis.get_operator(
    field_content, derivatives
)

# get the SU(2) tensors
su2_basis = operator.sun_tensors[su2]
# convert the tensors to FORM-readable output
su2_monom, su2_tensor = form_basis(op_type, su2_basis, su2)
# actual projection
P, G = tensor_projection(su2.N, su2_monom, su2_tensor, su2_example)
# P*G^(-1) gives the projection matrix
print(P)
print(G)


# TODO: Load Basis from autoeft

# a = eft_basis.load_basis(AUTOEFT_PATH / "eft/smeft", 5)  # 5/operators/4/psi(2) phi(2) + h.c.
# print(AUTOEFT_PATH / "eft/smeft/5/operators/4/psi(2) phi(2) + h.c./L L H H.yml")
# with open(AUTOEFT_PATH / "eft/smeft/5/operators/4/psi(2) phi(2) + h.c./L L H H.yml", "r") as file:
#     ops_dim5 = safe_load(file)
# b = eft_basis.load_operator(a[0], ops_dim5)

