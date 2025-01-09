from pymatgen.symmetry.groups import SpaceGroup
import json

spacegroup_sym_ops = dict()

symbols = []

for i in range(1, 231):
    spacegroup = SpaceGroup.from_int_number(i)
    symbols.append(spacegroup.symbol)

    sym_ops = []

    for sym_op in spacegroup.symmetry_ops:
        transform = sym_op.affine_matrix.transpose().tolist()
        sym_ops.append(transform)

    spacegroup_sym_ops[spacegroup.int_number] = sym_ops


file = open("assets/spacegroup_symmetry_operations.json", "w")
file.write(json.dumps(spacegroup_sym_ops))

file = open("assets/spacegroup_symbols.json", "w")
file.write(json.dumps(symbols))
