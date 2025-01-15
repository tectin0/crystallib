from pymatgen.symmetry.groups import SpaceGroup, SYMM_DATA
import numpy as np
import json

spacegroup_sym_ops = dict()
spacegroup_generators = dict()

symbols = []


for i in range(1, 231):
    spacegroup = SpaceGroup.from_int_number(i)

    int_symbol = spacegroup.symbol

    # copy & paste from pymatgen
    data = SpaceGroup.sg_encoding[int_symbol]
    # TODO: Support different origin choices.
    enc = list(data["enc"])
    inversion = int(enc.pop(0))
    n_gen = int(enc.pop(0))
    symm_ops = [np.eye(4).transpose().tolist()]
    if inversion:
        symm_ops.append(
            np.array([[-1, 0, 0, 0], [0, -1, 0, 0], [0, 0, -1, 0], [0, 0, 0, 1]])
            .transpose()
            .tolist()
        )
    for _ in range(n_gen):
        matrix = np.eye(4)
        matrix[:3, :3] = SpaceGroup.gen_matrices[enc.pop(0)]
        matrix[0, 3] = SpaceGroup.translations[enc.pop(0)]
        matrix[1, 3] = SpaceGroup.translations[enc.pop(0)]
        matrix[2, 3] = SpaceGroup.translations[enc.pop(0)]
        symm_ops.append(matrix.transpose().tolist())

    spacegroup_generators[spacegroup.int_number] = symm_ops

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

file = open("assets/spacegroup_generators.json", "w")
file.write(json.dumps(spacegroup_generators))
