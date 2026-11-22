"""Writes spacegroup_symbols.json: pymatgen's symbol of each space group, by number."""

from pymatgen.symmetry.groups import SpaceGroup
import json

symbols = [SpaceGroup.from_int_number(i).symbol for i in range(1, 231)]

file = open("assets/spacegroup_symbols.json", "w")
file.write(json.dumps(symbols))
