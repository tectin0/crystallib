"""Writes spacegroup_settings.tsv: spglib's 530 space group settings (one per Hall symbol), in
spglib's order, one per line: number, setting choice, short and full Hermann-Mauguin symbol, Hall
symbol, and the setting's symmetry operations (centring included) as CIFs write them, separated by
`;`. Needs spglib (BSD-3-Clause, see spacegroup_settings-NOTICE.md)."""

from fractions import Fraction

import spglib


def operation_as_xyz(rotation, translation):
    rows = []
    for row, shift in zip(rotation, translation):
        terms = ""
        for coefficient, axis in zip(row, "xyz"):
            if coefficient:
                terms += ("-" if coefficient < 0 else "+") + axis
        shift = Fraction(float(shift)).limit_denominator(12) % 1
        if shift:
            terms += f"+{shift}"
        rows.append(terms.removeprefix("+"))
    return ",".join(rows)


with open("assets/spacegroup_settings.tsv", "w") as file:
    for hall_number in range(1, 531):
        setting = spglib.get_spacegroup_type(hall_number)
        symmetry = spglib.get_symmetry_from_database(hall_number)
        operations = ";".join(
            operation_as_xyz(rotation, translation)
            for rotation, translation in zip(symmetry["rotations"], symmetry["translations"])
        )
        fields = [
            setting.number,
            setting.choice,
            setting.international_short,
            setting.international_full,
            setting.hall_symbol,
            operations,
        ]
        file.write("\t".join(map(str, fields)) + "\n")
