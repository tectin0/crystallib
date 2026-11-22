# crystallib

Crystal structures and the crystallographic data to work with them:

- a phase's unit cell, atoms and measured values with their standard uncertainties
  (`Phase`, `Cell`, `Atom`, `Measured`), and the oxidation state a type symbol names (`Fe3+`,
  `O2-`; `type_symbol`); atoms near one another, also across the cell's edges (`geometry`);
  the shape the atoms around a centre form, among the coordination polyhedra (`shapes`);
- the symmetry operations and generators of the 230 space groups, symbols and numbers, their 530
  settings and the setting a Hermann–Mauguin symbol names (`SpaceGroupSetting`), and the
  positions a site takes in the unit cell, and the atoms of a phase's cell with how often their
  places are occupied (`cell_atoms`) (feature `symmetry`);
- element data: van der Waals, covalent, metallic and Shannon's ionic radii, Pauling and Allen
  electronegativities, oxidation states and bond valence parameters, and which partly occupied
  atoms can be present together (`occupancy`) (feature `data`);
- structure analysis: whether two atoms are bonded and why (`bonding`, per element pair
  `bond_pairs`), oxidation states guessed from bond valence sums (`oxidation`), the kind of bond
  two elements make (`bond_character`), which atoms of a disordered structure are present
  together (`disorder`), and hydrogen bonds and short contacts (`contacts`) (features `data` and
  `symmetry`).

```rust
use crystallib::{Cell, SpaceGroupSymmetryOperations, IONIC_RADII};

let mut cell = Cell::default();
(cell.a, cell.b, cell.c) = (4.0.into(), 4.0.into(), 4.0.into());
(cell.alpha, cell.beta, cell.gamma) = (90.0.into(), 90.0.into(), 90.0.into());
cell.space_group_number = 221;
let positions = SpaceGroupSymmetryOperations::site_positions_in_cell(
    &cell,
    cgmath::Point3::new(0.5, 0.0, 0.0),
    None,
);
assert_eq!(positions.map(|positions| positions.len()), Some(3));

assert_eq!(IONIC_RADII.get("Na", 1, 6), Some(1.02));
```

## Features

| Feature | What it adds |
|---------|--------------|
| `symmetry` | Space group symmetry operations, generators, symbols and numbers. |
| `data` | Element data: radii, electronegativities, oxidation states, bond valence parameters. |
| `cif` | Reading a phase out of a CIF (`cif_phase`), with the `cif` crate. |
| `distance-distributions` | How long bonds and distances typically are (`distances`), from the Materials Project. |
| `bvparm2020` | Brown's table of bond valence parameters, for non-profit use only (see below). Implies `data`. |
| `serde` | `Serialize` and `Deserialize` for the structure types. |
| `get-size` | `GetSize` for `Measured` and `AdpType`. |

None is on by default.

## Data and licenses

crystallib's code is licensed under the GNU General Public License v3.0 only (`LICENSE`), with
an additional permission to combine it with code under the GSAS-II Open Source License or terms
derived from it (`LICENSE-ADDITIONAL-PERMISSION`). The data in `assets/` comes from the sources
each file names in its header:

- `spacegroup_generators.json`, `spacegroup_symbols.json`: generated from
  [pymatgen](https://github.com/materialsproject/pymatgen) by `assets/spacegroups.py`.
- `spacegroup_settings.tsv`: the 530 space group settings and their symmetry operations, from
  [spglib](https://github.com/spglib/spglib) (BSD-3-Clause), written by
  `assets/spacegroup_settings.py`; see `spacegroup_settings-NOTICE.md`.
- `data/shannon-radii.csv`: Shannon (1976), from pymatgen; see `data/shannon-radii-NOTICE.md`.
- `data/electronegativity.toml`, `data/metallic-radii.toml`, `data/oxidation-states.toml`,
  `data/bond-valence-1991.toml`: the cited sources, as collected in pymatgen.
- `data/covalent-radii.toml`, `data/van-der-Waals-radii.toml`: the cited publications.
- `data/electronegativity-allen.toml`: the cited publications (Allen; Mann, Meek, Allen et al.), as
  tabulated on Wikipedia.

pymatgen is licensed under the MIT license; its text is in `data/shannon-radii-NOTICE.md`.

`data/bvparm2020.cif` is I. David Brown's table of bond valence parameters, distributed by the
IUCr. It is not covered by crystallib's license: it may be used and distributed without fee for
non-profit purposes only, see `data/bvparm2020-NOTICE.md`. crystallib compiles it in only with the
`bvparm2020` feature.

`distances/distances.bin.zlib` holds distributions of bond lengths and distances derived from the
Materials Project's database, under CC BY 4.0; see `distances/NOTICE.md`. crystallib compiles it
in only with the `distance-distributions` feature.
