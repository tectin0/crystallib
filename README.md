# crystallib

Crystal structures and the crystallographic data to work with them:

- a phase's unit cell, atoms and measured values with their standard uncertainties
  (`Phase`, `Cell`, `Atom`, `Measured`);
- the symmetry operations and generators of the 230 space groups, symbols and numbers, and the
  positions a site takes in the unit cell (feature `symmetry`);
- element data: van der Waals, covalent, metallic and Shannon's ionic radii, Pauling and Allen
  electronegativities, oxidation states and bond valence parameters (feature `data`).

```rust
use crystallib::{Cell, SpaceGroupSymmetryOperations, IONIC_RADII};

let cell = Cell {
    a: 4.0.into(),
    b: 4.0.into(),
    c: 4.0.into(),
    alpha: 90.0.into(),
    beta: 90.0.into(),
    gamma: 90.0.into(),
    space_group_number: 221,
    ..Default::default()
};
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
| `bvparm2020` | Brown's table of bond valence parameters, for non-profit use only (see below). Implies `data`. |
| `serde` | `Serialize` and `Deserialize` for the structure types. |
| `get-size` | `GetSize` for `Measured` and `AdpType`. |

None is on by default.

## Data and licenses

crystallib's code is licensed under the GNU General Public License v3.0 only (`LICENSE`). The
data in `assets/` comes from the sources each file names in its header:

- `spacegroup_symmetry_operations.json`, `spacegroup_generators.json`, `spacegroup_symbols.json`:
  generated from [pymatgen](https://github.com/materialsproject/pymatgen) by
  `assets/spacegroups.py`.
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
