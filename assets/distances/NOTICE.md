# Distance distributions: source and licence

`distances.bin.zlib` holds distributions of bond lengths and interatomic distances derived from
the Materials Project's database release 2025-09-25, used under the Creative Commons Attribution
4.0 International licence (CC BY 4.0, https://creativecommons.org/licenses/by/4.0/).

Only materials whose `builder_meta.license` is `BY-C` (CC BY 4.0), that are matched to an
experimental structure (`theoretical: false`) and not deprecated were used: 49,283 structures.
They were changed: each structure's charges were guessed and its bonds found by structal, and its
distances were binned per kind of atom pair and summed over all structures (see
structal's `structal-lib/src/harnesses.rs`, the ignored test `generate`, which writes this
file).

Materials Project Data was accessed on 2026-11-19 from
https://registry.opendata.aws/materials-project.

A. Jain, S. P. Ong, G. Hautier, W. Chen, W. D. Richards, S. Dacek, S. Cholia, D. Gunter,
D. Skinner, G. Ceder and K. A. Persson, "Commentary: The Materials Project: A materials genome
approach to accelerating materials innovation", APL Materials 1, 011002 (2013),
https://doi.org/10.1063/1.4812323.
