//! How long bonds and distances typically are, from the experimentally matched structures of the
//! Materials Project (release 2025-09-25, CC BY 4.0, see `assets/distances/NOTICE.md`), to
//! compare a structure's bonds and distances with.
//!
//! A bond's distribution is that of its centre's species (element and charge, `Fe3+`), the
//! centre's coordination number and its partner's species, with charges from
//! [`crate::oxidation::with_guessed_oxidation_states`], bonds by the default
//! [`BondRule`](crate::bonding::BondRule). Where too few structures have one, [`bond`] falls
//! back to any coordination number, then to the elements alone. A distance between any two
//! atoms has the distribution of all distances between their species up to [`CUTOFF`] (see
//! [`pair`]).
//!
//! The structures are DFT relaxations, about half of them with GGA (or GGA+U), half with
//! r2SCAN. Against the same table built from the COD's measured structures (the ignored test
//! `compare`), a key's median runs 0.015 Å longer (the median over 6,356 keys both have);
//! cation–oxygen bonds 0.00–0.03 Å, K–O 0.05 Å.

use std::sync::LazyLock;

use std::collections::HashMap;

use crate::elements::element_of;

/// The width of a bond distribution's bins (Å).
pub const BOND_BIN: f32 = 0.01;

/// The width of the bins of a distribution of all distances (see [`pair`]) (Å).
pub const PAIR_BIN: f32 = 0.02;

/// The longest distance a distribution covers (Å).
pub const CUTOFF: f32 = 4.0;

/// The fewest structures a distribution comes from; a key with fewer has none.
pub const MIN_STRUCTURES: u32 = 10;

/// Where the distributions come from, as a host tells its users.
pub const SOURCE: &str = "From the Materials Project (materialsproject.org), release 2025-09-25, \
     under CC BY 4.0: its DFT relaxations (GGA or r2SCAN) of the 49,283 structures it matches to \
     experimental ones, each structure weighing the same. Their bonds run about 0.01–0.03 Å \
     longer than measured ones. A. Jain et al., APL Materials 1, 011002 (2013).";

/// The coordination number of a key that stands for any.
const ANY_COORDINATION: u8 = 0;

/// The table, zlib-compressed, as written by this module's generator (the ignored test
/// `generate`).
const DATA: &[u8] = include_bytes!("../assets/distances/distances.bin.zlib");

/// The table, read on first use.
static TABLE: LazyLock<Table> = LazyLock::new(|| {
    let table = miniz_oxide::inflate::decompress_to_vec_zlib(DATA)
        .map_err(|error| format!("it doesn't decompress: {error}"))
        .and_then(|bytes| Table::parse(&bytes));
    table.unwrap_or_else(|error| {
        log::error!("The bundled distance distributions can't be read: {error}");
        Table::default()
    })
});

/// How distances between two kinds of atoms are distributed.
#[derive(Debug, Clone, PartialEq)]
pub struct Distribution {
    /// How many structures it comes from.
    pub structures: u32,
    /// The distance (Å) at the start of the first bin.
    pub start: f32,
    /// The width of each bin (Å).
    pub bin: f32,
    /// Each bin's share, from `start` on in steps of `bin`; they sum to 1. Each structure
    /// weighs the same, however many distances it has.
    pub weights: Vec<f32>,
}

impl Distribution {
    /// The share of distances shorter than `distance`.
    pub fn fraction_below(&self, distance: f32) -> f32 {
        let position = (distance - self.start) / self.bin;
        let whole = position.floor().max(0.0) as usize;
        let below: f32 = self.weights.iter().take(whole).sum();
        let partial = self.weights.get(whole).map_or(0.0, |weight| {
            weight * (position - whole as f32).clamp(0.0, 1.0)
        });
        (below + partial).clamp(0.0, 1.0)
    }

    /// The distance (Å) below which the share `fraction` of the distances lie.
    pub fn quantile(&self, fraction: f32) -> f32 {
        let mut sum = 0.0;
        for (index, weight) in self.weights.iter().enumerate() {
            if *weight > 0.0 && sum + weight >= fraction {
                let within = (fraction - sum) / weight;
                return self.start + (index as f32 + within) * self.bin;
            }
            sum += weight;
        }
        self.start + self.weights.len() as f32 * self.bin
    }

    /// The distance (Å) at the middle of the bin with the largest share.
    pub fn mode(&self) -> f32 {
        let index = self
            .weights
            .iter()
            .enumerate()
            .max_by(|a, b| a.1.total_cmp(b.1))
            .map_or(0, |(index, _)| index);
        self.start + (index as f32 + 0.5) * self.bin
    }
}

/// The distribution of a bond's length, and the key it was found under.
#[derive(Debug, Clone, PartialEq)]
pub struct BondDistribution {
    /// The centre's species (`Fe3+`), or its element where the species had too few structures.
    pub centre: String,
    /// The centre's coordination number, or `None` where it stands for any.
    pub coordination: Option<usize>,
    /// The partner's species or element, like `centre`.
    pub partner: String,
    /// How its lengths are distributed.
    pub distribution: Distribution,
}

/// The distribution of the lengths of bonds from an atom of type `centre` (a type symbol,
/// `Fe3+`) with `coordination` bonds to atoms of type `partner`. Falls back to any coordination
/// number, then to the elements with and without it, where a key has fewer than
/// [`MIN_STRUCTURES`]; `None` where all of them do.
pub fn bond(centre: &str, coordination: usize, partner: &str) -> Option<BondDistribution> {
    let (centre_species, partner_species) = (species(centre), species(partner));
    let (centre_element, partner_element) = (element_of(centre), element_of(partner));
    let coordination = u8::try_from(coordination).unwrap_or(u8::MAX);

    [
        (&centre_species, coordination, &partner_species),
        (&centre_species, ANY_COORDINATION, &partner_species),
        (&centre_element, coordination, &partner_element),
        (&centre_element, ANY_COORDINATION, &partner_element),
    ]
    .into_iter()
    .find_map(|(centre, coordination, partner)| {
        let distribution = TABLE.get(Kind::Bond, centre, coordination, partner)?;
        Some(BondDistribution {
            centre: centre.clone(),
            coordination: (coordination != ANY_COORDINATION).then_some(usize::from(coordination)),
            partner: partner.clone(),
            distribution,
        })
    })
}

/// The distribution of the lengths of bonds between atoms of the elements `a` and `b`, whatever
/// their charges and coordination numbers; `None` where fewer than [`MIN_STRUCTURES`] have one.
pub fn element_bond(a: &str, b: &str) -> Option<BondDistribution> {
    let (a, b) = (element_of(a), element_of(b));
    [(&a, &b), (&b, &a)]
        .into_iter()
        .find_map(|(centre, partner)| {
            let distribution = TABLE.get(Kind::Bond, centre, ANY_COORDINATION, partner)?;
            Some(BondDistribution {
                centre: centre.clone(),
                coordination: None,
                partner: partner.clone(),
                distribution,
            })
        })
}

/// The distribution of all distances up to [`CUTOFF`] between atoms of the types `a` and `b`
/// (type symbols), by their species, else by their elements; `None` where both have fewer than
/// [`MIN_STRUCTURES`]. Returns the two keys it was found under with it.
pub fn pair(a: &str, b: &str) -> Option<(String, String, Distribution)> {
    [(species(a), species(b)), (element_of(a), element_of(b))]
        .into_iter()
        .find_map(|(a, b)| {
            let (first, second) = if a <= b { (&a, &b) } else { (&b, &a) };
            let distribution = TABLE.get(Kind::Pair, first, ANY_COORDINATION, second)?;
            Some((a.clone(), b.clone(), distribution))
        })
}

/// The species a type symbol names: its element with the charge it names, written the same
/// way whatever way the type symbol writes it (`Fe+3` and `Fe3+` are `Fe3+`, `Na+` is `Na1+`,
/// `Cu0+` is `Cu0`), or the element alone where it names none.
pub fn species(type_symbol: &str) -> String {
    let element = element_of(type_symbol);
    match crate::type_symbol::state_of(type_symbol) {
        None => element,
        Some(0.0) => format!("{element}0"),
        Some(state) => {
            let sign = if state > 0.0 { '+' } else { '-' };
            format!("{element}{}{sign}", state.abs())
        }
    }
}

/// What a record of the table describes.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
enum Kind {
    /// The lengths of bonds from a centre to a partner.
    Bond = 0,
    /// All distances between two species or elements.
    Pair = 1,
}

impl Kind {
    /// The width of the bins of a record of this kind (Å).
    fn bin(self) -> f32 {
        match self {
            Kind::Bond => BOND_BIN,
            Kind::Pair => PAIR_BIN,
        }
    }
}

/// The parsed table: each record's distribution by its key.
#[derive(Debug, Default)]
struct Table {
    /// Each key's index into `names`, by the name.
    indices: HashMap<String, u16>,
    /// Each record by its kind, its first key, its coordination number and its second key.
    records: HashMap<(Kind, u16, u8, u16), Distribution>,
}

impl Table {
    /// The table's magic number and format version.
    const MAGIC: &[u8; 4] = b"SDD2";

    /// The distribution of the record of `kind` with the keys `a`, `coordination` and `b`, where
    /// it has at least [`MIN_STRUCTURES`].
    fn get(&self, kind: Kind, a: &str, coordination: u8, b: &str) -> Option<Distribution> {
        let key = (
            kind,
            *self.indices.get(a)?,
            coordination,
            *self.indices.get(b)?,
        );
        self.records
            .get(&key)
            .filter(|distribution| distribution.structures >= MIN_STRUCTURES)
            .cloned()
    }

    /// Reads the table `bytes` hold (see [`write`] for its layout).
    fn parse(bytes: &[u8]) -> Result<Self, String> {
        let mut reader = Reader(bytes);
        if reader.take(4)? != Self::MAGIC {
            return Err("not a distance table of this version".to_string());
        }

        let names = (0..reader.u32()?)
            .map(|_| {
                let length = usize::from(reader.u8()?);
                String::from_utf8(reader.take(length)?.to_vec()).map_err(|error| error.to_string())
            })
            .collect::<Result<Vec<_>, _>>()?;
        let indices = names
            .iter()
            .enumerate()
            .map(|(index, name)| {
                Ok((
                    name.clone(),
                    u16::try_from(index).map_err(|e| e.to_string())?,
                ))
            })
            .collect::<Result<_, String>>()?;

        let mut records = HashMap::default();
        for _ in 0..reader.u32()? {
            let kind = match reader.u8()? {
                0 => Kind::Bond,
                1 => Kind::Pair,
                other => return Err(format!("unknown record kind {other}")),
            };
            let (a, coordination, b) = (reader.u16()?, reader.u8()?, reader.u16()?);
            let structures = reader.u32()?;
            let bin = kind.bin();
            let first_bin = reader.u16()?;
            let weights: Vec<f32> = (0..reader.u16()?)
                .map(|_| Ok(dequantize(reader.u8()?)))
                .collect::<Result<_, String>>()?;

            let sum: f32 = weights.iter().sum();
            let distribution = Distribution {
                structures,
                start: f32::from(first_bin) * bin,
                bin,
                weights: weights
                    .iter()
                    .map(|weight| weight / sum.max(f32::EPSILON))
                    .collect(),
            };
            records.insert((kind, a, coordination, b), distribution);
        }

        Ok(Self { indices, records })
    }
}

/// A bin's weight relative to the record's largest, from its byte in the table: the square root
/// of that share, times 255, so a share as small as 1/65025 still shows.
fn dequantize(byte: u8) -> f32 {
    (f32::from(byte) / 255.0).powi(2)
}

/// The byte [`dequantize`] reads back as about `share` (0 to 1).
#[cfg(test)]
fn quantize(share: f64) -> u8 {
    (share.sqrt() * 255.0).round() as u8
}

/// Reads little-endian numbers from the front of a byte slice.
struct Reader<'a>(&'a [u8]);

impl Reader<'_> {
    fn take(&mut self, count: usize) -> Result<&[u8], String> {
        if self.0.len() < count {
            return Err("the table ends early".to_string());
        }
        let (taken, rest) = self.0.split_at(count);
        self.0 = rest;
        Ok(taken)
    }

    fn u8(&mut self) -> Result<u8, String> {
        Ok(self.take(1)?[0])
    }

    fn u16(&mut self) -> Result<u16, String> {
        Ok(u16::from_le_bytes(
            self.take(2)?.try_into().unwrap_or_default(),
        ))
    }

    fn u32(&mut self) -> Result<u32, String> {
        Ok(u32::from_le_bytes(
            self.take(4)?.try_into().unwrap_or_default(),
        ))
    }
}

/// One record to write: its key, its structures and its weights per bin of [`BOND_BIN`] from 0
/// on (each structure's distances of the key summing to 1).
#[cfg(test)]
struct Record {
    kind: Kind,
    a: String,
    coordination: u8,
    b: String,
    structures: u32,
    weights: Vec<f64>,
}

/// The table of `records` as bytes: the magic number; the count of key names, then each as its
/// length (u8) and UTF-8 bytes; the count of records, then each as its kind (u8), first key
/// (u16, an index into the names), coordination number (u8, 0 for any), second key (u16), its
/// structures (u32), its first bin (u16, in bins of its kind's width, [`Kind::bin`]), its count
/// of bins (u16) and each bin's weight (u8, see [`dequantize`]), from the first bin with any
/// weight to the last. All little-endian.
#[cfg(test)]
fn write(records: &[Record]) -> Vec<u8> {
    let mut names: Vec<&str> = records
        .iter()
        .flat_map(|record| [record.a.as_str(), record.b.as_str()])
        .collect();
    names.sort_unstable();
    names.dedup();
    let index = |name: &str| names.binary_search(&name).unwrap() as u16;

    let mut bytes = Table::MAGIC.to_vec();
    bytes.extend((names.len() as u32).to_le_bytes());
    for name in &names {
        bytes.push(name.len() as u8);
        bytes.extend(name.as_bytes());
    }

    bytes.extend((records.len() as u32).to_le_bytes());
    for record in records {
        let merged = (record.kind.bin() / BOND_BIN).round() as usize;
        let weights: Vec<f64> = record
            .weights
            .chunks(merged)
            .map(|bins| bins.iter().sum())
            .collect();
        let largest = weights.iter().copied().fold(0.0, f64::max);
        let quantized: Vec<u8> = weights
            .iter()
            .map(|weight| quantize(weight / largest))
            .collect();
        let first = quantized.iter().position(|&weight| weight > 0).unwrap_or(0);
        let last = quantized
            .iter()
            .rposition(|&weight| weight > 0)
            .unwrap_or(0);

        bytes.push(record.kind as u8);
        bytes.extend(index(&record.a).to_le_bytes());
        bytes.push(record.coordination);
        bytes.extend(index(&record.b).to_le_bytes());
        bytes.extend(record.structures.to_le_bytes());
        bytes.extend((first as u16).to_le_bytes());
        bytes.extend(((last + 1 - first) as u16).to_le_bytes());
        bytes.extend(&quantized[first..=last]);
    }
    bytes
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_species_is_written_one_way() {
        assert_eq!(species("Fe3+"), "Fe3+");
        assert_eq!(species("Fe+3"), "Fe3+");
        assert_eq!(species("Na+"), "Na1+");
        assert_eq!(species("O2-"), "O2-");
        assert_eq!(species("Fe2.5+"), "Fe2.5+");
        assert_eq!(species("Cu0+"), "Cu0");
        assert_eq!(species("Fe"), "Fe");
        assert_eq!(species("Fe1"), "Fe");
    }

    #[test]
    fn a_written_table_reads_back() {
        let records = [
            Record {
                kind: Kind::Bond,
                a: "Si4+".to_string(),
                coordination: 4,
                b: "O2-".to_string(),
                structures: 12,
                weights: [vec![0.0; 160], vec![1.0, 3.0, 1.0], vec![0.0; 10]].concat(),
            },
            Record {
                kind: Kind::Bond,
                a: "Al3+".to_string(),
                coordination: 4,
                b: "O2-".to_string(),
                structures: 3,
                weights: [vec![0.0; 175], vec![1.0]].concat(),
            },
            Record {
                kind: Kind::Pair,
                a: "O2-".to_string(),
                coordination: ANY_COORDINATION,
                b: "Si4+".to_string(),
                structures: 10,
                weights: [vec![0.0; 161], vec![1.0, 1.0]].concat(),
            },
        ];
        let table = Table::parse(&write(&records)).unwrap();

        let bond = table.get(Kind::Bond, "Si4+", 4, "O2-").unwrap();
        assert_eq!(bond.structures, 12);
        assert!((bond.start - 1.60).abs() < 1e-6);
        assert_eq!(bond.weights.len(), 3);
        assert!((bond.weights[1] - 0.6).abs() < 1e-3, "{bond:?}");
        assert!((bond.mode() - 1.615).abs() < 1e-5);
        assert!((bond.fraction_below(1.61) - 0.2).abs() < 1e-3);
        assert!((bond.quantile(0.5) - 1.615).abs() < 1e-3);
        assert_eq!(table.get(Kind::Bond, "Si4+", 6, "O2-"), None);
        assert_eq!(
            table.get(Kind::Bond, "Al3+", 4, "O2-"),
            None,
            "too few structures"
        );

        let pair = table
            .get(Kind::Pair, "O2-", ANY_COORDINATION, "Si4+")
            .unwrap();
        assert!((pair.bin - PAIR_BIN).abs() < 1e-6);
        assert!((pair.start - 1.60).abs() < 1e-6);
        assert_eq!(
            pair.weights.len(),
            2,
            "two bins of 0.01 Å each in one of 0.02 Å"
        );
    }

    #[test]
    fn a_tail_bin_far_below_the_largest_still_shows() {
        assert!(dequantize(quantize(1e-4)) > 0.0);
        assert!((dequantize(quantize(0.25)) - 0.25).abs() < 0.005);
        assert_eq!(dequantize(quantize(1.0)), 1.0);
    }

    #[test]
    fn the_bundled_table_knows_silicate_and_rock_salt_bonds() {
        let silicon = bond("Si4+", 4, "O2-").unwrap();
        assert_eq!(
            (silicon.centre.as_str(), silicon.coordination),
            ("Si4+", Some(4))
        );
        assert!(
            (1.58..1.68).contains(&silicon.distribution.mode()),
            "{silicon:?}"
        );

        let sodium = bond("Na+", 6, "Cl-").unwrap();
        assert!(
            (2.75..2.90).contains(&sodium.distribution.mode()),
            "{sodium:?}"
        );

        let (_, _, distances) = pair("O2-", "Si4+").unwrap();
        assert!(distances.fraction_below(1.5) < 0.01);

        for (a, b) in [("O", "Ba"), ("Ba", "O")] {
            let barium = element_bond(a, b).unwrap();
            assert_eq!(barium.coordination, None);
            assert!(
                (2.70..3.00).contains(&barium.distribution.mode()),
                "{barium:?}"
            );
        }
    }

    /// Writes the distance table to `$DISTANCES_OUT` from either the Materials Project's summary
    /// collection in `$MP_SUMMARY` (its CC BY materials matched to an experimental structure,
    /// not deprecated) or the CIFs whose paths `$CIF_LIST` lists (fully occupied, at most 1000
    /// atoms in the cell). For each structure, charges are guessed and bonds found by the default
    /// [`BondRule`](crate::bonding::BondRule); each key's distances in a structure weigh 1 together. A bond is recorded under its
    /// centre's species, coordination number and partner's species, and again under any
    /// coordination number and under the elements; every distance up to [`CUTOFF`] under its
    /// pair of species and of elements.
    #[test]
    #[ignore = "needs MP_SUMMARY or CIF_LIST, and DISTANCES_OUT"]
    fn generate() {
        use std::sync::atomic::{AtomicUsize, Ordering};

        let out = std::env::var("DISTANCES_OUT").expect("DISTANCES_OUT");
        let sources = generator::sources();
        let progress = AtomicUsize::new(0);
        let counted = AtomicUsize::new(0);

        // A failing structure panics deep in the pipeline now and then; it is skipped quietly.
        std::panic::set_hook(Box::new(|_| {}));
        let threads = std::thread::available_parallelism().map_or(4, |n| n.get());
        let totals = std::thread::scope(|scope| {
            let workers: Vec<_> = (0..threads)
                .map(|thread| {
                    let (sources, progress, counted) = (&sources, &progress, &counted);
                    scope.spawn(move || {
                        let mut totals = generator::Totals::default();
                        for source in sources.iter().skip(thread).step_by(threads) {
                            for cif in source.cifs() {
                                let done = progress.fetch_add(1, Ordering::Relaxed) + 1;
                                if done % 10_000 == 0 {
                                    println!("{done} structures");
                                }

                                let result =
                                    std::panic::catch_unwind(|| generator::structure(&cif));
                                if let Ok(Some(keys)) = result {
                                    totals.add(keys);
                                    counted.fetch_add(1, Ordering::Relaxed);
                                }
                            }
                        }
                        totals
                    })
                })
                .collect();
            workers
                .into_iter()
                .map(|worker| worker.join().unwrap())
                .reduce(generator::Totals::merge)
                .unwrap_or_default()
        });
        drop(std::panic::take_hook());

        let records = totals.records();
        let bytes = write(&records);
        let compressed = miniz_oxide::deflate::compress_to_vec_zlib(&bytes, 10);
        std::fs::write(&out, &compressed).unwrap();
        println!(
            "{} of {} structures counted; {} records, {} bytes, {} compressed",
            counted.into_inner(),
            progress.into_inner(),
            records.len(),
            bytes.len(),
            compressed.len()
        );
    }

    /// Prints, for each bond key with a coordination number that both tables `$DISTANCES_A` and
    /// `$DISTANCES_B` have with at least [`MIN_STRUCTURES`], the medians and their difference,
    /// and the median and the 10 and 90 % quantiles of those differences.
    #[test]
    #[ignore = "needs DISTANCES_A and DISTANCES_B"]
    fn compare() {
        let read = |variable: &str| {
            let path = std::env::var(variable).expect(variable);
            let compressed = std::fs::read(path).unwrap();
            Table::parse(&miniz_oxide::inflate::decompress_to_vec_zlib(&compressed).unwrap())
                .unwrap()
        };
        let (a, b) = (read("DISTANCES_A"), read("DISTANCES_B"));
        let names = |table: &Table| {
            let mut names = vec![String::new(); table.indices.len()];
            for (name, &index) in &table.indices {
                names[usize::from(index)] = name.clone();
            }
            names
        };
        let a_names = names(&a);

        let mut differences = Vec::new();
        println!(
            "centre\tcoordination\tpartner\ta structures\tb structures\ta median\tb median\tdifference"
        );
        for (&(kind, centre, coordination, partner), first) in &a.records {
            let (centre, partner) = (
                &a_names[usize::from(centre)],
                &a_names[usize::from(partner)],
            );
            if kind != Kind::Bond
                || coordination == ANY_COORDINATION
                || first.structures < MIN_STRUCTURES
            {
                continue;
            }
            let Some(second) = b.get(Kind::Bond, centre, coordination, partner) else {
                continue;
            };
            let (median_a, median_b) = (first.quantile(0.5), second.quantile(0.5));
            println!(
                "{centre}\t{coordination}\t{partner}\t{}\t{}\t{median_a:.3}\t{median_b:.3}\t{:.3}",
                first.structures,
                second.structures,
                median_a - median_b
            );
            differences.push(median_a - median_b);
        }

        differences.sort_by(f32::total_cmp);
        let quantile = |q: f32| differences[((differences.len() - 1) as f32 * q) as usize];
        println!(
            "{} keys; differences a − b: median {:.3}, 10 % {:.3}, 90 % {:.3}",
            differences.len(),
            quantile(0.5),
            quantile(0.1),
            quantile(0.9)
        );
    }

    /// The generator's reading, counting and summing.
    mod generator {
        use std::io::BufRead;

        use super::super::*;
        use crate::bonding::BondRule;

        /// The most atoms a cell may hold to be counted.
        const MAX_ATOMS: usize = 1000;
        const BINS: usize = (CUTOFF / BOND_BIN) as usize;

        /// A file of structures.
        pub(super) enum Source {
            /// One gzipped JSON Lines file of the Materials Project's summary collection.
            Summary(std::path::PathBuf),
            /// One CIF.
            Cif(std::path::PathBuf),
        }

        /// The files `$MP_SUMMARY` or `$CIF_LIST` name.
        pub(super) fn sources() -> Vec<Source> {
            if let Ok(folder) = std::env::var("MP_SUMMARY") {
                let mut files: Vec<_> = std::fs::read_dir(folder)
                    .unwrap()
                    .map(|entry| entry.unwrap().path())
                    .filter(|path| path.to_string_lossy().ends_with(".jsonl.gz"))
                    .collect();
                files.sort();
                return files.into_iter().map(Source::Summary).collect();
            }
            let list = std::env::var("CIF_LIST").expect("MP_SUMMARY or CIF_LIST");
            std::fs::read_to_string(list)
                .unwrap()
                .lines()
                .filter(|line| !line.is_empty())
                .map(|line| Source::Cif(line.into()))
                .collect()
        }

        impl Source {
            /// The CIFs of the structures to count.
            pub(super) fn cifs(&self) -> Box<dyn Iterator<Item = Vec<u8>> + '_> {
                match self {
                    Source::Cif(path) => Box::new(std::fs::read(path).ok().into_iter()),
                    Source::Summary(path) => {
                        let file = std::fs::File::open(path).unwrap();
                        let lines = std::io::BufReader::new(flate2::read::GzDecoder::new(file))
                            .lines()
                            .map_while(Result::ok);
                        Box::new(lines.filter_map(|line| {
                            let record: serde_json::Value = serde_json::from_str(&line).ok()?;
                            let counted = record["theoretical"] == false
                                && record["deprecated"] == false
                                && record["builder_meta"]["license"] == "BY-C";
                            counted.then(|| p1_cif(&record["material_id"], &record["structure"]))?
                        }))
                    }
                }
            }
        }

        /// A P1 CIF of a pymatgen `Structure` as JSON, or `None` where it lacks a field.
        fn p1_cif(id: &serde_json::Value, structure: &serde_json::Value) -> Option<Vec<u8>> {
            use std::fmt::Write;

            let lattice = &structure["lattice"];
            let mut cif = format!(
                "data_{}\n_space_group_name_H-M_alt 'P 1'\n_space_group_IT_number 1\n",
                id.as_str()?
            );
            for (tag, key) in [
                ("length_a", "a"),
                ("length_b", "b"),
                ("length_c", "c"),
                ("angle_alpha", "alpha"),
                ("angle_beta", "beta"),
                ("angle_gamma", "gamma"),
                ("volume", "volume"),
            ] {
                writeln!(cif, "_cell_{tag} {:.8}", lattice[key].as_f64()?).ok()?;
            }
            cif.push_str("loop_\n_space_group_symop_operation_xyz\n'x, y, z'\n");
            cif.push_str(
                "loop_\n_atom_site_label\n_atom_site_type_symbol\n_atom_site_fract_x\n\
                 _atom_site_fract_y\n_atom_site_fract_z\n_atom_site_occupancy\n",
            );
            for (index, site) in structure["sites"].as_array()?.iter().enumerate() {
                let [x, y, z] = [0, 1, 2].map(|axis| site["abc"][axis].as_f64());
                for species in site["species"].as_array()? {
                    let element = species["element"].as_str()?;
                    let occupancy = species["occu"].as_f64()?;
                    writeln!(
                        cif,
                        "{element}{} {element} {:.8} {:.8} {:.8} {occupancy:.6}",
                        index + 1,
                        x?,
                        y?,
                        z?
                    )
                    .ok()?;
                }
            }
            Some(cif.into_bytes())
        }

        /// A record's key: its kind, first key, coordination number and second key.
        pub(super) type Key = (Kind, String, u8, String);

        /// One structure's distances per key, in bins of [`BOND_BIN`] from 0, or `None` where it isn't
        /// counted (unreadable, partly occupied, too large).
        pub(super) fn structure(cif: &[u8]) -> Option<HashMap<Key, Vec<u32>>> {
            let (_, phase) = crate::cif_phase::try_phase_from_cif_bytes(cif)?;
            if phase
                .atoms
                .iter()
                .any(|atom| (atom.occupancy.value - 1.0).abs() > 0.01)
            {
                return None;
            }
            let phase = crate::oxidation::with_guessed_oxidation_states(&phase).into_owned();
            let atoms = crate::cell_atoms::cell_atoms(&phase);
            if atoms.is_empty() || atoms.len() > MAX_ATOMS {
                return None;
            }
            let chemistry = crate::oxidation::site_chemistry(&phase);
            let disorder = crate::disorder::of_atoms(&phase, &atoms, &chemistry);
            let bond_atoms = crate::bonding::cell_bond_atoms(&phase, &atoms, &chemistry, &disorder);
            let fractional: Vec<_> = atoms.iter().map(|atom| atom.fractional).collect();

            let mut keys: HashMap<Key, Vec<u32>> = HashMap::default();
            let mut add = |key: Key, distance: f32| {
                let bin = (distance / BOND_BIN) as usize;
                if bin < BINS {
                    keys.entry(key).or_insert_with(|| vec![0; BINS])[bin] += 1;
                }
            };

            let species: Vec<_> = bond_atoms
                .iter()
                .map(|atom| species(atom.type_symbol))
                .collect();
            let elements: Vec<_> = bond_atoms
                .iter()
                .map(|atom| element_of(atom.type_symbol))
                .collect();

            // One search out to the longest bond limit or the cutoff, whichever is longer.
            let bond_rule = BondRule::default();
            let longest_limit = bond_atoms
                .iter()
                .flat_map(|a| bond_atoms.iter().map(move |b| (a, b)))
                .map(|(a, b)| bond_rule.limit(a.type_symbol, b.type_symbol).max)
                .fold(CUTOFF, f32::max);
            let contacts =
                crate::geometry::periodic_contacts(&phase.cell, &fractional, longest_limit)?;

            let mut bonds = vec![Vec::new(); atoms.len()];
            for &(from, to, distance) in &contacts {
                if bond_rule.are_bonded(bond_atoms[from], bond_atoms[to], distance) {
                    bonds[from].push((to, distance));
                }
            }
            for (centre, bonds) in bonds.iter().enumerate() {
                let coordination = u8::try_from(bonds.len()).unwrap_or(u8::MAX);
                for &(partner, length) in bonds {
                    for (a, b) in [(&species, &species), (&elements, &elements)] {
                        for count in [coordination, ANY_COORDINATION] {
                            let key = (Kind::Bond, a[centre].clone(), count, b[partner].clone());
                            add(key, length);
                        }
                    }
                }
            }

            for &(from, to, distance) in contacts.iter().filter(|contact| contact.2 < CUTOFF) {
                for names in [&species, &elements] {
                    let (a, b) = (&names[from], &names[to]);
                    let (a, b) = if a <= b { (a, b) } else { (b, a) };
                    add(
                        (Kind::Pair, a.clone(), ANY_COORDINATION, b.clone()),
                        distance,
                    );
                }
            }

            Some(keys)
        }

        /// Each key's summed weights and structures.
        #[derive(Default)]
        pub(super) struct Totals(HashMap<Key, (u32, Vec<f64>)>);

        impl Totals {
            /// Adds one structure's distances per key, each key's weighing 1 together.
            pub(super) fn add(&mut self, keys: HashMap<Key, Vec<u32>>) {
                for (key, counts) in keys {
                    let sum: u32 = counts.iter().sum();
                    let (structures, weights) =
                        self.0.entry(key).or_insert_with(|| (0, vec![0.0; BINS]));
                    *structures += 1;
                    for (weight, count) in weights.iter_mut().zip(counts) {
                        *weight += f64::from(count) / f64::from(sum);
                    }
                }
            }

            pub(super) fn merge(mut self, other: Self) -> Self {
                for (key, (structures, weights)) in other.0 {
                    let (total, sums) = self.0.entry(key).or_insert_with(|| (0, vec![0.0; BINS]));
                    *total += structures;
                    for (sum, weight) in sums.iter_mut().zip(weights) {
                        *sum += weight;
                    }
                }
                self
            }

            /// The records of the keys with at least [`MIN_STRUCTURES`], sorted by key.
            pub(super) fn records(self) -> Vec<Record> {
                let mut records: Vec<_> = self
                    .0
                    .into_iter()
                    .filter(|(_, (structures, _))| *structures >= MIN_STRUCTURES)
                    .map(
                        |((kind, a, coordination, b), (structures, weights))| Record {
                            kind,
                            a,
                            coordination,
                            b,
                            structures,
                            weights,
                        },
                    )
                    .collect();
                records.sort_by(|x, y| {
                    (x.kind as u8, &x.a, x.coordination, &x.b).cmp(&(
                        y.kind as u8,
                        &y.a,
                        y.coordination,
                        &y.b,
                    ))
                });
                records
            }
        }
    }
}
