//! Reading a [`Phase`] out of a CIF document.
//!
//! Tags are matched by their current names (see [`cif::dictionary::canonical`]), so an older tag
//! (`_symmetry_Int_Tables_number`) or one written the DDLm way (`_atom_site.label`) is read like
//! the current one.

use crate::type_symbol::{Charge, split, with_state};
use crate::{
    AdpType, Atom, Atoms, Cell, IntoSpaceGroupNumber, IntoSpaceGroupSymbol, Measured, Phase,
};
use cif::dictionary::{B_PER_U, canonical, type_symbol_from_label};
use cif::{DataBlock, Item, Loop, Value};

/// The name and phase of the first data block of the CIF `bytes` that describes a phase, or
/// `None` if the bytes aren't a CIF or none of its blocks describes one.
pub fn try_phase_from_cif_bytes(bytes: &[u8]) -> Option<(String, Phase)> {
    let document = match cif::read(&String::from_utf8_lossy(bytes)) {
        Ok(document) => document,
        Err(error) => {
            log::debug!("Not a CIF: {error}");
            return None;
        }
    };

    document
        .blocks
        .iter()
        .find_map(|block| match phase_of(block) {
            Ok(phase) => Some((block.name.clone(), phase)),
            Err(reason) => {
                log::debug!("Data block {} describes no phase: {reason}", block.name);
                None
            }
        })
}

/// The phase `block` describes: its cell, space group and atom sites.
pub fn phase_of(block: &DataBlock) -> Result<Phase, String> {
    Ok(Phase {
        cell: cell_of(block)?,
        atoms: atoms_of(block)?,
    })
}

fn cell_of(block: &DataBlock) -> Result<Cell, String> {
    let number = |tag: &str| {
        value(block, tag)
            .and_then(measured)
            .ok_or_else(|| format!("no number for {tag}"))
    };

    let symbol = value(block, "_space_group_name_H-M_alt")
        .and_then(Value::as_str)
        .map(str::to_string);
    let space_group_number = value(block, "_space_group_IT_number")
        .and_then(Value::as_str)
        .and_then(|text| text.parse::<u8>().ok());

    let (space_group, space_group_number) = match (symbol, space_group_number) {
        (Some(symbol), Some(number)) => (symbol, number),
        (Some(symbol), None) => {
            let number = symbol
                .into_space_group_number()
                .ok_or_else(|| format!("unknown space group symbol {symbol}"))?;
            (symbol, number)
        }
        (None, Some(number)) => {
            let symbol = number
                .into_space_group_symbol()
                .ok_or_else(|| format!("unknown space group number {number}"))?;
            (symbol.to_string(), number)
        }
        (None, None) => return Err("no space group symbol or number".to_string()),
    };

    Ok(Cell {
        a: number("_cell_length_a")?,
        b: number("_cell_length_b")?,
        c: number("_cell_length_c")?,
        alpha: number("_cell_angle_alpha")?,
        beta: number("_cell_angle_beta")?,
        gamma: number("_cell_angle_gamma")?,
        volume: number("_cell_volume")?,
        space_group,
        space_group_number,
        symmetry_operations: symmetry_operations(block),
    })
}

/// The symmetry operations `block` lists (`_space_group_symop_operation_xyz`, or the older
/// `_symmetry_equiv_pos_as_xyz`), as written; none if it lists none.
fn symmetry_operations(block: &DataBlock) -> Vec<String> {
    const TAG: &str = "_space_group_symop_operation_xyz";

    match loop_with(block, TAG) {
        Some(operations) => {
            let Some(column) = column(operations, TAG) else {
                return Vec::new();
            };
            operations
                .rows
                .iter()
                .filter_map(|row| Some(row.get(column)?.as_str()?.to_string()))
                .collect()
        }
        None => value(block, TAG)
            .and_then(Value::as_str)
            .map(|operation| vec![operation.to_string()])
            .unwrap_or_default(),
    }
}

/// The atom sites of `block`'s `_atom_site_label` loop, each with its anisotropic displacement
/// parameters from the row of the `_atom_site_aniso_label` loop that names it.
fn atoms_of(block: &DataBlock) -> Result<Atoms, String> {
    let sites = loop_with(block, "_atom_site_label").ok_or("no atom sites")?;
    let aniso = loop_with(block, "_atom_site_aniso_label");
    let oxidation_numbers = oxidation_numbers(block);

    let atoms = sites
        .rows
        .iter()
        .map(|row| {
            let site_value = |tag: &str| column(sites, tag).and_then(|c| row.get(c));
            let required = |tag: &str| {
                site_value(tag)
                    .and_then(measured)
                    .ok_or_else(|| format!("an atom site has no number for {tag}"))
            };

            let label = site_value("_atom_site_label")
                .and_then(Value::as_str)
                .ok_or("an atom site has no label")?
                .to_string();

            let type_ = site_value("_atom_site_type_symbol")
                .and_then(Value::as_str)
                .map_or_else(|| type_symbol_from_label(&label), str::to_string);
            let type_ = with_oxidation_number(type_, &oxidation_numbers);

            let u_iso_or_equiv = site_value("_atom_site_U_iso_or_equiv")
                .and_then(measured)
                .or_else(|| {
                    site_value("_atom_site_B_iso_or_equiv")
                        .and_then(measured)
                        .map(|b| b.scaled(1.0 / B_PER_U))
                })
                .unwrap_or_default();

            let adp_type = site_value("_atom_site_adp_type")
                .and_then(Value::as_str)
                .and_then(|text| text.parse::<AdpType>().ok())
                .unwrap_or_default();

            let aniso_row = aniso.and_then(|l| {
                let label_column = column(l, "_atom_site_aniso_label")?;
                let row = l
                    .rows
                    .iter()
                    .find(|row| row.get(label_column).and_then(Value::as_str) == Some(&label))?;
                Some((l, row))
            });
            let u = |tag: &str| {
                aniso_row
                    .and_then(|(l, row)| column(l, tag).and_then(|c| row.get(c)))
                    .and_then(measured)
                    .unwrap_or_default()
            };

            Ok(Atom {
                x: required("_atom_site_fract_x")?,
                y: required("_atom_site_fract_y")?,
                z: required("_atom_site_fract_z")?,
                // The CIF dictionary's default: a site the file gives no occupancy is full.
                occupancy: site_value("_atom_site_occupancy")
                    .and_then(measured)
                    .unwrap_or(Measured::new(1.0)),
                multiplicity: site_value("_atom_site_site_symmetry_multiplicity")
                    .and_then(Value::number)
                    .map(|(multiplicity, _)| multiplicity),
                disorder_assembly: site_value("_atom_site_disorder_assembly")
                    .and_then(Value::as_str)
                    .map(str::to_string),
                disorder_group: site_value("_atom_site_disorder_group")
                    .and_then(Value::as_str)
                    .map(str::to_string),
                adp_type,
                u_iso_or_equiv,
                u11: u("_atom_site_aniso_U_11"),
                u22: u("_atom_site_aniso_U_22"),
                u33: u("_atom_site_aniso_U_33"),
                u12: u("_atom_site_aniso_U_12"),
                u13: u("_atom_site_aniso_U_13"),
                u23: u("_atom_site_aniso_U_23"),
                label,
                type_,
            })
        })
        .collect::<Result<Vec<_>, String>>()?;

    Ok(Atoms(atoms))
}

/// The oxidation number `_atom_type_oxidation_number` gives each `_atom_type_symbol`, written as a
/// loop or, for a single type, as two pairs.
pub fn oxidation_numbers(block: &DataBlock) -> Vec<(String, f64)> {
    let entry = |symbol: Option<&Value>, number: Option<&Value>| {
        Some((symbol?.as_str()?.to_string(), number?.number()?.0))
    };

    match loop_with(block, "_atom_type_symbol") {
        Some(types) => {
            let (Some(symbol), Some(number)) = (
                column(types, "_atom_type_symbol"),
                column(types, "_atom_type_oxidation_number"),
            ) else {
                return Vec::new();
            };
            types
                .rows
                .iter()
                .filter_map(|row| entry(row.get(symbol), row.get(number)))
                .collect()
        }
        None => entry(
            value(block, "_atom_type_symbol"),
            value(block, "_atom_type_oxidation_number"),
        )
        .into_iter()
        .collect(),
    }
}

/// `type_symbol` with the charge its oxidation number gives it (`Ba` with 2 is `Ba2+`, `O` with
/// -2 is `O2-`, `B0+` with 0.83 is `B0.83+`, see [`with_state`]). A type symbol that already
/// names a charge other than zero, or whose oxidation number is zero or not given, is returned as
/// it is.
fn with_oxidation_number(type_symbol: String, oxidation_numbers: &[(String, f64)]) -> String {
    if Charge::of(&type_symbol).is_some() {
        return type_symbol;
    }

    match oxidation_numbers
        .iter()
        .find(|(symbol, _)| *symbol == type_symbol)
    {
        Some(&(_, number)) if number != 0.0 => with_state(split(&type_symbol).0, number),
        _ => type_symbol,
    }
}

fn measured(value: &Value) -> Option<Measured> {
    value.number().map(|(value, esd)| Measured { value, esd })
}

/// The value of `tag`: of its pair, or of its loop if that has a single row.
fn value<'a>(block: &'a DataBlock, tag: &str) -> Option<&'a Value> {
    let tag = canonical(tag);

    block.items.iter().find_map(|item| match item {
        Item::Pair { tag: t, value } if canonical(t) == tag => Some(value),
        Item::Loop(l) if l.rows.len() == 1 => column(l, &tag).and_then(|c| l.rows[0].get(c)),
        _ => None,
    })
}

fn loop_with<'a>(block: &'a DataBlock, tag: &str) -> Option<&'a Loop> {
    block.items.iter().find_map(|item| match item {
        Item::Loop(l) if column(l, tag).is_some() => Some(l),
        _ => None,
    })
}

fn column(l: &Loop, tag: &str) -> Option<usize> {
    let tag = canonical(tag);
    l.tags.iter().position(|t| canonical(t) == tag)
}

#[cfg(test)]
mod tests {
    use super::*;

    const SITES: &str = "
data_test
_cell_length_a 4
_cell_length_b 4
_cell_length_c 4
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
_cell_volume 64
";

    fn phase(text: &str) -> Phase {
        try_phase_from_cif_bytes(text.as_bytes())
            .expect("the CIF describes a phase")
            .1
    }

    #[test]
    fn older_and_ddlm_tags_read_like_current_ones() {
        let phase = phase(&format!(
            "{SITES}
_symmetry_Int_Tables_number 225
loop_
_atom_site.label
_atom_site.fract_x
_atom_site.fract_y
_atom_site.fract_z
_atom_site.occupancy
_atom_site_thermal_displace_type
Ni1 0 0 0 1 Uiso
"
        ));

        assert_eq!(phase.cell.space_group_number, 225);
        assert!(!phase.cell.space_group.is_empty());

        let ni = &phase.atoms.0[0];
        assert_eq!(ni.label, "Ni1");
        assert_eq!(ni.type_, "Ni", "the type symbol comes from the label");
        assert_eq!(ni.adp_type, AdpType::Uiso);
    }

    #[test]
    fn b_iso_is_read_as_u_iso() {
        let phase = phase(&format!(
            "{SITES}
_space_group_IT_number 225
loop_
_atom_site_label
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
_atom_site_B_iso_or_equiv
Ni1 0 0 0 1 0.414(8)
"
        ));

        let u = phase.atoms.0[0].u_iso_or_equiv;
        assert!((u.value - 0.414 / B_PER_U).abs() < 1e-12);
        assert!((u.esd.unwrap() - 0.008 / B_PER_U).abs() < 1e-12);
    }

    #[test]
    fn anisotropic_parameters_belong_to_the_site_their_row_names() {
        let phase = phase(&format!(
            "{SITES}
_space_group_IT_number 1
loop_
_atom_site_label
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
_atom_site_adp_type
O1 0 0 0 1 Uiso
Ni1 0.5 0.5 0.5 1 Uani
loop_
_atom_site_aniso_label
_atom_site_aniso_U_11
_atom_site_aniso_U_22
_atom_site_aniso_U_33
Ni1 0.01 0.02 0.03
"
        ));

        let [o, ni] = &phase.atoms.0[..] else {
            panic!("two sites");
        };
        assert_eq!(o.u11, Measured::default());
        assert_eq!(ni.adp_type, AdpType::Uani);
        assert_eq!(ni.u11.value, 0.01);
        assert_eq!(ni.u33.value, 0.03);
    }

    #[test]
    fn a_block_without_a_space_group_describes_no_phase() {
        let text = format!(
            "{SITES}
loop_
_atom_site_label
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
Ni1 0 0 0 1
"
        );

        assert!(try_phase_from_cif_bytes(text.as_bytes()).is_none());
        assert!(try_phase_from_cif_bytes(b"this is not cif").is_none());
    }

    fn types(text: &str) -> Vec<String> {
        phase(text)
            .atoms
            .0
            .into_iter()
            .map(|atom| atom.type_)
            .collect()
    }

    #[test]
    fn a_type_takes_the_charge_its_oxidation_number_gives() {
        let types = types(&format!(
            "{SITES}
_space_group_IT_number 221
loop_
_atom_type_symbol
_atom_type_oxidation_number
Ba 2
O -2
Ti 4
Fe 2.5
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
Ba1 Ba 0 0 0 1
Ti1 Ti 0.5 0.5 0.5 1
O1 O 0.5 0.5 0 1
Fe1 Fe 0.25 0.25 0.25 1
"
        ));

        assert_eq!(types, ["Ba2+", "Ti4+", "O2-", "Fe2.5+"]);
    }

    #[test]
    fn a_type_naming_a_charge_or_with_oxidation_number_zero_is_kept() {
        let types = types(&format!(
            "{SITES}
_space_group_IT_number 1
loop_
_atom_type_symbol
_atom_type_oxidation_number
C0+ 0
Na 0
Fe3+ 3
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
C1 C0+ 0 0 0 1
Na1 Na 0.5 0 0 1
Fe1 Fe3+ 0 0.5 0 1
Cl1 Cl 0 0 0.5 1
"
        ));

        assert_eq!(types, ["C0+", "Na", "Fe3+", "Cl"]);
    }

    #[test]
    fn a_type_naming_zero_takes_its_oxidation_number() {
        // As the ICSD writes boron in a closo-borate: type symbol B0+, oxidation number 0.83.
        let types = types(&format!(
            "{SITES}
_space_group_IT_number 1
loop_
_atom_type_symbol
_atom_type_oxidation_number
B0+ 0.83
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
B1 B0+ 0 0 0 1
"
        ));

        assert_eq!(types, ["B0.83+"]);
    }

    #[test]
    fn a_single_type_written_as_pairs_gives_its_oxidation_number() {
        let types = types(&format!(
            "{SITES}
_space_group_IT_number 225
_atom_type_symbol Ni
_atom_type_oxidation_number 2
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
Ni1 Ni 0 0 0 1
"
        ));

        assert_eq!(types, ["Ni2+"]);
    }

    #[test]
    fn a_block_keeps_its_symmetry_operations_as_written() {
        let operations = |tag: &str| {
            phase(&format!(
                "{SITES}
_space_group_IT_number 2
loop_
_space_group_symop_id
{tag}
1 'x, y, z'
2 '-x, -y, -z'
loop_
_atom_site_label
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
Ni1 0 0 0 1
"
            ))
            .cell
            .symmetry_operations
        };

        assert_eq!(
            operations("_space_group_symop_operation_xyz"),
            ["x, y, z", "-x, -y, -z"]
        );
        assert_eq!(
            operations("_symmetry_equiv_pos_as_xyz"),
            ["x, y, z", "-x, -y, -z"]
        );
    }

    #[test]
    fn a_block_without_symmetry_operations_has_none() {
        let phase = phase(&format!(
            "{SITES}
_space_group_IT_number 225
loop_
_atom_site_label
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
Ni1 0 0 0 1
"
        ));

        assert!(phase.cell.symmetry_operations.is_empty());
    }
}
