//! How a type symbol names a charge: its element (or the CIF's own type) followed by the
//! magnitude and sign of the oxidation state (`Fe3+`, `O2-`, `Fe2.5+`, `Na+`, also `Fe+3`).
//! This is the one place that reads and writes that format.

/// The sign of an ion's charge.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Deserialize, serde::Serialize))]
pub enum Charge {
    /// A cation.
    Positive,
    /// An anion.
    Negative,
}

impl Charge {
    /// The sign of the charge a type symbol names, if it names one that isn't zero: `C0+` names
    /// none (see [`state_of`]).
    pub fn of(type_symbol: &str) -> Option<Self> {
        let state = state_of(type_symbol)?;
        if state > 0.0 {
            Some(Self::Positive)
        } else if state < 0.0 {
            Some(Self::Negative)
        } else {
            None
        }
    }
}

/// The oxidation state a type symbol names at its end (`Fe3+` 3, `O2-` −2, `Na+` 1, `Fe+3` 3,
/// `Fe2.5+` 2.5, `C0+` 0), or `None` if it names none.
pub fn state_of(type_symbol: &str) -> Option<f64> {
    split(type_symbol).1
}

/// A type symbol without the state it names at its end (`Fe3+` is `Fe`, `B0+` is `B`), and that
/// state; a type symbol naming none is returned whole.
pub fn split(type_symbol: &str) -> (&str, Option<f64>) {
    let symbol = type_symbol.trim_end();
    let is_number = |c: char| c.is_ascii_digit() || c == '.';
    let whole = (type_symbol, None);

    // `Fe3+` (digits, then the sign) or `Fe+3` (the sign, then digits).
    let (base, sign, magnitude) = match symbol.strip_suffix(['+', '-']) {
        Some(rest) => {
            let base = rest.trim_end_matches(is_number);
            (base, symbol.chars().last(), &rest[base.len()..])
        }
        None => {
            let rest = symbol.trim_end_matches(is_number);
            match rest.strip_suffix(['+', '-']) {
                Some(base) => (base, rest.chars().last(), &symbol[rest.len()..]),
                None => return whole,
            }
        }
    };

    let magnitude = if magnitude.is_empty() {
        1.0
    } else {
        match magnitude.parse::<f64>() {
            Ok(magnitude) => magnitude,
            Err(_) => return whole,
        }
    };
    match sign {
        Some('+') => (base, Some(magnitude)),
        Some('-') => (base, Some(-magnitude)),
        _ => whole,
    }
}

/// `symbol` with the oxidation state `state` appended: `Fe` and 3 give `Fe3+`, `O` and −2
/// `O2-`, `Fe` and 2.5 `Fe2.5+`. A zero state leaves it as it is.
pub fn with_state(symbol: &str, state: f64) -> String {
    if state == 0.0 {
        return symbol.to_string();
    }
    let sign = if state > 0.0 { '+' } else { '-' };
    format!("{symbol}{}{sign}", state.abs())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_type_symbol_names_its_state_at_its_end() {
        assert_eq!(state_of("Fe3+"), Some(3.0));
        assert_eq!(state_of("O2-"), Some(-2.0));
        assert_eq!(state_of("Na+"), Some(1.0));
        assert_eq!(state_of("Fe+3"), Some(3.0));
        assert_eq!(state_of("Fe2.5+"), Some(2.5));
        assert_eq!(state_of("Cl1-"), Some(-1.0));
        assert_eq!(state_of("C0+"), Some(0.0));
        assert_eq!(state_of("Fe"), None);
        assert_eq!(state_of("O1"), None);
        assert_eq!(state_of(""), None);
    }

    #[test]
    fn a_type_symbol_names_its_charge_at_its_end() {
        assert_eq!(Charge::of("Fe3+"), Some(Charge::Positive));
        assert_eq!(Charge::of("O2-"), Some(Charge::Negative));
        assert_eq!(Charge::of("Na+"), Some(Charge::Positive));
        assert_eq!(Charge::of("Fe+3"), Some(Charge::Positive));
        assert_eq!(Charge::of("Fe2.5+"), Some(Charge::Positive));
        assert_eq!(Charge::of("Cl1-"), Some(Charge::Negative));
    }

    #[test]
    fn a_type_symbol_without_a_charge_or_with_zero_names_none() {
        assert_eq!(Charge::of("C0+"), None);
        assert_eq!(Charge::of("C+0"), None);
        assert_eq!(Charge::of("Fe"), None);
        assert_eq!(Charge::of("O1"), None);
        assert_eq!(Charge::of(""), None);
    }

    #[test]
    fn a_type_symbol_splits_into_its_base_and_state() {
        assert_eq!(split("Fe3+"), ("Fe", Some(3.0)));
        assert_eq!(split("Fe+3"), ("Fe", Some(3.0)));
        assert_eq!(split("B0+"), ("B", Some(0.0)));
        assert_eq!(split("Na+"), ("Na", Some(1.0)));
        assert_eq!(split("Fe"), ("Fe", None));
        assert_eq!(split("O1"), ("O1", None));
    }

    #[test]
    fn a_written_state_reads_back() {
        assert_eq!(with_state("Fe", 3.0), "Fe3+");
        assert_eq!(with_state("O", -2.0), "O2-");
        assert_eq!(with_state("Fe", 2.5), "Fe2.5+");
        assert_eq!(with_state("C", 0.0), "C");
        for state in [-3.0, -1.0, 1.0, 2.5, 6.0] {
            assert_eq!(state_of(&with_state("X", state)), Some(state));
        }
    }
}
