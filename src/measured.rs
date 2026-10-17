/// A measured or refined number with its standard uncertainty (esd), where one is known.
#[cfg_attr(feature = "serde", derive(serde::Deserialize, serde::Serialize))]
#[cfg_attr(feature = "serde", serde(default))]
#[cfg_attr(feature = "get-size", derive(get_size2::GetSize))]
#[derive(Debug, Default, Clone, Copy, PartialEq, PartialOrd)]
pub struct Measured {
    pub value: f64,
    pub esd: Option<f64>,
}

impl Measured {
    /// A value without an esd.
    pub const fn new(value: f64) -> Self {
        Self { value, esd: None }
    }

    /// Multiplies the value and its esd by `factor`.
    pub fn scaled(self, factor: f64) -> Self {
        Self {
            value: self.value * factor,
            esd: self.esd.map(|esd| (esd * factor).abs()),
        }
    }
}

impl From<f64> for Measured {
    fn from(value: f64) -> Self {
        Self::new(value)
    }
}

#[cfg(test)]
mod tests {
    use super::Measured;

    #[test]
    fn scaling_scales_the_esd_and_keeps_it_positive() {
        let measured = Measured {
            value: 2.0,
            esd: Some(0.1),
        };
        assert_eq!(
            measured.scaled(-3.0),
            Measured {
                value: -6.0,
                esd: Some(0.30000000000000004),
            }
        );
        assert_eq!(Measured::new(1.0).scaled(2.0), Measured::new(2.0));
    }
}
