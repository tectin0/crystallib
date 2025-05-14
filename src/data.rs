use std::{collections::HashMap, sync::LazyLock};

pub const VAN_DER_WAALS_RADII_RAW: &[u8] =
    include_bytes!("../assets/data/van-der-Waals-radii.toml");

pub struct VanDerWaalsRadii(LazyLock<HashMap<String, f32>>);

impl VanDerWaalsRadii {
    pub const fn load() -> Self {
        Self(LazyLock::new(|| {
            let data = std::str::from_utf8(VAN_DER_WAALS_RADII_RAW)
                .expect("Failed to convert van_der_Waals data to str");
            toml::from_str(data).expect("Failed to parse van_der_Waals data form toml")
        }))
    }

    pub fn get(&self, atom_type: &str) -> &f32 {
        self.0.get(atom_type).unwrap_or_else(|| {
            log::warn!(
                "Unknown atom type: {}. Using default radius of 1.0.",
                atom_type
            );
            &1.0
        })
    }
}

pub static VAN_DER_WAALS_RADII: VanDerWaalsRadii = VanDerWaalsRadii::load();

#[cfg(test)]
#[test]
fn test_van_der_waals_radii() {
    assert_eq!(VAN_DER_WAALS_RADII.get("H"), &1.2);
}

pub const COVALENT_RADII_RAW: &[u8] = include_bytes!("../assets/data/covalent-radii.toml");

pub struct CovalentRadii(LazyLock<HashMap<String, f32>>);

impl CovalentRadii {
    pub const fn load() -> Self {
        Self(LazyLock::new(|| {
            let data = std::str::from_utf8(COVALENT_RADII_RAW)
                .expect("Failed to convert covalent radii data to str");
            toml::from_str(data).expect("Failed to parse covalent radii data form toml")
        }))
    }

    pub fn get(&self, atom_type: &str) -> &f32 {
        self.0.get(atom_type).unwrap_or_else(|| {
            log::warn!(
                "Unknown atom type: {}. Using default radius of 1.0.",
                atom_type
            );
            &1.0
        })
    }
}

pub static COVALENT_RADII: CovalentRadii = CovalentRadii::load();

#[cfg(test)]
#[test]
fn test_covalent_radii() {
    assert_eq!(COVALENT_RADII.get("H"), &0.31);
}
