#[test]
fn test() {
    const ATOM: &str = r#"{
"label": "Ti1",
"type": "Ti",
"x": {"value": 0.5},
"y": {"value": 0.5},
"z": {"value": 0.5},
"occupancy": {"value": 1.0},
"multiplicity": 1.0,
"adp_type": "Uiso",
"u_iso_or_equiv": {"value": 0.0087},
"U11": {"value": 0.0},
"U22": {"value": 0.0},
"U33": {"value": 0.0},
"U12": {"value": 0.0},
"U13": {"value": 0.0},
"U23": {"value": 0.0}}"#;

    let atom: crystallib::Atom = serde_json::from_str(ATOM).unwrap();

    const CELL: &str = r#"{
"a": {"value": 4.0094},
"b": {"value": 4.0094},
"c": {"value": 4.0094},
"alpha": {"value": 90.0},
"beta": {"value": 90.0},
"gamma": {"value": 90.0},
"volume": {"value": 64.45},
"space_group": "P m -3 m"
}"#;

    let cell: crystallib::Cell = serde_json::from_str(CELL).unwrap();

    let phase = crystallib::Phase {
        cell,
        atoms: crystallib::Atoms(vec![atom]),
    };

    assert_eq!(phase.atoms[0].label, "Ti1");
    assert_eq!(phase.atoms[0].type_, "Ti");
    assert_eq!(phase.atoms[0].x.value, 0.5);
    assert_eq!(phase.atoms[0].y.value, 0.5);
    assert_eq!(phase.atoms[0].z.value, 0.5);
    assert_eq!(phase.atoms[0].occupancy.value, 1.0);
    assert_eq!(phase.atoms[0].multiplicity, Some(1.0));
    assert_eq!(phase.atoms[0].adp_type, crystallib::AdpType::Uiso);
    assert_eq!(phase.atoms[0].u_iso_or_equiv.value, 0.0087);
    assert_eq!(phase.atoms[0].u11.value, 0.0);
    assert_eq!(phase.atoms[0].u22.value, 0.0);
    assert_eq!(phase.atoms[0].u33.value, 0.0);
    assert_eq!(phase.atoms[0].u12.value, 0.0);
    assert_eq!(phase.atoms[0].u13.value, 0.0);
    assert_eq!(phase.atoms[0].u23.value, 0.0);

    assert_eq!(phase.cell.a.value, 4.0094);
    assert_eq!(phase.cell.b.value, 4.0094);
    assert_eq!(phase.cell.c.value, 4.0094);
    assert_eq!(phase.cell.alpha.value, 90.0);
    assert_eq!(phase.cell.beta.value, 90.0);
    assert_eq!(phase.cell.gamma.value, 90.0);
    assert_eq!(phase.cell.volume.value, 64.45);
    assert_eq!(phase.cell.space_group, "P m -3 m");
}
