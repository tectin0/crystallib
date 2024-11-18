#[test]
fn test() {
    const ATOM: &str = r#"{
"label": "Ti1",
"type": "Ti",
"x": 0.5,
"y": 0.5,
"z": 0.5,
"occupancy": 1.0,
"multiplicity": 1.0,
"adp_type": "Uiso",
"u_iso_or_equiv": 0.0087,
"U11": 0.0,
"U22": 0.0,
"U33": 0.0,
"U12": 0.0,
"U13": 0.0,
"U23": 0.0}"#;

    let atom: crystallib::Atom = serde_json::from_str(ATOM).unwrap();

    const CELL: &str = r#"{
"a": 4.0094,
"b": 4.0094,
"c": 4.0094,
"alpha": 90.0,
"beta": 90.0,
"gamma": 90.0,
"volume": 64.45,
"space_group": "P m -3 m"
}"#;

    let cell: crystallib::Cell = serde_json::from_str(CELL).unwrap();

    let phase = crystallib::Phase {
        cell,
        atoms: crystallib::Atoms(vec![atom]),
    };

    assert_eq!(phase.atoms[0].label, "Ti1");
    assert_eq!(phase.atoms[0].type_, "Ti");
    assert_eq!(phase.atoms[0].x, 0.5);
    assert_eq!(phase.atoms[0].y, 0.5);
    assert_eq!(phase.atoms[0].z, 0.5);
    assert_eq!(phase.atoms[0].occupancy, 1.0);
    assert_eq!(phase.atoms[0].multiplicity, 1.0);
    assert_eq!(phase.atoms[0].adp_type, "Uiso");
    assert_eq!(phase.atoms[0].u_iso_or_equiv, 0.0087);
    assert_eq!(phase.atoms[0].u11, 0.0);
    assert_eq!(phase.atoms[0].u22, 0.0);
    assert_eq!(phase.atoms[0].u33, 0.0);
    assert_eq!(phase.atoms[0].u12, 0.0);
    assert_eq!(phase.atoms[0].u13, 0.0);
    assert_eq!(phase.atoms[0].u23, 0.0);

    assert_eq!(phase.cell.a, 4.0094);
    assert_eq!(phase.cell.b, 4.0094);
    assert_eq!(phase.cell.c, 4.0094);
    assert_eq!(phase.cell.alpha, 90.0);
    assert_eq!(phase.cell.beta, 90.0);
    assert_eq!(phase.cell.gamma, 90.0);
    assert_eq!(phase.cell.volume, 64.45);
    assert_eq!(phase.cell.space_group, "P m -3 m");
}
