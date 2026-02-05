import os
import textwrap
import pytest

from crunchflow.input.databases.aqueous import AqueousDatabase


def _write(tmp_path, name, text):
    p = tmp_path / name
    p.write_text(textwrap.dedent(text), encoding="utf-8")
    return p


def test_read_example_db_smoke():
    """
    Reads the provided sample database and asserts we parsed something sensible.
    """
    src = os.path.join("tests", "data", "AqueousKinetics.dbs")
    db = AqueousDatabase.from_file(src)

    # Basic smoke checks
    assert isinstance(db.reactions, dict)
    assert isinstance(db.kinetics, dict)

    # Expect 26 and 24 reactions and kinetics, respectively
    assert len(db.reactions) == 26
    assert len(db.kinetics) >= 24

    # Check FeIIoxidation reaction
    assert 'FeIIoxidation' in db.reactions.keys()
    assert db.reactions['FeIIoxidation'].name == 'FeIIoxidation'
    assert db.reactions['FeIIoxidation'].type == 'catabolic'
    assert db.reactions['FeIIoxidation'].logK == 8.4725
    assert db.reactions['FeIIoxidation'].reaction.stoich == {'Fe+++': 1.0, 'Fe++': -1.0, "'O2(aq)'": -0.25, 'H+': -1.0}


def test_interleaved_blocks_parse(tmp_path):
    """
    Interleaved &Aqueous and &AqueousKinetics blocks should parse correctly.
    """
    content = """
    &Aqueous
      name          = AceCrO42HCO3Cr3
      type          = catabolic
      stoichiometry = -1.0 Acetate -1.0 CrO4-- 1.0 HCO3- 1.0 Cr+++ 1.0 H2O
      keq           = 24.4994
    /
    &AqueousKinetics
      name          = AceCrO42HCO3Cr3
      label         = default
      type          = monod
      rate25C       = 2.0e3
      monod_terms   = 'tot_Acetate' 2.03E-5 'tot_CrO4--' 1.06E-5
    /
    &Aqueous
      name          = Sulfate_reduction
      type          = catabolic
      stoichiometry = -0.05625 SO4-- -0.125 Acetate 0.05625 'H2S(aq)'
      keq           = 500.0
    /
    """
    path = _write(tmp_path, "interleaved.dbs", content)
    db = AqueousDatabase.from_file(str(path))

    assert "AceCrO42HCO3Cr3" in db.reactions
    assert "Sulfate_reduction" in db.reactions
    assert "AceCrO42HCO3Cr3" in db.kinetics

    kin = db.kinetics["AceCrO42HCO3Cr3"]
    assert kin.reaction_type in ("monod", "MonodBiomass")
    assert kin.rate25C == pytest.approx(2.0e3)
    assert kin.monod_terms and len(kin.monod_terms) == 2


def test_edit_stoichiometry_and_monod_and_write(tmp_path):
    """
    Change stoichiometry (&Aqueous), rate25C and monod terms (&AqueousKinetics), then write.
    """
    src = os.path.join("tests", "data", "AqueousKinetics.dbs")
    db = AqueousDatabase.from_file(src)

    # Edit stoichiometry
    db.reactions["Sulfate_reduction"].reaction.stoich['SO4--'] = -0.06625
    db.reactions["Sulfate_reduction"].reaction.stoich["'H2S(aq)'"] = 0.06625

    # Edit kinetics
    db.kinetics["Sulfate_reduction"].rate25C = 1e5
    db.kinetics["Sulfate_reduction"].monod_terms['tot_Acetate'] = 0.0001

    # Write to file
    out = tmp_path / "edit-me.out.dbs"
    db.write(str(out))

    # Read back in and check
    new_db = AqueousDatabase.from_file(str(out))
    correct_stoich = {
        'SO4--': -0.06625,
        'Acetate': -0.125,
        'H+': -0.14,
        "'NH3(aq)'": -0.0275,
        'H2O': -0.07125,
        'C5H7O2NSO4': 0.0275,
        "'CO2(aq)'": 0.015,
        'HCO3-': 0.0975,
        "'H2S(aq)'": 0.06625}
    correct_monod_terms = {'tot_Acetate': 0.0001, 'tot_SO4--': 0.001}

    assert new_db.reactions["Sulfate_reduction"].reaction.stoich == correct_stoich
    assert new_db.kinetics["Sulfate_reduction"].rate25C == 1e5
    assert new_db.kinetics["Sulfate_reduction"].monod_terms == correct_monod_terms


def test_round_trip_preserves_key_fields(tmp_path):
    """
    Read, write, read: key kinetic fields should persist.
    """
    content = """
    &Aqueous
      name          = CH4_SO4
      type          = catabolic
      stoichiometry = -1.0 'Methane(aq)' -1.0 SO4-- 1.0 HCO3- 1.0 HS- 1.0 H2O
      keq           = 3.0
    /
    &AqueousKinetics
      name          = CH4_SO4
      label         = default
      type          = monod
      rate25C       = 1.8e-8
      monod_terms   = 'tot_SO4--' 1.6E-3 'tot_Methane(aq)' 1.0E-4
    /
    """
    p = _write(tmp_path, "rt.dbs", content)
    db1 = AqueousDatabase.from_file(str(p))

    tmp = tmp_path / "rt.out.dbs"
    db1.write(str(tmp))
    db2 = AqueousDatabase.from_file(str(tmp))

    # Compare essentials
    assert "CH4_SO4" in db2.reactions and "CH4_SO4" in db2.kinetics
    k1 = db1.kinetics["CH4_SO4"]; k2 = db2.kinetics["CH4_SO4"]
    assert k2.reaction_type == k1.reaction_type
    assert k2.rate25C == pytest.approx(k1.rate25C)
    assert k2.monod_terms == k1.monod_terms
