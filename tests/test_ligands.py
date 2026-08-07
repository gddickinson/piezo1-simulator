"""The curated modulators, and the claim the resource exists to protect.

No PIEZO structure with a bound small-molecule modulator has been deposited.
A pocket drawn on a structure looks exactly like one observed in it, so the
evidence level has to travel with the site — and the build has to *verify* the
absence rather than assert it.
"""

from __future__ import annotations

import json

import pytest

from piezo1.config import RESOURCE_DIR
from piezo1.core.ligands import SITE_EVIDENCE_ORDER, load_ligands


@pytest.fixture(scope="module")
def ligands():
    result = load_ligands()
    if not len(result):
        pytest.skip(result.note)
    return result


def test_the_expected_modulators_are_present(ligands):
    keys = {l.key for l in ligands.ligands}
    assert keys == {"yoda1", "yoda2", "dooku1", "jedi1", "jedi2", "gsmtx4"}
    assert len(ligands.by_role("activator")) == 4
    assert ligands.get("gsmtx4").kind == "peptide"


def test_no_site_claims_a_bound_structure(ligands):
    """The point of the resource. None exists, so none may be claimed."""
    assert not ligands.any_observed_site
    for ligand in ligands.ligands:
        assert ligand.site_evidence in SITE_EVIDENCE_ORDER
        assert ligand.site_evidence != "bound_structure", ligand.key
        assert not ligand.site_is_observed


def test_the_one_residue_level_site_says_it_is_inferred(ligands):
    """Yoda1's pocket is from MD, and the text has to say so where it is read."""
    sited = ligands.with_sites()
    assert [l.key for l in sited] == ["yoda1"]
    yoda1 = sited[0]
    assert yoda1.site_residues == (1718, 2075, 2078)
    assert yoda1.site_evidence == "docking_md"
    assert "INFERRED" in yoda1.site_text()
    assert "not from a bound structure" in yoda1.site_text()


def test_a_ligand_without_a_site_explains_why(ligands):
    """Silence would read as "not looked at" rather than "deliberately absent"."""
    for ligand in ligands.ligands:
        if not ligand.has_site:
            assert ligand.site_note, ligand.key
            assert "no residue-level site" in ligand.site_text()

    # GsMTx4 specifically: it acts on the bilayer, so residues would mislead.
    gsmtx4 = ligands.get("gsmtx4")
    assert "bilayer" in gsmtx4.site_note.lower()


def test_potencies_match_the_project_ground_truth(ligands):
    """The two numbers CLAUDE.md lists as scientific anchors."""
    yoda1 = ligands.get("yoda1").potency
    assert yoda1["measure"] == "EC50"
    assert yoda1["value"] == pytest.approx(26.6)
    assert yoda1["unit"] == "uM"
    assert yoda1["citation"] == "syeda2015"

    gsmtx4 = ligands.get("gsmtx4").potency
    assert gsmtx4["measure"] == "Kd"
    assert gsmtx4["value"] == pytest.approx(155.0)
    assert gsmtx4["unit"] == "nM"
    assert gsmtx4["citation"] == "bae2011"


def test_every_citation_resolves(ligands):
    references = {r["key"] for r in json.loads(
        (RESOURCE_DIR / "references.json").read_text())["references"]}
    for ligand in ligands.ligands:
        if ligand.potency:
            assert ligand.potency["citation"] in references, ligand.key
        if ligand.site_citation:
            assert ligand.site_citation in references, ligand.key


def test_chemistry_was_fetched_and_verified(ligands):
    """A wrong CID must not pass silently, so the build checks the InChIKey."""
    for ligand in ligands.ligands:
        if ligand.kind != "small_molecule":
            continue
        assert ligand.pubchem_cid, ligand.key
        assert ligand.formula, ligand.key
        assert ligand.inchikey and len(ligand.inchikey) == 27, ligand.key
    assert ligands.get("yoda1").formula == "C13H8Cl2N4S2"


def test_no_deposited_structure_contains_a_modulator():
    """The claim, verified against the structures rather than asserted.

    If a Yoda1-bound entry is ever deposited and downloaded, this fails and the
    resource is out of date — which is the correct outcome.
    """
    import sys
    from pathlib import Path

    sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
    from build_ligands import deposited_modulators

    unexpected = deposited_modulators()
    assert unexpected == [], (
        f"unexpected heteroatoms in deposited structures: {unexpected}; one "
        f"may be a bound modulator, in which case ligands.json is out of date")


def test_the_resource_states_the_limitation_at_the_top(ligands):
    assert "INFERRED" in ligands.note
    assert "no piezo structure" in ligands.note.lower()
    assert set(ligands.evidence_levels) <= set(SITE_EVIDENCE_ORDER)
