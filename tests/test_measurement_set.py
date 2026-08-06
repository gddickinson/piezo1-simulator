"""Click-to-measure interaction logic.

``MeasurementSet`` is deliberately Qt-free so the whole interaction — how many
atoms each kind needs, when a measurement completes, what gets exported — is
testable without a display.
"""

import numpy as np
import pytest

from piezo1.analysis.measure import (MEASUREMENT_KINDS, Measurement,
                                     MeasurementSet)


@pytest.fixture
def ms():
    return MeasurementSet()


def test_kinds_need_the_right_number_of_atoms():
    assert MEASUREMENT_KINDS == {"distance": 2, "angle": 3, "dihedral": 4}
    for kind, n in MEASUREMENT_KINDS.items():
        s = MeasurementSet(kind=kind)
        assert s.required == n


def test_distance_completes_on_the_second_pick(ms):
    assert ms.add_atom(0, [0, 0, 0], "A") is None
    m = ms.add_atom(1, [3, 4, 0], "B")
    assert m is not None
    assert m.kind == "distance"
    assert m.value == pytest.approx(5.0)
    assert m.units == "A"
    assert m.atoms == (0, 1)
    assert m.labels == ("A", "B")


def test_pending_clears_after_completion(ms):
    ms.add_atom(0, [0, 0, 0])
    ms.add_atom(1, [1, 0, 0])
    assert ms.pending == []
    assert len(ms.measurements) == 1


def test_angle_and_dihedral():
    s = MeasurementSet(kind="angle")
    s.add_atom(0, [1, 0, 0])
    s.add_atom(1, [0, 0, 0])
    m = s.add_atom(2, [0, 1, 0])
    assert m.value == pytest.approx(90.0)
    assert m.units == "deg"

    s.set_kind("dihedral")
    for i, p in enumerate(([1, 1, 0], [0, 1, 0], [0, 0, 0], [-1, 0, 0])):
        m = s.add_atom(i, p)
    assert abs(m.value) == pytest.approx(180.0)


def test_same_atom_twice_is_ignored(ms):
    """A double-click should not produce a zero-length distance."""
    ms.add_atom(7, [1, 2, 3], "X")
    assert ms.add_atom(7, [1, 2, 3], "X") is None
    assert ms.pending == [7]
    assert ms.measurements == []
    # A different atom still completes it.
    m = ms.add_atom(8, [1, 2, 4], "Y")
    assert m.value == pytest.approx(1.0)


def test_changing_kind_discards_a_partial_pick(ms):
    ms.add_atom(0, [0, 0, 0])
    assert ms.pending
    ms.set_kind("angle")
    assert ms.pending == []
    assert ms.kind == "angle"


def test_unknown_kind_rejected(ms):
    with pytest.raises(ValueError, match="unknown measurement kind"):
        ms.set_kind("volume")


def test_remove_and_clear(ms):
    for i in range(4):
        ms.add_atom(i, [i, 0, 0])
    assert len(ms.measurements) == 2
    ms.remove(0)
    assert len(ms.measurements) == 1
    ms.remove(99)                       # out of range must not raise
    ms.clear()
    assert ms.measurements == [] and ms.pending == []


def test_anchor_is_the_centroid_of_the_picks(ms):
    ms.add_atom(0, [0, 0, 0])
    m = ms.add_atom(1, [4, 0, 0])
    assert np.allclose(m.anchor, [2, 0, 0])


def test_anchor_of_an_empty_measurement_is_the_origin():
    assert np.allclose(Measurement("distance", 1.0, "A").anchor, [0, 0, 0])


def test_csv_export_round_trips(ms):
    ms.add_atom(0, [0, 0, 0], "CYS2411A.SG")
    ms.add_atom(1, [2.04, 0, 0], "CYS2415A.SG")
    csv_text = ms.to_csv()
    lines = csv_text.strip().splitlines()
    assert lines[0] == "kind,value,units,atoms,selection,note"
    assert "2.04" in lines[1]
    assert "CYS2411A.SG - CYS2415A.SG" in lines[1]

    import csv
    import io
    rows = list(csv.DictReader(io.StringIO(csv_text)))
    assert len(rows) == 1
    assert float(rows[0]["value"]) == pytest.approx(2.04)


def test_text_export_is_one_line_per_measurement(ms):
    for i in range(4):
        ms.add_atom(i, [i, 0, 0], f"a{i}")
    assert len(ms.to_text().splitlines()) == 2


def test_labels_default_to_the_atom_index(ms):
    ms.add_atom(11, [0, 0, 0])
    m = ms.add_atom(12, [1, 0, 0])
    assert m.labels == ("11", "12")


def test_measurement_on_real_geometry(human_structure):
    """The known C2411–C2415 disulfide must come back at ~2.04 Å.

    Independently detected at the same distance by
    :mod:`piezo1.analysis.interactions`, from a different code path.
    """
    st = human_structure
    s = MeasurementSet()
    picked = []
    for residue in (2411, 2415):
        idx = np.flatnonzero((st.chain == "A") & (st.res_seq == residue)
                             & (st.atom_name == "SG"))
        assert len(idx) == 1
        picked.append(int(idx[0]))
    s.add_atom(picked[0], st.xyz[picked[0]])
    m = s.add_atom(picked[1], st.xyz[picked[1]])
    assert m.value == pytest.approx(2.04, abs=0.05)
