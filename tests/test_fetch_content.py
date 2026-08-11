"""Round 77: a download that verifies *what* arrived, not only how much.

The size guard has been necessary and insufficient twice. Round 60 found an
Ensembl endpoint returning an HTML error because the content type had been sent
as a query parameter instead of a header. Round 65 found two 127-byte error
pages **stored as structures** — the files existed, had plausible names, and
every step downstream treated them as data.

The size check caught those two only because they happened to be tiny. A CDN
error page, a login redirect, a rate-limit notice and a maintenance page are
all comfortably over 200 bytes, so only looking at the content can catch them.

The load-bearing test here is `test_the_guard_rejects_what_the_size_check_lets
_through`: a check that only fired on things the old guard already caught would
be worth nothing.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from piezo1.io.fetch import CONTENT_CHECKS, FetchResult, _download

#: A realistic CDN error page. The point is the length: 554 bytes sails past a
#: 200-byte floor.
ERROR_PAGE = (b"<!DOCTYPE html><html><head><title>404 Not Found</title></head>"
              b"<body><h1>Not Found</h1>"
              + b"<p>The requested resource was not found on this server.</p>" * 8
              + b"</body></html>")

RATE_LIMIT = json.dumps({"error": "rate limit exceeded",
                         "detail": "x" * 400}).encode()


class _FakeResponse:
    def __init__(self, data: bytes):
        self._data = data

    def read(self) -> bytes:
        return self._data

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


@pytest.fixture
def serve(monkeypatch):
    """Make the next `_download` return exactly these bytes."""
    def _serve(payload: bytes):
        monkeypatch.setattr("urllib.request.urlopen",
                            lambda *a, **k: _FakeResponse(payload))
    return _serve


# --------------------------------------------- the calibration that matters

def test_the_error_page_is_large_enough_to_pass_the_size_check():
    """Otherwise the rest of this file proves nothing new."""
    assert len(ERROR_PAGE) > 200, "the planted page must beat the old guard"
    assert len(RATE_LIMIT) > 200


def test_the_guard_rejects_what_the_size_check_lets_through(serve, tmp_path):
    """The whole point of the round.

    An HTML error page, well over the size floor, arriving where an mmCIF was
    expected. The old guard wrote it to disk under a structure's name.
    """
    serve(ERROR_PAGE)
    dest = tmp_path / "8YEZ.cif"
    result = _download("https://example.invalid/8YEZ.cif", dest, kind="cif")

    assert not result.ok, "an HTML error page was accepted as a structure"
    assert "mmCIF" in result.error
    assert not dest.exists(), (
        "the rejected payload was written anyway; a later step would read it")


def test_without_a_kind_the_old_behaviour_is_unchanged(serve, tmp_path):
    """The guard is opt-in per fetcher, so this records what that costs."""
    serve(ERROR_PAGE)
    dest = tmp_path / "unchecked.bin"
    result = _download("https://example.invalid/x", dest)
    assert result.ok and dest.exists(), (
        "an unkinded download should still behave as it did; if this changes, "
        "every caller needs a kind")


def test_nothing_is_written_until_the_content_is_checked(serve, tmp_path):
    """A rejected download must not leave a file behind.

    This is what made Round 65's failure durable: the bad file persisted, and
    because `_download` skips anything that already exists, every later run
    returned the error page from cache without touching the network.
    """
    serve(ERROR_PAGE)
    dest = tmp_path / "9VMX.cif"
    _download("https://example.invalid/9VMX.cif", dest, kind="cif")
    assert not dest.exists()

    # And the cache-skip really would have kept serving it.
    dest.write_bytes(ERROR_PAGE)
    again = _download("https://example.invalid/9VMX.cif", dest, kind="cif")
    assert again.ok and not again.downloaded, (
        "an existing file is served from cache without being re-checked, which "
        "is why nothing may be written unverified in the first place")


# ------------------------------------------------------- each kind of payload

@pytest.mark.parametrize("kind", sorted(CONTENT_CHECKS))
def test_every_kind_rejects_an_html_error_page(kind):
    assert CONTENT_CHECKS[kind](ERROR_PAGE), (
        f"the {kind} check accepts an HTML error page")


def test_json_accepts_json_and_rejects_a_truncated_body():
    check = CONTENT_CHECKS["json"]
    assert not check(json.dumps({"a": [1, 2, 3]}).encode())
    assert check(b'{"a": [1, 2, ')
    # A rate-limit notice IS valid JSON, so the content check cannot catch it.
    # Recorded rather than papered over: the callers that parse it will.
    assert not check(RATE_LIMIT)


def test_cif_needs_coordinates_not_merely_a_header():
    check = CONTENT_CHECKS["cif"]
    assert check(b"data_8YEZ\n_entry.id 8YEZ\n" + b"# padding\n" * 40), \
        "a header-only mmCIF carries no atoms and must be rejected"
    assert not check(b"data_8YEZ\nloop_\n_atom_site.group_PDB\nATOM 1 N\n")


def test_fasta_needs_a_sequence_not_merely_a_header():
    check = CONTENT_CHECKS["fasta"]
    assert check(b">sp|Q92508|PIEZ1_HUMAN description only\n")
    assert not check(b">sp|Q92508|PIEZ1_HUMAN\nMEPHVLGAVLYWLLLPCALLAA\n")
    assert check(b"MEPHVLGAVLY\n"), "a bare sequence is not FASTA"


def test_sdf_needs_atoms_and_a_terminator():
    check = CONTENT_CHECKS["sdf"]
    assert check(b"Yoda1\n  no counts line here\n$$$$\n")
    assert check(b"Yoda1\n  4  3  0  0  0  0            999 V2000\n")
    assert not check(b"Yoda1\n  4  3  0  0  0  0            999 V2000\n$$$$\n")


# -------------------------------------------------- the real fetchers use it

def test_every_downloader_declares_what_it_expects():
    """A fetcher with no `kind` is one that can still store an error page."""
    import inspect

    from piezo1.io import fetch

    source = inspect.getsource(fetch)
    body = source[source.index("# Individual fetchers"):]
    calls = body.count("_download(")
    kinded = body.count("kind=")
    assert calls == kinded, (
        f"{calls - kinded} of {calls} download calls do not say what they "
        f"expect back; each one can still write an error page to disk")


def test_a_good_structure_is_still_accepted(serve, tmp_path):
    """The guard must not reject real data — checked on a real file.

    A checker that says no to everything is as useless as one that says yes,
    and this is the half that catches it.
    """
    from piezo1.config import STRUCTURE_DIR

    real = STRUCTURE_DIR / "8YEZ.cif"
    if not real.exists():
        pytest.skip("8YEZ not downloaded; run python -m piezo1.io.fetch")

    serve(real.read_bytes())
    dest = tmp_path / "copy.cif"
    result = _download("https://example.invalid/8YEZ.cif", dest, kind="cif")
    assert result.ok, f"a real mmCIF was rejected: {result.error}"
    assert dest.exists()


def test_a_truncated_structure_is_rejected(serve, tmp_path):
    """An unlooked-for benefit, found by getting this test wrong first.

    Writing the check above I served the first 200 kB of a real 8YEZ and it was
    refused — mmCIF metadata runs well past that before the first coordinate.
    The guard was right and the test was wrong, and the behaviour is worth
    keeping: a connection that drops mid-transfer leaves a file that opens,
    parses, and contains a fraction of the molecule.
    """
    from piezo1.config import STRUCTURE_DIR

    real = STRUCTURE_DIR / "8YEZ.cif"
    if not real.exists():
        pytest.skip("8YEZ not downloaded; run python -m piezo1.io.fetch")

    serve(real.read_bytes()[:200_000])
    dest = tmp_path / "truncated.cif"
    result = _download("https://example.invalid/8YEZ.cif", dest, kind="cif")
    assert not result.ok, "a structure truncated before its coordinates passed"
    assert not dest.exists()
