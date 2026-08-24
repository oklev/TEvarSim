"""
Tests that a library's per-element annotation survives into the TE pool.

Simulate reads only the pool FASTA that TErandom writes -- it never opens --consensus --
so anything the library says about an element has to be carried onto the pool records or
it is lost. The pool writer used to pass description="" and drop it.

The synthesised full-length LTR record is the awkward case. TErandom splices a -LTR and
an -I record into one NAME-FULL element, and it used to describe itself as "Merged
LTR-I-LTR sequence for X", which would overwrite whatever the library had said. It takes
the LTR part's annotation instead: a target site duplication is made by the integration
machinery cutting the host at the element's ends, so it is a property of the LTR, and the
-LTR record is where a library states it.

No external test runner is required::

    PYTHONPATH=. python tests/test_pool_headers.py

The test functions are also discoverable by pytest.
"""
import os
import sys
import tempfile

from Bio import SeqIO

# Import the working-tree package regardless of any installed copy.
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from TEvarSim import build_pool  # noqa: E402
from TEvarSim.utils import description_tail, parse_tsd_tag  # noqa: E402


class _Args:
    def __init__(self, **kw):
        self.__dict__.update(kw)


ELEMENT = "ACGTTGCA" * 40
INTERNAL = "GGTTAACC" * 60


def _build(d, records, nINS=8):
    """Run TEPoolBuilder over a library and hand back the pool records it wrote.

    Every mutation rate is zeroed so the IDs are predictable: with an empty modification
    suffix each pool ID is the library ID plus a trailing underscore, which is existing
    behaviour and not what these tests are about.
    """
    consensus = os.path.join(d, "library.fa")
    with open(consensus, "w") as f:
        for header, seq in records:
            f.write(f"{header}\n{seq}\n")
    args = _Args(consensus=consensus, nINS=nINS, outprefix=os.path.join(d, "pool"),
                 snp_rate=0.0, indel_rate=0.0, indel_ins=0.4, indel_geom_p=0.7,
                 truncated_ratio=0.0, truncated_max_length=0.5,
                 polyA_ratio=0.0, polyA_min=5, polyA_max=20, seed=1)
    build_pool.TEPoolBuilder(args)._run()
    return list(SeqIO.parse(os.path.join(d, "pool.fa"), "fasta"))


def test_library_annotation_is_carried_onto_pool_records():
    """The tag has to reach the pool, because the pool is all Simulate ever sees."""
    with tempfile.TemporaryDirectory() as d:
        pool = _build(d, [(">copia#LTR/Copia TSD=5", ELEMENT)])
        assert pool
        for rec in pool:
            assert rec.id.startswith("copia#LTR/Copia"), rec.id
            assert parse_tsd_tag(description_tail(rec), rec.id) == (5, 5), rec.description


def test_untagged_library_writes_bare_headers():
    """A library that says nothing must produce the same pool it always did."""
    with tempfile.TemporaryDirectory() as d:
        pool_path = os.path.join(d, "pool.fa")
        _build(d, [(">copia#LTR/Copia", ELEMENT)])
        with open(pool_path) as f:
            headers = [ln.rstrip("\n") for ln in f if ln.startswith(">")]
        assert headers
        for header in headers:
            assert " " not in header, header


def test_full_length_record_inherits_the_ltr_tag():
    """The synthesised element takes the LTR's annotation, not a description of itself."""
    with tempfile.TemporaryDirectory() as d:
        pool = _build(d, [(">TY1-LTR#LTR/Copia TSD=5", ELEMENT),
                          (">TY1-I#LTR/Copia", INTERNAL)])
        assert pool
        for rec in pool:
            assert rec.id.startswith("TY1-FULL#LTR/Copia"), rec.id
            assert parse_tsd_tag(description_tail(rec), rec.id) == (5, 5), rec.description
            assert "Merged LTR-I-LTR" not in rec.description, rec.description


def test_full_length_falls_back_to_the_internal_tag():
    """A library that tagged only the internal record still gets its value through."""
    with tempfile.TemporaryDirectory() as d:
        pool = _build(d, [(">TY1-LTR#LTR/Copia", ELEMENT),
                          (">TY1-I#LTR/Copia TSD=6", INTERNAL)])
        for rec in pool:
            assert parse_tsd_tag(description_tail(rec), rec.id) == (6, 6), rec.description


def test_unpaired_ltr_parts_keep_their_own_tags():
    """An -I with no matching -LTR is emitted as itself, annotation and all."""
    with tempfile.TemporaryDirectory() as d:
        pool = _build(d, [(">TY1-I#LTR/Copia TSD=7", INTERNAL)], nINS=4)
        for rec in pool:
            assert rec.id.startswith("TY1-I#LTR/Copia"), rec.id
            assert parse_tsd_tag(description_tail(rec), rec.id) == (7, 7), rec.description


def test_a_malformed_tag_fails_at_pool_build():
    """TErandom is the only command that reads the library, so it is where a typo must die."""
    with tempfile.TemporaryDirectory() as d:
        try:
            _build(d, [(">copia#LTR/Copia TSD=9-2", ELEMENT)])
        except ValueError as exc:
            assert "copia#LTR/Copia" in str(exc), exc
            return
        raise AssertionError("an inverted TSD tag should not have built a pool")


def test_a_trailing_space_is_not_a_tag():
    """Real libraries have trailing whitespace in headers; that is absence, not a typo."""
    with tempfile.TemporaryDirectory() as d:
        pool = _build(d, [(">TY1-LTR#LTR/Copia ", ELEMENT),
                          (">TY1-I#LTR/Copia ", INTERNAL)])
        for rec in pool:
            assert parse_tsd_tag(description_tail(rec), rec.id) is None, rec.description


if __name__ == "__main__":
    failures = 0
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            try:
                fn()
                print(f"PASS {name}")
            except AssertionError as exc:
                failures += 1
                print(f"FAIL {name}: {exc}")
    print("\nall passed" if not failures else f"\n{failures} failure(s)")
    sys.exit(1 if failures else 0)
