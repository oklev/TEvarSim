"""
Tests for the target site duplication Simulate writes into an insertion.

A TSD is the stretch of host sequence the integration machinery cut twice, so the copy
carried into the ALT allele must be the bases immediately 5' of the insertion point, and
INFO/TSD must be how many of them there are. Neither held.

The slice read seq[start-tsd_len+1 : start+1] -- one base too far right. It duplicated
tsd_len-1 of the left flank's bases and then re-copied seq[start], the first base of the
right flank, so every insertion carried a spurious 1 bp tandem repeat at its 3' junction.
The net inserted length was still tsd_len, which is why nothing caught it: no length check
could see it, and the VCF and the genome were built from the same wrong string, so they
agreed with each other. It was only wrong against the reference.

The same slice went silently empty for an insertion closer to the contig start than
tsd_len, because a negative slice start reads from the far end of the contig.

No external test runner is required::

    PYTHONPATH=. python tests/test_tsd.py

The test functions are also discoverable by pytest.
"""
import os
import sys
import tempfile

# Import the working-tree package regardless of any installed copy.
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from TEvarSim import simulate  # noqa: E402


class _Args:
    def __init__(self, **kw):
        self.__dict__.update(kw)


# A non-repetitive contig, so a duplicated run is locatable and unambiguous. "ACGT"*n would
# make every candidate TSD match every other one and the assertions below would pass on the
# broken slice too.
CHR_T = ("GATTACACCGTAAGCTTGGCAATCGTTAAGGCCTATGACTGCAAGTTCCGGATCAGTTAACGGCTTAAGCCT"
         "GTCAAGGTTCCAAGGTTACGCATGGCCTTAAGGCATTCCGGAATTCCGGTTAACCGGTTAACCGGAATTCC"
         "AGGCTTAACCGGATCGATCGATTCCGGAAGGCCTTAACCGGTTAACCGGAATTCCGGAATTCCGGAATTCC")
ELEMENT = "TTTTTGGGGGAAAAACCCCC"

# The insertion point used throughout. Chosen so that CHR_T[INS_POS-1] != CHR_T[INS_POS]:
# the old slice re-copied CHR_T[INS_POS] onto the end of the duplication, and that doubling
# is only detectable where the two colliding bases differ.
INS_POS = 99


def _fixture(d, bed_rows, pool_records=((">elem#LTR/Copia", ELEMENT),)):
    ref = os.path.join(d, "ref.fa")
    with open(ref, "w") as f:
        f.write(">chrT\n")
        for i in range(0, len(CHR_T), 60):
            f.write(CHR_T[i:i + 60] + "\n")
    pool = os.path.join(d, "pool.fa")
    with open(pool, "w") as f:
        for header, seq in pool_records:
            f.write(f"{header}\n{seq}\n")
    bed = os.path.join(d, "events.bed")
    with open(bed, "w") as f:
        for row in bed_rows:
            f.write("\t".join(str(field) for field in row) + "\n")
    return ref, pool, bed


def _args(d, bed_rows, pool_records=((">elem#LTR/Copia", ELEMENT),), **overrides):
    ref, pool, bed = _fixture(d, bed_rows, pool_records)
    kw = dict(ref=ref, pool=pool, bed=bed, outprefix=os.path.join(d, "Sim"),
              num=4, af_min=1.0, af_max=1.0, tsd_min=6, tsd_max=6,
              sense_strand_ratio=0.5, diverse=False, diverse_config=None, seed=1,
                 af_dist="uniform", af_mean=None, allow_zero_carriers=False)
    kw.update(overrides)
    return _Args(**kw)


def _run(d, bed_rows, pool_records=((">elem#LTR/Copia", ELEMENT),), **overrides):
    """Run a whole simulation and hand back the VCF records and genome 0's sequence."""
    args = _args(d, bed_rows, pool_records, **overrides)
    simulate.Simulator(args)._run()
    records = []
    with open(args.outprefix + ".vcf") as f:
        for line in f:
            if not line.startswith("#"):
                records.append(line.rstrip("\n").split("\t"))
    with open(f"{args.outprefix}_0.fa") as f:
        genome = "".join(line.strip() for line in f if not line.startswith(">"))
    return records, genome


def _info(record):
    return dict(field.partition("=")[::2] for field in record[7].split(";"))


def _ins(pos, te_id="elem#LTR/Copia"):
    """A point insertion, as TErandom writes one: start == end."""
    return ("chrT", pos, pos, te_id, ".", "+")


def test_realised_duplication_equals_the_requested_length():
    """The ALT allele ends with the tsd_len bases immediately 5' of the insertion point."""
    pos = INS_POS
    with tempfile.TemporaryDirectory() as d:
        records, genome = _run(d, [_ins(pos)], tsd_min=6, tsd_max=6)
        assert len(records) == 1, records
        alt, info = records[0][4], _info(records[0])
        assert info["TSD"] == "6", info
        assert alt == CHR_T[pos - 1] + ELEMENT + CHR_T[pos - 6:pos], alt
        assert genome == CHR_T[:pos] + ELEMENT + CHR_T[pos - 6:pos] + CHR_T[pos:], genome


def test_no_spurious_tandem_duplication_at_the_junction():
    """The regression guard: the base at the insertion point must not appear twice.

    The old slice re-copied seq[start] on the end of the duplication, so the genome read
    ...<TE><tsd-1 bases><seq[start]><seq[start]>... The fixture position is chosen so that
    the two bases that would collide differ, making the doubling detectable.
    """
    pos = INS_POS
    assert CHR_T[pos - 1] != CHR_T[pos], "fixture must not already repeat at the junction"
    with tempfile.TemporaryDirectory() as d:
        _, genome = _run(d, [_ins(pos)], tsd_min=6, tsd_max=6)
        junction = pos + len(ELEMENT) + 6
        assert genome[junction] != genome[junction - 1], genome[junction - 8:junction + 4]
        assert genome == CHR_T[:pos] + ELEMENT + CHR_T[pos - 6:pos] + CHR_T[pos:]


def test_zero_tsd_inserts_only_the_element():
    """tsd_len 0 means no duplication -- the case the bg-prefix override already relied on."""
    pos = INS_POS
    with tempfile.TemporaryDirectory() as d:
        records, genome = _run(d, [_ins(pos)], tsd_min=0, tsd_max=0)
        assert records[0][4] == CHR_T[pos - 1] + ELEMENT, records[0][4]
        assert genome == CHR_T[:pos] + ELEMENT + CHR_T[pos:]


def test_tsd_near_the_contig_start_is_clamped():
    """An insertion closer to the start than tsd_len duplicates what there is, not nothing.

    A negative slice start reads from the far end of the contig, so the unclamped form
    produced an empty duplication while INFO/TSD still claimed tsd_len bases.
    """
    pos = 3
    with tempfile.TemporaryDirectory() as d:
        records, genome = _run(d, [_ins(pos)], tsd_min=10, tsd_max=10)
        alt = records[0][4]
        assert alt == CHR_T[pos - 1] + ELEMENT + CHR_T[0:pos], alt
        assert genome == CHR_T[:pos] + ELEMENT + CHR_T[0:pos] + CHR_T[pos:]
        # Specifically not the tail of the contig, which is what the negative slice reached for.
        assert CHR_T[-10:] not in alt


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
