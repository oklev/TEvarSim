"""
Tests for LTR-LTR recombination (excision) events.

An excision converts a full-length LTR element (LTR-Internal-LTR) in the reference into a
solo LTR. In the simulated genome this is a substitution: genotype 0 keeps the full element,
genotype 1 leaves only the solo LTR.

No external test runner is required::

    PYTHONPATH=. python tests/test_excision.py

The test functions are also discoverable by pytest.
"""
import os
import sys
import tempfile

# Import the working-tree package regardless of any installed copy.
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from TEvarSim import simulate, TE_real  # noqa: E402


class _Args:
    def __init__(self, **kw):
        self.__dict__.update(kw)


# ---- shared fixture sequence -------------------------------------------------
PREFIX = "C" * 50
LEFT_LTR = "AGCTAGCTAGCTAGCTAGCT"   # 20 bp
INTERNAL = "G" * 30
RIGHT_LTR = LEFT_LTR                 # identical flanking LTR
SUFFIX = "T" * 50
REF_SEQ = PREFIX + LEFT_LTR + INTERNAL + RIGHT_LTR + SUFFIX
ELEM_START = len(PREFIX)                     # 50 (0-based BED start)
ELEM_END = ELEM_START + len(LEFT_LTR) + len(INTERNAL) + len(RIGHT_LTR)  # 120
LTR_LEN = len(LEFT_LTR)                       # 20


def _write_fixture(d):
    ref = os.path.join(d, "ref.fa")
    with open(ref, "w") as f:
        f.write(">chrT\n")
        for i in range(0, len(REF_SEQ), 60):
            f.write(REF_SEQ[i:i + 60] + "\n")
    pool = os.path.join(d, "pool.fa")
    with open(pool, "w") as f:
        f.write(">dummy#LTR/Gypsy\nACGT\n")
    bed = os.path.join(d, "exc.bed")
    with open(bed, "w") as f:
        f.write(
            f"chrT\t{ELEM_START}\t{ELEM_END}\t"
            f"EXC-chrT-{ELEM_START}-{ELEM_END}-{LTR_LEN}-LTR/Gypsy-TestElem\t.\t+\t{LTR_LEN}\n"
        )
    return ref, pool, bed


def test_excision_simulate_end_to_end():
    """Genotype 0 keeps the full element; genotype 1 keeps only the solo LTR."""
    with tempfile.TemporaryDirectory() as d:
        ref, pool, bed = _write_fixture(d)
        out = os.path.join(d, "Sim")
        args = _Args(ref=ref, pool=pool, bed=bed, outprefix=out, num=10,
                     af_min=0.4, af_max=0.6, tsd_min=5, tsd_max=20,
                     sense_strand_ratio=0.5, diverse=False, diverse_config=None, seed=1,
                     af_dist="uniform", af_mean=None, allow_zero_carriers=False)
        simulate.Simulator(args)._run()

        # ---- VCF ----
        exc_line = None
        with open(out + ".vcf") as f:
            for line in f:
                if line.startswith("chrT") and "TYPE=EXC" in line:
                    exc_line = line.rstrip("\n").split("\t")
        assert exc_line is not None, "no EXC record written to VCF"
        pos, ref_allele, alt_allele, info = exc_line[1], exc_line[3], exc_line[4], exc_line[7]
        assert pos == str(ELEM_START), f"POS should be element start, got {pos}"
        assert ref_allele == REF_SEQ[ELEM_START - 1:ELEM_END], "REF should be anchor + full element"
        assert alt_allele == REF_SEQ[ELEM_START - 1] + REF_SEQ[ELEM_START:ELEM_START + LTR_LEN], \
            "ALT should be anchor + solo LTR"
        assert f"LTRLEN={LTR_LEN}" in info
        assert f"SVLEN={len(alt_allele) - len(ref_allele)}" in info

        gts = exc_line[9:]

        # ---- genomes ----
        full = LEFT_LTR + INTERNAL + RIGHT_LTR
        solo_only_genome = PREFIX + LEFT_LTR + SUFFIX
        full_genome = REF_SEQ
        for i, g in enumerate(gts):
            seq = _read_single_fasta(f"{out}_{i}.fa")
            if g == "1":
                assert seq == solo_only_genome, f"genome {i} (GT=1) should be solo-LTR only"
                assert full not in seq
            else:
                assert seq == full_genome, f"genome {i} (GT=0) should keep the full element"


def test_vcf_declares_and_writes_eventtype_per_alt():
    """EVENTTYPE is Number=A: one event type per ALT allele, declared in the header."""
    with tempfile.TemporaryDirectory() as d:
        ref, pool, bed = _write_fixture(d)
        out = os.path.join(d, "Sim")
        args = _Args(ref=ref, pool=pool, bed=bed, outprefix=out, num=10,
                     af_min=0.4, af_max=0.6, tsd_min=5, tsd_max=20,
                     sense_strand_ratio=0.5, diverse=False, diverse_config=None, seed=1,
                     af_dist="uniform", af_mean=None, allow_zero_carriers=False)
        simulate.Simulator(args)._run()

        header, records = [], []
        with open(out + ".vcf") as f:
            for line in f:
                (header if line.startswith("##") else
                 records if not line.startswith("#CHROM") else []).append(line.rstrip("\n"))
        assert any(l.startswith('##INFO=<ID=EVENTTYPE,Number=A,Type=String,') for l in header), \
            "EVENTTYPE must be declared as Number=A"
        assert records, "no records written"
        for line in records:
            fields = line.split("\t")
            info = dict(f.partition("=")[::2] for f in fields[7].split(";"))
            assert "EVENTTYPE" in info, f"record without EVENTTYPE: {fields[:5]}"
            # a single simulated generation writes one event per record
            assert info["EVENTTYPE"].split(",") == [info["TYPE"]], info
            assert len(info["EVENTTYPE"].split(",")) == len(fields[4].split(",")), \
                f"EVENTTYPE must have one value per ALT: {info['EVENTTYPE']} vs {fields[4][:40]}"


def test_tereal_emits_excision_from_repeatmasker():
    """TEreal detects a full-length LTR element and emits an EXC BED line with the LTR length,
    and removes that locus from the deletion pool."""
    rm = os.path.join(os.path.dirname(__file__), os.pardir, "testData", "ltr_element.out")
    with tempfile.TemporaryDirectory() as d:
        ins = os.path.join(d, "ins.fa")
        with open(ins, "w") as f:
            f.write(">chr21-500000-len20-TestAlu#SINE/Alu-polyA5-strand+\n" + "A" * 20 + "\n")
        out = os.path.join(d, "TR")
        args = _Args(knownINS=ins, existingTEs=os.path.abspath(rm), TEtype=["Gypsy", "L1"],
                     CHR="chr21", DELlen=100, nMIN=0, nSV=0, outprefix=out,
                     nINS=0, nDEL=1, nEXC=1, seed=3)
        TE_real.RealTE(args)._run()

        exc, dels = [], []
        with open(out + ".bed") as f:
            for line in f:
                cols = line.rstrip("\n").split("\t")
                if cols[3].startswith("EXC"):
                    exc.append(cols)
                elif cols[3].startswith("DEL"):
                    dels.append(cols)
        assert len(exc) == 1, "expected exactly one excision"
        assert len(exc[0]) == 7, "EXC line must carry the LTR length as a 7th column"
        assert int(exc[0][6]) == 299, "5' LTR fragment length (100300-100001)"
        # the L1 is the only remaining deletion; the LTR element became an excision instead
        assert len(dels) == 1 and "L1" in dels[0][3]
        # no deletion overlaps the excised element
        es, ee = int(exc[0][1]), int(exc[0][2])
        for dcols in dels:
            ds, de = int(dcols[1]), int(dcols[2])
            assert not (ds < ee and de > es), "a deletion overlaps the excised element"


def _read_single_fasta(path):
    seq = []
    with open(path) as f:
        for line in f:
            if not line.startswith(">"):
                seq.append(line.strip())
    return "".join(seq)


if __name__ == "__main__":
    failures = 0
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            try:
                fn()
                print(f"PASS {name}")
            except AssertionError as e:
                failures += 1
                print(f"FAIL {name}: {e}")
    sys.exit(1 if failures else 0)
