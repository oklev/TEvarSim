"""
Tests for genotype comparison by allele content in ``tevarsim Compare``.

The truth and prediction VCFs are built independently and enumerate their ALT alleles in
different orders, so comparing raw GT indices reports a genotype error whenever one file
carries an extra allele at a site -- a solo LTR left behind by an excision, or a second
element stacked on the first. Compare therefore matches alleles by content (sequence, or
net length change within a tolerance) instead of by index.

No external test runner is required::

    PYTHONPATH=. python tests/test_compare_alleles.py

The test functions are also discoverable by pytest.
"""
import os
import sys
import tempfile

# Import the working-tree package regardless of any installed copy.
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from TEvarSim import compare_vcf  # noqa: E402


class _Args:
    def __init__(self, **kw):
        self.__dict__.update(kw)


# ---- shared fixture sequences ------------------------------------------------
ANCHOR = "A"
ELEMENT = "GATTACA" * 847          # 5929 bp, stands in for a full-length Ty1
ELEMENT_JITTERED = ELEMENT[:5923]  # same element seen 6 bp shorter across a TSD
SOLO_LTR = ELEMENT[:336]           # what an excision leaves behind
STACKED = ELEMENT + ELEMENT        # two elements at one site

VCF_HEADER = (
    "##fileformat=VCFv4.2\n"
    "##contig=<ID=chrT,length=1000000>\n"
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
)


def _write_vcf(path, records):
    """records: list of (pos, varID, ref, [alts], gt_string)"""
    with open(path, "w") as f:
        f.write(VCF_HEADER)
        for pos, varID, ref, alts, gt in records:
            f.write("\t".join([
                "chrT", str(pos), varID, ref, ",".join(alts), ".", "PASS", ".", "GT", gt,
            ]) + "\n")


def _compare(d, truth_records, pred_records, **kw):
    """Run one Compare and return (tp, fp, fn, genotype_accuracy)."""
    truth = os.path.join(d, "truth.vcf")
    pred = os.path.join(d, "pred.vcf")
    _write_vcf(truth, truth_records)
    _write_vcf(pred, pred_records)
    args = _Args(truth=truth, pred=pred, predType="VCF", TEtype="TY1-FULL", INSonly=True,
                 nHap=1, max_dist=100, gt_len_tol=50, truthID="S1", predID="S1",
                 outprefix=os.path.join(d, "out"))
    args.__dict__.update(kw)
    cmp = compare_vcf.CompareVCF(args)
    cmp._run()
    return cmp


TRUTH_INS = (1000, "TY1-FULL#LTR/Copia_1INDEL", ANCHOR, [ANCHOR + ELEMENT], "1")


def test_same_element_different_alt_index_matches():
    """The motivating case: identical element, but the prediction indexes it as ALT 2."""
    with tempfile.TemporaryDirectory() as d:
        cmp = _compare(
            d,
            [TRUTH_INS],
            # prediction carries a solo LTR as ALT1, pushing the element to ALT2
            [(1000, "1.1", ANCHOR, [ANCHOR + SOLO_LTR, ANCHOR + ELEMENT_JITTERED], "2")],
        )
        assert cmp.nMatch == 1, cmp.nMatch
        assert cmp.gt_accuracy == 1.0, cmp.gt_accuracy
    print("PASS test_same_element_different_alt_index_matches")


def test_different_element_same_alt_index_mismatches():
    """The converse: indices agree but the prediction calls a solo LTR, not the element."""
    with tempfile.TemporaryDirectory() as d:
        cmp = _compare(
            d,
            [TRUTH_INS],
            [(1000, "1.1", ANCHOR, [ANCHOR + SOLO_LTR], "1")],
        )
        assert cmp.nMatch == 1, cmp.nMatch
        assert cmp.gt_accuracy == 0.0, cmp.gt_accuracy
    print("PASS test_different_element_same_alt_index_mismatches")


def test_length_tolerance_absorbs_tsd_jitter_but_not_a_stack():
    """A few bp of breakpoint jitter is the same allele; a doubled element is not."""
    with tempfile.TemporaryDirectory() as d:
        cmp = _compare(d, [TRUTH_INS], [(1000, "1.1", ANCHOR, [ANCHOR + ELEMENT_JITTERED], "1")])
        assert cmp.gt_accuracy == 1.0, cmp.gt_accuracy
    with tempfile.TemporaryDirectory() as d:
        cmp = _compare(d, [TRUTH_INS], [(1000, "1.1", ANCHOR, [ANCHOR + STACKED], "1")])
        assert cmp.gt_accuracy == 0.0, cmp.gt_accuracy
    print("PASS test_length_tolerance_absorbs_tsd_jitter_but_not_a_stack")


def test_symbolic_allele_falls_back_to_index():
    """<*> carries no sequence, so only an index comparison is possible."""
    ref = compare_vcf.resolve_alleles(ANCHOR, (ANCHOR + ELEMENT,), (1,))
    sym = compare_vcf.resolve_alleles(ANCHOR, ("<*>", ANCHOR + ELEMENT), (1,))
    assert sym[0].symbolic and sym[0].size is None
    # truth allele 1 vs a symbolic prediction allele 1: indices agree, so treated as a match
    assert compare_vcf.genotype_agrees(ref, sym, 50) is True
    # against symbolic allele index 2 the indices disagree
    sym2 = compare_vcf.resolve_alleles(ANCHOR, (ANCHOR + ELEMENT, "<*>"), (2,))
    assert compare_vcf.genotype_agrees(ref, sym2, 50) is False
    print("PASS test_symbolic_allele_falls_back_to_index")


def test_bed_prediction_falls_back_to_index_comparison():
    """A BED prediction has no allele sequences; genotype_agrees defers to the caller."""
    ref = compare_vcf.resolve_alleles(ANCHOR, (ANCHOR + ELEMENT,), (1,))
    assert compare_vcf.genotype_agrees(ref, None, 50) is None
    assert compare_vcf.genotype_agrees(None, ref, 50) is None
    print("PASS test_bed_prediction_falls_back_to_index_comparison")


def test_diploid_allele_multiset_is_order_insensitive():
    """1/2 and 2/1 describing the same two alleles agree regardless of ALT ordering."""
    a = compare_vcf.resolve_alleles(ANCHOR, (ANCHOR + ELEMENT, ANCHOR + SOLO_LTR), (1, 2))
    b = compare_vcf.resolve_alleles(ANCHOR, (ANCHOR + SOLO_LTR, ANCHOR + ELEMENT), (2, 1))
    assert compare_vcf.genotype_agrees(a, b, 50) is True
    # same indices, but the second file's ALT2 is a stacked pair rather than a solo LTR
    c = compare_vcf.resolve_alleles(ANCHOR, (ANCHOR + SOLO_LTR, ANCHOR + STACKED), (2, 1))
    assert compare_vcf.genotype_agrees(a, c, 50) is False
    print("PASS test_diploid_allele_multiset_is_order_insensitive")


def test_parse_te_ids_returns_family_not_superfamily():
    """--TEtype filters on the family, which is what separates TY1 from TY2/TY4/TY5."""
    assert compare_vcf.parse_te_ids("TY1-FULL#LTR/Copia_15INDEL") == ("TY1-FULL", "Copia")
    assert compare_vcf.parse_te_ids("TY1-FULL#LTR/Copia") == ("TY1-FULL", "Copia")
    assert compare_vcf.parse_te_ids(
        "EXC-chrXIV_RagTag_2-611713-617630-335-LTR/Copia-TY1-FULL") == ("TY1-FULL", "Copia")
    assert compare_vcf.parse_te_ids(
        "EXC-chrVII_RagTag_2_0-412357-417740-354-LTR/Gypsy-TY3_1p-FULL") == ("TY3_1p-FULL", "Gypsy")
    assert compare_vcf.parse_te_ids("DEL-chrI-100-200-LTR/Copia-TY2-LTR") == ("TY2-LTR", "Copia")
    # a DEL ID built without a trailing family name degrades to the superfamily
    assert compare_vcf.parse_te_ids("DEL-chrI-100-200-LTR/Copia") == ("Copia", "Copia")
    print("PASS test_parse_te_ids_returns_family_not_superfamily")


def test_parse_te_ids_reads_the_family_out_of_a_callset_insertion_id():
    """
    The MEI callset names an insertion after the locus it was lifted from and annotates the
    polyA tail and strand: chr21-211282-len320-AluY#SINE/Alu-polyA23-strand+. Neither part
    says which family the element belongs to, so reading the ID literally made every
    insertion its own family and left "--TEtype AluY" matching nothing.
    """
    assert compare_vcf.parse_te_ids(
        "chr21-211282-len320-AluY#SINE/Alu-polyA23-strand+") == ("AluY", "Alu")
    assert compare_vcf.parse_te_ids(
        "chr1-683234-len2747-SVA_F#Retroposon/SVA-polyA49-strand-") == ("SVA_F", "SVA")
    assert compare_vcf.parse_te_ids(
        "chr1-4031386-len674-L1Ambig#LINE/L1-polyA52-strand-") == ("L1Ambig", "L1")
    # the README's undecorated CHR-POS-ID form
    assert compare_vcf.parse_te_ids("chr1-683234-AluSp#SINE/Alu") == ("AluSp", "Alu")
    print("PASS test_parse_te_ids_reads_the_family_out_of_a_callset_insertion_id")


def test_parse_te_ids_leaves_an_undecorated_family_name_alone():
    """A family name that merely contains dashes is not a CHR-POS-ID: its second
    dash-separated field is not a position."""
    assert compare_vcf.parse_te_ids("TY1-FULL#LTR/Copia") == ("TY1-FULL", "Copia")
    assert compare_vcf.parse_te_ids("TY3_1p-FULL#LTR/Gypsy") == ("TY3_1p-FULL", "Gypsy")
    assert compare_vcf.parse_te_ids("AluY#SINE/Alu") == ("AluY", "Alu")
    assert compare_vcf.parse_te_ids("L1HS#LINE/L1") == ("L1HS", "L1")
    print("PASS test_parse_te_ids_leaves_an_undecorated_family_name_alone")


def test_parse_te_ids_keeps_a_hyphenated_superfamily_whole():
    """hAT-Charlie and TcMar-Tigger are single superfamilies; only the callset's own
    trailing annotations are stripped, not everything after the first dash."""
    assert compare_vcf.parse_te_ids("Charlie1#DNA/hAT-Charlie") == ("Charlie1", "hAT-Charlie")
    assert compare_vcf.parse_te_ids(
        "chr1-100-len300-Tigger1#DNA/TcMar-Tigger-polyA20-strand+") == ("Tigger1", "TcMar-Tigger")
    print("PASS test_parse_te_ids_keeps_a_hyphenated_superfamily_whole")


def test_callset_ids_collapse_into_one_filterable_family():
    """Two AluY insertions lifted from different loci are one family, so --TEtype AluY
    benchmarks both rather than raising on an unknown family."""
    a = (1000, "chr21-211282-len320-AluY#SINE/Alu-polyA23-strand+", ANCHOR,
         [ANCHOR + ELEMENT], "1")
    b = (5000, "chr21-980614-len334-AluY#SINE/Alu-polyA37-strand-", ANCHOR,
         [ANCHOR + ELEMENT], "1")
    with tempfile.TemporaryDirectory() as d:
        cmp = _compare(d, [a, b],
                       [(1000, "1.1", ANCHOR, [ANCHOR + ELEMENT_JITTERED], "1"),
                        (5000, "1.2", ANCHOR, [ANCHOR + ELEMENT_JITTERED], "1")],
                       TEtype="AluY")
        assert (cmp.nMatch, cmp.recall, cmp.precision) == (2, 1.0, 1.0), \
            (cmp.nMatch, cmp.recall, cmp.precision)
    print("PASS test_callset_ids_collapse_into_one_filterable_family")


def test_family_filter_separates_same_superfamily_families():
    """TY1 and TY2 share the Copia superfamily; only the TY1 record must be compared."""
    ty2 = (2000, "TY2-FULL#LTR/Copia_3INDEL", ANCHOR, [ANCHOR + ELEMENT], "1")
    with tempfile.TemporaryDirectory() as d:
        cmp = _compare(d, [TRUTH_INS, ty2],
                       [(1000, "1.1", ANCHOR, [ANCHOR + ELEMENT_JITTERED], "1")],
                       TEtype="TY1-FULL")
        # only the TY1 truth record is loaded, so the single prediction is a clean TP
        assert (cmp.nMatch, cmp.recall, cmp.precision) == (1, 1.0, 1.0), \
            (cmp.nMatch, cmp.recall, cmp.precision)
    with tempfile.TemporaryDirectory() as d:
        cmp = _compare(d, [TRUTH_INS, ty2],
                       [(1000, "1.1", ANCHOR, [ANCHOR + ELEMENT_JITTERED], "1")],
                       TEtype="TY2-FULL")
        # filtering to TY2 leaves the TY1 prediction unmatched: 1 FP and 1 FN
        assert (cmp.nMatch, cmp.recall, cmp.precision) == (0, 0, 0), \
            (cmp.nMatch, cmp.recall, cmp.precision)
    print("PASS test_family_filter_separates_same_superfamily_families")


def test_family_filter_is_case_insensitive():
    with tempfile.TemporaryDirectory() as d:
        cmp = _compare(d, [TRUTH_INS],
                       [(1000, "1.1", ANCHOR, [ANCHOR + ELEMENT_JITTERED], "1")],
                       TEtype="ty1-full")
        assert cmp.nMatch == 1, cmp.nMatch
    print("PASS test_family_filter_is_case_insensitive")


def test_unknown_family_raises_with_the_available_families():
    """The old superfamily value must fail loudly, not silently compare against nothing."""
    with tempfile.TemporaryDirectory() as d:
        try:
            _compare(d, [TRUTH_INS], [(1000, "1.1", ANCHOR, [ANCHOR + ELEMENT], "1")],
                     TEtype="Copia")
        except ValueError as e:
            assert "TY1-FULL" in str(e), str(e)
            assert "Copia" in str(e), str(e)
        else:
            raise AssertionError("expected ValueError for an unknown --TEtype")
    print("PASS test_unknown_family_raises_with_the_available_families")


def test_haploid_comparison_is_reported_as_haplotype_accuracy():
    """With one allele per sample there is no genotype beyond which haplotype is carried."""
    import contextlib
    import io
    with tempfile.TemporaryDirectory() as d:
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            cmp = _compare(d, [TRUTH_INS],
                           [(1000, "1.1", ANCHOR, [ANCHOR + ELEMENT_JITTERED], "1")])
        assert cmp.accuracy_label == "haplotype", cmp.accuracy_label
        assert "haplotype accuracy" in buf.getvalue(), buf.getvalue()
        assert "genotype accuracy" not in buf.getvalue()
    print("PASS test_haploid_comparison_is_reported_as_haplotype_accuracy")


def test_diploid_comparison_is_still_reported_as_genotype_accuracy():
    import contextlib
    import io
    diploid_truth = (1000, "TY1-FULL#LTR/Copia_1INDEL", ANCHOR, [ANCHOR + ELEMENT], "1/1")
    with tempfile.TemporaryDirectory() as d:
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            cmp = _compare(d, [diploid_truth],
                           [(1000, "1.1", ANCHOR, [ANCHOR + ELEMENT_JITTERED], "1/1")])
        assert cmp.accuracy_label == "genotype", cmp.accuracy_label
        assert "genotype accuracy" in buf.getvalue(), buf.getvalue()
    print("PASS test_diploid_comparison_is_still_reported_as_genotype_accuracy")


def test_label_follows_the_data_not_the_nHap_flag():
    """A diploid file compared with --nHap 1 is still a genotype comparison."""
    diploid = (1000, "TY1-FULL#LTR/Copia_1INDEL", ANCHOR, [ANCHOR + ELEMENT], "0/1")
    with tempfile.TemporaryDirectory() as d:
        cmp = _compare(d, [diploid],
                       [(1000, "1.1", ANCHOR, [ANCHOR + ELEMENT_JITTERED], "0/1")], nHap=1)
        assert cmp.accuracy_label == "genotype", cmp.accuracy_label
    print("PASS test_label_follows_the_data_not_the_nHap_flag")


def test_match_file_is_parseable_and_flags_mismatches():
    """The match file must be real CSV with unix line endings so awk/cut work on it."""
    import csv as _csv
    with tempfile.TemporaryDirectory() as d:
        _compare(d, [TRUTH_INS], [(1000, "1.1", ANCHOR, [ANCHOR + SOLO_LTR], "1")])
        raw = open(os.path.join(d, "out.csv"), "rb").read()
        assert b"\r" not in raw, "match file must not use \\r\\n line endings"
        rows = list(_csv.DictReader(open(os.path.join(d, "out.csv"))))
        assert len(rows) == 1, rows
        assert rows[0]["gt_match"] == "MISMATCH", rows[0]
        assert rows[0]["truth_allele_bp"] == str(len(ELEMENT)), rows[0]
        assert rows[0]["pred_allele_bp"] == str(len(SOLO_LTR)), rows[0]
    print("PASS test_match_file_is_parseable_and_flags_mismatches")


if __name__ == "__main__":
    test_same_element_different_alt_index_matches()
    test_different_element_same_alt_index_mismatches()
    test_length_tolerance_absorbs_tsd_jitter_but_not_a_stack()
    test_symbolic_allele_falls_back_to_index()
    test_bed_prediction_falls_back_to_index_comparison()
    test_diploid_allele_multiset_is_order_insensitive()
    test_parse_te_ids_returns_family_not_superfamily()
    test_parse_te_ids_reads_the_family_out_of_a_callset_insertion_id()
    test_parse_te_ids_leaves_an_undecorated_family_name_alone()
    test_parse_te_ids_keeps_a_hyphenated_superfamily_whole()
    test_callset_ids_collapse_into_one_filterable_family()
    test_family_filter_separates_same_superfamily_families()
    test_family_filter_is_case_insensitive()
    test_unknown_family_raises_with_the_available_families()
    test_haploid_comparison_is_reported_as_haplotype_accuracy()
    test_diploid_comparison_is_still_reported_as_genotype_accuracy()
    test_label_follows_the_data_not_the_nHap_flag()
    test_match_file_is_parseable_and_flags_mismatches()
