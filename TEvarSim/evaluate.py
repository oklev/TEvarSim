'''
``tevarsim Evaluate`` -- score a prediction against every simulated event at once.

``Compare`` answers "how well was *this* genome called?". It takes one truth sample and one
prediction sample and reports site-level precision/recall plus a genotype accuracy for that
pair, so benchmarking a whole simulation means looping it over every genome and stitching
the per-genome outputs back together by hand.

``Evaluate`` turns the question around. The unit of analysis is the simulated *event*, not
the genome: each record in the truth VCF is one event, and each event is scored once,
against all of the simulated genomes at the same time. Truth and prediction samples are
paired by name (the pairing ``Compare`` is normally driven with, ``-I $sample -J $sample``),
so no genome has to be named on the command line. Each event then carries

  * whether the prediction recovered the locus at all, and how far off its breakpoint and
    allele length were,
  * which genomes were called as carriers versus which genomes actually carry it, and
  * whether the genotype of each paired genome agrees by allele content,

and those per-event records are aggregated into overall metrics plus breakdowns by event
type, TE family, event size and allele frequency.

Allele comparison, family parsing and the confusion-count arithmetic are shared with
``Compare`` so the two subcommands agree on what "the same allele" means.
'''
import json
import sys
from collections import OrderedDict

import pysam

from .compare_vcf import (
    calculate_metrics,
    convert_to_ploidy,
    genotype_agrees,
    parse_te_ids,
    resolve_alleles,
)

# Size strata, in bp of net length change. The defaults are chosen around the shapes an LTR
# simulation produces: a solo LTR (~340 bp), a full-length element (~6 kb) and a stacked
# pair (~12 kb) each land in their own bin.
DEFAULT_SIZE_BINS = (100, 500, 1000, 5000, 10000)
# Allele-frequency strata. Rare events are the hard ones, so the low end is finer.
DEFAULT_AF_BINS = (0.1, 0.25, 0.5, 0.75)
# Above this many simulated genomes an exact-carrier-count table is more noise than signal.
MAX_GENOMES_FOR_CARRIER_TABLE = 20


class Event:
    '''One simulated event: a single record of the truth VCF, with all of its genomes.'''

    def __init__(self, chrom, pos, varID, ref, alts, info, family, superfamily,
                 gts, descs, alt_descs):
        self.chrom = chrom
        self.pos = pos
        self.id = varID
        self.ref = ref
        self.alts = alts
        self.info = info
        self.family = family
        self.superfamily = superfamily
        self.gts = gts              # truth sample -> GT tuple
        self.descs = descs          # truth sample -> resolved Allele tuple (or None)
        self.alt_descs = alt_descs  # every non-reference allele of the record
        self.match = None           # the prediction record paired with this event


class PredRecord:
    '''One record of the prediction, kept with every sample it genotypes.'''

    def __init__(self, chrom, pos, varID, ref, alts, gts, descs, alt_descs, nonvariant):
        self.chrom = chrom
        self.pos = pos
        self.id = varID
        self.ref = ref
        self.alts = alts
        self.gts = gts
        self.descs = descs
        self.alt_descs = alt_descs
        self.nonvariant = nonvariant  # allele indices that mean "no variant here"
        self.matched = False


# ---- loading -----------------------------------------------------------------


def alt_descriptors(ref, alts, skip=()):
    '''
    Resolve every ALT of a record into an Allele, independent of any sample's genotype.

    ``skip`` names ALT strings that do not describe a variant, such as ``<*>``, so they
    neither pad the allele list nor offer themselves as a match.
    '''
    if not alts:
        return ()
    descs = resolve_alleles(ref, alts, tuple(range(1, len(alts) + 1)))
    if descs is None:
        return ()
    return tuple(d for d, a in zip(descs, alts) if a not in skip)


def event_type(record, info):
    '''
    The event class of a truth record: ``INS``, ``DEL``, ``EXC`` or ``OTHER``.

    ``INFO/TYPE`` is authoritative when present (every VCF ``Simulate`` writes carries it);
    a hand-built or third-party truth file without it falls back to the REF/ALT shapes.
    '''
    declared = info.get("TYPE")
    if declared:
        return str(declared)
    alt = record.alts[0] if record.alts else ""
    if len(record.ref) == 1 and len(alt) > 1:
        return "INS"
    if len(record.ref) > 1 and len(alt) == 1:
        return "DEL"
    return "OTHER"


def event_history(record, info):
    '''
    What happened to this record's element over the whole simulation, e.g. ``INS`` for a
    plain insertion, ``INS,EXC`` for one that was later excised to a solo LTR, and
    ``nested INS`` for the inner element of a stacked pair.

    ``INFO/EVENTTYPE`` carries one value per ALT allele and ``INFO/MEPRESENT`` says which
    of those alleles actually hold this record's element, so only the alleles that hold it
    say anything about its history. A truth VCF straight out of ``Simulate`` writes one
    event per record, so the history is just the type; a backtracked lineage VCF merges an
    element's later history onto the same record and the history is where that shows up.
    '''
    n_alts = len(record.alts or ())
    raw_types = info.get("EVENTTYPE")
    if not raw_types:
        return event_type(record, info)
    types = [str(t) for t in _as_tuple(raw_types)]
    types += [types[-1] if types else "?"] * (n_alts - len(types))
    raw_present = info.get("MEPRESENT")
    if raw_present:
        present = [int(v) == 1 for v in _as_tuple(raw_present)]
        present += [True] * (n_alts - len(present))
    else:
        present = [True] * n_alts
    history = [t for t, p in zip(types, present) if p]
    collapsed = [t for i, t in enumerate(history) if i == 0 or t != history[i - 1]]
    label = ",".join(collapsed) or "-"
    # The element of a stacked pair that is absent from the locus's first allele is the one
    # that landed inside the other.
    if present and not present[0]:
        label = f"nested {label}"
    return label


def _as_tuple(value):
    return value if isinstance(value, (tuple, list)) else (value,)


def event_size(alt_descs):
    '''
    The size of an event, as the net length change of its largest allele, keeping the sign
    so that an excision reads negative. Returns None when no allele carries a sequence.
    '''
    sizes = [d.size for d in alt_descs if not d.symbolic and d.size is not None]
    if not sizes:
        return None
    return max(sizes, key=abs)


def load_truth_events(vcf_file, INSonly, TEtype):
    '''
    Read every record of the truth VCF as one event, keeping all of its genomes.

    Mirrors ``Compare``'s ``--TEtype``/``--INSonly`` filtering, including its refusal to run
    against a family that is not in the file: filtering to a family that does not exist
    yields an empty benchmark whose 0% detection rate says nothing about the prediction.
    '''
    vcf = pysam.VariantFile(vcf_file)
    samples = list(vcf.header.samples)
    if not samples:
        raise ValueError(f"Truth VCF {vcf_file} has no sample columns to evaluate")
    wanted = TEtype.casefold() if TEtype is not None else None
    seen_families = set()
    events = []
    for record in vcf:
        varID = record.id or f"{record.chrom}:{record.pos}"
        info = dict(record.info)
        etype = event_type(record, info)
        if INSonly and etype != "INS":
            continue
        family, superfamily = parse_te_ids(varID)
        seen_families.add(family)
        if wanted is not None and family.casefold() != wanted:
            continue
        gts, descs = {}, {}
        for sample in samples:
            gt = record.samples[sample].get("GT")
            gt = tuple(gt) if gt is not None else (None,)
            gts[sample] = gt
            descs[sample] = resolve_alleles(record.ref, record.alts, gt)
        alt_descs = alt_descriptors(record.ref, record.alts)
        event = Event(record.chrom, record.pos, varID, record.ref, tuple(record.alts or ()),
                      info, family, superfamily, gts, descs, alt_descs)
        event.type = etype
        event.history = event_history(record, info)
        event.size = event_size(alt_descs)
        events.append(event)
    if wanted is not None and not any(f.casefold() == wanted for f in seen_families):
        raise ValueError(
            f"--TEtype '{TEtype}' matches no TE family in {vcf_file}. "
            f"Families present: {', '.join(sorted(seen_families))}. "
            "Note --TEtype matches the family (e.g. TY1-FULL), not the superfamily (Copia).")
    return events, samples


def load_pred_vcf(vcf_file):
    '''
    Read the prediction VCF, keeping every record that any sample calls as a variant.

    A record all of whose samples are reference (or ``<*>``, the "no variant in this
    genome" allele) describes nothing that could match a simulated event, and counting it
    among the predictions would deflate precision for free.
    '''
    vcf = pysam.VariantFile(vcf_file)
    samples = list(vcf.header.samples)
    records = []
    for record in vcf:
        nonvariant = {0}
        if record.alts and "<*>" in record.alts:
            nonvariant.add(record.alts.index("<*>") + 1)
        gts, descs = {}, {}
        called = False
        for sample in samples:
            gt = record.samples[sample].get("GT")
            gt = tuple(gt) if gt is not None else (None,)
            gts[sample] = gt
            descs[sample] = resolve_alleles(record.ref, record.alts, gt)
            if any(a is not None and a not in nonvariant for a in gt):
                called = True
        # A sample-less VCF (a sites-only call set) still describes predictions; keep it.
        if samples and not called:
            continue
        alt_descs = alt_descriptors(record.ref, record.alts, skip=("<*>",))
        records.append(PredRecord(record.chrom, record.pos, record.id or ".", record.ref,
                                  tuple(record.alts or ()), gts, descs, alt_descs,
                                  nonvariant))
    return records, samples


def load_pred_bed(bed_file):
    '''
    Read a BED prediction. A BED carries no allele sequences and no sample columns, so it
    can only be scored at the locus level: detection and breakpoint offset, no genotypes.
    '''
    records = []
    with open(bed_file) as fin:
        for line in fin:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.split()
            chrom, start = fields[0], int(fields[1])
            varID = fields[3] if len(fields) > 3 else "."
            records.append(PredRecord(chrom, start, varID, None, (), {}, {}, (), {0}))
    return records, []


# ---- sample pairing ----------------------------------------------------------


def read_sample_map(path):
    '''Read a two-column ``truth_sample<TAB>pred_sample`` mapping file.'''
    pairs = []
    with open(path) as fin:
        for lineno, line in enumerate(fin, 1):
            if line.startswith("#") or not line.strip():
                continue
            fields = line.split()
            if len(fields) < 2:
                raise ValueError(
                    f"{path}:{lineno}: expected 'truth_sample<TAB>pred_sample', got {line!r}")
            pairs.append((fields[0], fields[1]))
    return pairs


def pair_samples(truth_samples, pred_samples, sample_map=None):
    '''
    Decide which prediction sample stands for which simulated genome.

    An explicit ``--sample_map`` wins. Otherwise the samples are paired by name, which is
    the pairing a per-genome ``Compare`` loop is normally driven with; as a convenience a
    lone sample on each side is paired whatever it is called. Nothing is paired by column
    order: a silently wrong pairing would report confident, meaningless genotype accuracy.
    '''
    truth_set, pred_set = set(truth_samples), set(pred_samples)
    if sample_map:
        pairs = []
        for tname, pname in sample_map:
            if tname not in truth_set:
                raise ValueError(f"--sample_map names '{tname}', which is not a truth sample")
            if pname not in pred_set:
                raise ValueError(f"--sample_map names '{pname}', which is not a prediction sample")
            pairs.append((tname, pname))
        return pairs, "sample map"
    shared = [s for s in truth_samples if s in pred_set]
    if shared:
        return [(s, s) for s in shared], "matching names"
    if len(truth_samples) == 1 and len(pred_samples) == 1:
        return [(truth_samples[0], pred_samples[0])], "the only sample on each side"
    return [], "nothing"


# ---- matching ----------------------------------------------------------------


def alleles_overlap(truth_descs, pred_descs, tol):
    '''Does any non-reference allele of the event also appear in the prediction record?'''
    for t in truth_descs:
        for p in pred_descs:
            if t.symbolic or p.symbolic:
                if t.symbolic and p.symbolic and t.seq == p.seq:
                    return True
                continue
            if t.seq == p.seq or abs(t.size - p.size) <= tol:
                return True
    return False


def best_length_error(truth_descs, pred_descs):
    '''
    The smallest length disagreement between any pair of truth and predicted alleles, i.e.
    how far the prediction is from the simulated element at its closest. None when either
    side carries no comparable sequence.
    '''
    errors = [p.size - t.size
              for t in truth_descs if not t.symbolic and t.size is not None
              for p in pred_descs if not p.symbolic and p.size is not None]
    if not errors:
        return None
    return min(errors, key=abs)


def match_events(events, records, max_dist, tol):
    '''
    Pair each simulated event with at most one prediction record, and vice versa.

    Every candidate pair within ``max_dist`` is scored, and the pairs are then taken in
    order of preference: allele agreement first, breakpoint proximity second. The one-to-one
    constraint is what keeps a stacked pair of elements honest -- two events at one locus
    need two prediction records to both count as recovered, rather than both claiming the
    single record that happens to sit there.
    '''
    by_chrom = {}
    for record in records:
        by_chrom.setdefault(record.chrom, []).append(record)
    candidates = []
    for i, event in enumerate(events):
        for j, record in enumerate(by_chrom.get(event.chrom, [])):
            distance = abs(record.pos - event.pos)
            if distance > max_dist:
                continue
            agree = alleles_overlap(event.alt_descs, record.alt_descs, tol)
            candidates.append((not agree, distance, i, event.chrom, j))
    candidates.sort()
    taken_events, taken_records = set(), set()
    for _disagree, _distance, i, chrom, j in candidates:
        if i in taken_events or (chrom, j) in taken_records:
            continue
        taken_events.add(i)
        taken_records.add((chrom, j))
        events[i].match = by_chrom[chrom][j]
        by_chrom[chrom][j].matched = True


# ---- per-event scoring -------------------------------------------------------


def is_carrier(gt, nonvariant):
    '''Does this genotype hold the variant at all? Missing calls count as non-carriers.'''
    return any(a is not None and a not in nonvariant for a in gt)


def match_ploidy(t_gt, t_descs, p_gt, p_descs):
    '''
    Replicate the shorter genotype so a haploid truth can be compared against a diploid
    call (``1`` versus ``1/1``), the same widening ``Compare`` applies. Genotypes whose
    lengths are not multiples of one another are left alone and simply fail to agree.
    '''
    lt, lp = len(t_gt), len(p_gt)
    if lt == lp or lt == 0 or lp == 0:
        return t_gt, t_descs, p_gt, p_descs
    if lp > lt and lp % lt == 0:
        k = lp // lt
        return t_gt * k, (t_descs * k if t_descs else t_descs), p_gt, p_descs
    if lt > lp and lt % lp == 0:
        k = lt // lp
        return t_gt, t_descs, p_gt * k, (p_descs * k if p_descs else p_descs)
    return t_gt, t_descs, p_gt, p_descs


def score_event(event, pairs, tol):
    '''Turn one simulated event and its matched prediction into a scored record.'''
    match = event.match
    truth_carriers, pred_carriers = [], []
    tp = fp = fn = 0
    compared = concordant = 0
    for tsample, psample in pairs:
        t_gt = event.gts.get(tsample, (None,))
        t_carrier = is_carrier(t_gt, {0})
        if t_carrier:
            truth_carriers.append(tsample)
        if match is None:
            p_carrier = False
        else:
            p_gt = match.gts.get(psample, (None,))
            p_carrier = is_carrier(p_gt, match.nonvariant)
            if p_carrier:
                pred_carriers.append(psample)
            # A genotype is only comparable where a prediction record describes the locus.
            wt_gt, wt_descs, wp_gt, wp_descs = match_ploidy(
                t_gt, event.descs.get(tsample), p_gt, match.descs.get(psample))
            same = genotype_agrees(wt_descs, wp_descs, tol)
            if same is None:
                same = sorted(a for a in wt_gt if a is not None) == \
                       sorted(a for a in wp_gt if a is not None)
            compared += 1
            concordant += bool(same)
        tp += t_carrier and p_carrier
        fp += (not t_carrier) and p_carrier
        fn += t_carrier and (not p_carrier)

    # Allele frequency is a property of the simulation, so it is read off every genome the
    # truth VCF holds -- not only the ones that could be paired with a prediction.
    alleles = [a for gt in event.gts.values() for a in gt if a is not None]
    an = len(alleles)
    ac = sum(1 for a in alleles if a != 0)
    n_truth_carriers = sum(1 for gt in event.gts.values() if is_carrier(gt, {0}))

    scored = OrderedDict()
    scored["chrom"] = event.chrom
    scored["pos"] = event.pos
    scored["id"] = event.id
    scored["type"] = event.type
    scored["history"] = event.history
    scored["family"] = event.family
    scored["superfamily"] = event.superfamily
    scored["size_bp"] = event.size
    scored["allele_frequency"] = round(ac / an, 4) if an else None
    scored["allele_count"] = ac
    scored["allele_number"] = an
    scored["n_carrier_genomes"] = n_truth_carriers
    scored["detected"] = match is not None
    if match is None:
        scored["match"] = None
    else:
        length_error = best_length_error(event.alt_descs, match.alt_descs)
        scored["match"] = OrderedDict([
            ("id", match.id),
            ("pos", match.pos),
            ("pos_offset", match.pos - event.pos),
            ("allele_bp", event_size(match.alt_descs)),
            ("length_error", length_error),
            ("allele_match", alleles_overlap(event.alt_descs, match.alt_descs, tol)
                             if match.alt_descs else None),
        ])
    scored["carriers"] = OrderedDict([
        ("truth", truth_carriers),
        ("predicted", pred_carriers),
        ("tp", tp), ("fp", fp), ("fn", fn),
    ])
    scored["genotypes"] = OrderedDict([("compared", compared), ("concordant", concordant)])
    return scored


# ---- aggregation -------------------------------------------------------------


def _percentile(values, q):
    if not values:
        return None
    ordered = sorted(values)
    k = (len(ordered) - 1) * q
    low = int(k)
    high = min(low + 1, len(ordered) - 1)
    return round(ordered[low] + (ordered[high] - ordered[low]) * (k - low), 2)


def _spread(values):
    '''Median / mean / 90th percentile / max of a list of absolute errors.'''
    if not values:
        return OrderedDict([("n", 0), ("median", None), ("mean", None),
                            ("p90", None), ("max", None)])
    return OrderedDict([
        ("n", len(values)),
        ("median", _percentile(values, 0.5)),
        ("mean", round(sum(values) / len(values), 2)),
        ("p90", _percentile(values, 0.9)),
        ("max", max(values)),
    ])


def aggregate(scored_events):
    '''Roll a list of scored events up into one block of metrics.'''
    n = len(scored_events)
    detected = [e for e in scored_events if e["detected"]]
    allele_ok = sum(1 for e in detected if e["match"]["allele_match"])
    allele_comparable = sum(1 for e in detected if e["match"]["allele_match"] is not None)
    tp = sum(e["carriers"]["tp"] for e in scored_events)
    fp = sum(e["carriers"]["fp"] for e in scored_events)
    fn = sum(e["carriers"]["fn"] for e in scored_events)
    compared = sum(e["genotypes"]["compared"] for e in scored_events)
    concordant = sum(e["genotypes"]["concordant"] for e in scored_events)
    f1, precision, recall = calculate_metrics(tp, fp, fn)
    offsets = [abs(e["match"]["pos_offset"]) for e in detected]
    errors = [abs(e["match"]["length_error"]) for e in detected
              if e["match"]["length_error"] is not None]
    return OrderedDict([
        ("n_events", n),
        ("n_detected", len(detected)),
        ("detection_rate", round(len(detected) / n, 4) if n else None),
        ("n_allele_concordant", allele_ok),
        ("allele_concordance_rate",
         round(allele_ok / allele_comparable, 4) if allele_comparable else None),
        ("breakpoint_offset_bp", _spread(offsets)),
        ("allele_length_error_bp", _spread(errors)),
        ("carriers", OrderedDict([
            ("tp", tp), ("fp", fp), ("fn", fn),
            ("recall", round(recall, 4)), ("precision", round(precision, 4)),
            ("f1", round(f1, 4)),
        ])),
        ("genotypes", OrderedDict([
            ("compared", compared),
            ("concordant", concordant),
            ("concordance", round(concordant / compared, 4) if compared else None),
        ])),
    ])


def stratify(scored_events, key):
    '''Group events by ``key(event)`` and aggregate each group. Sorted by group size.'''
    groups = OrderedDict()
    for event in scored_events:
        groups.setdefault(key(event), []).append(event)
    ordered = sorted(groups.items(), key=lambda kv: (-len(kv[1]), str(kv[0])))
    return OrderedDict((str(label), aggregate(members)) for label, members in ordered)


def size_bin(size, edges):
    '''Label the size stratum an event falls in, by the magnitude of its length change.'''
    if size is None:
        return "unknown"
    magnitude = abs(size)
    labels = _size_bin_labels(edges)
    for edge, label in zip(edges, labels):
        if magnitude < edge:
            return label
    return labels[-1]


def _size_bin_labels(edges):
    labels = [f"<{_bp(edges[0])}"]
    for low, high in zip(edges, edges[1:]):
        labels.append(f"{_bp(low)}-{_bp(high)}")
    labels.append(f">={_bp(edges[-1])}")
    return labels


def _bp(value):
    if value >= 1000 and value % 1000 == 0:
        return f"{value // 1000}kb"
    return f"{value}bp"


def af_bin(af, edges):
    '''Label the allele-frequency stratum an event falls in.'''
    if af is None:
        return "unknown"
    previous = 0.0
    for edge in edges:
        if af <= edge:
            return f"{previous:g}-{edge:g}"
        previous = edge
    return f"{previous:g}-1"


# ---- reporting ---------------------------------------------------------------


def _pct(value):
    return "n/a" if value is None else f"{value * 100:.2f}%"


def _num(value):
    return "-" if value is None else f"{value:g}" if isinstance(value, float) else str(value)


def _table(headers, rows):
    '''Render a fixed-width table. Numbers right-aligned, the first column left-aligned.'''
    if not rows:
        return ["  (no events)"]
    widths = [max(len(str(h)), max(len(str(r[i])) for r in rows))
              for i, h in enumerate(headers)]
    def render(cells):
        out = [str(cells[0]).ljust(widths[0])]
        out += [str(c).rjust(w) for c, w in zip(cells[1:], widths[1:])]
        return "  " + "  ".join(out)
    lines = [render(headers), "  " + "  ".join("-" * w for w in widths)]
    lines += [render(r) for r in rows]
    return lines


def _stratum_rows(strata):
    rows = []
    for label, stats in strata.items():
        rows.append([
            label,
            stats["n_events"],
            stats["n_detected"],
            _pct(stats["detection_rate"]),
            _pct(stats["allele_concordance_rate"]),
            _num(stats["breakpoint_offset_bp"]["median"]),
            _num(stats["carriers"]["f1"]),
            _pct(stats["genotypes"]["concordance"]),
        ])
    return rows


STRATUM_HEADERS = ["stratum", "events", "found", "detection",
                   "allele ok", "med |off|", "carrier F1", "genotype"]


def format_report(summary, meta):
    '''The printed summary: the headline metrics, then one table per breakdown.'''
    lines = []
    lines.append("")
    lines.append("tevarsim Evaluate")
    lines.append(f"  truth : {meta['truth']}")
    lines.append(f"          {meta['n_events']} simulated events over "
                 f"{meta['n_truth_samples']} genomes")
    lines.append(f"  pred  : {meta['pred']}")
    lines.append(f"          {meta['n_pred_records']} records over "
                 f"{meta['n_pred_samples']} samples")
    lines.append(f"  paired: {len(meta['sample_pairs'])} genome(s) by {meta['pairing']}")
    if meta["sample_pairs"]:
        shown = ", ".join(f"{t}={p}" if t != p else t for t, p in meta["sample_pairs"][:6])
        if len(meta["sample_pairs"]) > 6:
            shown += f", ... (+{len(meta['sample_pairs']) - 6})"
        lines.append(f"          {shown}")
    else:
        lines.append("          no genomes could be paired -- carrier and genotype "
                     "statistics are unavailable.")
        lines.append("          Pass --sample_map truth_sample<TAB>pred_sample to pair them.")
    if meta["filters"]:
        lines.append(f"  filter: {meta['filters']}")
        lines.append("          the truth is filtered but the prediction is not, so "
                     "unmatched predictions include events the filter removed.")

    overall = summary["overall"]
    lines.append("")
    lines.append(f"Event detection ({overall['n_events']} simulated events)")
    lines.append(f"  recovered            : {overall['n_detected']} / {overall['n_events']}"
                 f"  ({_pct(overall['detection_rate'])})")
    lines.append(f"  allele concordant    : {overall['n_allele_concordant']} of the "
                 f"recovered  ({_pct(overall['allele_concordance_rate'])})")
    lines.append(f"  unmatched predictions: {summary['predictions']['unmatched']} of "
                 f"{summary['predictions']['total']}")
    lines.append(f"  event precision      : {_pct(summary['predictions']['precision'])}"
                 "   (prediction records matching a simulated event)")

    offsets = overall["breakpoint_offset_bp"]
    errors = overall["allele_length_error_bp"]
    lines.append("")
    lines.append("Accuracy of the recovered events")
    lines.append(f"  breakpoint |offset| : median {_num(offsets['median'])} bp, "
                 f"mean {_num(offsets['mean'])} bp, p90 {_num(offsets['p90'])} bp, "
                 f"max {_num(offsets['max'])} bp   (n={offsets['n']})")
    lines.append(f"  allele |length err| : median {_num(errors['median'])} bp, "
                 f"mean {_num(errors['mean'])} bp, p90 {_num(errors['p90'])} bp, "
                 f"max {_num(errors['max'])} bp   (n={errors['n']})")

    carriers = overall["carriers"]
    genotypes = overall["genotypes"]
    lines.append("")
    if meta["sample_pairs"]:
        lines.append(f"Carrier accuracy over {len(meta['sample_pairs'])} paired genomes "
                     f"({carriers['tp'] + carriers['fn']} simulated carrier genomes)")
        lines.append(f"  TP {carriers['tp']}, FP {carriers['fp']}, FN {carriers['fn']}")
        lines.append(f"  recall {carriers['recall']:.4f}, "
                     f"precision {carriers['precision']:.4f}, F1 {carriers['f1']:.4f}")
        lines.append(f"  {meta['genotype_label']} concordance: "
                     f"{genotypes['concordant']} / {genotypes['compared']} "
                     f"({_pct(genotypes['concordance'])})")
    else:
        lines.append("Carrier accuracy: not available (no genomes paired)")

    for title, key in (("By event type", "by_event_type"),
                       ("By event history", "by_event_history"),
                       ("By TE family", "by_family"),
                       ("By TE superfamily", "by_superfamily"),
                       ("By event size", "by_size"),
                       ("By allele frequency", "by_allele_frequency"),
                       ("By number of carrier genomes", "by_carrier_count")):
        strata = summary.get(key)
        if not strata:
            continue
        lines.append("")
        lines.append(title)
        lines.extend(_table(STRATUM_HEADERS, _stratum_rows(strata)))
    lines.append("")
    return "\n".join(lines)


# ---- driver ------------------------------------------------------------------


class Evaluator:
    def __init__(self, args):
        self.truth_file = args.truth
        self.pred_file = args.pred
        self.predType = getattr(args, "predType", "VCF")
        self.outprefix = args.outprefix
        self.TEtype = getattr(args, "TEtype", None)
        self.INSonly = getattr(args, "INSonly", False)
        self.max_dist = getattr(args, "max_dist", 100)
        self.gt_len_tol = getattr(args, "gt_len_tol", 50)
        self.nHap = getattr(args, "nHap", 1)
        self.sample_map_file = getattr(args, "sample_map", None)
        self.size_bins = tuple(getattr(args, "size_bins", None) or DEFAULT_SIZE_BINS)
        self.af_bins = tuple(getattr(args, "af_bins", None) or DEFAULT_AF_BINS)

    def _run(self):
        # Ground truth. --nHap merges consecutive haplotype columns into one individual,
        # exactly as Compare does, before anything is paired with the prediction.
        truth_file = self.truth_file
        if self.nHap > 1:
            truth_file = f"{self.outprefix}.polyhap.vcf"
            convert_to_ploidy(self.truth_file, self.nHap, truth_file)
        events, truth_samples = load_truth_events(truth_file, self.INSonly, self.TEtype)

        if self.predType == "VCF":
            records, pred_samples = load_pred_vcf(self.pred_file)
        else:
            records, pred_samples = load_pred_bed(self.pred_file)

        sample_map = read_sample_map(self.sample_map_file) if self.sample_map_file else None
        pairs, pairing = pair_samples(truth_samples, pred_samples, sample_map)

        match_events(events, records, self.max_dist, self.gt_len_tol)
        scored = [score_event(event, pairs, self.gt_len_tol)
                  for event in sorted(events, key=lambda e: (e.chrom, e.pos))]

        matched = sum(1 for r in records if r.matched)
        unmatched = len(records) - matched
        filters = ", ".join(
            part for part in (
                f"--TEtype {self.TEtype}" if self.TEtype else "",
                "--INSonly" if self.INSonly else "") if part)

        summary = OrderedDict()
        summary["overall"] = aggregate(scored)
        summary["predictions"] = OrderedDict([
            ("total", len(records)),
            ("matched", matched),
            ("unmatched", unmatched),
            ("precision", round(matched / len(records), 4) if records else None),
        ])
        summary["by_event_type"] = stratify(scored, lambda e: e["type"])
        if any(e["history"] != e["type"] for e in scored):
            summary["by_event_history"] = stratify(scored, lambda e: e["history"])
        summary["by_family"] = stratify(scored, lambda e: e["family"])
        # The family is parsed exactly as Compare's --TEtype parses it, so the labels here
        # are the values that flag accepts. For a pool whose IDs name the source locus
        # (the MEI callset's chr1-683234-AluSp#SINE/Alu) that makes almost every event its
        # own family, and the superfamily is the level that actually groups them -- so add
        # that table, but only where it groups anything the family table did not.
        superfamilies = {e["superfamily"] for e in scored if e["superfamily"]}
        if superfamilies and len(superfamilies) < len(summary["by_family"]):
            summary["by_superfamily"] = stratify(scored, lambda e: e["superfamily"] or "-")
        summary["by_size"] = _order_bins(
            stratify(scored, lambda e: size_bin(e["size_bp"], self.size_bins)),
            _size_bin_labels(self.size_bins))
        summary["by_allele_frequency"] = _order_bins(
            stratify(scored, lambda e: af_bin(e["allele_frequency"], self.af_bins)),
            _af_bin_labels(self.af_bins))
        if len(truth_samples) <= MAX_GENOMES_FOR_CARRIER_TABLE:
            summary["by_carrier_count"] = OrderedDict(sorted(
                stratify(scored, lambda e: e["n_carrier_genomes"]).items(),
                key=lambda kv: int(kv[0])))

        # With one allele per genome there is no genotype to get right beyond which
        # haplotype is carried, so the concordance that is reported says "haplotype".
        haploid = all(len(gt) == 1 for e in events for gt in e.gts.values()) and \
                  all(len(gt) == 1 for r in records for gt in r.gts.values())
        meta = OrderedDict([
            ("truth", self.truth_file),
            ("pred", self.pred_file),
            ("predType", self.predType),
            ("n_events", len(events)),
            ("n_truth_samples", len(truth_samples)),
            ("n_pred_records", len(records)),
            ("n_pred_samples", len(pred_samples)),
            ("sample_pairs", pairs),
            ("pairing", pairing),
            ("genotype_label", "haplotype" if haploid else "genotype"),
            ("filters", filters),
            ("max_dist", self.max_dist),
            ("gt_len_tol", self.gt_len_tol),
            ("nHap", self.nHap),
            ("size_bins", list(self.size_bins)),
            ("af_bins", list(self.af_bins)),
        ])

        report = format_report(summary, meta)
        print(report)

        out_path = f"{self.outprefix}.json"
        with open(out_path, "w") as fo:
            json.dump(OrderedDict([("meta", meta), ("summary", summary),
                                   ("events", scored)]), fo, indent=2)
            fo.write("\n")
        # Flush first: the report goes to a block-buffered stdout when it is piped, and the
        # unbuffered stderr note would otherwise land above a report it comes after.
        sys.stdout.flush()
        print(f"[INFO] per-event results written to {out_path}", file=sys.stderr)

        self.meta = meta
        self.summary = summary
        self.events = scored
        return self


def _order_bins(strata, labels):
    '''Put binned strata back in bin order; ``stratify`` sorts by size, which scrambles it.'''
    ordered = OrderedDict((label, strata[label]) for label in labels if label in strata)
    for label, stats in strata.items():
        if label not in ordered:
            ordered[label] = stats
    return ordered


def _af_bin_labels(edges):
    labels, previous = [], 0.0
    for edge in edges:
        labels.append(f"{previous:g}-{edge:g}")
        previous = edge
    labels.append(f"{previous:g}-1")
    return labels


def run(args):
    Evaluator(args)._run()
