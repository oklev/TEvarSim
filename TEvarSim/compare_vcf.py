import csv
from collections import namedtuple

import pysam

# One allele of a genotype, resolved against its own record's REF/ALT fields so that
# truth and prediction can be compared by allele *content* rather than by ALT index.
# The two VCFs are built independently and enumerate their ALTs in different orders --
# a site where the prediction has an extra allele (a solo LTR left by an excision, or a
# second element stacked on the first) shifts every later index by one, which made an
# index-equality check report a genotype error even when both files described the same
# haplotype.
#   idx      allele index as written in the GT field (0 = REF)
#   size     net length change of this allele vs REF, in bp (None if symbolic)
#   seq      the allele sequence, or the symbolic ID (e.g. "<*>")
#   symbolic True for symbolic alleles, whose content cannot be compared
Allele = namedtuple("Allele", "idx size seq symbolic")


def resolve_alleles(ref, alts, gt):
    '''
    Resolve a GT tuple into per-allele descriptors using the record's own REF/ALT.
    Returns None if any allele cannot be resolved, so the caller can fall back.
    '''
    descs = []
    for a in gt:
        if a is None:
            return None
        if a == 0:
            descs.append(Allele(0, 0, ref, False))
            continue
        if alts is None or a > len(alts):
            return None
        alt = alts[a - 1]
        symbolic = alt.startswith("<")
        descs.append(Allele(a, None if symbolic else len(alt) - len(ref), alt, symbolic))
    return tuple(descs)


def alleles_agree(t, p, tol):
    '''
    Do a truth allele and a predicted allele describe the same thing?
    Exact sequence identity first; otherwise net length change within `tol` bp, which
    absorbs the TSD and breakpoint jitter that makes the same element differ by a few
    bp between the two files while still separating a ~340bp solo LTR from a ~5.9kb
    full-length element from a ~11.9kb stacked pair.
    '''
    if t is None or p is None:
        return False
    # Symbolic alleles carry no comparable sequence; only an index check is possible.
    if t.symbolic or p.symbolic:
        if t.symbolic and p.symbolic:
            return t.seq == p.seq
        return t.idx == p.idx
    if t.seq == p.seq:
        return True
    return abs(t.size - p.size) <= tol


def all_haploid(*variant_sets):
    '''
    Do all the loaded variants carry a single allele?

    Read off the genotypes themselves rather than --nHap, so the answer follows the data:
    convert_to_ploidy merges haplotype columns before anything is loaded, and a file can
    disagree with the ploidy it was invoked with. With one allele per sample there is no
    genotype to get right beyond which haplotype is carried, so the accuracy that is
    reported is a haplotype accuracy and says so.
    '''
    for variants in variant_sets:
        for records in variants.values():
            for _pos, gt, _descs in records:
                if len(gt) != 1:
                    return False
    return True


def genotype_agrees(tdescs, pdescs, tol):
    '''
    Compare two resolved genotypes as unordered allele multisets. Falls back to raw
    index equality when either side has no resolvable alleles (e.g. a BED prediction).
    '''
    if tdescs is None or pdescs is None:
        return None
    if len(tdescs) != len(pdescs):
        return False
    # Pair the alleles by size so the comparison does not depend on ALT ordering,
    # generalising the old sorted(gt) == sorted(gt) check to allele content.
    key = lambda d: (d.size is None, d.size, d.idx)
    return all(alleles_agree(a, b, tol)
               for a, b in zip(sorted(tdescs, key=key), sorted(pdescs, key=key)))


def parse_te_ids(varID):
    '''
    Split a truth VCF variant ID into (family, superfamily).

    Filtering on the superfamily cannot separate the Ty families: TY1, TY2, TY4, TY5 and
    TSU4 are all LTR/Copia, so "--TEtype Copia" would pull all of them into one comparison.
    The family name is the discriminating level, so that is what --TEtype matches.

    Two ID shapes occur:
      INS      <family>#<class>/<superfamily>[_<modifications>]
               e.g. TY1-FULL#LTR/Copia_15INDEL          -> (TY1-FULL, Copia)
      DEL/EXC  <type>-<chrom>-<start>-<end>-...-<class>/<superfamily>-<family>
               e.g. EXC-chrXIV-611713-617630-335-LTR/Copia-TY1-FULL -> (TY1-FULL, Copia)
    Some DEL IDs are built without a trailing family name; those fall back to the
    superfamily so they stay filterable rather than silently dropping out.
    '''
    if varID.startswith("DEL") or varID.startswith("EXC"):
        tail = varID.rsplit("/", 1)[-1] if "/" in varID else varID.rsplit("-", 1)[-1]
        superfamily, _, family = tail.partition("-")
        return (family or superfamily), superfamily
    family = varID.split("#")[0]
    superfamily = varID.split("/")[1].split("_")[0] if "/" in varID else ""
    return family, superfamily


def _fmt_gt(gt):
    '''Render a GT tuple for the match file, e.g. (1, 2) -> "1/2".'''
    return "/".join("." if a is None else str(a) for a in gt)


def _fmt_sizes(descs):
    '''Render the net length change of each allele, e.g. "5929" or "<*>".'''
    if descs is None:
        return ""
    return "/".join(d.seq if d.symbolic else str(d.size) for d in descs)


def load_truth_vcf(vcf_file, sampleID, INSonly, TEtype):
    variants = {}
    vcf = pysam.VariantFile(vcf_file)

    if sampleID not in vcf.header.samples:
        raise ValueError(f"Sample ID '{sampleID}' not found in VCF header")
    
    wanted = TEtype.casefold() if TEtype is not None else None
    seen_families = set()
    for record in vcf:
        varID = record.id
        if INSonly and (varID.startswith("DEL") or varID.startswith("EXC")):
            continue
        TEfamily, _ = parse_te_ids(varID)
        seen_families.add(TEfamily)
        if wanted is not None and TEfamily.casefold() != wanted:
            continue
        sample = record.samples[sampleID]
        gt = sample.get('GT')
        # filter positions withou any variants
        if any(allele is None for allele in gt) or all(allele == 0 for allele in gt):
            continue
        descs = resolve_alleles(record.ref, record.alts, gt)
        variants.setdefault(record.chrom, []).append((record.pos, gt, descs))
    # A --TEtype naming no family in the truth VCF yields an empty benchmark and a
    # meaningless comparison (recall 0 against nothing), so say so rather than report it.
    if wanted is not None and not any(f.casefold() == wanted for f in seen_families):
        raise ValueError(
            f"--TEtype '{TEtype}' matches no TE family in {vcf_file}. "
            f"Families present: {', '.join(sorted(seen_families))}. "
            "Note --TEtype matches the family (e.g. TY1-FULL), not the superfamily (Copia)."
        )
    return variants

def load_vcf(vcf_file, sampleID):
    variants = {}
    vcf = pysam.VariantFile(vcf_file)

    if sampleID not in vcf.header.samples:
        raise ValueError(f"Sample ID '{sampleID}' not found in VCF header")
    
    for record in vcf.fetch():
        sample = record.samples[sampleID]
        non_var_gts = {0}
        if record.alts and "<*>" in record.alts:
            non_var_gts.add(record.alts.index("<*>")+1)
        gt = sample.get('GT')
        # filt positions withou any variants
        if any(allele is None for allele in gt) or all(allele in non_var_gts for allele in gt):
            continue
        descs = resolve_alleles(record.ref, record.alts, gt)
        variants.setdefault(record.chrom, []).append((record.pos, gt, descs))
    return variants

def load_bed(bed_file):
    variants = {}
    with open(bed_file, 'r') as fin:
        for line in fin:
            chrom, start, _, gt, *_ = line.strip().split()
            # tuple of genotypes
            if "/" in gt and "." not in gt:
                gt = tuple(map(int, gt.split("/")))
            else:
                gt = (None,)
            start = int(start)
            # A BED prediction carries no allele sequences, so no descriptors: the
            # genotype check falls back to comparing raw allele indices.
            variants.setdefault(chrom, []).append((start, gt, None))
    return variants

def calculate_metrics(tp, fp, fn):
    '''
    Calculate F1, precision, recall from confusion counts.
    '''
    # Parameter check
    assert tp >= 0 and fp >= 0 and fn >= 0, "tp, fp, fn must be non-negative integers"
    # Calculate Recall and Precision
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0
    # Calculate F1
    if (recall + precision) > 0:
        F1 = 2 * recall * precision / (recall + precision)
    else:
        F1 = 0
    return F1, precision, recall


class CompareVCF:
    def __init__(self, args):
        self.nHap = args.nHap
        self.truth_file = args.truth
        self.compare_file = args.pred
        self.predType = args.predType
        self.TEtype = args.TEtype
        self.INSonly = args.INSonly
        self.max_dist = args.max_dist
        # bp of slack when matching two alleles by length; getattr keeps older callers working
        self.gt_len_tol = getattr(args, "gt_len_tol", 50)
        self.truthID = args.truthID
        self.predID = args.predID
        self.MatchFile = args.outprefix
        # inner variables
        self.nMatch = 0
    
    def _run(self):
        # ground truth data
        if self.nHap > 1:
            self.convert_to_ploidy()
            bench_var = load_truth_vcf("polyhap.vcf", self.truthID, self.INSonly, self.TEtype)
        else:
            bench_var = load_truth_vcf(self.truth_file, self.truthID, self.INSonly, self.TEtype)
        # prediction data
        if self.predType == "VCF":
            pred_var = load_vcf(self.compare_file, self.predID)
        else:
            self.predID = "sample" if self.predID is None else self.predID
            pred_var = load_bed(self.compare_file)
        # metrics
        tp, fp, fn, gtDiff = 0, 0, 0, 0
        # different chromosomes
        bench_chr = set(bench_var.keys())
        pred_chr = set(pred_var.keys())
        fn_chr = bench_chr - pred_chr
        if len(fn_chr) > 0:
            fn += sum(len(bench_var[chrom]) for chrom in fn_chr)
        fp_chr = pred_chr - bench_chr
        if len(fp_chr) > 0:
            fp += sum(len(pred_var[chrom]) for chrom in fp_chr)
        # common chromosomes
        common_chr = pred_chr & bench_chr
        print(f"Common chromosomes: {common_chr}")
        if self.TEtype is None:
            tetype = "TEs"
        else:
            tetype = self.TEtype
        print(f"Nmber of {tetype} in truth VCF: {sum(len(bench_var[chrom]) for chrom in bench_chr)}")
        print(f"Number of {tetype} in prediction VCF: {sum(len(pred_var[chrom]) for chrom in pred_chr)}")
        # Truncate the match file once per run and give it a header. It is appended to
        # once per chromosome below, so it cannot be opened for writing in that loop.
        with open(self.MatchFile + ".csv", "w", newline="") as fo:
            csv.writer(fo, lineterminator="\n").writerow(
                ["chrom", "truth_pos", "truth_gt", "truth_allele_bp",
                 "pred_pos", "pred_gt", "pred_allele_bp", "gt_match"])
        for chrom in common_chr:
            bench_pos, bench_gt, bench_desc = zip(*bench_var[chrom])
            compare_pos, compare_gt, compare_desc = zip(*pred_var[chrom])
            tp_tmp, fp_tmp, fn_tmp, gtDiff_tmp = self.calculate_confusion_counts(
                bench_pos, bench_gt, bench_desc, compare_pos, compare_gt, compare_desc, chrom)
            tp += tp_tmp
            fp += fp_tmp
            fn += fn_tmp
            gtDiff += gtDiff_tmp

        self.F1, self.precision, self.recall = calculate_metrics(tp, fp, fn)
        print(f"For {tetype} occurrence sites:")
        print(f"True Positive: {tp}, False Positive: {fp}, False negative: {fn}")
        print(f"Recall: {self.recall:.4f}, Precision: {self.precision:.4f}, F1: {self.F1:.4f}")
        self.gt_accuracy = 1 - gtDiff/self.nMatch if self.nMatch else None
        self.accuracy_label = "haplotype" if all_haploid(bench_var, pred_var) else "genotype"
        if self.nMatch:
            print(f"The {self.accuracy_label} accuracy for matched {tetype}: {self.gt_accuracy}")
        else:
            print(f"The {self.accuracy_label} accuracy for matched {tetype}: N/A")
    
    def convert_to_ploidy(self):
        # self.nHap
        # self.truth_file
        with open(self.truth_file, 'r') as fin, open("polyhap.vcf", 'w') as fout:
            for line in fin:
                if line.startswith('##'):
                    fout.write(line)
                elif line.startswith('#CHROM'):
                    fields = line.strip().split('\t')
                    fixed_cols = fields[:9]
                    samples = fields[9:]
                    if len(samples) % self.nHap != 0:
                        raise ValueError(f"Error: number of haplotypes ({len(samples)}) is not divisible by ploidy ({self.nHap})")
                    merged_samples = []
                    for i in range(0, len(samples), self.nHap):
                        group = samples[i:i+self.nHap]
                        merged_name = "_".join(group)
                        merged_samples.append(merged_name)
                    fout.write('\t'.join(fixed_cols + merged_samples) + '\n')
                else:
                    fields = line.strip().split('\t')
                    fixed_cols = fields[:9]
                    genotypes = fields[9:]
                    merged_gts = []
                    for i in range(0, len(genotypes), self.nHap):
                        group = genotypes[i:i+self.nHap]
                        merged_gt = '|'.join(group)
                        merged_gts.append(merged_gt)
                    fout.write('\t'.join(fixed_cols + merged_gts) + '\n')

    def calculate_confusion_counts(self, bench_pos, bench_gt, bench_desc,
                                   compare_pos, compare_gt, compare_desc, chrom):
        # print("Calculating confusion counts...")
        # get matched data
        i, j = 0, 0
        b_len = len(bench_pos)
        c_len = len(compare_pos)
        match_idxb = []
        match_idxc = []
        while i < b_len and j < c_len:
            pos_a = bench_pos[i]
            pos_b = compare_pos[j]
            if pos_b > pos_a + self.max_dist:
                i += 1
            elif pos_a - self.max_dist <= pos_b <= pos_a + self.max_dist:
                match_idxb.append(i)
                match_idxc.append(j)
                i += 1
                j += 1
            else:
                j += 1

        # calculate confusion counts
        tp = len(match_idxb)
        self.nMatch += tp
        #print(f"Number of matched TEs: {self.nMatch}")
        # FP: only in comparison VCF
        fp = c_len - tp
        # FN: only in bench vcf
        fn = b_len - tp
        bgt = [bench_gt[i] for i in match_idxb]
        cgt = [compare_gt[i] for i in match_idxc]
        bdesc = [bench_desc[i] for i in match_idxb]
        cdesc = [compare_desc[i] for i in match_idxc]
        # Guard on a non-empty match set: a chromosome can carry truth and prediction
        # variants yet match none of them (tp == 0), leaving bgt/cgt empty -- indexing
        # bgt[0] then raised IndexError and killed the whole comparison.
        if bgt and len(bgt[0]) != len(cgt[0]):
            # print("Warning: the number of allele is different between truth and prediction files")
            # print("Warning: filling the missing alleles in truth VCF, e.g., 1 -> 1/1")
            for i in range(len(bgt)):
                bgt[i] = bgt[i] + bgt[i]
                if bdesc[i] is not None:
                    bdesc[i] = bdesc[i] + bdesc[i]

        # Genotype agreement, by allele content where both files provide sequences and by
        # raw allele index otherwise (BED predictions, symbolic alleles).
        genotype_same = 0
        with open(self.MatchFile + ".csv", "a", newline="") as fo:
            writer = csv.writer(fo, lineterminator="\n")
            for k in range(tp):
                same = genotype_agrees(bdesc[k], cdesc[k], self.gt_len_tol)
                if same is None:
                    same = sorted(bgt[k]) == sorted(cgt[k])
                genotype_same += same
                writer.writerow([chrom,
                                 bench_pos[match_idxb[k]], _fmt_gt(bgt[k]), _fmt_sizes(bdesc[k]),
                                 compare_pos[match_idxc[k]], _fmt_gt(cgt[k]), _fmt_sizes(cdesc[k]),
                                 "match" if same else "MISMATCH"])
        gt_diff = tp - genotype_same

        return tp, fp, fn, gt_diff

def run(args):
    CompareVCF(args)._run()

