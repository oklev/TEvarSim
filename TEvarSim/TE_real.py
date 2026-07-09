import random
import os
import logging
from .utils import bgSV, make_min_TE

def check_output_file(output_path):
    if os.path.exists(output_path):
        raise FileExistsError(f"Output file '{output_path}' already exists. Please choose a different name.")

def CHRnorm(chrA,chrB):
    ### normalize chrA according the style if chrB
    if chrA.isdigit() and chrB.isdigit():
         return chrA
    elif chrA.isdigit() and not chrB.isdigit():
        if not chrB[3:].isdigit():
            raise ValueError(f"Your Chromosome name is {chrB}, it must be either plain numbers or prefixed with chr")
        else:
            chrA = chrB[:3] + chrA
    elif not chrA.isdigit() and chrB.isdigit():
            chrA = chrA[3:]
    else:
        if chrA[:3].upper() != chrB[:3].upper():
            raise ValueError(f"Your Chromosome name is {chrB}, it must be either plain numbers or prefixed with chr")
        else:
            chrA = chrB[:3] + chrA[3:]
    return chrA

class RealTE:
    def __init__(self, args):
        # self.reference = args.reference
        self.INSfile = args.knownINS
        self.DELfile = args.existingTEs
        self.TEtype = args.TEtype
        self.DELlen = args.DELlen
        self.CHR = args.CHR
        self.nINS = args.nINS
        self.nDEL = args.nDEL
        self.nEXC = args.nEXC  # number of LTR-LTR recombinations (excisions)
        self.nMIN = args.nMIN
        self.nSV = args.nSV
        self.outprefix = args.outprefix
        self.random_seed = args.seed
        if self.random_seed is not None:
            random.seed(self.random_seed)

    def _run(self):
        """
        Entry point for the `generate-bed` subcommand.
        Dispatches to appropriate parser depending on input file type.
        """
        # 1. 解析文件
        self.parse_INS_file()
        _, ext = os.path.splitext(self.DELfile)
        if ext.lower() == ".txt":
            self.parse_DEL_ucsc()
        elif ext.lower() == ".out":
            self.parse_DEL_repeatmasker()
        else:
            raise ValueError("DEL file name must end with .txt (UCSC) or .out (RepeatMasker).")

        # 2. Remove overlapping regions between INS and DEL
        self.remove_dup_DEL()

        # 3. Generate BED file
        self.build_bed()
        logging.info(f"Generated TE BED file: {self.outprefix}.bed")

        # 4. Add background SVs if specified
        if self.nSV > 0:
            bedin = self.outprefix + ".bed"
            bedout = self.outprefix + ".bgSV.bed"
            fastaout = os.path.splitext(self.INSfile)[0] + ".bgSV.fa"
            bg_ins_ratio = self.nINS / (self.nINS + self.nDEL) if (self.nINS + self.nDEL) > 0 else 0.5
            bgSV(bedin, bedout, self.nSV, bg_ins_ratio, self.INSfile, fastaout)
            logging.info(f"Generated BED file with background SVs: {bedout}")

    def build_bed(self):
        # 1. select excisions (full-length LTR elements) first
        if self.nEXC > len(self.EXC):
            raise ValueError(
                f"Requested {self.nEXC} excisions but only {len(self.EXC)} full-length LTR elements "
                f"available in --existingTEs (check --TEtype: the LTR family must be included)."
            )
        EXClist = random.sample(self.EXC, self.nEXC)
        exc_spans = [(e[0], e[1]) for e in EXClist]
        # drop deletion candidates overlapping a chosen excision element
        del_pool = [
            d for d in self.DEL
            if not any(d[0] < ee and d[1] > es for (es, ee) in exc_spans)
        ]
        # 2. select insertions and deletions
        if self.nINS > len(self.INS):
            raise ValueError(f"Requested {self.nINS} insertions but only {len(self.INS)} available in --knownINS.")
        if self.nDEL > len(del_pool):
            raise ValueError(
                f"Requested {self.nDEL} deletions but only {len(del_pool)} available in --existingTEs "
                f"after removing excision loci."
            )
        if self.nMIN > 0:
            INSlist = make_min_TE(self.INS, self.nMIN, self.nINS, self.TEtype)
            DELlist = make_min_TE(del_pool, self.nMIN, self.nDEL, self.TEtype)
        else:
            INSlist = random.sample(self.INS, self.nINS)
            DELlist = random.sample(del_pool, self.nDEL)
        merged = INSlist + DELlist + EXClist

        merged.sort()
        # bedfile output
        bed_name = self.outprefix + ".bed"
        check_output_file(bed_name)
        with open(bed_name, "w") as f:
            for start, end, teID, *rest in merged:
                cols = [self.CHR, str(start), str(end), teID]
                # EXC records carry (family, "EXC", ltr_len); emit strand + ltr_len columns
                if len(rest) >= 3 and rest[1] == "EXC":
                    cols += [".", "+", str(rest[2])]
                f.write("\t".join(cols) + "\n")

    
    def parse_DEL_ucsc(self):
        """
        Parse UCSC repeat annotation file (.txt) and return a list of parsed records.
        self.DEL
        """
        self.DEL = []
        self.EXC = []  # UCSC .txt lacks LTR fragment structure; excision is unsupported here
        if not self.TEtype:
            self.TEtype = {'Alu', 'L1', 'SVA'}
        ## nomalize the chr ID
        #with open(self.DELfile, "r") as f:
        #    first_line = f.readline().strip()
        #    chrom = first_line.split('\t')[5]
        #    print(f"chromosome from DEL file: {chrom}")
        #    self.CHR = CHRnorm(self.CHR, chrom)
        # process file
        with open(self.DELfile) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                fields = line.strip().split('\t')
                chrom = fields[5]
                # nomalize the chr ID
                self.CHR = CHRnorm(self.CHR, chrom)
                # repClass = fields[12][:3]
                repClass = fields[12]
                if chrom == self.CHR and repClass in self.TEtype:
                    start = int(fields[6])
                    end = int(fields[7])
                    if end - start > self.DELlen:
                        teID = f"DEL-{chrom}-{start}-{end}-{fields[12]}-{fields[10]}"
                        self.DEL.append((start, end, teID, repClass, "DEL"))
                    # self.TEpool[teID] = fasta.fetch(self.CHR, start, end)

    def parse_DEL_repeatmasker(self):
        """
        Parse repeatmasker output file (.out).
        self.DEL : per-line TE deletions (existing behaviour).
        self.EXC : full-length LTR elements (LTR-I-LTR), eligible for excision. Detected by
                   grouping fragments that share a RepeatMasker ID (field 15).
        """
        self.DEL = []
        self.EXC = []
        if not self.TEtype:
            self.TEtype = {'Alu', 'L1', 'SVA'}
        groups = {}
        # process file
        with open(self.DELfile) as f:
            for line in f:
                if not line.strip() or not line.strip()[0].isdigit():
                    continue
                fields = line.strip().split()
                chrom = fields[4]
                # nomalize the chr ID
                self.CHR = CHRnorm(self.CHR, chrom)
                if chrom != self.CHR:
                    continue
                start = int(fields[5])
                end = int(fields[6])
                if "/" in fields[10]:
                    repClass = fields[10].split("/")[1]
                    if repClass in self.TEtype and end - start > self.DELlen:
                        teID = f"DEL-{chrom}-{start}-{end}-{fields[10]}-{fields[9]}"
                        self.DEL.append((start, end, teID, repClass, "DEL"))
                # accumulate fragments by RepeatMasker ID for full-length LTR detection
                groups.setdefault(fields[14], []).append(fields)
        for frags in groups.values():
            self._add_exc_candidate(frags)

    def _add_exc_candidate(self, frags):
        """Detect a full-length LTR element (LTR-I-LTR) and record it as an excision candidate."""
        names = [fr[9] for fr in frags]
        if len(frags) < 3:
            return
        # flanked by LTR fragments with an internal (-I) region between them
        if not (names[0].endswith("-LTR") and names[-1].endswith("-LTR")):
            return
        if not any(n.endswith("-I") for n in names):
            return
        chrom = frags[0][4]
        start = min(int(fr[5]) for fr in frags)
        end = max(int(fr[6]) for fr in frags)
        class_fam = frags[0][10]
        family = class_fam.split("/")[1] if "/" in class_fam else class_fam
        if family not in self.TEtype:
            return
        ltr_len = int(frags[0][6]) - int(frags[0][5])  # 5' LTR fragment length
        if not (0 < ltr_len < end - start):
            return
        name = frags[0][9].rsplit("-", 1)[0] + "-FULL"  # strip -LTR, mark full-length
        exc_id = f"EXC-{chrom}-{start}-{end}-{ltr_len}-{class_fam}-{name}"
        self.EXC.append((start, end, exc_id, family, "EXC", ltr_len))
        
    
    def parse_INS_file(self):
        """
        Parse inbuilt TE insertion file to bed file
        self.INS
        """
        # nomalize the chr ID
        self.INS = []
        with open(self.INSfile, "r") as f:
            chrom = f.readline().strip().lstrip(">").split("-")[0]
            self.CHR = CHRnorm(self.CHR, chrom)
        # process file
        with open(self.INSfile) as fin:
            for line in fin:
                if line.startswith(">"):
                    teID = line.strip().lstrip(">")
                    info = teID.split("-")
                    chrom = info[0]
                    if chrom == self.CHR:
                        # bed file is 0-based position
                        start = int(info[1]) - 1 
                        end = int(info[1])
                        repClass = teID.split("-")[3].split("/")[1]  
                        self.INS.append((start, end, teID, repClass, "INS"))


    def remove_dup_DEL(self):
        """
        Remove DELs that overlap with INSs.
        """
        i, j = 0, 0
        rm_idx = []

        while i < len(self.INS) and j < len(self.DEL):
            a_start, a_end, *_ = self.INS[i]
            b_start, b_end, *_ = self.DEL[j]

            # If B[j] overlaps with A[i]
            if a_start < b_end and a_end > b_start:
                rm_idx.append(j)

            if a_end < b_end:
                i += 1
            else:
                j += 1
        # remove
        delete_set = set(rm_idx)
        self.DEL = [item for i, item in enumerate(self.DEL) if i not in delete_set]

def run(args):
    RealTE(args)._run()
