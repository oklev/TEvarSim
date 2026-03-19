import numpy as np
import random
import logging
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from .utils import sample_TEins, bgSV, make_min_TE, pick_stranded

def is_dna(seq: Seq) -> bool:
    return set(str(seq)) <= set("ATGCN")

def check_TEid(seq: str) -> bool:
    return "#" in seq and "/" in seq

def check_output_file(output_path):
    if os.path.exists(output_path):
        raise FileExistsError(f"Output file '{output_path}' already exists. Please choose a different name.")

def CHRnorm(chr,chr_list):
    ### normalize chr according the style of chrs in fasta index
    if chr in chr_list:
        return chr
    chr_permutations = {
        chr.replace("chr",""),
        chr.replace("chr","CHR"),
        chr.replace("CHR","chr"),
        chr.replace("CHR",""),
        f"chr{chr}",
        f"CHR{chr}"
    }
    existing_chrs = set(iter(chr_list))
    matches = chr_permutations & existing_chrs
    if len(matches) == 1:
        return matches.pop()
    raise ValueError(f"Chromosome {chr} not found in chromosome list: {','.join(iter(chr_list))}")

class TEtype(list):
    def __init__(self,te_list):
        if te_list is None:
            te_list = []
        super().__init__(te_list)
    def __contains__(self, other):
        if not self:
            return True
        else:
            return super().__contains__(other)

class RandomTE:
    def __init__(self, args):
        self.pool_fasta = args.outprefix + ".fa"
        self.out_fasta = args.outprefix + ".bgSV.fa"
        self.DELfile = args.knownDEL
        self.sense_strand_ratio = args.sense_strand_ratio
        self.CHR = {}
        if not os.path.isfile(f"{args.ref}.fai"):
            raise FileNotFoundError(f"Reference fasta index required: {args.ref}.fai not found.")
        with open(f"{args.ref}.fai") as fin:
            for line in fin:
                line = line.strip().split("\t")
                if line:
                    self.CHR[line[0]] = int(line[1])
        if args.CHR:
            chosen_chrs = [CHRnorm(chr,self.CHR) for chr in args.CHR.split(",")]
            self.CHR = {chr:chr_len for chr,chr_len in self.CHR.items() if chr in chosen_chrs}
        if args.regions:
            try:
                with open(args.regions) as fin:
                    regions = fin.read().strip().split("\n")
            except FileNotFoundError:
                raise FileNotFoundError(f"Regions file {args.regions} not found.")
            self.regions = []
            for region in regions:
                region = region.split("\t")
                if region[0] in self.CHR:
                    region[1] = int(region[1])
                    region[2] = int(region[2])
                    if len(region) > 5:
                        self.regions.append(region)
                    else:
                        self.regions.append(region[:3]+[".",".","+"])
                        self.regions.append(region[:3]+[".",".","-"])
        else:
            self.regions = [[chr, 0, chr_len - 1, ".", ".", strand] 
                for strand in ["+", "-"] 
                for chr, chr_len in self.CHR.items()]
        if args.exclude:
            try:
                with open(args.exclude) as fin:
                    excluded_regions = fin.read().strip().split("\n")
            except FileNotFoundError:
                raise FileNotFoundError(f"Excluded regions file {args.exclude} not found.")
            for excluded_region in excluded_regions:
                excluded_region = excluded_region.split("\t")
                excluded_region[1] = int(excluded_region[1])
                excluded_region[2] = int(excluded_region[2])
                new_regions = []
                for region in self.regions:
                    if region[0] == excluded_region[0] and len(excluded_region) < 6 or region[5] == excluded_region[5]:
                        if excluded_region[1] < region[2] and excluded_region[2] > region[1]:
                            if region[1] < excluded_region[1]:
                                new_regions.append((region[0],region[1],excluded_region[1],".",".",region[5]))
                            if region[2] > excluded_region[2]:
                                new_regions.append((region[0],excluded_region[2],region[2],".",".",region[5]))
                            continue
                    new_regions.append(region)
                self.regions = new_regions

        self.nTE = args.nTE
        self.ins_ratio = args.ins_ratio
        self.TEtype = TEtype(args.TEtype)
        self.prefix = args.outprefix
        self.DELlen = args.DELlen
        self.TEdistance = args.TEdistance
        self.nSV = args.nSV # number of background SVs
        self.nMIN = args.nMIN # minimum number of TEs for each family
        self.random_seed = args.seed
        if self.random_seed is not None:
            np.random.seed(self.random_seed)
            random.seed(self.random_seed)


        if args.sense_strand_ratio is not None:
            n_sense = round(self.nTE*args.sense_strand_ratio)
            self.target_strands = (["+"] * n_sense) + (["-"] * (self.nTE - n_sense))
            np.random.shuffle(self.target_strands)
        else:
            self.target_strands = [None] * self.nTE

    def _run(self):
        # 1. parse DEL file
        self.parse_DEL()
        # 2. parse TEpool
        self.nDEL = min(int(round(self.nTE* (1-self.ins_ratio))),len(self.DEL))
        logging.info(f"Generating {self.nTE - self.nDEL} INS and {self.nDEL} DEL for chromosome(s) {','.join(self.CHR.keys())}")
        if self.nMIN >= self.nTE:
            raise ValueError(f"minumum number of a TE family ({self.nMIN}) should be less than nTE ({self.nTE})")
        if self.nMIN > 0:
            self.DEL = make_min_TE(self.DEL, self.nMIN, self.nDEL, self.TEtype, self.target_strands[:self.nDEL])
        else:
            self.DEL = pick_stranded(self.DEL, self.nDEL, self.target_strands[:self.nDEL])
        self.nINS = self.nTE - len(self.DEL)
        self.parse_TEpool()
        # 3. Generate BED file
        self.build_bed()
        logging.info(f"Generated TE BED file: {self.prefix}.bed")
        # 4. Add background SVs if specified
        if self.nSV > 0:
            bedin = self.prefix + ".bed"
            bedout = self.prefix + ".bgSV.bed"
            bgSV(bedin, bedout, self.nSV, self.ins_ratio, self.pool_fasta, self.out_fasta)
            logging.info(f"Generated BED file and fasta file with background SVs are: {bedout} and {self.out_fasta}")

    def build_bed(self):
        # merge INS and DEL
        merged = self.INS + self.DEL
        merged.sort()
        # bedfile output
        bed_name = self.prefix + ".bed"
        check_output_file(bed_name)
        with open(bed_name, "w") as f:
            for chrom, start, end, teID, __, __, strand, *__ in merged:
                f.write("\t".join([chrom, str(start), str(end), teID, ".", strand]) + "\n")
        logging.info(f"Generated BED file: {bed_name}")
    
    def parse_DEL(self):
        self.DEL = []
        # process file
        ext = self.DELfile.split(".")[-1]
        if ext == "out":
            # Parse RepeatMasker input
            self.parse_DEL_repeatmasker()
            return
        elif ext == "txt":
            # Parse UCSC input
            def readline(line):
                if line.startswith("#"):
                    return
                fields = line.strip().split('\t')
                #       chrom,      repClass    start           end             strand      name        class_fam
                return  fields[5],  fields[12], int(fields[6]), int(fields[7]), fields[9],  fields[10], fields[12]
            return
        elif ext == "bed":
            # Parse bed input
            def readline(line):
                if not line.strip():
                    return
                fields = line.strip().split("\t")
                name, class_fam = fields[3].split("#")
                #       chrom,      repClass                start           end             strand     name    class_fam
                return fields[0],   class_fam.split("/"),  int(fields[1]), int(fields[2]), fields[5],  name,   class_fam
        else:
            raise ValueError(f"Input file format not recognized for --knownDEL: {ext}")
        with open(self.DELfile) as f:
            for line in f:
                try:
                    chrom, repClass, start, end, strand, name, class_fam = readline(line)
                except ValueError as e:
                    if "not enough values to unpack" in str(e):
                        continue
                    raise
                if end - start < self.DELlen:
                    continue
                if repClass in self.TEtype:
                    teID = f"DEL-{chrom}-{start}-{end}-{class_fam}-{name}"
                    self.DEL.append((chrom, start, end, teID, repClass,"DEL",strand))
    
    def parse_DEL_repeatmasker(self):
        repeatmasker_records = {}
        with open(self.DELfile) as f:
            for line in f:
                if not line.strip() or not line.strip()[0].isdigit(): 
                    continue
                line = line.strip().split()
                if line[8] == "C":
                    line[8] = "-"
                for i in (5,6,14):
                    line[i] = int(line[i])
                if line[14] not in repeatmasker_records:
                    repeatmasker_records[line[14]] = []
                repeatmasker_records[line[14]].append(line)
        for te_info in repeatmasker_records.values():
            chrom = te_info[0][4]
            start = min(record[5] for record in te_info)
            end = max(record[6] for record in te_info)
            if end - start < self.DELlen:
                continue
            match_name = {}
            match_class_fam = {}
            match_strand = {"+":0,"-":0}
            for record in te_info:
                if record[9] not in match_name:
                    match_name[record[9]] = 0
                if record[10] not in match_class_fam:
                    match_class_fam[record[10]] = 0
                record_len = record[6] - record[5]
                match_name[record[9]] += record_len
                match_class_fam[record[10]] += record_len
                match_strand[record[8]] += record_len
            if match_strand["+"] > match_strand["-"]:
                strand = "+"
            else:
                strand = "-"
            if len(match_name) == 1:
                name = next(iter(match_name))
                class_fam = next(iter(match_class_fam))
            else:
                ty_i = {n[:-2]:n_len for n,n_len in match_name.items() if n[-2:] == "-I"}
                if ty_i:
                    match_name = ty_i
                    print(te_info)
                match_name = dict(sorted(match_name.items(),key=lambda item: item[1]))
                match_class_fam = dict(sorted(match_class_fam.items(),key=lambda item: item[1]))
                match_max = next(iter(match_name.values()))
                matches = [n for n in match_name if match_name[n] >= 0.2*match_max]
                name = "|".join(matches)
                if ty_i and te_info[0][9][-4:] == "-LTR" and te_info[-1][9][-4:] == "-LTR":
                    name += "-FULL"
                    print(name)
                class_fam = next(iter(match_class_fam))
            repClass = class_fam.split("/")[0]
            if repClass in self.TEtype:
                teID = f"DEL-{chrom}-{start}-{end}-{class_fam}-{name}"
                self.DEL.append((chrom, start, end, teID, repClass,"DEL",strand))
    
    def parse_TEpool(self):
        self.INS = []
        records = list(SeqIO.parse(self.pool_fasta, "fasta"))
        INSpos = sample_TEins(self.regions, self.DEL, self.nINS, TEdistance=self.TEdistance, target_strands=self.target_strands[-self.nINS:])
        for record, (chrom, pos, strand) in zip(records,INSpos):
            self.INS.append((chrom, pos, pos, record.id, record.id.split("/")[1], "INS", strand))

class TEPoolBuilder:
    def __init__(self, args):
        # base
        self.consensus = args.consensus
        self.nTE = args.nTE
        self.ins_ratio = args.ins_ratio
        self.nINS = round(self.nTE * self.ins_ratio)
        self.prefix = args.outprefix
        # SNP and INDEL
        self.snp_rate = args.snp_rate
        self.indel_rate = args.indel_rate
        self.indel_ins = args.indel_ins
        self.indel_geom_p = args.indel_geom_p
        # truncate
        self.truncated_ratio = args.truncated_ratio
        self.truncated_max_length = args.truncated_max_length
        # polyA
        self.polyA_ratio = args.polyA_ratio
        self.polyA_min = args.polyA_min
        self.polyA_max = args.polyA_max
        # other
        self.random_seed = args.seed
        if self.random_seed is not None:
            np.random.seed(self.random_seed)
        logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
        
        # data stored
        self.BASES = ['A', 'C', 'G', 'T']
        self.CHOICES_DICT = {b: [c for c in self.BASES if c != b] for b in self.BASES}
        self.current_seq = []
        self.seqID_suffix = ""

    def _run(self):
        records = list(SeqIO.parse(self.consensus, "fasta"))
        if not records:
            raise ValueError(f"No sequences found in consensus FASTA: {self.consensus}")
        logging.info(f"Found {len(records)} sequences in consensus FASTA.")
        # check if all sequence are legal
        parsed_records = []
        potential_LTRs = {}
        for record in records:
            record.seq = record.seq.upper()
            if not is_dna(record.seq):
                raise ValueError(f"Invalid characters found in {record.id}")
            if not check_TEid(record.id):
                raise ValueError(f"Invalid format in TE ID: {record.id}. Please follow the format '>TEname#class/superfamily', e.g., '>AluY#SINE/Alu'")
            if record.id.split("#")[0][-2:] == "-I" or record.id.split("#")[0][-4:] == "-LTR":
                id, class_fam = record.id.split("#")
                ltr_name, ltr_part = id.rsplit("-",1)
                if ltr_name not in potential_LTRs:
                    potential_LTRs[ltr_name] = {"class_fam":class_fam}
                potential_LTRs[ltr_name][ltr_part] = record
            else:
                parsed_records.append(record)
        if potential_LTRs:
            for ltr_name, ltr_parts in potential_LTRs.items():
                if "I" in ltr_parts and "LTR" in ltr_parts:
                    ltr_rec = ltr_parts["LTR"]
                    internal_rec = ltr_parts["I"]
                    
                    new_id = f"{ltr_name}-FULL#{ltr_parts['class_fam']}"
                    new_record = SeqRecord(
                        ltr_rec.seq + internal_rec.seq + ltr_rec.seq,
                        id=new_id,
                        description=f"Merged LTR-I-LTR sequence for {ltr_name}"
                    )
                    
                    parsed_records.append(new_record)
                else:
                    for part_rec in ltr_parts.values():
                        parsed_records.append(part_rec)
            records = parsed_records
                    

        # generate random masks for truncation and polyA
        truncated_num = np.random.random(self.nINS)
        polyA_num = np.random.random(self.nINS)
        trunc_mask = truncated_num < self.truncated_ratio
        polyA_mask = polyA_num < self.polyA_ratio
        mutate = np.select(
            [trunc_mask & polyA_mask, trunc_mask, polyA_mask],
            [3, 2, 1],
            default=0
        )

        # map numbers to corresponding functions
        func_map = {
            0: lambda: self.INDEL_mutate().SNP_mutate(),
            1: lambda: self.INDEL_mutate().SNP_mutate().apply_polyA(),
            2: lambda: self.INDEL_mutate().SNP_mutate().apply_truncate(),
            3: lambda: self.INDEL_mutate().SNP_mutate().apply_truncate().apply_polyA()
        }

        out_records = []
        for j in mutate:
            idx = np.random.randint(0, len(records))
            record = records[idx]
            self.current_seq = list(str(record.seq))
            self.seqID_suffix = ""  # reset suffix each loop
            func_map[j]()
            new_id = f"{record.id}_{self.seqID_suffix}"
            new_record = SeqRecord(Seq("".join(self.current_seq)), id=new_id, description="")
            out_records.append(new_record)

        # output
        output_fasta = self.prefix + ".fa"
        check_output_file(output_fasta)
        SeqIO.write(out_records, output_fasta, "fasta")
        logging.info(f"Generated {len(out_records)} sequences -> {output_fasta}")
    

    def apply_truncate(self):
        max_trunc = int(len(self.current_seq) * self.truncated_max_length)
        if max_trunc < 1:
            return self
        trunc_len = np.random.randint(1, max_trunc + 1)  # inclusive upper bound
        self.current_seq = self.current_seq[trunc_len:]
        self.seqID_suffix += f"{trunc_len}truncate"
        return self

    def apply_polyA(self):
        polyA_len = np.random.randint(self.polyA_min, self.polyA_max + 1)
        self.current_seq.extend(["A"] * polyA_len)
        self.seqID_suffix += f"{polyA_len}polyA"
        return self
    
    def SNP_mutate(self):
        L = len(self.current_seq)
        n_snp = np.random.poisson(self.snp_rate * L)
        if n_snp < 1:
            return self
        self.seqID_suffix += f"{n_snp}SNP"
        # SNP positions, ensure not exceeding sequence length
        snp_positions = np.random.choice(L, size=min(n_snp, L), replace=False)
        for pos in snp_positions:
            current_base = self.current_seq[pos]
            if current_base in self.CHOICES_DICT:
                self.current_seq[pos] = np.random.choice(self.CHOICES_DICT[current_base])
        return self
    
    def INDEL_mutate(self):
        # total INDEL number
        L = len(self.current_seq)
        n_indel = np.random.poisson(self.indel_rate * L)
        if n_indel < 1:
            return self
        self.seqID_suffix += f"{n_indel}INDEL"
        if n_indel == 0:
            return self
        # generate INDEL lengths
        indel_len_list = np.random.geometric(self.indel_geom_p, n_indel)
        indel_len_list = np.clip(indel_len_list, a_min=None, a_max=30)
        # select positions for INDELs
        indel_positions = np.random.choice(L, size=min(n_indel, L), replace=False)
        for pos, indel_len in zip(sorted(indel_positions, reverse=True), indel_len_list):
            if np.random.random() < self.indel_ins:
                ins_seq = [np.random.choice(self.BASES) for _ in range(indel_len)]
                self.current_seq[pos:pos] = ins_seq
            else:
                del_end = min(pos + indel_len, len(self.current_seq))
                del self.current_seq[pos:del_end]
        return self

def run(args):
    # generate TE pool
    TEPoolBuilder(args)._run()
    # generate bed file with random TE insertions
    RandomTE(args)._run()


