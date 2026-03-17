import numpy as np
import random
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# ---------------- sampling TE insertions with min distance ----------------
def sample_TEins(regions, deletions, n: int, TEdistance: int, target_strands: list[str|None]):
    current_regions = np.array([[r[0], r[1], r[2], r[5]] for r in regions], dtype=object)
    
    def keep_TEdistance(regs, chrom, excl_start, excl_end):
        if len(regs) == 0:
            return regs
        
        reg_starts = regs[:, 1].astype(int)
        reg_ends = regs[:, 2].astype(int)
        
        mask = (regs[:, 0] != chrom) | (reg_ends <= excl_start) | (reg_starts >= excl_end)
        keep_regs = regs[mask]
        overlap_regs = regs[~mask]
        split_segments = []
        
        for r in overlap_regs:
            r_chrom, r_start, r_end, r_strand = r[0], int(r[1]), int(r[2]), r[3]
            if r_start < excl_start:
                split_segments.append([r_chrom, r_start, excl_start, r_strand])
            if r_end > excl_end:
                split_segments.append([r_chrom, excl_end, r_end, r_strand])
        
        if not split_segments:
            return keep_regs
        return np.vstack([keep_regs, np.array(split_segments, dtype=object)])

    # Subtract Deletions
    for d in (deletions or []):
        exclusion_start = d[1] - max(TEdistance, 1)
        exclusion_end = (d[2] + TEdistance) if TEdistance else (d[1] + 1)
        current_regions = keep_TEdistance(current_regions, d[0], exclusion_start, exclusion_end)

    # Sample n Positions
    positions = []
    for target_s in target_strands:        
        # Filter by strand if a target is specified
        if target_s is not None:
            eligible_regs = current_regions[current_regions[:, 3] == target_s]
            if len(eligible_regs) == 0:
                eligible_regs = current_regions 
        else:
            eligible_regs = current_regions

        lengths = (eligible_regs[:, 2].astype(int) - eligible_regs[:, 1].astype(int))
        total_len = np.sum(lengths)

        if total_len <= 0:
            raise ValueError("TEdistance too large or genomic space exhausted.")

        # Select region
        r_val = np.random.randint(0, total_len)
        cum_lengths = np.cumsum(lengths)
        idx = np.searchsorted(cum_lengths, r_val, side='right')

        # Calculate coordinates
        prev_cum = cum_lengths[idx-1] if idx > 0 else 0
        offset = r_val - prev_cum
        chosen = eligible_regs[idx]
        actual_pos = int(chosen[1]) + offset
        
        positions.append((chosen[0], actual_pos, chosen[3]))

        # Update regions to enforce TEdistance for the next pos
        current_regions = keep_TEdistance(
            current_regions, chosen[0], 
            actual_pos - TEdistance, 
            actual_pos + TEdistance
        )

    return positions

# ---------------- adding background SVs ----------------
def bgSV(bedin:str, bedout:str, nSV:int, ins_ratio:float, fasta_in:str, fasta_out:str):
    """
    Add background SVs (INS/DEL) to a BED file.
    Args:
        bedin: input BED file (only TE insertions/deletions)
        bedout: output BED file with added SVs
        nSV: number of SVs to add
        ins_ratio: ratio of insertions among the background SVs
        fasta_file: FASTA file of sequences with newly added insertions (INS)
    """
    # parse the input bed file
    TEs = []
    with open(bedin, "r") as fin:
        for line in fin:
            fields = line.strip().split("\t")
            chrom, start, end, teID = fields[0], int(fields[1]), int(fields[2]), fields[3]
            TEs.append((chrom, start, end, teID))

    # number of background SVs
    nINS = int(nSV * ins_ratio)
    nDEL = nSV - nINS
    
    # generating DEL and INS randomly
    SVmin, SVmax = 30, 300
    # DELlens = np.random.uniform(low=SVmin, high=SVmax, size=nDEL)
    INSlens = np.random.randint(low=SVmin, high=SVmax, size=nINS)
    # INS sequences
    bgINS_seqs = []
    for idx, ilen in enumerate(INSlens):
        seq = ''.join(random.choices('ATGC', k=ilen))
        bgINS_seqs.append((f"bgINS_{idx}_{ilen}", seq))
    
    # output all INS sequences to new fasta file
    records = list(SeqIO.parse(fasta_in, "fasta"))
    new_records = [SeqRecord(Seq(i), id=j) for j, i in bgINS_seqs]
    records.extend(new_records)
    SeqIO.write(records, fasta_out, "fasta")
    
    # intervals for background SV
    existing_intervals = [(i[1], i[2]) for i in TEs]
    empty_intervals = []
    prev_end = existing_intervals[0][1]
    for start, end in existing_intervals:
        if start > prev_end:
            empty_intervals.append((start - prev_end, prev_end, start))
        prev_end = max(prev_end, end)
    
    # sample DEL and INS positions
    # not allow for too much background SVs
    max_DEL = len(empty_intervals) // 2
    if nDEL > max_DEL:
        raise ValueError(f"Too much background deletions. We recommend to reduce the number of background deletions < {max_DEL}.")
    # select longer empty intervals for DEL
    empty_intervals.sort(reverse=True)
    candidate_DEL_intervals = empty_intervals[:nDEL]
    candidate_INS_intervals = empty_intervals[nDEL:]
    # background SV positions
    # bgSV_positions = []
    # sample DEL positions
    SVmin = min(empty_intervals[nDEL][0], SVmin)
    for idx, i in enumerate(candidate_DEL_intervals):
        length, e_start, e_end = i
        del_len = random.randint(SVmin, min(SVmax, length))
        del_start = random.randint(e_start, e_end - del_len)
        TEs.append((chrom, del_start, del_start + del_len, f"bgDEL_{idx}_{del_len}"))
    # sample INS positions
    sampled_points = set()
    length_eachINS_intervals = [i[0] for i in candidate_INS_intervals]
    while len(sampled_points) < nINS:
        interval = random.choices(candidate_INS_intervals, weights=length_eachINS_intervals, k=1)[0]
        _, start, end = interval
        point = random.randint(start, end)
        sampled_points.add(point)
    for i, j in zip(sampled_points, bgINS_seqs):
        TEs.append((chrom, i, i+1, j[0]))
    
    # output the new bed file
    with open(bedout, "w") as fout:
        TEs.sort()
        for te in TEs:
            fout.write(f"{te[0]}\t{te[1]}\t{te[2]}\t{te[3]}\n")

# ---------------- select TEs with the restrict of minimum number ----------------
def make_min_TE(TE_list: list, nMIN: int, nTE: int, TEtype: set, target_strands: list):
    # Organize TEs by family
    te_by_family = {}
    for entry in TE_list:
        family = entry[4]
        if family in TEtype:
            te_by_family.setdefault(family, []).append(entry)

    # Feasibility Check
    if len(te_by_family) * nMIN > nTE:
        raise ValueError(f"nMIN error: nMIN ({nMIN}) * num_families ({len(te_by_family)}) exceeds number of deletions requested ({nTE}).")

    selected_te = []
    
    # Satisfy nMIN for each family
    for family, members in te_by_family.items():
        if len(members) < nMIN:
            raise ValueError(f"Family {family} has fewer than {nMIN} members.")
        
        # Pick nMIN randomly to start
        selected_te.extend(random.sample(members, nMIN))
    
    nTE -= len(selected_te)
    if nTE == 0:
        return selected_te
    
    for te in selected_te:
        try:
            target_strands.remove(te[6])
        except ValueError:
            try:
                target_strands.remove(None)
            except ValueError:
                target_strands.pop()

    return selected_te + pick_stranded(
        [te for te in TE_list if te not in selected_te],
        nTE,
        target_strands
    )

def pick_stranded(TE_list: list, nTE: int, target_strands: list):
    selected_te = []
    n_either = target_strands.count(None)
    n_sense = target_strands.count("+")
    n_antisense = nTE - n_either - n_sense

    if n_sense:
        sense_pool = [x for x in TE_list if x[6] == "+"]
        if n_sense > len(sense_pool):
            selected_te += sense_pool
            n_either += n_sense - len(sense_pool)
        else:
            selected_te.extend(random.sample(sense_pool,n_sense))
    if n_antisense:
        antisense_pool = [x for x in TE_list if x[6] == "-"]
        if n_antisense > len(antisense_pool):
            selected_te += antisense_pool
            n_either += n_antisense - len(antisense_pool)
        else:
            selected_te.extend(random.sample(antisense_pool,n_antisense))
    if n_sense or n_antisense and n_either:
        TE_list = [x for x in TE_list if x not in selected_te]
    selected_te.extend(random.sample(TE_list,n_either))

    return selected_te