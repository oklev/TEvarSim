# TEvarSim

**TEvarSim** is a versatile genome simulation tool for generating polymorphic transposable element (TE) variants.   

Key features:
- Supports TE insertions, deletions, and LTR-LTR recombinations (excisions of full-length LTR elements into solo LTRs)
-	Simulates both real and random TE variants
-	Simulates both short- and long-read sequencing data
-	Supports population-scale genome simulation
-	Generates VCF files
-	Includes scripts to compare predicted vs. simulated variants, per genome (`Compare`) or across every simulated event at once (`Evaluate`)
---

## Installation

The version of TEvarSim available through pip is outdated; install it from the github repo instead.

```bash
# gfatools and repeatmasker are used for *TEpan*.
# mason and pbsim3 are used for short-reads simulation and long-reads simulation, respectively.
# You may skip installing the software if you do not use the corresponding functionality.
conda create -n tevarsim -c bioconda gfatools repeatmasker mason pbsim3
conda activate tevarsim
#git clone TEvarSim
cd TEvarSim
pip install .
```
## Quick start
Example data can be found in the **testData** directory   

The number of each event class is specified explicitly with `--nINS` (insertions), `--nDEL`
(deletions), and `--nEXC` (LTR-LTR recombinations / excisions). Set at least one to a value > 0.

**1. Simulate 4 insertions and 2 deletions from known TEs**
```bash
tevarsim TEreal --knownINS MEI.fa --existingTEs rptmsk.out --CHR 21 --nINS 4 --nDEL 2
```
- `MEI.fa` is known pTE insertion, from paper [Logsdon, G.A. et al. Nature, 2025](https://www.nature.com/articles/s41586-025-09140-6)  
- `rptmsk.out` is the repeatmasker output for chr21_tiny.fa.

**1b. Add LTR-LTR recombinations (excisions)**
```bash
tevarsim TEreal --knownINS MEI.fa --existingTEs ltr_element.out --CHR 21 \
    --TEtype Gypsy --nDEL 1 --nEXC 1
```
- An excision converts a full-length LTR element (LTR-Internal-LTR) into a solo LTR.
- Full-length elements are identified from the fragment structure in a RepeatMasker `.out`
  file, so `--nEXC` requires `.out` input (not UCSC `.txt`). Include the relevant LTR family
  in `--TEtype`. `testData/ltr_element.out` is a small demonstration fixture.

**2. Simulate deletions plus random TE insertions**
```bash
tevarsim TErandom --consensus human_TE.fa --existingTEs rptmsk.out --ref chr21_tiny.fa --nINS 4 --nDEL 2 --CHR chr21,chr22
```
- `TEconsensus.fa` is human TE consensus sequences from Dfam
- Add `--nEXC N` (with a RepeatMasker `.out` `--existingTEs` file) to also simulate excisions.
- You can try adding `--CHR chr21,chr22` to restrict the analysis to just two 'chromosomes.'
- You can try adding `--regions regions.bed` or `--exclude regions.bed` to test out the functionality to restrict the simulation to certain regions or exclude them from the simulation.

**3. Simulate pTE from pangenome graph**
```bash
# Fetch pangenome graph from HPRC
curl https://human-pangenomics.s3.amazonaws.com/pangenomes/freeze/freeze1/minigraph/hprc-v1.0-minigraph-grch38.gfa.gz > hprc-v1.0-minigraph-grch38.gfa.gz
tevarsim TEpan --gfa hprc-v1.0-minigraph-grch38.gfa.gz --lib Homo_sapiens_DFAM.fa  --CHR chr21 --nINS 4 --nDEL 2
```
- `hprc-v1.0-minigraph-grch38.gfa.gz` is downloaded from [HPRC](https://data.humanpangenome.org/alignments)
- LTR-LTR recombination (`--nEXC`) is not yet supported for TEpan.

**4. Simulate 10 genomes with 6 pTE**  
```bash
tevarsim Simulate --ref chr21_tiny.fa --bed TEreal.bed --num 10 --pool MEI.fa
# if you want to generate sequence vairiations of the same TE between genomes, run below commonds
tevarsim Simulate --ref chr21_tiny.fa --bed TEreal.bed --num 10 --pool MEI.fa --diverse --diverse_config diverse.config
```
- `chr21_tiny.fa` is the reference sequence
- `real.bed` is the position of pTE positions that generated from `tevarsim TEreal`
- `diverse` : Introduce sequence diversity among individuals for the same TE event (which is suitable for evaluating methods that require a TE panel as input)
- `diverse_config` : A configuration file of parameters for introducing sequence diversity among individuals for the same TE event (optional)

**5. Generate sequencing reads from simulated genome** 
```bash
tevarsim Readsim --type short --genome Sim.fa --depth 10 
tevarsim Readsim --type long --genome Sim.fa --depth 10
```
- `type` : short reads or long reads

**6. Benchmark a caller against the simulation**
```bash
# one genome at a time: how well was this genome called?
tevarsim Compare --truth Sim.vcf --pred calls.vcf --truthID Sim_0 --predID Sim_0 -O compare_Sim_0
# every simulated event at once, over all genomes: which events were recovered, and where does it fail?
tevarsim Evaluate --truth Sim.vcf --pred calls.vcf -O evaluate
```
- `Compare` needs a genome named on both sides; `Evaluate` pairs the genomes by name itself and
  reports per-event results aggregated across all of them

**7. The complete workflow** 
- We offered the complete workflow for the short/long-read based tool benchmarking. Please see the file `workflow.sh`

## Flowchart
![flowchart](https://github.com/JanMiao/TEvarSim/blob/main/chart.png)  
- The known TE deletion information can be obtained from [UCSC annotaion file (.txt)](https://genome.ucsc.edu/cgi-bin/hgTables) or [repeatmasker annotation (.out)  ](https://www.repeatmasker.org/genomicDatasets/RMGenomicDatasets.html)
- The known TE insertion position can be obtained from our pre-built dataset (data/MEI_Callset_GRCh38.ALL.20241211.fasta). Any TE insertion sequence is acceptable , as long as the sequence ID follows the naming format **CHR-POS-ID**, e.g., **chr1-683234-AluSp#SINE/Alu**

## Usage
TEvarSim provides seven main command-line subcommands:
```bash
tevarsim <subcommand> [options]
```

### 1. TErandom
Generate pTE position from known deletion sites and random TE insertion.

**Required arguments:**
- `consensus` : Path to the TE consensus FASTA file. The sequenceIDs in the FASTA header should be >TEname#class/superfamily, e.g., >AluY#SINE/Alu
- `ref` : Reference genome FASTA
  
**Optional arguments:**
- `existingTEs` : Known TE annotation file (RepeatMasker .out or UCSC .txt); required when `--nDEL` or `--nEXC` > 0
- `CHR` : Chromosome to simulate TE insertions on (e.g., chr21 or 21)
- `nINS` : Number of random TE insertions to simulate (default: 0)
- `nDEL` : Number of TE deletions to simulate from `--existingTEs` (default: 0)
- `nEXC` : Number of LTR-LTR recombinations (excisions of full-length LTR elements) to simulate; requires a RepeatMasker `.out` `--existingTEs` file (default: 0)
- `outprefix` : Output prefix for TE pool FASTA (default: TErandom)
- `TEtype` : Which TE super families to be extracted from the TE deletion file (default: Alu, L1, and SVA). Specify the TE type by `--TEtype Alu --TEtype L1`
- `DELlen` : A minimum length of known TE deletions to be considered for simulating pTE deletions (default: 100 bp)
- `nMIN` : A minimum number of TE deletions for each TE super family to be simulated (default: 0)
- `TEdistance` : A minimum length of distance between two TE insertions (default: 500 bp)
- `nSV` : Number of background structural variants to simulate (default: 0)
- `regions` : Bed file with regions to include in the simulation.
- `exclude` : Bed file with regions to exclude from the simulation.
- `snp-rate` : SNP mutation rate per base (default: 0.02)
- `indel-rate` : INDEL mutation rate per base (default: 0.005)
- `indel-ins` : Proportion of indels that are insertions (default: 0.4)
- `indel-geom-p` : Geometric distribution parameter for indel lengths (default: 0.7)
- `truncated-ratio` : Proportion of sequences to truncate (default: 0.3)
- `truncated-max-length` : Maximum proportion of sequence to truncate (default: 0.5)
- `polyA-ratio` : Proportion of sequences to add polyA tail (default: 0.8)
- `polyA-min` : Minimum polyA length (default: 5)
- `polyA-max` : Maximum polyA length (default: 20)
- `seed` : Random seed (default: None)


### 2. TEreal
Automatically generate pTE positions from RepeatMasker or UCSC repeat annotations.

**Required arguments:**  
- `knownINS` : Known TE insertion file (The sequence ID should follow the naming format CHR-POS-ID, e.g., chr1-683234-AluSp#SINE/Alu)  
- `existingTEs` : Known TE annotation file (from RepeatMasker `.out` or UCSC `.txt`)  
- `CHR` : Chromosome used to simulate pTE  

**Optional arguments:**  
- `nINS` : Number of TE insertions to simulate from `--knownINS` (default: 0)
- `nDEL` : Number of TE deletions to simulate from `--existingTEs` (default: 0)
- `nEXC` : Number of LTR-LTR recombinations (excisions of full-length LTR elements into solo LTRs); requires a RepeatMasker `.out` `--existingTEs` file, and the LTR family must be included in `--TEtype` (default: 0)
- `DELlen` : A minimum length of known TE deletions to be considered for simulating pTE deletions (default: 100 bp)
- `nMIN` : A minimum number of TE deletions for each TE super family to be simulated (default: 0)
- `nSV` : Number of background structural variants to simulate (default: 0)
- `outprefix` : Output prefix for BED file (default: real)  
- `TEtype` : Which TE super families to be extracted from the TE deletion file (default: Alu, L1, and SVA). Specify the TE type by `--TEtype Alu --TEtype L1`
- `seed` : Random seed (default: None)  

### 3. TEpan
Generate pTE position from Pangenome graph.

**Required arguments:**  
- `gfa` : GFA file of the pangenome graph
- `lib` : RepeatMasker library file [pre-formatted libraries](https://www.repeatmasker.org/~cgoubert/GraffiTE_libraries/).
- `CHR` : Chromosome used to simulate pTE  

**Optional arguments:**  
- `outprefix` : Output prefix for BED file (default: TEpan)
- `nINS` : Number of TE insertions to simulate from the pangenome graph (default: 0)
- `nDEL` : Number of TE deletions to simulate from the pangenome graph (default: 0)
- `minLen` : Minimum length of structural variants to consider (default: 250)  
- `cov` : Minimum TE coverage to consider a structural variant as TE (0-1, default: 0.5)
- `tmpDir` : Temporary directory for intermediate files (default: tmp_TEpan)
- `TEtype` : TEs to be extracted from the TE deletion file, with the default set as LINE, SINE, LTR, and RC. Specify the TE type by `--TEtype LINE --TEtype SINE`
- `seed` : Random seed (default: None)  

### 4. Simulate
Simulate pTE insertions/deletions and generate VCF and modified genome FASTA.

**Required arguments:**
- `ref` : Reference genome FASTA  
- `te-pool` : TE pool FASTA  
- `bed` : BED file of TE positions (can be generated by `ppte TEreal`)  
- `num` : Number of simulated genomes

**Optional arguments:**  
- `diverse` : Introduce sequence diversity among individuals for the same TE event, which is suitable for evaluating methods that require a TE panel as input(default: False)
- `diverse_config` : Path to a configuration file for introducing sequence diversity among individuals for the same TE event (default: None; requires --diverse)
- `outprefix` : Output prefix (default: Sim)  
- `af-min / --af-max` : Min/max allele frequency (default: 0.1/0.9)  
- `tsd-min / --tsd-max` : Min/max TSD length (default: 5/20)  
- `sense-strand-ratio` : Proportion of sense-strand insertions (default: 0.5)  
- `seed` : Random seed (default: None)  

**Event types in the output VCF** (`INFO/TYPE`):
- `INS` : a TE insertion (point event; the inserted sequence comes from the TE pool).
- `DEL` : a TE deletion (a reference span is removed).
- `EXC` : an LTR-LTR recombination. `POS` is the 5' start of the full-length element, `REF` is
  the full `LTR-Internal-LTR` element, and `ALT` is the solo LTR left behind. `INFO/LTRLEN`
  gives the solo LTR length and `INFO/SVLEN` the (negative) net length change. Genotype 0
  keeps the full element; genotype 1 carries the solo LTR. Excision events are encoded in the
  BED with an `EXC-` ID prefix and a 7th column carrying the LTR length.


### 5. Readsim
Generate short or long reads from the simulated genome.

**Required arguments:**
- `type` : Type of reads to simulate (short or long)  
- `genome` : Reference genome file (FASTA) where reads will be simulated from  
- `depth` : Depth of simulated reads    

**Optional arguments(long reads only):**  
- `Lerror` : Sequencing error rate for long reads (default: 0.15)  
- `Lmean` : Average read length for long reads (default: 9000)   
- `Lstd` : Standard deviation of read length for long reads (default: 7000)    

**Optional arguments(short reads only):**  
- `length` : Read length (default: 150)  
- `Fmean` : Average fragment length (default: 300)   
- `Fstd` : Fragment length standard deviation (default: 30)    

**Optional arguments(general):**  
- `seed` : Random seed for reproducibility (default: None)

### 6. Compare
Compare predicted VCF to the simulated VCF.

**Required arguments:**  
- `truth` : Ground truth VCF (generated by `ppte simulate`)  
- `pred` : Predicted file
- `predType` : Type of the predicted file (VCF or BED, default: VCF)
- `outprefix` : the Matched variants in two VCF file  
- `truthID` : Sample ID in truth VCF  
- `predID` : Sample ID in predicted VCF  

**Optional arguments:**  
- `TEtype`: TE family in truth VCF to consider in the comparison. This is the family, not the
  superfamily — TY1, TY2, TY4 and TY5 are all LTR/Copia. The family is read out of the variant ID,
  ignoring the source locus and the polyA/strand annotation that a known-insertion pool decorates
  it with, so `--TEtype AluY` matches every insertion of `chr*-*-len*-AluY#SINE/Alu-polyA*-strand*`
- `INSonly`: Only compare insertions in truth VCF
- `nHap` : Ploidy level, e.g., 2 for Humans  (default: 2)  
- `max_dist` : Maximum allowed distance for variant matching (default: 100 bp)  

**Example**:
```bash
tevarsim compare --truth sim.vcf --pred variants.vcf --truthID Hap1_Hap2 --predID Sample
```
In simulation files, genomes are named Hap1, Hap2, etc.; for polyploids, combine haplotype IDs with `_` for one individual.  

### 7. Evaluate
Score a prediction against **every simulated event at once**, aggregated over all genomes.

`Compare` answers "how well was *this* genome called?" — it takes one truth sample and one
prediction sample, so benchmarking a whole simulation means looping it over every genome and
stitching the per-genome outputs back together by hand. `Evaluate` turns the question around:
the unit of analysis is the simulated **event**, not the genome. Every record in the truth VCF
is scored once, using all the simulated genomes at the same time, and the per-event results are
aggregated into overall metrics plus breakdowns by event type, TE family, event size and allele
frequency. No genome is named on the command line.

**Required arguments:**
- `truth` : Ground truth VCF (generated by `tevarsim Simulate`)
- `pred` : Predicted file

**Optional arguments:**
- `predType` : Type of the predicted file (VCF or BED, default: VCF). A BED prediction carries no
  allele sequences and no sample columns, so only detection and breakpoint accuracy are reported
- `outprefix` : Output prefix for the per-event JSON (default: `evaluate`)
- `sample_map` : Two-column file pairing truth and prediction samples (`truth_sample<TAB>pred_sample`).
  Default: pair samples by name
- `TEtype` : TE family in truth VCF to consider, e.g. `TY1-FULL` (case-insensitive). This is the
  family, not the superfamily — TY1, TY2, TY4 and TY5 are all LTR/Copia. Default: all families
- `INSonly` : Only evaluate insertions in truth VCF
- `nHap` : Number of haplotype columns per individual in the truth VCF; >1 merges them before
  pairing with the prediction (default: 1)
- `max_dist` : Maximum allowed distance (bp) to consider a prediction and a simulated event the
  same locus (default: 100)
- `gt_len_tol` : Maximum allowed difference (bp) in allele length to consider two alleles the same
  (default: 50)
- `size_bins` : Upper edges (bp) of the event-size strata; repeat the flag per edge
  (default: 100 500 1000 5000 10000)
- `af_bins` : Upper edges of the allele-frequency strata; repeat the flag per edge
  (default: 0.1 0.25 0.5 0.75)

**Example**:
```bash
tevarsim Evaluate --truth sim.vcf --pred variants.vcf --INSonly --TEtype TY1-FULL -O eval
```

**How genomes are paired.** Truth and prediction samples are paired by name, which is the pairing
a per-genome `Compare` loop is normally driven with (`-I $sample -J $sample`); as a convenience a
lone sample on each side is paired whatever it is called. Nothing is ever paired by column order —
a silently wrong pairing would report confident, meaningless genotype accuracy. Use `--sample_map`
when the two files name their genomes differently. If nothing can be paired, detection and
breakpoint accuracy are still reported and the carrier statistics are marked unavailable.

**What is reported.** For each simulated event:
- whether the prediction recovered the locus at all, and whether the allele it called is the same
  one that was simulated (by sequence, or by net length change within `--gt_len_tol`);
- the breakpoint offset and allele length error of the matched call;
- which genomes were called as carriers versus which genomes actually carry it (carrier TP/FP/FN);
- whether the genotype of each paired genome agrees, by allele content.

Events and prediction records are matched one-to-one, preferring allele agreement over breakpoint
proximity. That is what keeps a stacked pair of elements honest: two elements at one locus need two
prediction records for both to count as recovered, rather than both claiming the single record that
happens to sit there.

**Output.** `--outprefix` writes two things:

- `<outprefix>.json` — the run's parameters, the aggregate summary, and one object per simulated
  event, so results can be reloaded for plotting or regression checks.
- `<outprefix>_loci/` — one self-contained JSON per locus, named `<chrom>_<pos>.json` after where
  the event was simulated. Every simulated event gets a file, including the ones the prediction
  never called, and so does every prediction that matched no simulated event (`"kind":
  "unmatched_prediction"`, recording its nearest simulated event so a call landing just outside
  `--max_dist` is distinguishable from one with nothing simulated near it). Two elements stacked at
  one position take a `-2`, `-3` suffix. Each file repeats enough of the run to be read on its own;
  strip its `run` and `kind` keys and what is left is exactly that locus's entry in the summary's
  `events` list, which also carries a `locus_file` pointing back at it. Files left by an earlier run
  of the same prefix are cleared first, so the directory only ever describes the current run.

**Note on the TE family.** The family is parsed exactly as `Compare --TEtype` parses it, so the
labels in the "By TE family" table are the values that flag accepts: `AluY` for an insertion whose
pool sequence ID is `chr21-211282-len320-AluY#SINE/Alu-polyA23-strand+`, since neither the source
locus nor the polyA/strand annotation says which family the element belongs to. A "By TE
superfamily" table is added whenever the superfamily groups events that the family table does not
(the AluY / AluYa5 / AluYb8 of a human simulation all being `Alu`).

