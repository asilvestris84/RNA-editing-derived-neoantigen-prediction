# RNA-editing-derived-neoantigen-prediction
Automated and fully reproducible pipeline for the systematic prediction and prioritization of RNA editing–derived neoantigens, based on integrated comparative analysis of RNA-seq and matched WGS/WES data, enabling robust detection of RNA-editing events, recoded peptide generation, MHC binding prediction, and downstream statistical analyses.
![Pipeline overview](/images/fcell-08-00728-g001.jpg)
*Figure adapted from Han X-J, Ma X, Yang L, Wei Y, Peng Y and Wei X (2020) Progress in Neoantigen Targeted Cancer Immunotherapies. Front. Cell Dev. Biol. 8:728. doi: 10.3389/fcell.2020.00728*

## Tools cited

This pipeline uses **Neoantimon** R package (DOI: 10.1093/bioinformatics/btaa901) for neoantigen prediction from annotated variants.  
Neoantimon performs per-sample peptide generation and MHC binding affinity prediction using **netMHCpan** and **netMHCIIpan**.

## Scope of the pipeline

The workflow supports:

- Detection of RNA editing events from RNA-seq data
- Comparative filtering using matched Whole Genome Sequencing (WGS) or Whole Exome Sequencing (WES)
- Identification of recoding RNA editing events
- Generation of edited and unedited peptide sequences
- MHC class I binding prediction
- Systematic prioritization of candidate neoantigens
- Downstream statistical analysis and visualization

## Script descriptions

### `last_Neoantimon_improved.R`
### Neoantigen prediction from RNA editing using Neoantimon (MHC class I and II)

### Function
This script is the core execution step of the pipeline. It runs Neoantimon on ANNOVAR-annotated RNA-seq variants to generate candidate neoantigens derived from RNA editing events. For each sample, it integrates:
- RNA editing–derived variants
- matched RNA expression
- sample-specific HLA typing
- MHC class I and class II binding prediction

### Input requirements
- One or more `hg38_multianno.txt` files (ANNOVAR output) located in the working directory  
- A text file containing one path per line to HLA typing results  
- A text file containing one path per line to RNA expression files  

### User edits required
The following hardcoded paths **must be edited** before execution:
```
hla_paths <- "path/to/hla_paths.txt"
gene_expression_path <- "path/to/gene_expression_paths.txt"
```
---
### `alle_gen.R`
### Add RNA-editing frequency (REDItools) to Neoantimon `.ALL` output tables

**Purpose**  
`alle_gen.R` post-processes Neoantimon per-sample `*ALL*` tables by **adding the RNA-editing frequency** for each editing event, **directly taken from REDItools output tables**.  
For each sample, it matches Neoantimon events to REDItools sites and writes an updated `.ALL` table including a new `Frequency` column.

---

### How it works
For each sample name provided:
1. Finds Neoantimon tables matching `*ALL*` inside a user-provided folder (supports wildcards/globs).
2. Reads the corresponding REDItools outTable for that sample.
3. Builds a REDItools site key as:
   - `Region_Position` (concatenation of `Region` and `Position` from REDItools)
4. Matches Neoantimon `Mutation_Position` against that key and assigns:
   - `Frequency = ot$Frequency[match(...)]`
5. Inserts `Frequency` immediately after `Mutation_Position` and writes the updated table.

---

### Input requirements
You must provide **three arguments** to the script:

1) **Sample name list file** (one sample ID per line)  
This is required because the script uses these sample IDs to select the correct files.

2) **Folder (or wildcard path) containing Neoantimon outputs**, where Neoantimon files include `*ALL*` in the name  
Example wildcard: `"/path/to/result.ID.SNV*"`

3) **A text file listing paths to REDItools outTables** (one per line)  
The script assumes each outTable contains at least these columns:
- `Region`
- `Position`
- `Frequency`

**Neoantimon `.ALL` tables must contain:**
- `Mutation_Position` (used for matching)

---

### User edits required
None. This script is fully driven by command-line arguments.

---

### How to run
```bash
Rscript alle_gen.R <sample_names.txt> "<neoantimon_output_folder_or_glob>" <reditools_outTable_paths.txt>
```
---

### `calc_priority_score2.R`
### PriorityScore computation for RNA editing–derived neoantigens (expression × frequency + binding gain)

**Purpose**  
`calc_priority_score2.R` computes a **hierarchical PriorityScore** to systematically rank candidate neoantigens by integrating:  
1) **RNA expression support** (`Total_RNA`)  
2) **Editing frequency** (`Frequency`)  
3) **Mutant vs wild-type binding improvement** (via `Mut_Rank` and `Wt_Rank`)

This script is designed to prioritize candidates that are **supported by RNA evidence** (expressed and frequent) while still accounting for the **gain in predicted binding**.

---

### Input requirements
Provide one or more **tab-delimited** files (e.g. `.txt`, `.tsv`, `.tab`) containing at least these columns:

- `Mut_Rank`  (mutant peptide binding rank/percentile)
- `Wt_Rank`   (wild-type peptide binding rank/percentile)
- `Frequency` (RNA editing frequency nG/(nA+nG)
- `Total_RNA` (RNA abundance TPM)

If any of these columns are missing, the script stops with an error.

---

### How the PriorityScore is computed

#### 1) Binding improvement (ΔRank)
The script computes the mutant-to-wild-type ratio:
\[
\Delta Rank = \frac{Mut\_Rank}{Wt\_Rank}
\]
Lower values (< 1) indicate that the mutant is predicted to bind better than the wild-type.

#### 2) Normalization
The script applies min–max normalization:

**Positive normalization** (higher = better):
\[
x_{norm} = \frac{x - \min(x)}{\max(x) - \min(x)}
\]
applied to:
- `Total_RNA` → `RNA_norm`
- `Frequency` → `Freq_norm`

**Inverse normalization** (lower = better, then inverted):
\[
x_{invnorm} = 1 - \frac{x - \min(x)}{\max(x) - \min(x)}
\]
applied to:
- `DeltaRank` → `DeltaRank_norm`

**Edge case handling**  
If a variable has no range (max = min) or non-finite values, the script assigns `0.5` to all entries for that normalized feature.

#### 3) Final PriorityScore formula
The final score is:
\[
PriorityScore = (RNA_{norm} \times Freq_{norm}) + 0.25 \times \Delta Rank_{norm}
\]

- `RNA_norm × Freq_norm` enforces that top candidates are **both expressed and frequent**
- `0.25 × DeltaRank_norm` adds a weighted contribution from **binding improvement**
- The weight **α = 0.25** is explicitly hardcoded in the script.

---

### How to run
Run the script by passing one or more input files as arguments:

```bash
Rscript calc_priority_score2.R file1.tsv file2.tsv
```
---

### `comm_plotgra.R`
### Neoantigen co-occurrence analysis and non-redundant strong binder export

### Purpose
This script performs a downstream aggregation of Neoantimon output tables to:
1. Compute **sample-by-sample co-occurrence** of RNA editing–derived neoantigens.
2. Generate a **non-redundant table of strong-binding neoantigens** shared across samples.

The analysis is based exclusively on **strong binders**, defined using the Neoantimon binding rank.

---

### Input requirements
- A directory containing Neoantimon output tables with filenames matching: filtered_*.ALL.txt
Each `.ALL.txt` file must contain at least the following columns:
- `Mut_Rank`
- `Evaluated_Mutant_Peptide`
- `Gene`
- `NM_ID`
- `Change`
- `Mutation_Position`

---

### Strong binder definition
A neoantigen is considered a **strong binder** if: Mut_Rank ≤ 0.5 for MHC I; Mut_Rank ≤ 1 for MHC II
The script explicitly converts `Mut_Rank` to numeric before applying this threshold.

---

### What the script does

#### 1. Neoantigen collection per sample
For each Neoantimon `.ALL.txt` file:
- Strong binders are filtered based on `Mut_Rank`.
- Unique mutant peptide sequences (`Evaluated_Mutant_Peptide`) are stored per sample.

#### 2. Co-occurrence matrix computation
- A binary presence/absence matrix (peptide × sample) is constructed.
- A sample-by-sample co-occurrence matrix is computed as:

C = t(P) %*% P

where `P` is the binary presence matrix.
- The diagonal and upper triangle are masked to avoid redundancy.
- Sample identifiers are extracted from filenames (e.g. `ICGC_XXXX`).
#### 3. Non-redundant strong binder table
- Strong binders from all samples are pooled.
- Entries are collapsed by `Evaluated_Mutant_Peptide`.
- For each unique peptide, the following fields are reported:
- `Gene` (numeric prefixes removed)
- `NM_ID`
- `Change`
- `Mutation_Position`
- `Peptide_Sequence`
- `Samples` (comma-separated list)
- `N_Samples` (number of distinct samples)

---

### User edits required
Before execution, the Neoantimon output directory must be edited at the top of the script:

```r
data_folder <- "path/to/neoantimon_output_directory/"
```
The directory must contain the filtered_*.ALL.txt files.
---
---
### `make_PS_heatmap.R`
### Heatmap of neoantigen PriorityScore across samples

### Purpose
This script generates a **peptide-by-sample heatmap** summarizing the **PriorityScore** of RNA editing–derived neoantigens across samples.  
Each row represents a unique neoantigen (gene + mutant peptide), while each column represents a sample.

The heatmap is intended to visualize the distribution and recurrence of high-priority neoantigens across the cohort.

---

### Input requirements
- A **text file** containing **one full path per line** to PriorityScore tables produced by `calc_priority_score2.R`.

Each PriorityScore table must contain at least the following columns:
- `Gene`
- `Evaluated_Mutant_Peptide`
- `PriorityScore`

---

### User edits required
Before execution, edit the following parameter at the top of the script:

```r
file_list <- "path/to/priority_score_file_list.txt"
```
---

### `bubble_plot.R`
### Bubble plot of RNA editing–derived neoantigens across GTEx tissues

### Purpose
This script generates a **bubble plot** summarizing the distribution of RNA editing–derived neoantigens across **GTEx healthy tissues**.  
The visualization integrates, for each RNA editing site:

- tissue specificity (GTEx tissues)
- median RNA editing frequency
- prevalence across healthy samples

The plot is designed to contextualize candidate neoantigens with respect to their presence in normal tissues.

---

### Input requirements
A **tab-delimited file** containing GTEx annotations for RNA editing sites, with at least the following columns:

- `Gene`
- `Mutation_Position_Hg38`
- `Healthy_Tissues`  
- `Median_Editing_Frequency_.`
- `SRR_Samples_.`

The script assumes that:
- `Healthy_Tissues` may contain **comma-separated values**
- `Median_Editing_Frequency_.` and `SRR_Samples_.` may also be comma-separated and are coerced to numeric

---

### User edits required
Before execution, the **input file path must be edited** in the script:

```r
data <- read.delim("path/to/common_strong_binders_with_GTEx_data.tsv", stringsAsFactors = FALSE)
```
## Input data requirements

The pipeline assumes availability of:

- RNA-seq–derived RNA editing calls
- Matched WGS or WES variant data for genomic variant exclusion
- Gene and transcript annotations
- HLA typing information
- Expression quantification data (gene or transcript level)

**Raw sequencing data are not included in this repository.**

---

## Reproducibility

All scripts provided in this repository were used to generate the results reported in the associated manuscript.  
Exact parameter choices, filtering criteria, and thresholds are documented within the scripts.

---

## Code availability

All scripts used to perform the analyses presented in this study are publicly available in this GitHub repository.

---

## Citation

If you use this pipeline or code, please cite the associated manuscript and this repository.  
A DOI will be provided upon publication or repository archiving.

---

## License

This repository is distributed for academic and research use.  
Please refer to the LICENSE file for details.

---

## Contact

For questions or issues related to the pipeline, please open an issue in this repository or contact me at domenico.silvestris@uniba.it.
