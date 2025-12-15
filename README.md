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
Provide one or more **tab-delimited** files (e.g. `.txt`, `.tsv`, `.tab`) containing at least these columns: :contentReference[oaicite:1]{index=1}

- `Mut_Rank`  (mutant peptide binding rank/percentile)
- `Wt_Rank`   (wild-type peptide binding rank/percentile)
- `Frequency` (RNA editing frequency nG/(nA+nG)
- `Total_RNA` (RNA abundance TPM)

If any of these columns are missing, the script stops with an error. :contentReference[oaicite:2]{index=2}

---

### How the PriorityScore is computed

#### 1) Binding improvement (ΔRank)
The script computes the mutant-to-wild-type ratio:
\[
\Delta Rank = \frac{Mut\_Rank}{Wt\_Rank}
\]
Lower values (< 1) indicate that the mutant is predicted to bind better than the wild-type. :contentReference[oaicite:3]{index=3}

#### 2) Normalization
The script applies min–max normalization:

**Positive normalization** (higher = better):
\[
x_{norm} = \frac{x - \min(x)}{\max(x) - \min(x)}
\]
applied to:
- `Total_RNA` → `RNA_norm`
- `Frequency` → `Freq_norm` :contentReference[oaicite:4]{index=4}

**Inverse normalization** (lower = better, then inverted):
\[
x_{invnorm} = 1 - \frac{x - \min(x)}{\max(x) - \min(x)}
\]
applied to:
- `DeltaRank` → `DeltaRank_norm` :contentReference[oaicite:5]{index=5}

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
Run the script by passing one or more input files as arguments: :contentReference[oaicite:8]{index=8}

```bash
Rscript calc_priority_score2.R file1.tsv file2.tsv
```
---

### `comm_plotgra.R`
Produces comparative graphical summaries across experimental conditions or sample groups.

---

### `make_PS_heatmap.R`
Generates heatmaps summarizing prioritization scores across samples or conditions.

---

### `bubble_plot.R`
Generates bubble plots to visualize prioritized neoantigens, typically integrating binding affinity, expression level, and prioritization score.

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
