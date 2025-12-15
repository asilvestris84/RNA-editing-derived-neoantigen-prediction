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
For each sample, it matches Neoantimon events to REDItools sites and writes an updated `.ALL` table including a new `Frequency` column. :contentReference[oaicite:0]{index=0}

---

### How it works (what it actually does)
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
This is required because the script uses these sample IDs to select the correct files. :contentReference[oaicite:1]{index=1}

2) **Folder (or wildcard path) containing Neoantimon outputs**, where Neoantimon files include `*ALL*` in the name  
Example wildcard: `"/path/to/result.ID.SNV*"` :contentReference[oaicite:2]{index=2}

3) **A text file listing paths to REDItools outTables** (one per line)  
The script assumes each outTable contains at least these columns:
- `Region`
- `Position`
- `Frequency` :contentReference[oaicite:3]{index=3}

**Neoantimon `.ALL` tables must contain:**
- `Mutation_Position` (used for matching) :contentReference[oaicite:4]{index=4}

---

### User edits required
None. This script is fully driven by command-line arguments. :contentReference[oaicite:5]{index=5}

---

### How to run
```bash
Rscript alle_gen.R <sample_names.txt> "<neoantimon_output_folder_or_glob>" <reditools_outTable_paths.txt>
```
---

### `calc_priority_score2.R`
Computes a composite prioritization score for candidate neoantigens based on multiple features, including predicted MHC binding, expression-related metrics, and recoding properties.

---

### `bubble_plot.R`
Generates bubble plots to visualize prioritized neoantigens, typically integrating binding affinity, expression level, and prioritization score.

---

### `comm_plotgra.R`
Produces comparative graphical summaries across experimental conditions or sample groups.

---

### `make_PS_heatmap.R`
Generates heatmaps summarizing prioritization scores across samples or conditions.

---

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
