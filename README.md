# RNA-editing-derived-neoantigen-prediction
Automated and fully reproducible pipeline for the systematic prediction and prioritization of RNA editing–derived neoantigens, based on integrated comparative analysis of RNA-seq and matched WGS/WES data, enabling robust detection of RNA-editing events, recoded peptide generation, MHC binding prediction, and downstream statistical analyses.
![Pipeline overview](/images/fcell-08-00728-g001.jpg)
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
Main analysis script implementing the core pipeline for RNA editing–derived neoantigen prediction.  
It integrates RNA editing calls, DNA-based variant filtering, peptide generation, HLA binding prediction, and output formatting for downstream analyses.

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

### `new_recoding_heatmap_matrix2.R`
Builds heatmap-ready matrices focusing on recoding RNA editing events and their derived neoantigens.

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
