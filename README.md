# Pertussis genomic transmission: public reproducibility package

This repository accompanies the manuscript **“Local persistence, reseeding, and lineage growth in post-pandemic pertussis resurgence.”** Release `v1.0.0` is the frozen public package used for submission to the *Journal of Infection*.

## Scope

The package contains harmonised public-data derivatives, frozen model inputs, the final core-SNP alignment and trees, lineage and phylogeographic results, posterior summaries, figure source data, rendering scripts, analysis code, tests, and provenance records.

All source genome sequences and surveillance series were publicly available. Sequence accession identifiers are retained, and public source URLs are recorded in the registry or source tables where a single URL is applicable, so that source records can be retrieved from their original repositories. No controlled-access or directly identifying participant data are included.

The Australian specimen-level sequencing-process table was obtained from Supplementary Appendix 2 of Fong et al. (2026). That third-party table is not redistributed here because its independent redistribution terms were not established. Cohort and tree tables retain public BioSample and run accessions, but source-specific sample aliases, exact Ct values, specimen types, profile-success fields, and article-supplied marker annotations are blank for those records. The repository contains only the aggregate calibration outputs and derived summaries used in the manuscript. Consult the source article and its appendix for the original records.

## Repository layout

- `data/derived/`: harmonised public-data derivatives and source registry
- `data/model_inputs/`: frozen inputs to the joint model
- `results/`: final phylogeny, lineage, phylogeographic, cgMLST, and model summaries
- `figures/`: publication figures and their tabular source data
- `scripts/model/`: Stan model
- `scripts/pipeline/`: data-processing and analysis programs
- `scripts/figures/`: figure-rendering programs
- `scripts/qa/`: package manifest builder
- `tests/`: checks of final model tables and reported outputs
- `provenance/`: checksums, citation audit, and software-session information

See `DATA_DICTIONARY.md` for a directory- and table-level description.

## Reuse and reproducibility

Large fitted R objects and compiled Stan objects are excluded because of file size. Their sizes and SHA-256 checksums are recorded in `provenance/EXTERNAL_LARGE_OBJECTS.tsv`; the repository instead includes the Stan source, frozen inputs, posterior summaries, diagnostics, and figure source data.

Some early pipeline scripts document the original high-performance-computing workflow and require the corresponding bioinformatics command-line tools and source genomes downloaded from the listed public accessions. Portable model, summary, figure, and validation workflows can be inspected or rerun from the released inputs. `provenance/R_SESSION_INFO.txt` records the principal R environment. Because no containerised analysis environment is supplied, this release is a results-audit package rather than a one-command reconstruction of every upstream computation.

No repository-level licence has yet been assigned. Copyright and reuse terms therefore remain with the respective rights holders, and third-party source data retain their original terms. A licence should be selected before encouraging code reuse.

## Citation

Please cite the accompanying manuscript. Machine-readable citation metadata are provided in `CITATION.cff`.
