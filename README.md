# Sampled ancestry and MT28-associated lineage expansion in pertussis resurgence

Public data, code, figures, and audit records for the *Journal of Infection*
Letter **“Sampled ancestry and MT28-associated lineage expansion in
post-pandemic pertussis resurgence.”**

This `main` branch uses a Letter-first design. The original long-form public
package remains permanently available at tag [`v1.0.0`](https://github.com/xmusphlkg/pertussis-genomic-transmission-public/tree/v1.0.0).

## Start here

| Goal | Entry point |
|---|---|
| Read the two-figure result | [Figure 1](figures/letter/Figure_1_observation_structure.pdf) and [Figure 2](figures/letter/Figure_2_growth_robustness.pdf) |
| Trace each Letter claim to evidence | [`LETTER_EVIDENCE_MAP.md`](LETTER_EVIDENCE_MAP.md) |
| Inspect panel-level source data | [`figures/source_data/`](figures/source_data/) |
| Inspect frozen model inputs | [`data/model_inputs/`](data/model_inputs/) |
| Understand released files | [`DATA_DICTIONARY.md`](DATA_DICTIONARY.md) |
| Re-render the figures | [`scripts/figures/`](scripts/figures/) |
| Validate the public package | `python3 scripts/qa/validate_public_package.py` |

## Result in one paragraph

A uniformly quality-controlled 989-genome core-SNP tree identified one
monophyletic 288-tip population containing all 99 tips with a source-published
MT28 sublineage label. In the project-adjusted model, this MT28-associated
genomic lineage had a relative net-growth multiplier of 1.112 (95% credible
interval 1.078–1.156) and exceeded every comparator lineage. China-only and
Japan-only fits remained above the lineage reference; the sparse Australian
fit was uninformative. The signal persisted across project specifications,
country and dominant-project omissions, deterministic cohort-cap weighting,
and six model-input refits.

## Main figures

### Figure 1 — lineage definition and sampled-observation context

[![Figure 1](figures/letter/Figure_1_observation_structure.png)](figures/letter/Figure_1_observation_structure.pdf)

The phylogeny defines the MT28-associated genomic lineage. The accompanying
panels show reported cases, public-genome timing, deterministic cohort-cap
weighting, and Ct-dependent profile recovery. Sampled-tree shares are not
national prevalence estimates.

### Figure 2 — relative growth and robustness

[![Figure 2](figures/letter/Figure_2_growth_robustness.png)](figures/letter/Figure_2_growth_robustness.pdf)

The figure combines the shared estimate, direct pairwise comparisons,
country-only fits, country/project/input perturbations, conditional anchored
scenarios, and recovery diagnostics. Relative net-growth multipliers are model
summaries, not biological-fitness or transmissibility parameters.

Three compact supplementary figures are in [`figures/letter/`](figures/letter/):
genomic validation, posterior-predictive checks, and ancestry/input
sensitivity.

## Interpretation boundaries

This package supports a bounded conclusion: archive composition alone is
unlikely to explain the observed MT28-associated expansion signal in China and
Japan. It does **not** establish:

- national lineage prevalence;
- individual transmission links or counts of introductions;
- a causal contribution of the lineage to the notification rebound;
- biological fitness, transmissibility, or an attributable fraction;
- a molecular mechanism for expansion; or
- an identifiable import scale.

Project adjustment addresses measurable lineage enrichment among included
projects. It cannot correct specimen acquisition, sequencing failure, or public
archive inclusion without denominators.

## Repository design

The release has two evidence layers.

### Letter evidence layer

- `figures/letter/`: two main and three supplementary Letter figures
- `figures/source_data/`: compact panel-level source tables
- `data/model_inputs/`: frozen model inputs
- `data/model_sensitivity/`: frozen inputs for cohort, prior, threshold, timing,
  and root perturbations
- `results/lineage_nomenclature/`: MT28 label-to-clade crosswalk
- `results/model_main/` and `results/model_no_project/`: primary summaries
- `results/model_growth_robustness/`: country, project, and omission analyses
- `results/model_input_sensitivity_summary/`: six input-refit summaries
- `results/cohort_gate_sensitivity/`: pre-rescue versus final-cohort audit
- `results/molecular_characterisation/`: annotation coverage and limits

### Extended audit layer

`data/derived/`, `results/phylogeny/`, `results/phylogeography*/`,
`results/cgmlst/`, `results/lineages/`, `scripts/pipeline/`, and `tests/`
preserve the upstream data roles, tree construction, genomic validation, model
implementation, and sensitivity checks. These files support audit and reuse;
they are not additional Letter figures.

## Reproduction

From the repository root:

```bash
Rscript scripts/figures/render_letter_figures.R
Rscript scripts/figures/render_letter_supplementary_figures.R
python3 scripts/qa/validate_public_package.py
python3 -m pytest -q
```

The R scripts require the packages recorded in
[`provenance/R_SESSION_INFO.txt`](provenance/R_SESSION_INFO.txt). Large fitted
Stan objects are excluded; frozen inputs, Stan source, posterior summaries,
diagnostics, and checksums are provided instead. Consequently, this is a
results-audit and figure-reproduction package, not a one-command reconstruction
of every upstream bioinformatics step.

## Public-data and redistribution boundary

All source genome sequences and surveillance series were publicly available;
accessions and source URLs are retained in released tables. No controlled-access
or directly identifying participant data are included.

The Australian specimen-level sequencing-process table originated in
Supplementary Appendix 2 of Fong et al. (2026) and is not redistributed because
independent redistribution terms were not established. Only aggregate
calibration outputs used by the Letter are included. Public BioSample/run
accessions used in the tree remain in released tables, but exact Ct values,
source-specific aliases, specimen types, profile-success fields, and
article-supplied marker annotations are removed from per-record public files.
Consult the source article for the original records.

No repository-level licence has yet been assigned. Code and data reuse therefore
remain subject to the rights of their respective holders and the original
third-party sources.

## Citation and versions

Machine-readable metadata are in [`CITATION.cff`](CITATION.cff). See
[`CHANGELOG.md`](CHANGELOG.md) for the transition from the long-form `v1.0.0`
package to the Letter-first v2 design.
