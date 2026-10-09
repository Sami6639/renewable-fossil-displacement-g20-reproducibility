# Renewable Addition and Fossil Displacement in G20 Power Systems

[![Reproducibility](https://img.shields.io/badge/reproducibility-verified-188977)](verification/verification_report.json)
[![Python](https://img.shields.io/badge/Python-3.12-3776AB)](environment.yml)
[![License: MIT](https://img.shields.io/badge/code%20license-MIT-yellow.svg)](LICENSE)

This repository contains the data-processing code, derived analytical panel,
country-year regime classifications, transition outputs, robustness diagnostics,
and verification information for:

> Küçükoğlu, S. *Renewable Addition and Fossil Displacement in G20 Power
> Systems: Regime Persistence and Transition Dynamics*.

The study follows all 19 country members of the G20 from 2000 to 2024 and
classifies annual joint changes in renewable and fossil electricity generation
into four observed regimes: fossil displacement, additive expansion, fossil
resurgence, and joint contraction. The archive is designed to reproduce every
reported country-year state, transition probability, temporal trend,
material-change sensitivity, country-cluster bootstrap result, leave-one-country-
out diagnostic, and secondary conditional model.

## Reproduce the study

### Windows and Anaconda

Double-click `RUN_REPRODUCIBILITY.bat`. The script creates or updates the Conda
environment and runs the complete workflow.

### Command line

```bash
conda env create -f environment.yml
conda run -n g20-power-repro python code/run_all.py
```

If the environment already exists:

```bash
conda env update -n g20-power-repro -f environment.yml --prune
conda run -n g20-power-repro python code/run_all.py
```

The first run retrieves the two public OWID source files and accepts them only
when both SHA-256 checksums match the frozen 1 September 2026 analytical source.
The verified raw files remain local and are excluded from Git; all derived data
and reported outputs are included in the repository.

## Repository contents

| Path | Contents |
| --- | --- |
| `data/derived/` | Balanced 2000–2024 analytical panel (475 country-years) |
| `classifications/` | Country-year regimes, transition pairs, regime spells, and secondary-model samples |
| `results/` | Frequencies, transition matrices, temporal trends, country profiles, spell summaries, and conditional estimates |
| `robustness/` | 20,000 country-cluster bootstrap replications and all LOCO diagnostics |
| `figures/` | Publication-quality PNG and editable SVG figures |
| `code/` | Numbered analysis scripts and the shared audited pipeline |
| `verification/` | Source manifest, output checksums, and manuscript-alignment report |

The complete workflow is:

1. `code/01_prepare_data.py`
2. `code/02_construct_regimes.py`
3. `code/03_transition_analysis.py`
4. `code/04_robustness.py`
5. `code/05_verify.py`

`code/run_all.py` executes all five stages in order.

## Frozen source and verification

The source data are the publicly available [Our World in Data Energy
dataset](https://github.com/owid/energy-data) and its codebook. The exact frozen
files are identified by these SHA-256 checksums:

```text
77b3db513f02f5fffb69fe02832907ce70b01d3906fc2c5dd40fa47e3ee7d0f3  owid-energy-data.csv
3cc9b7db0d921496e2988568ce3aee5ed41f50431dd234a0663b5f0a4b2e32bb  owid-energy-codebook.csv
```

The automated verification gate confirms, among other checks:

- 19 countries and 475 country-years;
- 13 exact-boundary observations under the primary definition;
- 442 consecutive directional transition pairs;
- primary counts of 178 additive, 152 displacement, 111 resurgence, and 21
  joint-contraction observations;
- the four principal transition probabilities reported in the manuscript; and
- the expected samples of 330, 156, and 145 observations for the temporal and
  secondary models.

See [`verification/source_manifest.md`](verification/source_manifest.md) for
source lineage and [`verification/verification_report.json`](verification/verification_report.json)
for the machine-readable audit.

## Scope and interpretation

The regimes describe observed joint movements in absolute generation. “Fossil
displacement” does not assert that renewable generation caused the fossil
decline. The secondary models are associational and do not identify policy,
technology, or demand mechanisms.

No restricted, confidential, or personal data are included.

## Citation

Please cite the manuscript and this repository. Machine-readable citation
metadata are available in [`CITATION.cff`](CITATION.cff).

## Licensing and attribution

The analysis code is released under the [MIT License](LICENSE). Source and
derived data retain the attribution and reuse conditions of their upstream
providers; see [`LICENSES.md`](LICENSES.md). The repository does not alter or
replace the authoritative OWID codebook.
