# Reviewer-revision reproducibility companion

This companion was prepared on 9 October 2026 from the public repository at https://github.com/Sami6639/renewable-fossil-displacement-g20-reproducibility, commit e518a95082ba6f2b4df3b3eea433fb6752f35aac. The public repository was read and cloned; no remote files were changed or published. The `repo/` directory preserves the supplied repository snapshot, excluding its Git history in the delivered ZIP.

## Run

Python 3 with NumPy, pandas, and SciPy is required. The verified local versions were Python 3.12.14, NumPy 2.3.5, pandas 2.2.3, and SciPy 1.17.0. Execute `python audit_duration.py` from any working directory. Paths resolve relative to the script. No network access is needed. The original repository pipeline remains untouched; this companion does not require downloading the large raw OWID dataset.

Outputs are `audit_results.json`, `reconstructed_displacement_episodes.csv`, `matched_horizon_episode_results.csv`, and `nested_horizon_summary.csv`. `run_log.txt` records the validated execution. Run `python independent_validation.py` after the main script for an independent reconstruction, country-demeaned trend fit, and exact interval checks. All tabular rates use percentages; differences between nested rates use percentage points.

## Analysis protocol and scope

The additional analyses address the reviewer's request for substantive insight beyond counts and trends. The nested same-cohort protocol was specified before calculation in the revision workflow; this is not a prospectively registered study or a causal design. All requested nested comparisons and bounded sensitivities are retained, rather than selecting on statistical significance.

Start with the original sign-based displacement definition: annual renewable generation increases and fossil generation decreases. A displacement episode begins after a non-displacement or boundary year and continues across consecutive displacement years. No primary episode starts in 2000, so every episode's pre-onset baseline is observed directly in the supplied panel.

For each onset s, compare three nested outcomes on the same horizon-H cohort:

1. Joint continuity: every year s through s+H−1 is displacement (SD2 when H=2; SD3 when H=3).
2. Fossil-only continuity: fossil generation decreases in every year of the same window, whether or not renewables increase.
3. Endpoint retention: fossil output at s+H−1 is strictly below fossil output at s−1.

For H=2 use starts through 2023 (86 episodes); for H=3 use starts through 2022 (80 episodes). These are complete-follow-up cohorts, not an imputed survival outcome. Baselines are read directly from observed generation levels. Retention does not mean the initial decline was fully retained, that the path was monotonic, that renewables caused it, or that a transition was complete. It means output remains below the pre-onset level at the specified endpoint. The windows include H annual changes and end H−1 calendar years after the onset year.

The conditional retention rate restricts each same-horizon cohort to episodes that fail joint continuity. Paired gaps subtract nested binary indicators within each episode before aggregation and country resampling. Descriptive changes in fossil output are calculated relative to pre-onset fossil generation; medians and interquartile ranges avoid treating episodes of different size as comparable TWh units. Episode-level TWh sums are retained in the machine-readable audit only: overlapping episode windows and different baselines make them inappropriate as aggregate G20 reductions.

Uncertainty uses 20,000 country-cluster bootstrap draws with NumPy default_rng seed 20260901. Countries are alphabetically ordered, all 19 are eligible for resampling, and countries are drawn with replacement as intact collections of eligible episodes. Percentile intervals describe sensitivity to country composition and are not causal confidence statements. LOCO recomputes each pooled statistic after omitting each country. The finite population of 19 G20 country members and the small number of country clusters limit the inferential interpretation.

Sensitivity checks retain the same primary cohorts and (a) exclude 2020 episode starts or (b) require the endpoint reduction to exceed 0.10% or 0.25% of pre-onset total electricity generation. These endpoint thresholds do not redefine the primary displacement state or its original annual threshold sensitivity. Their denominators remain the original horizon-specific cohorts except for the explicit 2020 exclusion.

## Reproduction and important corrections

The companion independently reconstructs 475 country-years, 19 countries, 152 displacement years, 87 displacement episodes, 7 right-censored episodes, fixed-horizon SD2 37/86 and SD3 14/80, primary and threshold Kaplan–Meier point estimates, 80 demand-conditioned observations in 64 runs, original fixed-horizon country-FE trend estimates, and the exact transition counts. The repository's committed checksum list passes for its supplied files.

Leaving displacement is not equivalent to rising fossil generation. Among 145 directional transitions from displacement, 65 remain in displacement, 36 enter additive expansion, 33 enter fossil resurgence, and 11 enter joint contraction. Thus 55.2% leave the joint state, but 47.6% enter rising-fossil states and 7.6% continue fossil contraction with falling renewables.

Original SD2 trend inference uses normal-reference p=0.0450 with a CR1 country-cluster standard error. Keeping the same coefficient and standard error but using a t reference with G−1 degrees of freedom gives p=0.0603 (19 clusters). This is a sensitivity, not a replacement selected for its significance. Both are retained. SD3 is imprecise under either reference. Omitting starts in 2020 leaves 18 contributing countries for these trend models.

The public repository contains the earlier baseline code and outputs, but no sustained-displacement extension script. The manuscript's duration point estimates can nevertheless be reconstructed from the existing data. This companion's explicitly specified 20,000-draw bootstrap gives SD2 33.8–52.9% and SD3 10.1–26.9%, differing slightly from the manuscript's 33.7–53.0% and 10.0–27.0%. The original extension's precise random-draw ordering is unavailable, so its exact bootstrap quantiles have not been verified. Do not mix the new and original intervals without identifying their provenance.

This is a reproduction from the committed derived panel. The frozen raw OWID files and field-level codebook were not downloaded or independently validated in this revision. The original repository documents their checksums and lineage; reconstructing the complete source-to-panel processing still requires those matching raw files. A URL for a mutable current source is not an immutable archive of the frozen source.

## Files, licensing, and redistribution

The original repository's `LICENSE` covers its code and `LICENSES.md` preserves the upstream attribution/reuse conditions for the OWID-derived measurements. This companion does not change those conditions or claim ownership of third-party data. Original source/provider attribution remains in `repo/verification/source_manifest.md`. The companion is delivered privately for review; its creation does not publish the extension in the public repository.
