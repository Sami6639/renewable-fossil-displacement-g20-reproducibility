"""Reproducibility pipeline for the G20 renewable–fossil regime study.

The module intentionally uses only NumPy, pandas, and matplotlib so the full
analysis can be run in a small, auditable Conda environment.  Paths are
resolved relative to the repository root, never the current working directory.
"""

from __future__ import annotations

import hashlib
import json
import math
import urllib.request
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
RAW = DATA / "raw"
DERIVED = DATA / "derived"
CLASSIFICATIONS = ROOT / "classifications"
RESULTS = ROOT / "results"
ROBUSTNESS = ROOT / "robustness"
FIGURES = ROOT / "figures"
VERIFICATION = ROOT / "verification"

SOURCE_URL = "https://owid-public.owid.io/data/energy/owid-energy-data.csv"
CODEBOOK_URL = "https://raw.githubusercontent.com/owid/energy-data/master/owid-energy-codebook.csv"
SOURCE_SHA256 = "77b3db513f02f5fffb69fe02832907ce70b01d3906fc2c5dd40fa47e3ee7d0f3"
CODEBOOK_SHA256 = "3cc9b7db0d921496e2988568ce3aee5ed41f50431dd234a0663b5f0a4b2e32bb"
BOOTSTRAP_SEED = 20260901
BOOTSTRAP_REPLICATIONS = 20_000

COUNTRIES = {
    "ARG": "Argentina",
    "AUS": "Australia",
    "BRA": "Brazil",
    "CAN": "Canada",
    "CHN": "China",
    "FRA": "France",
    "DEU": "Germany",
    "IND": "India",
    "IDN": "Indonesia",
    "ITA": "Italy",
    "JPN": "Japan",
    "KOR": "South Korea",
    "MEX": "Mexico",
    "RUS": "Russia",
    "SAU": "Saudi Arabia",
    "ZAF": "South Africa",
    "TUR": "Türkiye",
    "GBR": "United Kingdom",
    "USA": "United States",
}

STATE_ORDER = [
    "Additive expansion",
    "Fossil displacement",
    "Fossil resurgence",
    "Joint contraction",
]
ALL_STATES = STATE_ORDER + ["Boundary"]
SPECIFICATIONS = {
    "primary": ("Primary sign-based", 0.0),
    "threshold_010": ("0.10% material-change threshold", 0.001),
    "threshold_025": ("0.25% material-change threshold", 0.0025),
}

KEEP_COLUMNS = [
    "country",
    "year",
    "iso_code",
    "renewables_electricity",
    "fossil_electricity",
    "electricity_generation",
    "electricity_demand",
    "fossil_share_elec",
    "net_elec_imports_share_demand",
    "wind_electricity",
    "solar_electricity",
    "nuclear_electricity",
]


def ensure_directories() -> None:
    for path in [RAW, DERIVED, CLASSIFICATIONS, RESULTS, ROBUSTNESS, FIGURES, VERIFICATION]:
        path.mkdir(parents=True, exist_ok=True)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def download_verified(url: str, destination: Path, expected_sha256: str) -> None:
    if destination.exists() and sha256(destination) == expected_sha256:
        return
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix(destination.suffix + ".download")
    with urllib.request.urlopen(url, timeout=120) as response, temporary.open("wb") as out:
        while True:
            block = response.read(1024 * 1024)
            if not block:
                break
            out.write(block)
    observed = sha256(temporary)
    if observed != expected_sha256:
        temporary.unlink(missing_ok=True)
        raise RuntimeError(
            f"Checksum mismatch for {url}. Expected {expected_sha256}, observed {observed}. "
            "The public source may have been updated; use the frozen source identified in verification/source_manifest.md."
        )
    temporary.replace(destination)


def prepare_data() -> pd.DataFrame:
    """Download, verify, filter, and derive the balanced 2000–2024 panel."""
    ensure_directories()
    source = RAW / "owid-energy-data.csv"
    codebook = RAW / "owid-energy-codebook.csv"
    download_verified(SOURCE_URL, source, SOURCE_SHA256)
    download_verified(CODEBOOK_URL, codebook, CODEBOOK_SHA256)

    full = pd.read_csv(source, usecols=KEEP_COLUMNS)
    lagged = full.loc[
        full["iso_code"].isin(COUNTRIES) & full["year"].between(1999, 2024), KEEP_COLUMNS
    ].copy()
    lagged["country_source"] = lagged["country"]
    lagged["country"] = lagged["iso_code"].map(COUNTRIES)
    lagged = lagged.sort_values(["iso_code", "year"]).reset_index(drop=True)

    required = ["renewables_electricity", "fossil_electricity", "electricity_generation"]
    if lagged.groupby("iso_code")["year"].nunique().ne(26).any():
        raise RuntimeError("The 1999–2024 lag-buffer panel is not balanced for all 19 countries.")
    if lagged[required].isna().any().any():
        raise RuntimeError("Primary electricity variables contain missing observations.")

    groups = lagged.groupby("iso_code", sort=False)
    lagged["renewables_change_twh"] = groups["renewables_electricity"].diff()
    lagged["fossil_change_twh"] = groups["fossil_electricity"].diff()
    lagged["prior_generation_twh"] = groups["electricity_generation"].shift(1)
    lagged["electricity_demand_growth_pct"] = groups["electricity_demand"].pct_change(fill_method=None) * 100
    lagged["wind_solar_share_elec"] = (
        (lagged["wind_electricity"].fillna(0) + lagged["solar_electricity"].fillna(0))
        / lagged["electricity_generation"]
        * 100
    )
    lagged["nuclear_share_elec"] = (
        lagged["nuclear_electricity"].fillna(0) / lagged["electricity_generation"] * 100
    )

    panel = lagged.loc[lagged["year"].between(2000, 2024)].copy()
    if len(panel) != 475:
        raise RuntimeError(f"Expected 475 main-period observations, found {len(panel)}.")
    panel.to_csv(DERIVED / "analytical_panel_2000_2024.csv", index=False, float_format="%.10g")
    return panel


def _classify(d_renewable: float, d_fossil: float, prior_generation: float, threshold: float) -> str:
    if pd.isna(d_renewable) or pd.isna(d_fossil):
        return "Boundary"
    if threshold > 0 and (
        abs(d_renewable) < threshold * prior_generation
        or abs(d_fossil) < threshold * prior_generation
    ):
        return "Boundary"
    if d_renewable == 0 or d_fossil == 0:
        return "Boundary"
    if d_renewable > 0 and d_fossil > 0:
        return "Additive expansion"
    if d_renewable > 0 and d_fossil < 0:
        return "Fossil displacement"
    if d_renewable < 0 and d_fossil > 0:
        return "Fossil resurgence"
    return "Joint contraction"


def construct_regimes(panel: pd.DataFrame | None = None) -> pd.DataFrame:
    if panel is None:
        panel = pd.read_csv(DERIVED / "analytical_panel_2000_2024.csv")
    regimes = panel.copy()
    for key, (_, threshold) in SPECIFICATIONS.items():
        regimes[f"state_{key}"] = [
            _classify(dr, df, g, threshold)
            for dr, df, g in zip(
                regimes["renewables_change_twh"],
                regimes["fossil_change_twh"],
                regimes["prior_generation_twh"],
            )
        ]
    regimes.to_csv(CLASSIFICATIONS / "country_year_regimes.csv", index=False, float_format="%.10g")

    pairs = []
    spells = []
    for spec_key, (spec_label, _) in SPECIFICATIONS.items():
        state_col = f"state_{spec_key}"
        for iso, group in regimes.groupby("iso_code", sort=False):
            group = group.sort_values("year").reset_index(drop=True)
            for idx in range(len(group) - 1):
                current, following = group.iloc[idx], group.iloc[idx + 1]
                if int(following["year"]) != int(current["year"]) + 1:
                    continue
                pairs.append(
                    {
                        "specification": spec_label,
                        "specification_key": spec_key,
                        "country": current["country"],
                        "iso_code": iso,
                        "year_t": int(current["year"]),
                        "year_t1": int(following["year"]),
                        "state_t": current[state_col],
                        "state_t1": following[state_col],
                        "directional_pair": current[state_col] in STATE_ORDER and following[state_col] in STATE_ORDER,
                    }
                )

            run_state = None
            run_start = None
            run_length = 0
            prior_year = None
            for _, row in group.iterrows():
                state = row[state_col]
                year = int(row["year"])
                continues = state in STATE_ORDER and state == run_state and prior_year is not None and year == prior_year + 1
                if continues:
                    run_length += 1
                else:
                    if run_state in STATE_ORDER:
                        spells.append(
                            {
                                "specification": spec_label,
                                "specification_key": spec_key,
                                "country": group.iloc[0]["country"],
                                "iso_code": iso,
                                "state": run_state,
                                "start_year": run_start,
                                "end_year": prior_year,
                                "duration_years": run_length,
                            }
                        )
                    run_state = state if state in STATE_ORDER else None
                    run_start = year if state in STATE_ORDER else None
                    run_length = 1 if state in STATE_ORDER else 0
                prior_year = year
            if run_state in STATE_ORDER:
                spells.append(
                    {
                        "specification": spec_label,
                        "specification_key": spec_key,
                        "country": group.iloc[0]["country"],
                        "iso_code": iso,
                        "state": run_state,
                        "start_year": run_start,
                        "end_year": prior_year,
                        "duration_years": run_length,
                    }
                )

    pd.DataFrame(pairs).to_csv(CLASSIFICATIONS / "transition_pairs.csv", index=False)
    pd.DataFrame(spells).to_csv(CLASSIFICATIONS / "regime_spells.csv", index=False)
    return regimes


def transition_matrix(pairs: pd.DataFrame, spec_key: str) -> pd.DataFrame:
    directional = pairs.loc[(pairs["specification_key"] == spec_key) & pairs["directional_pair"]].copy()
    counts = pd.crosstab(directional["state_t"], directional["state_t1"]).reindex(
        index=STATE_ORDER, columns=STATE_ORDER, fill_value=0
    )
    rows = []
    label = SPECIFICATIONS[spec_key][0]
    for origin in STATE_ORDER:
        denominator = int(counts.loc[origin].sum())
        for destination in STATE_ORDER:
            number = int(counts.loc[origin, destination])
            rows.append(
                {
                    "specification": label,
                    "specification_key": spec_key,
                    "current_regime": origin,
                    "next_year_regime": destination,
                    "n": number,
                    "origin_state_n": denominator,
                    "probability": number / denominator if denominator else np.nan,
                    "probability_pct": 100 * number / denominator if denominator else np.nan,
                }
            )
    return pd.DataFrame(rows)


def _normal_cdf(value: float) -> float:
    return 0.5 * (1.0 + math.erf(value / math.sqrt(2.0)))


def clustered_ols(data: pd.DataFrame, outcome: str, predictors: list[str], cluster: str) -> dict:
    """OLS with country fixed effects and finite-sample-corrected cluster covariance."""
    model = data[[outcome, cluster, *predictors]].dropna().copy()
    dummies = pd.get_dummies(model[cluster], prefix="country", drop_first=True, dtype=float)
    x_named = pd.concat(
        [pd.Series(1.0, index=model.index, name="intercept"), model[predictors].astype(float), dummies],
        axis=1,
    )
    x = x_named.to_numpy(float)
    y = model[outcome].to_numpy(float)
    xtx_inv = np.linalg.pinv(x.T @ x)
    beta = xtx_inv @ x.T @ y
    residual = y - x @ beta
    meat = np.zeros((x.shape[1], x.shape[1]))
    groups = model[cluster].to_numpy()
    unique_groups = pd.unique(groups)
    for group in unique_groups:
        mask = groups == group
        score = x[mask].T @ residual[mask]
        meat += np.outer(score, score)
    n, k, g = len(y), np.linalg.matrix_rank(x), len(unique_groups)
    correction = (g / (g - 1)) * ((n - 1) / (n - k))
    covariance = correction * xtx_inv @ meat @ xtx_inv
    standard_errors = np.sqrt(np.maximum(np.diag(covariance), 0))
    estimates = dict(zip(x_named.columns, beta))
    ses = dict(zip(x_named.columns, standard_errors))
    output = {}
    for predictor in predictors:
        estimate = float(estimates[predictor])
        se = float(ses[predictor])
        z = estimate / se if se else np.nan
        p = 2 * (1 - _normal_cdf(abs(z))) if np.isfinite(z) else np.nan
        output[predictor] = {
            "estimate": estimate,
            "clustered_se": se,
            "p_value": p,
            "ci_low": estimate - 1.96 * se,
            "ci_high": estimate + 1.96 * se,
        }
    return {"n": n, "events": int(y.sum()), "coefficients": output}


def _trend_model(regimes: pd.DataFrame, state_col: str) -> dict:
    sample = regimes.loc[regimes[state_col].isin(["Additive expansion", "Fossil displacement"])].copy()
    sample["displacement"] = (sample[state_col] == "Fossil displacement").astype(int)
    sample["year_trend"] = sample["year"] - 2000
    result = clustered_ols(sample, "displacement", ["year_trend"], "iso_code")
    result["sample"] = sample
    return result


def _secondary_sample(regimes: pd.DataFrame, origin: str, destination: str) -> pd.DataFrame:
    predictor_cols = [
        "electricity_demand_growth_pct",
        "fossil_share_elec",
        "net_elec_imports_share_demand",
        "wind_solar_share_elec",
        "nuclear_share_elec",
    ]
    frames = []
    for _, group in regimes.groupby("iso_code", sort=False):
        group = group.sort_values("year").copy()
        group["next_state"] = group["state_primary"].shift(-1)
        group["next_year"] = group["year"].shift(-1)
        sample = group.loc[
            (group["state_primary"] == origin)
            & (group["next_year"] == group["year"] + 1)
            & group["next_state"].isin(STATE_ORDER)
        ].copy()
        sample["outcome"] = (sample["next_state"] == destination).astype(int)
        frames.append(sample)
    sample = pd.concat(frames, ignore_index=True)
    sample = sample.dropna(subset=predictor_cols).copy()
    for col in predictor_cols:
        # Sample-standard-deviation scaling matches the prespecified models
        # reported in Table S3 of the manuscript.
        sample[f"z_{col}"] = (sample[col] - sample[col].mean()) / sample[col].std(ddof=1)
    sample["year_trend"] = sample["year"] - 2000
    return sample


def analyze_results(regimes: pd.DataFrame | None = None) -> None:
    if regimes is None:
        regimes = pd.read_csv(CLASSIFICATIONS / "country_year_regimes.csv")
    pairs = pd.read_csv(CLASSIFICATIONS / "transition_pairs.csv")
    pairs["directional_pair"] = pairs["directional_pair"].astype(str).str.lower().eq("true")
    spells = pd.read_csv(CLASSIFICATIONS / "regime_spells.csv")

    frequency_rows = []
    matrices = []
    for key, (label, _) in SPECIFICATIONS.items():
        state_col = f"state_{key}"
        counts = regimes[state_col].value_counts().reindex(ALL_STATES, fill_value=0)
        for state, count in counts.items():
            frequency_rows.append(
                {
                    "specification": label,
                    "specification_key": key,
                    "regime": state,
                    "n": int(count),
                    "share_all": count / len(regimes),
                    "share_all_pct": 100 * count / len(regimes),
                }
            )
        matrices.append(transition_matrix(pairs, key))
    frequencies = pd.DataFrame(frequency_rows)
    full_matrix = pd.concat(matrices, ignore_index=True)
    frequencies.to_csv(RESULTS / "regime_frequencies.csv", index=False, float_format="%.10g")
    full_matrix.to_csv(RESULTS / "transition_matrix.csv", index=False, float_format="%.10g")

    primary = regimes.copy()
    primary["period"] = pd.cut(
        primary["year"],
        bins=[1999, 2004, 2009, 2014, 2019, 2024],
        labels=["2000–2004", "2005–2009", "2010–2014", "2015–2019", "2020–2024"],
    )
    temporal = (
        primary.groupby(["period", "state_primary"], observed=False)
        .size()
        .rename("n")
        .reset_index()
        .rename(columns={"state_primary": "regime"})
    )
    temporal["period_n"] = temporal.groupby("period", observed=False)["n"].transform("sum")
    temporal["share"] = temporal["n"] / temporal["period_n"]
    temporal["share_pct"] = 100 * temporal["share"]
    temporal.to_csv(RESULTS / "temporal_regime_composition.csv", index=False, float_format="%.10g")

    country_profiles = (
        primary.groupby(["country", "iso_code", "state_primary"]).size().rename("n").reset_index()
    )
    country_profiles["country_n"] = country_profiles.groupby("iso_code")["n"].transform("sum")
    country_profiles["share_pct"] = 100 * country_profiles["n"] / country_profiles["country_n"]
    country_profiles.to_csv(RESULTS / "country_regime_profiles.csv", index=False, float_format="%.10g")

    spell_summary = (
        spells.loc[spells["specification_key"] == "primary"]
        .groupby("state")["duration_years"]
        .agg(spell_count="size", mean_years="mean", median_years="median", maximum_years="max")
        .reset_index()
    )
    spell_summary.to_csv(RESULTS / "regime_spell_summary.csv", index=False, float_format="%.10g")

    trend_rows = []
    for key, (label, _) in SPECIFICATIONS.items():
        trend = _trend_model(regimes, f"state_{key}")
        coef = trend["coefficients"]["year_trend"]
        trend_rows.append(
            {
                "specification": label,
                "specification_key": key,
                "n": trend["n"],
                "annual_change": coef["estimate"],
                "annual_change_pp": 100 * coef["estimate"],
                "clustered_se": coef["clustered_se"],
                "clustered_se_pp": 100 * coef["clustered_se"],
                "p_value": coef["p_value"],
                "ci_low": coef["ci_low"],
                "ci_high": coef["ci_high"],
                "ci_low_pp": 100 * coef["ci_low"],
                "ci_high_pp": 100 * coef["ci_high"],
            }
        )
    pd.DataFrame(trend_rows).to_csv(RESULTS / "temporal_trend_models.csv", index=False, float_format="%.10g")

    predictor_cols = [
        "electricity_demand_growth_pct",
        "fossil_share_elec",
        "net_elec_imports_share_demand",
        "wind_solar_share_elec",
        "nuclear_share_elec",
    ]
    z_predictors = [f"z_{col}" for col in predictor_cols]
    model_definitions = [
        ("Additive expansion → fossil displacement", "Additive expansion", "Fossil displacement"),
        ("Fossil displacement → fossil displacement", "Fossil displacement", "Fossil displacement"),
    ]
    secondary_rows = []
    for outcome_label, origin, destination in model_definitions:
        sample = _secondary_sample(regimes, origin, destination)
        sample.to_csv(
            CLASSIFICATIONS / ("secondary_sample_addition_to_displacement.csv" if origin.startswith("Additive") else "secondary_sample_displacement_persistence.csv"),
            index=False,
            float_format="%.10g",
        )
        model = clustered_ols(sample, "outcome", ["year_trend", *z_predictors], "iso_code")
        for source_name, z_name in zip(predictor_cols, z_predictors):
            coef = model["coefficients"][z_name]
            secondary_rows.append(
                {
                    "outcome": outcome_label,
                    "predictor": source_name,
                    "n": model["n"],
                    "events": model["events"],
                    **coef,
                }
            )
    pd.DataFrame(secondary_rows).to_csv(RESULTS / "secondary_models.csv", index=False, float_format="%.10g")

    _make_figures(temporal, full_matrix, country_profiles, pd.DataFrame(secondary_rows))


def _bootstrap_transition_diagnostics(pairs: pd.DataFrame) -> pd.DataFrame:
    directional = pairs.loc[(pairs["specification_key"] == "primary") & pairs["directional_pair"]].copy()
    countries = list(COUNTRIES)
    count_arrays = {}
    for origin, destination, name in [
        ("Additive expansion", "Fossil displacement", "p_AD"),
        ("Additive expansion", "Additive expansion", "p_AA"),
        ("Fossil displacement", "Fossil displacement", "p_DD"),
        ("Fossil displacement", "Additive expansion", "p_DA"),
    ]:
        numerator = np.array(
            [((directional["iso_code"] == iso) & (directional["state_t"] == origin) & (directional["state_t1"] == destination)).sum() for iso in countries],
            dtype=float,
        )
        denominator = np.array(
            [((directional["iso_code"] == iso) & (directional["state_t"] == origin)).sum() for iso in countries],
            dtype=float,
        )
        count_arrays[name] = (numerator, denominator)

    rng = np.random.RandomState(BOOTSTRAP_SEED)
    draws = rng.multinomial(len(countries), [1 / len(countries)] * len(countries), size=BOOTSTRAP_REPLICATIONS)
    out = {"replication": np.arange(1, BOOTSTRAP_REPLICATIONS + 1)}
    for name, (numerator, denominator) in count_arrays.items():
        out[name] = (draws @ numerator) / (draws @ denominator)
    out["p_DD_minus_p_AA"] = out["p_DD"] - out["p_AA"]
    return pd.DataFrame(out)


def robustness_analysis(regimes: pd.DataFrame | None = None) -> None:
    if regimes is None:
        regimes = pd.read_csv(CLASSIFICATIONS / "country_year_regimes.csv")
    pairs = pd.read_csv(CLASSIFICATIONS / "transition_pairs.csv")
    pairs["directional_pair"] = pairs["directional_pair"].astype(str).str.lower().eq("true")

    bootstrap = _bootstrap_transition_diagnostics(pairs)
    bootstrap.to_csv(ROBUSTNESS / "bootstrap_20000.csv", index=False, float_format="%.10g")
    interval_rows = []
    for metric in ["p_AD", "p_AA", "p_DD", "p_DA", "p_DD_minus_p_AA"]:
        values = bootstrap[metric]
        interval_rows.append(
            {
                "metric": metric,
                "replications": BOOTSTRAP_REPLICATIONS,
                "seed": BOOTSTRAP_SEED,
                "ci_low": values.quantile(0.025),
                "ci_high": values.quantile(0.975),
                "ci_low_pct": 100 * values.quantile(0.025),
                "ci_high_pct": 100 * values.quantile(0.975),
            }
        )
    pd.DataFrame(interval_rows).to_csv(ROBUSTNESS / "bootstrap_intervals.csv", index=False, float_format="%.10g")

    loco_probability_rows = []
    directional = pairs.loc[(pairs["specification_key"] == "primary") & pairs["directional_pair"]].copy()
    for omitted_iso, omitted_country in COUNTRIES.items():
        sample = directional.loc[directional["iso_code"] != omitted_iso]
        for origin, destination, metric in [
            ("Additive expansion", "Fossil displacement", "p_AD"),
            ("Additive expansion", "Additive expansion", "p_AA"),
            ("Fossil displacement", "Fossil displacement", "p_DD"),
            ("Fossil displacement", "Additive expansion", "p_DA"),
        ]:
            denominator = (sample["state_t"] == origin).sum()
            numerator = ((sample["state_t"] == origin) & (sample["state_t1"] == destination)).sum()
            loco_probability_rows.append(
                {
                    "omitted_country": omitted_country,
                    "omitted_iso_code": omitted_iso,
                    "metric": metric,
                    "numerator": int(numerator),
                    "denominator": int(denominator),
                    "probability": numerator / denominator,
                    "probability_pct": 100 * numerator / denominator,
                }
            )
    pd.DataFrame(loco_probability_rows).to_csv(ROBUSTNESS / "leave_one_country_out_probabilities.csv", index=False, float_format="%.10g")

    trend_loco = []
    for omitted_iso, omitted_country in COUNTRIES.items():
        sample = regimes.loc[regimes["iso_code"] != omitted_iso].copy()
        trend = _trend_model(sample, "state_primary")
        coef = trend["coefficients"]["year_trend"]
        trend_loco.append(
            {
                "omitted_country": omitted_country,
                "omitted_iso_code": omitted_iso,
                "n": trend["n"],
                "annual_change": coef["estimate"],
                "annual_change_pp": 100 * coef["estimate"],
                "clustered_se": coef["clustered_se"],
                "p_value": coef["p_value"],
            }
        )
    pd.DataFrame(trend_loco).to_csv(ROBUSTNESS / "leave_one_country_out_temporal_trend.csv", index=False, float_format="%.10g")

    predictor_cols = [
        "electricity_demand_growth_pct",
        "fossil_share_elec",
        "net_elec_imports_share_demand",
        "wind_solar_share_elec",
        "nuclear_share_elec",
    ]
    z_predictors = [f"z_{col}" for col in predictor_cols]
    conditional_loco = []
    for outcome_label, origin, destination in [
        ("Additive expansion → fossil displacement", "Additive expansion", "Fossil displacement"),
        ("Fossil displacement → fossil displacement", "Fossil displacement", "Fossil displacement"),
    ]:
        for omitted_iso, omitted_country in COUNTRIES.items():
            reduced = regimes.loc[regimes["iso_code"] != omitted_iso]
            sample = _secondary_sample(reduced, origin, destination)
            model = clustered_ols(sample, "outcome", ["year_trend", *z_predictors], "iso_code")
            for source_name, z_name in zip(predictor_cols, z_predictors):
                coef = model["coefficients"][z_name]
                conditional_loco.append(
                    {
                        "outcome": outcome_label,
                        "omitted_country": omitted_country,
                        "omitted_iso_code": omitted_iso,
                        "predictor": source_name,
                        "n": model["n"],
                        **coef,
                    }
                )
    pd.DataFrame(conditional_loco).to_csv(ROBUSTNESS / "leave_one_country_out_secondary_models.csv", index=False, float_format="%.10g")


def _make_figures(temporal: pd.DataFrame, matrix: pd.DataFrame, country_profiles: pd.DataFrame, secondary: pd.DataFrame) -> None:
    colors = {
        "Additive expansion": "#D98E04",
        "Fossil displacement": "#188977",
        "Fossil resurgence": "#B64A55",
        "Joint contraction": "#5B6F8F",
        "Boundary": "#B8B8B8",
    }
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.titlesize": 12})

    pivot = temporal.pivot(index="period", columns="regime", values="share_pct").fillna(0).reindex(columns=ALL_STATES)
    ax = pivot.plot(kind="bar", stacked=True, figsize=(9, 5.3), color=[colors[x] for x in ALL_STATES], width=0.72)
    ax.set_xlabel("")
    ax.set_ylabel("Share of country-years (%)")
    ax.set_ylim(0, 100)
    ax.tick_params(axis="x", rotation=0)
    ax.legend(title="Regime", bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    plt.tight_layout()
    _save_figure("figure_1_regime_evolution")

    primary = matrix.loc[matrix["specification_key"] == "primary"].pivot(
        index="current_regime", columns="next_year_regime", values="probability_pct"
    ).reindex(index=STATE_ORDER, columns=STATE_ORDER)
    fig, ax = plt.subplots(figsize=(7.2, 5.6))
    image = ax.imshow(primary, cmap="YlGnBu", vmin=0, vmax=60)
    for row in range(len(STATE_ORDER)):
        for col in range(len(STATE_ORDER)):
            ax.text(col, row, f"{primary.iloc[row, col]:.1f}%", ha="center", va="center", color="black")
    ax.set_xticks(range(4), ["Addition", "Displacement", "Resurgence", "Contraction"], rotation=25, ha="right")
    ax.set_yticks(range(4), ["Addition", "Displacement", "Resurgence", "Contraction"])
    ax.set_xlabel("Next-year regime")
    ax.set_ylabel("Current regime")
    fig.colorbar(image, ax=ax, label="Transition probability (%)", shrink=0.84)
    plt.tight_layout()
    _save_figure("figure_2_transition_matrix")

    selected = country_profiles.loc[country_profiles["state_primary"].isin(["Fossil displacement", "Additive expansion"])]
    country_pivot = selected.pivot(index="country", columns="state_primary", values="share_pct").fillna(0)
    country_pivot = country_pivot.sort_values("Fossil displacement")
    ax = country_pivot[["Fossil displacement", "Additive expansion"]].plot(
        kind="barh", figsize=(8.5, 7.2), color=[colors["Fossil displacement"], colors["Additive expansion"]]
    )
    ax.set_xlabel("Share of annual observations (%)")
    ax.set_ylabel("")
    ax.legend(frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    plt.tight_layout()
    _save_figure("figure_3_country_profiles")

    threshold_rows = matrix.loc[
        matrix["current_regime"].isin(["Additive expansion", "Fossil displacement"])
        & matrix["next_year_regime"].isin(["Fossil displacement"])
    ].copy()
    threshold_rows["metric"] = np.where(
        threshold_rows["current_regime"] == "Additive expansion", "Addition → displacement", "Displacement persistence"
    )
    threshold_pivot = threshold_rows.pivot(index="specification", columns="metric", values="probability_pct").reindex(
        [v[0] for v in SPECIFICATIONS.values()]
    )
    ax = threshold_pivot.plot(kind="bar", figsize=(8.2, 5), color=[colors["Fossil displacement"], "#4B6FAE"])
    ax.set_xlabel("")
    ax.set_ylabel("Probability (%)")
    ax.set_xticklabels(["Primary", "0.10%", "0.25%"], rotation=0)
    ax.legend(frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    plt.tight_layout()
    _save_figure("figure_4_threshold_robustness")

    if not secondary.empty:
        fig, axes = plt.subplots(1, 2, figsize=(10, 5), sharey=True)
        labels = {
            "electricity_demand_growth_pct": "Demand growth",
            "fossil_share_elec": "Fossil share",
            "net_elec_imports_share_demand": "Net-import share",
            "wind_solar_share_elec": "Wind + solar share",
            "nuclear_share_elec": "Nuclear share",
        }
        for ax, (outcome, group) in zip(axes, secondary.groupby("outcome", sort=False)):
            group = group.copy()
            y = np.arange(len(group))
            ax.errorbar(
                group["estimate"], y,
                xerr=[group["estimate"] - group["ci_low"], group["ci_high"] - group["estimate"]],
                fmt="o", color="#235789", capsize=3,
            )
            ax.axvline(0, color="#777777", linewidth=1)
            ax.set_yticks(y, [labels[x] for x in group["predictor"]])
            ax.set_title(outcome.replace(" → ", "\n→ "))
            ax.set_xlabel("Coefficient (95% CI)")
            ax.spines[["top", "right"]].set_visible(False)
        plt.tight_layout()
        _save_figure("figure_s1_secondary_coefficients")


def _save_figure(stem: str) -> None:
    plt.savefig(FIGURES / f"{stem}.png", dpi=300, bbox_inches="tight")
    plt.savefig(FIGURES / f"{stem}.svg", bbox_inches="tight")
    plt.close()


def verify_outputs() -> dict:
    """Run manuscript-alignment checks and write a SHA-256 manifest."""
    regimes = pd.read_csv(CLASSIFICATIONS / "country_year_regimes.csv")
    pairs = pd.read_csv(CLASSIFICATIONS / "transition_pairs.csv")
    pairs["directional_pair"] = pairs["directional_pair"].astype(str).str.lower().eq("true")
    matrix = pd.read_csv(RESULTS / "transition_matrix.csv")
    trends = pd.read_csv(RESULTS / "temporal_trend_models.csv")
    secondary = pd.read_csv(RESULTS / "secondary_models.csv")

    checks = {
        "panel_rows_475": len(regimes) == 475,
        "countries_19": regimes["iso_code"].nunique() == 19,
        "primary_boundary_13": (regimes["state_primary"] == "Boundary").sum() == 13,
        "threshold_010_boundary_56": (regimes["state_threshold_010"] == "Boundary").sum() == 56,
        "threshold_025_boundary_109": (regimes["state_threshold_025"] == "Boundary").sum() == 109,
        "primary_directional_pairs_442": ((pairs["specification_key"] == "primary") & pairs["directional_pair"]).sum() == 442,
        "primary_addition_178": (regimes["state_primary"] == "Additive expansion").sum() == 178,
        "primary_displacement_152": (regimes["state_primary"] == "Fossil displacement").sum() == 152,
        "primary_resurgence_111": (regimes["state_primary"] == "Fossil resurgence").sum() == 111,
        "primary_contraction_21": (regimes["state_primary"] == "Joint contraction").sum() == 21,
    }
    key = matrix.loc[matrix["specification_key"] == "primary"].set_index(["current_regime", "next_year_regime"])
    checks.update(
        {
            "p_AD_26_5pct": abs(key.loc[("Additive expansion", "Fossil displacement"), "probability_pct"] - 26.5060) < 0.01,
            "p_AA_48_8pct": abs(key.loc[("Additive expansion", "Additive expansion"), "probability_pct"] - 48.7952) < 0.01,
            "p_DD_44_8pct": abs(key.loc[("Fossil displacement", "Fossil displacement"), "probability_pct"] - 44.8276) < 0.01,
            "p_DA_24_8pct": abs(key.loc[("Fossil displacement", "Additive expansion"), "probability_pct"] - 24.8276) < 0.01,
            "trend_n_330": int(trends.loc[trends["specification_key"] == "primary", "n"].iloc[0]) == 330,
            "secondary_n_156_145": sorted(secondary.groupby("outcome")["n"].first().astype(int).tolist()) == [145, 156],
        }
    )
    checks = {name: bool(passed) for name, passed in checks.items()}
    report = {
        "source_sha256": sha256(RAW / "owid-energy-data.csv"),
        "codebook_sha256": sha256(RAW / "owid-energy-codebook.csv"),
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }
    (VERIFICATION / "verification_report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    if not report["all_checks_pass"]:
        failed = [name for name, passed in checks.items() if not passed]
        raise RuntimeError(f"Manuscript-alignment verification failed: {failed}")

    files = []
    for path in ROOT.rglob("*"):
        relative = path.relative_to(ROOT)
        if not path.is_file():
            continue
        if ".git" in relative.parts or "__pycache__" in relative.parts:
            continue
        if relative.as_posix() in {
            "data/raw/owid-energy-data.csv",
            "data/raw/owid-energy-codebook.csv",
            "verification/SHA256SUMS.txt",
        }:
            continue
        files.append(path)
    lines = [f"{sha256(path)}  {path.relative_to(ROOT).as_posix()}" for path in sorted(files)]
    (VERIFICATION / "SHA256SUMS.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    return report


def run_all() -> dict:
    panel = prepare_data()
    regimes = construct_regimes(panel)
    analyze_results(regimes)
    robustness_analysis(regimes)
    return verify_outputs()


if __name__ == "__main__":
    result = run_all()
    print(json.dumps(result, indent=2))
