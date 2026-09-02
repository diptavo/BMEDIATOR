from __future__ import annotations

import csv
import json
import math
import os
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np


SCENARIOS = ("M0", "M1", "M2", "M3", "M4", "M5")
SCENARIO_INDEX = {name: idx for idx, name in enumerate(SCENARIOS)}
ALLELES = np.array([["A", "G"], ["C", "T"], ["G", "A"], ["T", "C"]], dtype=object)


@dataclass
class BenchmarkFiles:
    rf: Path
    pqtl: Path
    cancer: Path
    protein_info: Path
    truth: Path


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def load_config(path: Path) -> Dict:
    return json.loads(path.read_text())


def benchmark_scenarios(config: Dict, benchmark: str) -> Tuple[str, ...]:
    bench_cfg = config.get(benchmark, {})
    scenarios = bench_cfg.get("scenarios")
    if not scenarios:
        return SCENARIOS
    return tuple(str(s) for s in scenarios)


def scenario_alias_map(config: Dict, benchmark: str) -> Dict[str, str]:
    bench_cfg = config.get(benchmark, {})
    mapping = {str(k): str(v) for k, v in bench_cfg.get("scenario_aliases", {}).items()}
    return mapping


def effective_scenario(name: str, alias_map: Dict[str, str]) -> str:
    return alias_map.get(name, name)


def norm_sf(z: np.ndarray) -> np.ndarray:
    return 0.5 * np.vectorize(math.erfc)(np.abs(z) / math.sqrt(2.0))


def choose_alleles(rng: np.random.Generator, n: int) -> np.ndarray:
    idx = rng.integers(0, len(ALLELES), size=n)
    return ALLELES[idx]


def write_table(path: Path, rows: Sequence[Dict], fieldnames: Sequence[str], delimiter: str) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter=delimiter)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def read_table(path: Path, delimiter: str) -> List[Dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def scenario_parameters(cell: Dict, scenario: str) -> Dict[str, float]:
    params = {
        "beta1": float(cell["beta1"]) if scenario in {"M1", "M2", "M5"} else 0.0,
        "beta2": float(cell["beta2"]) if scenario in {"M1", "M4"} else 0.0,
        "beta3": float(cell.get("beta3", 0.0)) if scenario in {"M1", "M3", "M5"} else 0.0,
        "delta_sd": float(cell["delta_sd"]),
        "psi_sd": float(cell["psi_sd"]),
        "phi_sd": float(cell["phi_sd"]),
        "rho_m5": float(cell.get("rho_m5", cell.get("rho_m3", 0.0))),
        "direction_flip_rate": float(cell.get("direction_flip_rate", 0.0)),
    }
    overrides = cell.get("scenario_overrides", {}).get(scenario, {})
    for key, value in overrides.items():
        params[key] = float(value)
    return params


def _bounded_corr(value: float) -> float:
    return max(-0.99, min(0.99, float(value)))


def sampling_error_correlations(cell: Dict) -> Tuple[float, float, float]:
    spec = cell.get("sample_overlap", cell.get("sampling_error_correlation", {}))
    if isinstance(spec, (int, float)):
        corr = _bounded_corr(float(spec))
        return corr, corr, corr
    if not isinstance(spec, dict):
        return 0.0, 0.0, 0.0
    return (
        _bounded_corr(float(spec.get("rf_pqtl", 0.0))),
        _bounded_corr(float(spec.get("rf_outcome", spec.get("rf_cancer", 0.0)))),
        _bounded_corr(float(spec.get("pqtl_outcome", spec.get("pqtl_cancer", 0.0)))),
    )


def sample_observed_effects(
    rng: np.random.Generator,
    cell: Dict,
    set_name: str,
    rf_true: float,
    pqtl_true: float,
    outcome_true: float,
    rf_se: float,
    pqtl_se: float,
    outcome_se: float,
) -> Tuple[float, float, float]:
    rf_obs = float(rf_true)
    pqtl_obs = float(pqtl_true)
    outcome_obs = float(outcome_true)

    ld_r = float(cell.get("ld_r", 0.0))
    if set_name in set(str(x) for x in cell.get("ld_leakage_sets", ["B", "C"])):
        outcome_obs += ld_r * float(cell.get("linked_pqtl_outcome_leakage", 0.0)) * pqtl_true
        outcome_obs += ld_r * float(cell.get("linked_rf_outcome_leakage", 0.0)) * rf_true

    if rng.random() < float(cell.get("pqtl_contamination_rate", 0.0)):
        pqtl_obs += float(rng.normal(0.0, float(cell.get("pqtl_contamination_sd", 0.0))))
    if rng.random() < float(cell.get("outcome_contamination_rate", 0.0)):
        outcome_obs += float(rng.normal(0.0, float(cell.get("outcome_contamination_sd", 0.0))))
    if rng.random() < float(cell.get("rf_contamination_rate", 0.0)):
        rf_obs += float(rng.normal(0.0, float(cell.get("rf_contamination_sd", 0.0))))

    if bool(cell.get("sampling_noise", False)):
        rf_pqtl, rf_outcome, pqtl_outcome = sampling_error_correlations(cell)
        corr = np.array(
            [
                [1.0, rf_pqtl, rf_outcome],
                [rf_pqtl, 1.0, pqtl_outcome],
                [rf_outcome, pqtl_outcome, 1.0],
            ],
            dtype=float,
        )
        scales = np.array([rf_se, pqtl_se, outcome_se], dtype=float)
        cov = corr * np.outer(scales, scales)
        try:
            noise = rng.multivariate_normal(np.zeros(3), cov)
        except np.linalg.LinAlgError:
            cov = cov + np.eye(3) * 1e-10
            noise = rng.multivariate_normal(np.zeros(3), cov)
        rf_obs += float(noise[0])
        pqtl_obs += float(noise[1])
        outcome_obs += float(noise[2])

    return rf_obs, pqtl_obs, outcome_obs


def sample_set_counts(rng: np.random.Generator, benchmark: str, cell: Dict, scenario: str) -> Tuple[int, int, int]:
    scenario_counts = cell.get("scenario_set_counts", {}).get(scenario)
    if scenario_counts is not None:
        return int(scenario_counts["A"]), int(scenario_counts["B"]), int(scenario_counts["C"])
    scenario_dist = cell.get("scenario_set_count_distribution", {}).get(scenario)
    if scenario_dist is not None:
        return int(rng.choice(scenario_dist["A"])), int(rng.choice(scenario_dist["B"])), int(rng.choice(scenario_dist["C"]))
    if benchmark == "classification" or "set_counts" in cell:
        sets = cell["set_counts"]
        return int(sets["A"]), int(sets["B"]), int(sets["C"])
    dist = cell["set_count_distribution"]
    return int(rng.choice(dist["A"])), int(rng.choice(dist["B"])), int(rng.choice(dist["C"]))


def simulate_effects(
    rng: np.random.Generator,
    scenario: str,
    cell: Dict,
    n_a: int,
    n_b: int,
    n_c: int,
) -> Dict[str, np.ndarray]:
    params = scenario_parameters(cell, scenario)
    beta1 = params["beta1"]
    beta2 = params["beta2"]
    beta3 = params["beta3"]

    delta_sd = params["delta_sd"]
    psi_sd = params["psi_sd"]
    phi_sd = params["phi_sd"]
    rho = params["rho_m5"]
    direction_flipped = scenario == "M1" and rng.random() < params["direction_flip_rate"]
    if direction_flipped:
        beta2 = -beta2
    gamma_lo, gamma_hi = cell["gamma_range"]
    cis_lo, cis_hi = cell["cis_effect_range"]

    gamma_a = rng.uniform(gamma_lo, gamma_hi, size=n_a) if n_a else np.array([])
    gamma_c = rng.uniform(gamma_lo, gamma_hi, size=n_c) if n_c else np.array([])
    cis_b = rng.uniform(cis_lo, cis_hi, size=n_b) if n_b else np.array([])
    cis_c = rng.uniform(cis_lo, cis_hi, size=n_c) if n_c else np.array([])

    if scenario == "M5":
        cov = np.array(
            [
                [delta_sd ** 2, rho * delta_sd * psi_sd],
                [rho * delta_sd * psi_sd, psi_sd ** 2],
            ]
        )
        ac_joint = rng.multivariate_normal([0.0, 0.0], cov, size=n_a + n_c) if (n_a + n_c) else np.empty((0, 2))
        delta_ac = ac_joint[:, 0] if len(ac_joint) else np.array([])
        psi_ac = ac_joint[:, 1] if len(ac_joint) else np.array([])
    else:
        delta_ac = rng.normal(0.0, delta_sd, size=n_a + n_c) if (n_a + n_c) else np.array([])
        psi_ac = rng.normal(0.0, psi_sd, size=n_a + n_c) if (n_a + n_c) else np.array([])

    phi_b = rng.normal(0.0, phi_sd, size=n_b) if n_b else np.array([])
    delta_a = delta_ac[:n_a]
    delta_c = delta_ac[n_a:]
    psi_a = psi_ac[:n_a]
    psi_c = psi_ac[n_a:]

    alpha_a_true = beta1 * gamma_a + delta_a
    alpha_b_true = cis_b
    alpha_c_true = beta1 * gamma_c + cis_c + delta_c
    gamma_b_true = rng.normal(0.0, 0.005, size=n_b) if n_b else np.array([])

    gamma_out_a = (beta2 * beta1 + beta3) * gamma_a + beta2 * delta_a + psi_a
    gamma_out_b = beta2 * cis_b + phi_b
    gamma_out_c = (beta2 * beta1 + beta3) * gamma_c + beta2 * (cis_c + delta_c) + psi_c

    if scenario == "M0":
        alpha_a_true[:] = 0.0
        alpha_b_true[:] = 0.0
        alpha_c_true[:] = 0.0
        gamma_out_a[:] = 0.0
        gamma_out_b[:] = 0.0
        gamma_out_c[:] = 0.0
    elif scenario == "M2":
        gamma_out_a = psi_a
        gamma_out_b = phi_b
        gamma_out_c = psi_c
    elif scenario == "M3":
        alpha_a_true[:] = 0.0
        alpha_b_true[:] = 0.0
        alpha_c_true[:] = 0.0
        gamma_out_a = beta3 * gamma_a + psi_a
        gamma_out_b = phi_b
        gamma_out_c = beta3 * gamma_c + psi_c
    elif scenario == "M4":
        alpha_a_true[:] = 0.0
        alpha_c_true[:] = 0.0
        gamma_out_a = psi_a
        gamma_out_b = beta2 * cis_b + phi_b
        gamma_out_c = beta2 * cis_c + psi_c
    elif scenario == "M5":
        gamma_out_a = beta3 * gamma_a + psi_a
        gamma_out_b = phi_b
        gamma_out_c = beta3 * gamma_c + psi_c

    return {
        "true_beta1": beta1,
        "true_beta2": beta2,
        "true_beta3": beta3,
        "direction_flipped": direction_flipped,
        "gamma_a": gamma_a,
        "gamma_b": gamma_b_true,
        "gamma_c": gamma_c,
        "alpha_a": alpha_a_true,
        "alpha_b": alpha_b_true,
        "alpha_c": alpha_c_true,
        "Gamma_a": gamma_out_a,
        "Gamma_b": gamma_out_b,
        "Gamma_c": gamma_out_c,
    }


def write_dataset(
    outdir: Path,
    benchmark: str,
    cell_name: str,
    replicate: int,
    proteins: List[Dict],
) -> BenchmarkFiles:
    rep_dir = ensure_dir(outdir / benchmark / cell_name / f"rep_{replicate:04d}")
    rf_rows: List[Dict] = []
    pqtl_rows: List[Dict] = []
    cancer_rows: List[Dict] = []
    info_rows: List[Dict] = []
    truth_rows: List[Dict] = []

    for protein in proteins:
        pid = protein["protein_id"]
        gene = protein["gene_name"]
        chr_num = protein["chr"]
        gene_start = protein["gene_start"]
        gene_end = protein["gene_end"]
        rf_se = protein["rf_se"]
        pqtl_se = protein["pqtl_se"]
        outcome_se = protein["outcome_se"]
        set_counts = protein["set_counts"]

        info_rows.append({"PROTEIN": pid, "GENE": gene, "CHR": chr_num, "START": gene_start, "END": gene_end})
        truth_rows.append(
            {
                "protein_id": pid,
                "gene_name": gene,
                "true_scenario": protein["scenario"],
                "nA_true": set_counts["A"],
                "nB_true": set_counts["B"],
                "nC_true": set_counts["C"],
                "true_beta1": protein["true_beta1"],
                "true_beta2": protein["true_beta2"],
                "true_beta3": protein["true_beta3"],
                "true_mediated_effect": protein["true_beta1"] * protein["true_beta2"],
                "direction_flipped": int(bool(protein.get("direction_flipped", False))),
            }
        )

        for set_name in ("A", "B", "C"):
            for rec in protein["snps"][set_name]:
                rf_rows.append(
                    {
                        "SNP": rec["rsid"],
                        "A1": rec["a1"],
                        "A2": rec["a2"],
                        "FREQ": rec["freq"],
                        "BETA": rec["rf_beta"],
                        "SE": rf_se,
                        "P": rec["rf_p"],
                        "CHR": chr_num,
                        "BP": rec["bp"],
                    }
                )
                pqtl_rows.append(
                    {
                        "PROTEIN": pid,
                        "SNP": rec["rsid"],
                        "A1": rec["a1"],
                        "A2": rec["a2"],
                        "FREQ": rec["freq"],
                        "BETA": rec["pqtl_beta"],
                        "SE": pqtl_se,
                        "P": rec["pqtl_p"],
                        "CHR": chr_num,
                        "BP": rec["bp"],
                    }
                )
                cancer_rows.append(
                    {
                        "SNP": rec["rsid"],
                        "A1": rec["a1"],
                        "A2": rec["a2"],
                        "FREQ": rec["freq"],
                        "BETA": rec["outcome_beta"],
                        "SE": outcome_se,
                        "P": rec["outcome_p"],
                        "CHR": chr_num,
                        "BP": rec["bp"],
                    }
                )

    rf = rep_dir / "rf_sumstat.txt"
    pqtl = rep_dir / "pqtl_sumstat.txt"
    cancer = rep_dir / "cancer_sumstat.txt"
    protein_info = rep_dir / "protein_info.txt"
    truth = rep_dir / "truth.tsv"

    write_table(rf, rf_rows, ["SNP", "A1", "A2", "FREQ", "BETA", "SE", "P", "CHR", "BP"], " ")
    write_table(pqtl, pqtl_rows, ["PROTEIN", "SNP", "A1", "A2", "FREQ", "BETA", "SE", "P", "CHR", "BP"], " ")
    write_table(cancer, cancer_rows, ["SNP", "A1", "A2", "FREQ", "BETA", "SE", "P", "CHR", "BP"], " ")
    write_table(protein_info, info_rows, ["PROTEIN", "GENE", "CHR", "START", "END"], " ")
    write_table(
        truth,
        truth_rows,
        ["protein_id", "gene_name", "true_scenario", "nA_true", "nB_true", "nC_true", "true_beta1", "true_beta2", "true_beta3", "true_mediated_effect", "direction_flipped"],
        "\t",
    )
    return BenchmarkFiles(rf=rf, pqtl=pqtl, cancer=cancer, protein_info=protein_info, truth=truth)


def generate_proteins(
    rng: np.random.Generator,
    benchmark: str,
    cell: Dict,
    scenario_sequence: Iterable[str],
) -> List[Dict]:
    proteins: List[Dict] = []
    for idx, scenario in enumerate(scenario_sequence, start=1):
        n_a, n_b, n_c = sample_set_counts(rng, benchmark, cell, scenario)
        effects = simulate_effects(rng, scenario, cell, n_a, n_b, n_c)
        chr_num = int(rng.integers(1, 23))
        base = 10_000_000 + idx * 50_000
        gene_start = base
        gene_end = base + 20_000
        alleles = choose_alleles(rng, n_a + n_b + n_c)
        freqs = rng.uniform(0.08, 0.45, size=n_a + n_b + n_c)

        snps = {"A": [], "B": [], "C": []}
        allele_cursor = 0
        for set_name, count in (("A", n_a), ("B", n_b), ("C", n_c)):
            for j in range(count):
                a1, a2 = alleles[allele_cursor]
                freq = round(float(freqs[allele_cursor]), 4)
                allele_cursor += 1
                if set_name == "A":
                    rf_true = float(effects["gamma_a"][j])
                    pqtl_true = float(effects["alpha_a"][j])
                    outcome_true = float(effects["Gamma_a"][j])
                    bp = gene_start - 2_000_000 - j * 2500
                elif set_name == "B":
                    rf_true = float(effects["gamma_b"][j])
                    pqtl_true = float(effects["alpha_b"][j])
                    outcome_true = float(effects["Gamma_b"][j])
                    bp = gene_start + 2000 + j * 1200
                else:
                    rf_true = float(effects["gamma_c"][j])
                    pqtl_true = float(effects["alpha_c"][j])
                    outcome_true = float(effects["Gamma_c"][j])
                    bp = gene_start + 9000 + j * 1400

                rf_se = float(cell["rf_se"])
                pqtl_se = float(cell["pqtl_se"])
                outcome_se = float(cell["outcome_se"])
                rf_beta, pqtl_beta, outcome_beta = sample_observed_effects(
                    rng,
                    cell,
                    set_name,
                    rf_true,
                    pqtl_true,
                    outcome_true,
                    rf_se,
                    pqtl_se,
                    outcome_se,
                )
                rf_p = float(norm_sf(np.array([rf_beta / rf_se]))[0] * 2.0)
                pqtl_p = float(norm_sf(np.array([pqtl_beta / pqtl_se]))[0] * 2.0)
                outcome_p = float(norm_sf(np.array([outcome_beta / outcome_se]))[0] * 2.0)
                snps[set_name].append(
                    {
                        "rsid": f"rs{scenario.lower()}_{idx}_{set_name.lower()}{j + 1}",
                        "a1": str(a1),
                        "a2": str(a2),
                        "freq": freq,
                        "rf_beta": rf_beta,
                        "pqtl_beta": pqtl_beta,
                        "outcome_beta": outcome_beta,
                        "rf_p": rf_p,
                        "pqtl_p": pqtl_p,
                        "outcome_p": outcome_p,
                        "bp": int(bp),
                    }
                )

        params = scenario_parameters(cell, scenario)
        proteins.append(
            {
                "protein_id": f"{scenario}_P{idx:04d}",
                "gene_name": f"GENE_{scenario}_{idx:04d}",
                "scenario": scenario,
                "chr": chr_num,
                "gene_start": gene_start,
                "gene_end": gene_end,
                "rf_se": float(cell["rf_se"]),
                "pqtl_se": float(cell["pqtl_se"]),
                "outcome_se": float(cell["outcome_se"]),
                "true_beta1": effects["true_beta1"],
                "true_beta2": effects["true_beta2"],
                "true_beta3": effects["true_beta3"],
                "direction_flipped": effects["direction_flipped"],
                "set_counts": {"A": n_a, "B": n_b, "C": n_c},
                "snps": snps,
            }
        )
    return proteins


def classification_sequence(proteins_per_scenario: int, scenarios: Sequence[str] = SCENARIOS) -> List[str]:
    seq: List[str] = []
    for scenario in scenarios:
        seq.extend([scenario] * proteins_per_scenario)
    return seq


def calibration_sequence(
    rng: np.random.Generator,
    proteins_per_replicate: int,
    scenario_mix: Dict[str, float],
    scenarios: Sequence[str] = SCENARIOS,
) -> List[str]:
    probs = np.array([scenario_mix[s] for s in scenarios], dtype=float)
    probs = probs / probs.sum()
    return list(rng.choice(np.array(list(scenarios), dtype=object), size=proteins_per_replicate, p=probs))


def run_bmediator(binary: Path, files: BenchmarkFiles, out_prefix: Path, global_cfg: Dict) -> Tuple[Path, Path]:
    binary = binary.resolve()
    cmd = [
        str(binary),
        "--rf-sumstat",
        str(files.rf),
        "--pqtl-sumstat",
        str(files.pqtl),
        "--cancer-sumstat",
        str(files.cancer),
        "--protein-info",
        str(files.protein_info),
        "--out",
        str(out_prefix),
        "--p-thresh-rf",
        str(global_cfg["rf_p_threshold"]),
        "--p-thresh-cis",
        str(global_cfg["cis_p_threshold"]),
        "--max-eb-iter",
        str(global_cfg["max_eb_iter"]),
        "--max-cavi-iter",
        str(global_cfg["max_cavi_iter"]),
        "--threads",
        str(global_cfg["threads"]),
    ]
    optional_args = {
        "eb_tol": "--eb-tol",
        "elbo_tol": "--elbo-tol",
        "prior_p0": "--prior-p0",
        "prior_p1": "--prior-p1",
        "prior_p2": "--prior-p2",
        "prior_p3": "--prior-p3",
        "prior_p4": "--prior-p4",
        "prior_p5": "--prior-p5",
        "sigma2_beta1": "--sigma2-beta1",
        "sigma2_beta2": "--sigma2-beta2",
        "sigma2_beta3": "--sigma2-beta3",
        "direction_mode": "--direction-mode",
        "direction_weight": "--direction-weight",
        "direction_min_prob": "--direction-min-prob",
        "m1_min_cis_only": "--m1-min-cis-only",
        "m1_min_first_stage_z": "--m1-min-first-stage-z",
        "m1_min_second_stage_z": "--m1-min-second-stage-z",
        "m1_resid_corr_threshold": "--m1-resid-corr-threshold",
        "m1_resid_corr_penalty": "--m1-resid-corr-penalty",
    }
    if "cis_window_bp" in global_cfg:
        cmd.extend(["--cis-window", str(int(float(global_cfg["cis_window_bp"]) / 1000))])
    if bool(global_cfg.get("fixed_priors", False)):
        cmd.append("--fixed-priors")
    for key, flag in optional_args.items():
        if key in global_cfg:
            cmd.extend([flag, str(global_cfg[key])])
    extra = global_cfg.get("binary_options", {})
    if isinstance(extra, list):
        cmd.extend(str(value) for value in extra)
    elif isinstance(extra, dict):
        for key, value in extra.items():
            flag = str(key) if str(key).startswith("--") else "--" + str(key).replace("_", "-")
            if value is True:
                cmd.append(flag)
            elif value is not False and value is not None:
                cmd.extend([flag, str(value)])
    env = os.environ.copy()
    env.setdefault("MPLCONFIGDIR", str(out_prefix.parent / ".mplconfig"))
    subprocess.run(cmd, check=True, env=env)
    return Path(f"{out_prefix}.mediation"), Path(f"{out_prefix}.hyp")


def _resolve_column(rows: Sequence[Dict], candidates: Iterable[str]) -> str:
    if not rows:
        raise KeyError("empty table")
    columns = set(rows[0].keys())
    for cand in candidates:
        if cand in columns:
            return cand
    raise KeyError(f"could not resolve any of columns: {list(candidates)}")


def read_truth(path: Path) -> List[Dict]:
    rows = read_table(path, "\t")
    for row in rows:
        for key in ("nA_true", "nB_true", "nC_true"):
            row[key] = int(row[key])
        for key in ("true_beta1", "true_beta2", "true_beta3", "true_mediated_effect"):
            row[key] = float(row[key])
        if "direction_flipped" in row:
            row["direction_flipped"] = int(row["direction_flipped"])
    return rows


def read_results(path: Path) -> List[Dict]:
    rows = read_table(path, "\t")
    for row in rows:
        for key, value in list(row.items()):
            if value in ("YES", "NO", "", None):
                continue
            try:
                if any(ch in value for ch in (".", "e", "E")):
                    row[key] = float(value)
                else:
                    row[key] = int(value)
            except Exception:
                pass
    return rows


def attach_truth(results: List[Dict], truth: List[Dict]) -> List[Dict]:
    protein_col = _resolve_column(results, ["Protein", "protein_id"])
    result_map = {str(row[protein_col]): row for row in results}
    merged = []
    for row in truth:
        merged_row = dict(row)
        merged_row.update(result_map.get(str(row["protein_id"]), {}))
        merged.append(merged_row)
    return merged


def compute_benchmark_metrics(rows: List[Dict], alias_map: Dict[str, str] | None = None) -> List[Dict]:
    alias_map = alias_map or {}
    prob_cols = {
        "M0": _resolve_column(rows, ["P_M0", "prob_M0"]),
        "M1": _resolve_column(rows, ["P_M1", "prob_M1"]),
        "M2": _resolve_column(rows, ["P_M2", "prob_M2"]),
        "M3": _resolve_column(rows, ["P_M3", "prob_M3"]),
        "M4": _resolve_column(rows, ["P_M4", "prob_M4"]),
        "M5": _resolve_column(rows, ["P_M5", "prob_M5"]),
    }
    scored = []
    for row in rows:
        probs = []
        for scenario in SCENARIOS:
            value = row.get(prob_cols[scenario])
            try:
                probs.append(float(value))
            except Exception:
                probs.append(float("nan"))
        raw_pred = None if all(np.isnan(probs)) else SCENARIOS[int(np.nanargmax(probs))]
        true_scenario = effective_scenario(str(row["true_scenario"]), alias_map)
        pred = effective_scenario(raw_pred, alias_map) if raw_pred is not None else None
        enriched = dict(row)
        enriched["raw_true_scenario"] = row["true_scenario"]
        enriched["raw_pred_scenario"] = raw_pred
        enriched["true_scenario"] = true_scenario
        enriched["pred_scenario"] = pred
        enriched["analyzed"] = int(raw_pred is not None)
        enriched["correct_class"] = int(pred == true_scenario) if raw_pred is not None else None
        try:
            enriched["p_m1"] = float(row[prob_cols["M1"]])
        except Exception:
            enriched["p_m1"] = None
        raw_true = str(row["true_scenario"])
        try:
            enriched["true_state_prob"] = float(row[prob_cols[raw_true]]) if raw_true in prob_cols else None
        except Exception:
            enriched["true_state_prob"] = None
        enriched["is_true_m1"] = int(true_scenario == "M1")
        scored.append(enriched)

    scored.sort(key=lambda row: (row["p_m1"] is not None, row["p_m1"] or -float("inf")), reverse=True)
    cumulative = 0.0
    rank = 0
    for row in scored:
        if row["p_m1"] is None:
            row["rank"] = None
            row["estimated_bfdr"] = None
            continue
        rank += 1
        cumulative += 1.0 - row["p_m1"]
        row["rank"] = rank
        row["estimated_bfdr"] = cumulative / rank
    return scored


def summarize_fdr_power(rows: List[Dict], thresholds: Sequence[float] = (0.01, 0.05, 0.10)) -> List[Dict]:
    rows = [row for row in rows if row["estimated_bfdr"] is not None]
    total_true_m1 = sum(1 for row in rows if row["is_true_m1"] == 1)
    out = []
    for threshold in thresholds:
        selected = [row for row in rows if row["estimated_bfdr"] <= threshold]
        true_positives = sum(1 for row in selected if row["is_true_m1"] == 1)
        false_positives = len(selected) - true_positives
        out.append(
            {
                "bfdr_threshold": threshold,
                "n_selected": len(selected),
                "true_positives": true_positives,
                "false_positives": false_positives,
                "empirical_fdr": (false_positives / len(selected)) if selected else None,
                "power": (true_positives / total_true_m1) if total_true_m1 else None,
            }
        )
    return out


def summarize_classification(rows: List[Dict], scenarios: Sequence[str] | None = None) -> Tuple[List[Dict], List[Dict]]:
    scenarios = tuple(scenarios or tuple(dict.fromkeys(row["true_scenario"] for row in rows)))
    by_scenario = []
    confusion = []
    for true_scenario in scenarios:
        subset = [row for row in rows if row["true_scenario"] == true_scenario]
        if subset:
            analyzed = [row for row in subset if row["analyzed"] == 1]
            by_scenario.append(
                {
                    "true_scenario": true_scenario,
                    "n": len(subset),
                    "n_analyzed": len(analyzed),
                    "accuracy": float(np.mean([row["correct_class"] for row in analyzed])) if analyzed else None,
                    "mean_p_m1": float(np.mean([row["p_m1"] for row in analyzed])) if analyzed else None,
                    "mean_true_state_prob": float(np.mean([row["true_state_prob"] for row in analyzed])) if analyzed else None,
                }
            )
        for pred in scenarios:
            confusion.append(
                {
                    "true_scenario": true_scenario,
                    "pred_scenario": pred,
                    "count": sum(1 for row in rows if row["true_scenario"] == true_scenario and row["pred_scenario"] == pred),
                }
            )
    return by_scenario, confusion


def summarize_calibration(rows: List[Dict], n_bins: int = 10) -> Tuple[List[Dict], List[Dict]]:
    rows = [row for row in rows if row["p_m1"] is not None]
    bins = np.linspace(0.0, 1.0, n_bins + 1)
    cal_rows = []
    for idx in range(n_bins):
        lo = bins[idx]
        hi = bins[idx + 1]
        if idx == n_bins - 1:
            subset = [row for row in rows if lo <= row["p_m1"] <= hi]
        else:
            subset = [row for row in rows if lo <= row["p_m1"] < hi]
        if not subset:
            continue
        cal_rows.append(
            {
                "calibration_bin": idx,
                "n": len(subset),
                "mean_pred": float(np.mean([row["p_m1"] for row in subset])),
                "observed_m1": float(np.mean([row["is_true_m1"] for row in subset])),
                "bin_mid": float(0.5 * (lo + hi)),
            }
        )

    rank_breaks = [10, 25, 50, 100, len(rows)]
    start = 0
    bfdr_rows = []
    for stop in rank_breaks:
        if stop <= start:
            continue
        subset = rows[start: min(stop, len(rows))]
        if not subset:
            continue
        bfdr_rows.append(
            {
                "rank_bin": f"{start + 1}-{min(stop, len(rows))}",
                "n": len(subset),
                "mean_estimated_bfdr": float(np.mean([row["estimated_bfdr"] for row in subset])),
                "empirical_fdr": float(1.0 - np.mean([row["is_true_m1"] for row in subset])),
            }
        )
        start = stop
    return cal_rows, bfdr_rows
