#!/usr/bin/env python3
"""Evaluate fixed-prior factorized Bayesian FDR from completed task outputs."""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path

from summarize_factorized_calibration import as_float, mean, quantile


def posterior_from_log_bf(log_bf: float, prior: float) -> float:
    log_odds = log_bf + math.log(prior) - math.log1p(-prior)
    if log_odds >= 0.0:
        return 1.0 / (1.0 + math.exp(-log_odds))
    exp_value = math.exp(log_odds)
    return exp_value / (1.0 + exp_value)


def logsumexp(values: list[float]) -> float:
    maximum = max(values)
    return maximum + math.log(sum(math.exp(value - maximum) for value in values))


def slope_posterior(row: dict[str, str], leg: str, prior_slope: float,
                    prior_directional: float) -> float:
    log_bf_slope = as_float(row.get(f"factor_log_BF_slope_only_{leg}"))
    log_bf_directional = as_float(row.get(f"factor_log_BF_directional_only_{leg}"))
    log_bf_both = as_float(row.get(f"factor_log_BF_slope_directional_{leg}"))
    if not all(math.isfinite(value) for value in (
        log_bf_slope, log_bf_directional, log_bf_both
    )):
        return math.nan
    weights = [
        math.log1p(-prior_slope) + math.log1p(-prior_directional),
        log_bf_slope + math.log(prior_slope) + math.log1p(-prior_directional),
        log_bf_directional + math.log1p(-prior_slope) + math.log(prior_directional),
        log_bf_both + math.log(prior_slope) + math.log(prior_directional),
    ]
    return math.exp(logsumexp([weights[1], weights[3]]) - logsumexp(weights))


def task_metrics(
    rows: list[dict[str, str]], prior_xm: float, prior_my: float,
    prior_directional: float, alpha: float
) -> dict[str, float]:
    tested = []
    for row in rows:
        if row.get("identification_class") == "nonidentifiable":
            continue
        posterior_xm = slope_posterior(row, "XM", prior_xm, prior_directional)
        posterior_my = slope_posterior(row, "MY", prior_my, prior_directional)
        if not math.isfinite(posterior_xm) or not math.isfinite(posterior_my):
            continue
        posterior_mediation = posterior_xm * posterior_my
        tested.append((posterior_mediation, row.get("true_scenario") == "M1"))

    tested.sort(reverse=True, key=lambda value: value[0])
    selected = []
    cumulative_false_probability = 0.0
    for rank, (posterior, truth) in enumerate(tested, start=1):
        cumulative_false_probability += 1.0 - posterior
        bayes_q = cumulative_false_probability / rank
        if bayes_q <= alpha:
            selected.append((posterior, truth))

    true_total = sum(truth for _, truth in tested)
    true_selected = sum(truth for _, truth in selected)
    false_selected = len(selected) - true_selected
    return {
        "n_tested": float(len(tested)),
        "n_selected": float(len(selected)),
        "fdr": false_selected / len(selected) if selected else 0.0,
        "power": true_selected / true_total if true_total else math.nan,
        "mean_posterior_true": mean([posterior for posterior, truth in tested if truth]),
        "mean_posterior_null": mean([posterior for posterior, truth in tested if not truth]),
        "brier": mean([(posterior - float(truth)) ** 2 for posterior, truth in tested]),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--prior-xm", type=float, nargs="+", default=(0.05, 0.10, 0.25, 0.50))
    parser.add_argument("--prior-my", type=float, nargs="+", default=(0.01, 0.05, 0.10, 0.25))
    parser.add_argument("--prior-directional", type=float, nargs="+", default=(0.05, 0.10, 0.25))
    parser.add_argument("--alpha", type=float, default=0.05)
    args = parser.parse_args()
    if not 0.0 < args.alpha < 1.0:
        raise SystemExit("--alpha must be in (0,1)")
    if any(not 0.0 < prior < 1.0 for prior in
           args.prior_xm + args.prior_my + args.prior_directional):
        raise SystemExit("all prior probabilities must be in (0,1)")

    grouped: dict[tuple[str, str, float, float, float], list[dict[str, float]]] = defaultdict(list)
    paths = sorted(args.input.glob("*/*/rep_*/task_metrics.tsv"))
    if not paths:
        raise SystemExit(f"No task outputs under {args.input}")
    for path in paths:
        with path.open(newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        relative = path.relative_to(args.input)
        benchmark, cell = relative.parts[0], relative.parts[1]
        for prior_xm in args.prior_xm:
            for prior_my in args.prior_my:
                for prior_directional in args.prior_directional:
                    grouped[(benchmark, cell, prior_xm, prior_my,
                             prior_directional)].append(
                        task_metrics(rows, prior_xm, prior_my,
                                     prior_directional, args.alpha)
                    )

    fields = [
        "benchmark", "cell", "prior_XM", "prior_MY", "prior_directional", "replicates",
        "mean_tested", "mean_selected", "mean_fdr", "fdr_q95", "mean_power",
        "power_q05", "mean_posterior_true", "mean_posterior_null", "mean_brier",
    ]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for (benchmark, cell, prior_xm, prior_my, prior_directional), records in sorted(grouped.items()):
            writer.writerow({
                "benchmark": benchmark,
                "cell": cell,
                "prior_XM": prior_xm,
                "prior_MY": prior_my,
                "prior_directional": prior_directional,
                "replicates": len(records),
                "mean_tested": mean([record["n_tested"] for record in records]),
                "mean_selected": mean([record["n_selected"] for record in records]),
                "mean_fdr": mean([record["fdr"] for record in records]),
                "fdr_q95": quantile([record["fdr"] for record in records], 0.95),
                "mean_power": mean([record["power"] for record in records]),
                "power_q05": quantile([record["power"] for record in records], 0.05),
                "mean_posterior_true": mean([record["mean_posterior_true"] for record in records]),
                "mean_posterior_null": mean([record["mean_posterior_null"] for record in records]),
                "mean_brier": mean([record["brier"] for record in records]),
            })
    print(f"Wrote Bayesian prior sensitivity for {len(paths)} task families to {args.output}")


if __name__ == "__main__":
    main()
