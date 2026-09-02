from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
from typing import Dict, List, Optional, Sequence


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compute truth-free sensitivity-bounded mediation p-values and multiplicity-adjusted q-values "
            "from a BMEDIATOR .mediation file."
        )
    )
    parser.add_argument("--input", type=Path, required=True, help="BMEDIATOR .mediation TSV.")
    parser.add_argument("--output", type=Path, required=True, help="Output TSV.")
    parser.add_argument(
        "--bias-bound",
        type=float,
        default=None,
        help="Global absolute bound Delta on protein-to-outcome bias.",
    )
    parser.add_argument(
        "--bounds-file",
        type=Path,
        help="Optional TSV with Protein and bias_bound columns; overrides the global bound by protein.",
    )
    parser.add_argument(
        "--adjustment",
        choices=("BH", "BY"),
        default="BY",
        help="BY is valid under arbitrary dependence; BH is less conservative under PRDS assumptions.",
    )
    return parser.parse_args()


def finite_float(value, default: float = float("nan")) -> float:
    try:
        result = float(value)
        if math.isfinite(result):
            return result
    except Exception:
        pass
    return default


def normal_two_sided_p(z_value: float) -> float:
    if not math.isfinite(z_value):
        return 1.0
    return min(1.0, math.erfc(abs(z_value) / math.sqrt(2.0)))


def interval_null_p(beta: float, se: float, bias_bound: float) -> float:
    """Conservative p-value for H0: |beta| <= bias_bound."""
    if not math.isfinite(beta) or not math.isfinite(se) or se <= 0.0:
        return 1.0
    beyond = (abs(beta) - bias_bound) / se
    return normal_two_sided_p(beyond) if beyond > 0.0 else 1.0


def adjust_pvalues(p_values: Sequence[float], method: str) -> List[float]:
    count = len(p_values)
    if count == 0:
        return []
    harmonic = sum(1.0 / rank for rank in range(1, count + 1)) if method == "BY" else 1.0
    order = sorted(range(count), key=lambda idx: p_values[idx])
    adjusted = [1.0] * count
    previous = 1.0
    for reverse_rank, idx in enumerate(reversed(order), start=1):
        rank = count - reverse_rank + 1
        value = min(previous, p_values[idx] * count * harmonic / rank)
        adjusted[idx] = min(1.0, value)
        previous = value
    return adjusted


def read_bounds(path: Optional[Path]) -> Dict[str, float]:
    if path is None:
        return {}
    result: Dict[str, float] = {}
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "Protein" not in reader.fieldnames or "bias_bound" not in reader.fieldnames:
            raise ValueError("--bounds-file must contain Protein and bias_bound columns")
        for row in reader:
            bound = finite_float(row.get("bias_bound"))
            if not math.isfinite(bound) or bound < 0.0:
                raise ValueError(f"Invalid bias_bound for {row.get('Protein', '')}: {row.get('bias_bound', '')}")
            result[str(row["Protein"])] = bound
    return result


def main() -> None:
    args = parse_args()
    if args.bias_bound is None and args.bounds_file is None:
        raise ValueError("Specify --bias-bound, --bounds-file, or both")
    if args.bias_bound is not None and args.bias_bound < 0.0:
        raise ValueError("--bias-bound must be non-negative")

    per_protein_bounds = read_bounds(args.bounds_file)
    with args.input.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError("Input contains no result rows")

    required = {"Protein", "ivw_rf_to_pp_p", "ivw_pp_to_outcome_beta", "ivw_pp_to_outcome_se"}
    missing = required.difference(rows[0])
    if missing:
        raise ValueError(f"Input is missing required columns: {', '.join(sorted(missing))}")

    p_values: List[float] = []
    bounds: List[float] = []
    p_alpha_values: List[float] = []
    p_beta_values: List[float] = []
    for row in rows:
        protein = str(row["Protein"])
        if protein in per_protein_bounds:
            bound = per_protein_bounds[protein]
        elif args.bias_bound is not None:
            bound = args.bias_bound
        else:
            raise ValueError(f"No bias bound supplied for protein {protein}")
        p_alpha = finite_float(row.get("ivw_rf_to_pp_p"), 1.0)
        p_alpha = min(1.0, max(0.0, p_alpha))
        p_beta = interval_null_p(
            finite_float(row.get("ivw_pp_to_outcome_beta")),
            finite_float(row.get("ivw_pp_to_outcome_se")),
            bound,
        )
        bounds.append(bound)
        p_alpha_values.append(p_alpha)
        p_beta_values.append(p_beta)
        p_values.append(max(p_alpha, p_beta))

    q_values = adjust_pvalues(p_values, args.adjustment)
    output_fields = list(rows[0].keys()) + [
        "analytic_bias_bound",
        "analytic_p_rf_to_protein",
        "analytic_p_protein_to_outcome_beyond_bound",
        "analytic_p_mediation",
        "analytic_adjustment",
        "analytic_q_mediation",
    ]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=output_fields, delimiter="\t")
        writer.writeheader()
        for idx, row in enumerate(rows):
            enriched = dict(row)
            enriched.update(
                {
                    "analytic_bias_bound": bounds[idx],
                    "analytic_p_rf_to_protein": p_alpha_values[idx],
                    "analytic_p_protein_to_outcome_beyond_bound": p_beta_values[idx],
                    "analytic_p_mediation": p_values[idx],
                    "analytic_adjustment": args.adjustment,
                    "analytic_q_mediation": q_values[idx],
                }
            )
            writer.writerow(enriched)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
