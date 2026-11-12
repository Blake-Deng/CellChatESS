#!/usr/bin/env python3
"""Summarize the balanced CellChat 2.2 / AutoZyme 0.3.1 queue."""

from __future__ import annotations

import csv
import math
import random
import statistics
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results" / "recomputed_autozyme"
RESULTS.mkdir(exist_ok=True)
INPUT = ROOT / "results" / "autozyme_192" / "all_successful_runs.tsv"
BOOTSTRAP_SEED = 20261027
BOOTSTRAP_REPLICATES = 100_000


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], columns: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def quantile(values: list[float], probability: float) -> float:
    ordered = sorted(values)
    if len(ordered) == 1:
        return ordered[0]
    position = (len(ordered) - 1) * probability
    low = math.floor(position)
    high = math.ceil(position)
    if low == high:
        return ordered[low]
    return ordered[low] + (position - low) * (ordered[high] - ordered[low])


def median(values: list[float]) -> float:
    return statistics.median(values)


def main() -> None:
    rows = read_tsv(INPUT)
    if len(rows) != 576:
        raise SystemExit(f"Expected 576 successful method rows, observed {len(rows)}")
    if any(r["status"] != "ok" or r["exit_status"] != "0" for r in rows):
        raise SystemExit("Queue contains failed or nonzero-exit rows")

    runs: dict[tuple[str, str, str], dict[str, dict[str, str]]] = defaultdict(dict)
    for row in rows:
        runs[(row["dataset_id"], row["scale_label"], row["repeat"])][row["method"]] = row
    if len(runs) != 192 or any(set(methods) != {"original", "autozyme", "ess"} for methods in runs.values()):
        raise SystemExit("Expected 192 complete paired three-method timing blocks")

    methods = ("original", "autozyme", "ess")
    method_rows = []
    for method in methods:
        method_data = [row for row in rows if row["method"] == method]
        for metric in ("elapsed_sec", "full_process_wall_sec", "peak_RSS_GiB"):
            values = [float(row[metric]) for row in method_data]
            method_rows.append({
                "method": method,
                "metric": metric,
                "n": len(values),
                "median": median(values),
                "q1": quantile(values, 0.25),
                "q3": quantile(values, 0.75),
                "minimum": min(values),
                "maximum": max(values),
            })
    write_tsv(
        RESULTS / "method_resource_summary.tsv",
        method_rows,
        ["method", "metric", "n", "median", "q1", "q3", "minimum", "maximum"],
    )

    comparisons = (("original", "ess"), ("autozyme", "ess"), ("original", "autozyme"))
    paired_rows = []
    dataset_rows = []
    per_dataset: dict[tuple[str, str], list[float]] = defaultdict(list)
    for (dataset, scale, repeat), method_map in runs.items():
        for numerator, denominator in comparisons:
            ratio = float(method_map[numerator]["elapsed_sec"]) / float(method_map[denominator]["elapsed_sec"])
            per_dataset[(dataset, f"{numerator}/{denominator}")].append(ratio)

    for numerator, denominator in comparisons:
        all_ratios = [
            float(method_map[numerator]["elapsed_sec"]) / float(method_map[denominator]["elapsed_sec"])
            for method_map in runs.values()
        ]
        paired_rows.append({
            "comparison": f"{numerator}/{denominator}",
            "timing_metric": "inference_function_elapsed_sec",
            "n_paired_blocks": len(all_ratios),
            "median_paired_ratio": median(all_ratios),
            "q1_paired_ratio": quantile(all_ratios, 0.25),
            "q3_paired_ratio": quantile(all_ratios, 0.75),
            "numerator_faster_count": sum(value < 1 for value in all_ratios),
            "denominator_faster_count": sum(value > 1 for value in all_ratios),
            "equal_count": sum(value == 1 for value in all_ratios),
            "interpretation": f"ratio > 1 means {numerator} took longer than {denominator}",
        })

        dataset_medians = {
            dataset: median(per_dataset[(dataset, f"{numerator}/{denominator}")])
            for dataset in sorted({key[0] for key in per_dataset})
        }
        rng = random.Random(BOOTSTRAP_SEED)
        dataset_names = list(dataset_medians)
        bootstrap = []
        for _ in range(BOOTSTRAP_REPLICATES):
            sample = [dataset_medians[rng.choice(dataset_names)] for _ in dataset_names]
            bootstrap.append(median(sample))
        paired_rows[-1].update({
            "median_of_dataset_medians": median(list(dataset_medians.values())),
            "dataset_cluster_bootstrap_ci95_lower": quantile(bootstrap, 0.025),
            "dataset_cluster_bootstrap_ci95_upper": quantile(bootstrap, 0.975),
            "bootstrap_unit": "dataset median; resample 12 datasets with replacement",
            "bootstrap_replicates": BOOTSTRAP_REPLICATES,
            "bootstrap_seed": BOOTSTRAP_SEED,
        })
        for dataset in sorted(dataset_medians):
            values = per_dataset[(dataset, f"{numerator}/{denominator}")]
            dataset_rows.append({
                "dataset_id": dataset,
                "comparison": f"{numerator}/{denominator}",
                "n_paired_blocks": len(values),
                "median_paired_ratio": median(values),
                "q1_paired_ratio": quantile(values, 0.25),
                "q3_paired_ratio": quantile(values, 0.75),
                "numerator_faster_count": sum(value < 1 for value in values),
                "denominator_faster_count": sum(value > 1 for value in values),
            })

    write_tsv(
        RESULTS / "paired_speedup_summary.tsv",
        paired_rows,
        ["comparison", "timing_metric", "n_paired_blocks", "median_paired_ratio", "q1_paired_ratio",
         "q3_paired_ratio", "numerator_faster_count", "denominator_faster_count", "equal_count",
         "median_of_dataset_medians", "dataset_cluster_bootstrap_ci95_lower",
         "dataset_cluster_bootstrap_ci95_upper", "bootstrap_unit", "bootstrap_replicates", "bootstrap_seed",
         "interpretation"],
    )
    write_tsv(
        RESULTS / "dataset_paired_speedups.tsv",
        dataset_rows,
        ["dataset_id", "comparison", "n_paired_blocks", "median_paired_ratio", "q1_paired_ratio",
         "q3_paired_ratio", "numerator_faster_count", "denominator_faster_count"],
    )

    scale_rows = []
    for scale in sorted({key[1] for key in runs}):
        blocks = [method_map for (_, block_scale, _), method_map in runs.items() if block_scale == scale]
        for numerator, denominator in comparisons:
            values = [float(block[numerator]["elapsed_sec"]) / float(block[denominator]["elapsed_sec"]) for block in blocks]
            scale_rows.append({
                "scale_label": scale,
                "comparison": f"{numerator}/{denominator}",
                "n_paired_blocks": len(values),
                "median_paired_ratio": median(values),
                "q1_paired_ratio": quantile(values, 0.25),
                "q3_paired_ratio": quantile(values, 0.75),
                "numerator_faster_count": sum(value < 1 for value in values),
                "denominator_faster_count": sum(value > 1 for value in values),
            })
    write_tsv(
        RESULTS / "scale_paired_speedups.tsv",
        scale_rows,
        ["scale_label", "comparison", "n_paired_blocks", "median_paired_ratio", "q1_paired_ratio",
         "q3_paired_ratio", "numerator_faster_count", "denominator_faster_count"],
    )

    output_rows = []
    for method in ("autozyme", "ess"):
        method_data = [row for row in rows if row["method"] == method]
        output_rows.append({
            "method": method,
            "n_runs": len(method_data),
            "pvalues_elementwise_exact_count": sum(row["pval_exact"] == "TRUE" for row in method_data),
            "pvalues_elementwise_exact_fraction": sum(row["pval_exact"] == "TRUE" for row in method_data) / len(method_data),
            "maximum_probability_absolute_difference": max(float(row["prob_max_abs_diff"]) for row in method_data),
            "maximum_pvalue_absolute_difference": max(float(row["pval_max_abs_diff"]) for row in method_data),
        })
    write_tsv(
        RESULTS / "numerical_equivalence_summary.tsv",
        output_rows,
        ["method", "n_runs", "pvalues_elementwise_exact_count", "pvalues_elementwise_exact_fraction",
         "maximum_probability_absolute_difference", "maximum_pvalue_absolute_difference"],
    )

    print(f"Validated {len(rows)} successful runs, {len(runs)} complete paired blocks, 12 datasets.")
    print(f"Wrote summary tables to {RESULTS}")
    for row in paired_rows:
        print(row)


if __name__ == "__main__":
    main()
