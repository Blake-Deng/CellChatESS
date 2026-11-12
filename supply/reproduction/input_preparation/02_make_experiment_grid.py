#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path


def parse_scales(text):
    out = []
    for item in text.split(","):
        item = item.strip()
        if not item:
            continue
        out.append(item if item == "all" else int(item))
    return out


def read_manifest(path):
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return list(reader)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", default="results/data_manifest.tsv")
    ap.add_argument("--out", default="results/experiment_grid.csv")
    ap.add_argument("--max-datasets", type=int, default=12)
    ap.add_argument("--repeats", type=int, default=3)
    ap.add_argument("--scales", default="1000,5000,10000,25000,50000,all")
    ap.add_argument("--engines", default="both,accelerated")
    ap.add_argument("--ablations", default="no_accel_kernel,no_accel_pathway,no_accel_aggregate")
    ap.add_argument("--accel-algorithms", default="dense")
    args = ap.parse_args()

    manifest = read_manifest(args.manifest)
    # Prefer RDS/Seurat objects because these can be run directly by CellChat.
    rds = [r for r in manifest if r.get("kind") == "rds"]
    # Keep a range of sizes while avoiding duplicate sample ids.
    seen = set()
    selected = []
    for r in sorted(rds, key=lambda x: -int(x["bytes"])):
        key = (r["dataset_root"], r["sample_id"])
        if key in seen:
            continue
        seen.add(key)
        selected.append(r)
        if len(selected) >= args.max_datasets:
            break

    scales = parse_scales(args.scales)
    engines = [x.strip() for x in args.engines.split(",") if x.strip()]
    ablations = [x.strip() for x in args.ablations.split(",") if x.strip()]
    accel_algorithms = [x.strip() for x in args.accel_algorithms.split(",") if x.strip()]

    rows = []
    for ds_i, ds in enumerate(selected, start=1):
        dataset_id = f"{ds['dataset_root']}__{ds['sample_id']}".replace("/", "_").replace(" ", "_")
        for n_cells in scales:
            for rep in range(1, args.repeats + 1):
                for engine in engines:
                    if engine == "accelerated":
                        ablation_list = ablations
                    elif engine == "both":
                        ablation_list = ["full"]
                    else:
                        ablation_list = ["none"]
                    algorithm_list = accel_algorithms if engine in {"accelerated", "both"} else ["dense"]
                    for ablation in ablation_list:
                        for accel_algorithm in algorithm_list:
                            alg_suffix = "" if accel_algorithm == "dense" else f"__algorithm-{accel_algorithm}"
                            exp_id = (
                                f"{dataset_id}__cells-{n_cells}__rep-{rep}"
                                f"__engine-{engine}__ablation-{ablation}{alg_suffix}"
                            )
                            rows.append(
                                {
                                    "experiment_id": exp_id,
                                    "dataset_id": dataset_id,
                                    "input_path": ds["path"],
                                    "dataset_root": ds["dataset_root"],
                                    "sample_id": ds["sample_id"],
                                    "input_kind": ds["kind"],
                                    "input_bytes": ds["bytes"],
                                    "n_cells": n_cells,
                                    "repeat": rep,
                                    "seed": 100000 + ds_i * 1000 + rep,
                                    "engine": engine,
                                    "ablation": ablation,
                                    "accel_algorithm": accel_algorithm,
                                    "label_col": "auto",
                                }
                            )

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as fh:
        fieldnames = [
            "experiment_id", "dataset_id", "input_path", "dataset_root", "sample_id",
            "input_kind", "input_bytes", "n_cells", "repeat", "seed",
            "engine", "ablation", "accel_algorithm", "label_col",
        ]
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"selected_datasets={len(selected)}")
    print(f"grid_rows={len(rows)} out={out}")


if __name__ == "__main__":
    main()

