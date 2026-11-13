# Reproducing the recorded revision

The package source is at the repository root. `frozen/` preserves the exact server
scripts and input manifests used to build the revision evidence. Those scripts use
the original handoff layout (`software/`, `upstream/`, `reproduction/`, `audit/`);
they must run in a separate workspace, not directly in this source repository.

```sh
python3 supply/verify_supply.py
python3 supply/reproduction/make_workspace.py --output /absolute/path/new-reproduction --fetch-upstream
```

The optional download fetches the three exact upstream archives and checks their
SHA256: CellChat 2.2.0.9001, AutoZyme 0.3.1 and the separately pinned legacy CellChat
1.6.1 reference. The committed ESS candidate archive is already supplied once.
Create the dependency environment from `../environment/server/conda-explicit.txt`,
then follow `frozen/reproduction/README_SERVER.md` and `EXECUTION_PLAN_CN.md` within
the reconstructed workspace. The explicit Conda lock is the recorded Linux
platform environment; it is not a portable macOS installer.

1. Build the pinned packages in fresh separate R libraries and edit
   `reproduction/config/runtime.example.json` with the actual Rscript/library paths.
2. Follow [INPUTS.md](../../INPUTS.md) to request the frozen prepared objects,
   verify their SHA256 values and configure a local path-prefix map. The same
   guide identifies the included PBMC3k example, public source-study links and
   the external source objects required by the archived preparation scripts.
3. Run version preflight, input/structure tests and small three-method smoke tests.
4. `run_unified.py primary216` runs 216 two-method pairs (432 method processes).
   `run_unified.py autozyme192` runs 32 inputs in six method orders (192 blocks,
   576 method processes). These commands default to plan-only; execution requires
   `--execute`. Use a fresh `--output` directory, `--jobs 1` initially, and inspect
   the written execution plan and memory requirements before increasing concurrency.
5. Controlled A1/A3/A4/TRACE/PROFILE derivatives are under `variants/`; their patches
   and source identities distinguish them from the final performance build.
   All primary timings are single-threaded. OpenMP agreement is separate validation.
6. Use `summarize_and_plot.py` for completed fresh primary/AutoZyme campaigns.
   The full-627 paired Xenium scripts are separately preserved in `../xenium/scripts/`;
   their absolute frozen server paths must be adapted explicitly before re-execution.
   `../xenium/run_manifest.json` defines the actual paired scope. Earlier
   `run_server_xenium.py` / `run_xenium_profiled.R` are auxiliary PROFILE tooling and
   do not generate the final paired resource table.

This compact deposit contains all reported numerical comparison records and
measurements, not the multi-gigabyte prepared inputs or full output RDS objects.
Original server paths in records are provenance, not claims that those paths exist
on another machine. Preserve the manifest verification before and after transfer.
