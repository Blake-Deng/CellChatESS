# Five-minute reviewer checklist

1. From the repository root, run
   `Rscript examples/minimal_reproducible_example/run_review_demo.R`.
2. Confirm the final line starts with `PASS`.
3. Open `results/reviewer_demo_checks.tsv`; every logical check should be
   `TRUE`.
4. Open `R/accel_cellchat.R`; inspect the `threads = 1L`
   interface and recorded `threads_requested/threads_used`.
5. Open `src/cellchat_prob_cpp.cpp`; search for
   `if (threads_used == 1)` and `#pragma omp parallel`.

The serial and OpenMP branches call the same `process_bootstrap` worker. The
review demo checks exact equality between those branches and comparison against
an independently generated original CellChat reference.
