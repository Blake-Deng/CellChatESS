# Xenium aggregation on saved inference outputs

Each row in `runs.tsv` records one call to `CellChat::aggregateNet(obj, return.object = FALSE)` from CellChat 2.2.0.9001. The method column identifies the producer of the saved probability and p-value tensors. Both tensors had dimensions 2,357 × 2,357 × 627. All six count matrices matched the archived counts exactly.

We loaded the same prepared Xenium object and the corresponding saved tensors before timing. Each call ran serially in a fresh R process pinned to the same logical CPU, with sequential future execution and one BLAS/OpenMP thread. We completed a full garbage collection before starting the timer. We ran CellChat then CellChatESS in rounds 1 and 3 and reversed the order in round 2. The elapsed-time medians were 341.935 s for CellChat 2.2.0.9001 outputs and 341.823 s for CellChatESS 0.2.2 outputs. `Xenium_aggregateNet_3x.tsv` lists the three paired differences and the median elapsed time for each method. The difference column is blank in the median row.

To reproduce one row on Linux, use the pinned R library and provide the prepared input and saved result files listed in the Xenium input and output manifests. Use the same CPU and single-thread settings for every invocation. The six arguments identify the prepared RDS, result RDS, output producer, pair, execution position and output TSV.

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 RCPP_PARALLEL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 R_ENABLE_JIT=3 taskset -c "$CPU" Rscript --vanilla replay.R "$PREPARED" "$RESULT" CellChat_2.2.0.9001 1 1 run1.tsv
```

The parent Xenium folder documents the separate complete-workflow runs used for inference time, whole-process wall time and memory. Published copies of workflow logs and measurement extracts omit the superseded standalone aggregation field. The stage timestamps and resource records retain their measured values.
