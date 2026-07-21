# Minimal reproducible example for CellChatESS v0.2.1

This is the repository-integrated version of the 21 July 2026 reviewer demo. It
tests the central conditional-OpenMP claim on a prepared PBMC3k CellChat object:

> `threads = 1L` uses the serial fast path; `threads = 4L` requests OpenMP;
> both paths produce identical CellChat-compatible outputs.

## One-command validation

Prerequisites are R (>= 4.1), Rcpp and CellChat. From the repository root, run:

```sh
Rscript examples/minimal_reproducible_example/run_review_demo.R
```

The script installs the current repository source into the private
`.example_lib/` directory, leaving the user's normal R library unchanged. It
runs the one- and four-thread configurations with `nboot = 100` and
`seed.use = 1L`, compares both with the bundled original CellChat reference,
and writes machine-readable results under `results/`.

Expected final message:

```text
PASS: serial and OpenMP results are identical; both match original CellChat.
```

On a build without OpenMP, the four-thread request safely falls back to one
thread and the numerical checks still run.

## Included files

- `run_review_demo.R`: complete install, execution and validation workflow.
- `data/`: prepared PBMC3k CellChat object and original CellChat reference.
- `expected/acceptance.tsv`: explicit numerical acceptance criteria.
- `results/`: validated v0.2.1 outputs and session information.
- `REVIEWER_CHECKLIST.md`: short code-review route.

The manuscript's principal runtime and memory comparisons remain controlled
one-thread comparisons. This example's four-thread execution is a
resource-sensitivity validation of the optional v0.2.1 implementation.
