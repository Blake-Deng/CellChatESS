# CellChatESS 0.2.2

CellChatESS implements Efficient Signal Streaming for CellChat inference. It retains sparse expression data until the required genes have been selected, reuses bootstrap expression summaries and evaluates sender–receiver pairs supported by ligand and receptor expression. The supported workflow is single-dataset, non-spatial RNA inference with `type = "triMean"` and `population.size = FALSE`.

## Installation

Install CellChat and Rcpp first. The revision comparisons use CellChat 2.2.0.9001, pinned to upstream commit `75253cd0c9e68410e6e721a6d3a0419a1d7e358f`. Exact source archives and environment specifications are described in the [reproduction guide](supply/reproduction/README.md).

Install the released version from GitHub with

```r
install.packages("remotes")
remotes::install_github("Blake-Deng/CellChatESS", ref = "v0.2.2")
```

A downloaded R source package can also be installed with

```r
install.packages("CellChatESS_0.2.2.tar.gz", repos = NULL, type = "source")
```

Run inference on a prepared CellChat object with

```r
library(CellChat)
library(CellChatESS)

result <- computeCommunProbESS(
    object,
    type = "triMean",
    raw.use = TRUE,
    population.size = FALSE,
    nboot = 100L,
    seed.use = 1L,
    Kh = 0.5,
    n = 1,
    algorithm = "ess",
    threads = 1L
)
```

The input signaling matrix must be numeric, finite and nonnegative, with at least one positive value. Dense matrices and `dgCMatrix` inputs are accepted. The [supported workflows](SUPPORTED_WORKFLOWS.md) and [methods](METHODS.md) describe the input checks and calculation rules.

`algorithm = "dense"` selects the public dense Rcpp comparison path used for ablation stage A2. Optional OpenMP bootstrap parallelism is available through `threads`. All main performance comparisons use one thread. The output stores probabilities and p-values for every ligand–receptor interaction and sender–receiver group pair, requiring O(n_LR K²) space per array.

## PBMC3k example

The [prepared PBMC3k example](examples/minimal_reproducible_example/README.md) includes a small CellChat object and a script that compares original CellChat with one-thread and four-thread ESS runs.

```sh
Rscript examples/minimal_reproducible_example/run_review_demo.R
```

The script reports the installed CellChat version and actual ESS thread count. A compiler with OpenMP support is needed to exercise four native threads.

## Validation and reproducibility

The revision compares original CellChat 2.2.0.9001, AutoZyme 0.3.1 and CellChatESS 0.2.2 under the same input, parameter, thread and timing settings. The recorded CellChat 1.6.1 / ESS 0.1.1 experiments document the original development environment.

The [supply directory](supply/README.md) contains the original manuscript measurements, input hashes, source identities, environment specifications, numerical checks, plotting scripts and figure data. [Supplementary Data S7b](supply/data/Supplementary_Data_S7b.tsv) provides all 216 benchmark records. [Figure 1](supply/figures/Figure1.pdf) and [Figure S1](supply/figures/FigureS1.pdf) correspond to those manuscript records.

A complete independent replay finished on 12 November 2026. Its [report](supply/independent_reproduction/20261111_12/results/crosscheck/REPORT.md), method measurements, regenerated figures and input reconstruction checks are supplied alongside the manuscript records. All 216 main pairs, 192 three-method blocks, the full 627-LR Xenium comparison and 216 legacy pairs completed. The 129 distinct inputs were regenerated and matched in their ordered cells, groups, expression matrices, databases and LR tables. Resaving these objects produced different file hashes, which are documented in the reconstruction ledger.

In the routine and full-Xenium comparisons, empirical p-values and aggregate counts matched the original implementation element by element. Maximum absolute probability differences were 1.11 × 10⁻¹⁶ and 2.22 × 10⁻¹⁶, respectively. Measurements from the two runs are listed together in the independent replay report, including their different observed runtimes and peak memory.

The tested computational source is recorded by commit `139a61c33d4c09568f48f4afd0a16dc57bb22f4f` and numeric source ID `bb3bb6310b4af2d5ec000ace97d61339680b0e2a940508da3e072bb553f8417b`. Its [source manifest](supply/provenance/candidate_source_manifest.tsv), [original source archive](supply/provenance/archives/CellChatESS_0.2.2_source.tar.gz) and [Git bundle](supply/provenance/archives/CellChatESS_0.2.2_local.bundle) identify that snapshot. Release preparation updates the documentation and adds the independent replay evidence.

Verify the deposited files and the computational source with

```sh
python3 supply/verify_supply.py
```

Large raw inputs, prepared expression objects and full output tensors require separate retrieval. The [reproduction guide](supply/reproduction/README.md) records their source locations and expected hashes and explains how to reconstruct a working directory.

## License

CellChatESS is released under GPL-3. See [LICENSE](LICENSE).
