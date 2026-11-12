Historical preparation of the 216 frozen CellChat inputs

The July server snapshot of cellchat_benchmark.R is included with its run_grid.R launcher. The grid builder records the selected source object, sample size and seed, with automatic group-label selection. See cellchat_benchmark.R lines 42–81 and 279–293 for the order of sampling and assignment, both label-selection branches, and the saved checkpoint. The original benchmark environment inventory lists the Seurat version. It is a project-level record.

The audit joins the original successful baseline run metrics and historical checkpoint grid to the revised results using experiment_id (the parent directory of the saved checkpoint). Source object paths are relative to the original data share. All 216 historical group counts agree with the revised frozen-object counts. Historical and revised checkpoint identifiers agree for all 216 rows. The SHA256 column comes from the revised run manifest. The two Group datasets each have one source Seurat object in the original run metrics. For Group1 and Group2, the first two 25,000-cell repeats had 25 and 21 groups, respectively, compared with 24 and 20 at full size. Per-job output metrics retain the resulting group count. The original launcher did not separately log which metadata column was selected or record a cluster-branch indicator.

The archived runner accepts a source Seurat object supplied by the user. The original large Seurat inputs and prepared RDS checkpoints are external to this directory. The revision timed inference on the frozen checkpoints.

SHA256 of the archival tables used for the audit
all_metrics.tsv	435710e37049cbda8c2e133fdd9c0267bfae18be2d142ecf7436bd9fbd8a7f37
sparse_stream_checkpoint_grid.csv	8d230bea0f9a26b9317f66c5501d03333253070413e492363565b820302fa1e8
dataset_characteristics_216.tsv	82cc6ae5539588e08f9f179e64bf98fefd33f3894aebfdb0a945328392fcbe28
