# CellChatESS 0.2.0

- Adds an optional `threads` argument to `computeCommunProbESS()`.
- Uses a true serial fast path when `threads = 1L`, without entering an
  OpenMP parallel region.
- Uses deterministic OpenMP parallelism across bootstrap permutations when
  `threads > 1L`.
- Records `threads_requested` and `threads_used` in
  `object@options$CellChatESS`.
- Preserves exact probability, p-value, pathway, weight and count outputs
  between one- and four-thread execution across the 216-object validation.

# CellChatESS 0.1.1

- Refines the ESS expansion to Efficient Signal Streaming.
- Avoids broad generality wording and keeps the method framed as computationally efficient within the dense CellChat output boundary.

# CellChatESS 0.1.0

- Introduces Efficient Signal Streaming (ESS).
- Projects signaling data to ligand-receptor, complex and cofactor genes before dense conversion.
- Uses adaptive triMean quantile evaluation.
- Evaluates expression-supported sender-receiver combinations.
- Includes a self-contained PBMC3k example and large-scale validation summaries.
