# Conditional OpenMP changes from CellChatESS v0.1.1

## Hand-edited implementation files

- `R/accel_cellchat.R`
  - adds `threads = 1L` to `computeCommunProbESS()`;
  - validates the requested thread count;
  - passes it to the ESS C++ kernel;
  - records `threads_requested` and `threads_used`.
- `src/cellchat_prob_cpp.cpp`
  - includes OpenMP conditionally;
  - keeps `threads = 1` on a true serial fast path;
  - uses OpenMP across bootstrap permutations only when `threads > 1`;
  - gives each worker private scratch and rejection buffers;
  - deterministically reduces integer rejection counts;
  - reuses serial scratch buffers across ligand-receptor interactions.
- `src/Makevars` and `src/Makevars.win`
  - add R's portable OpenMP compile and link flags.

## Generated interface files

- `R/RcppExports.R`
- `src/RcppExports.cpp`

These were regenerated after adding the C++ `threads` argument.

## Documentation and metadata

- `DESCRIPTION`
- `README.md`
- `NEWS.md`
- `man/computeCommunProbESS.Rd`

## Validated behavior

- `threads = 1L`: no OpenMP parallel region; performance-neutral relative to
  v0.1.1.
- `threads > 1L`: OpenMP bootstrap parallelism.
- Builds without OpenMP fall back to one thread.
- Across 216 prepared objects, one- and four-thread CellChatESS results were
  exactly identical for probability, p-value, pathway probability, aggregate
  weight and aggregate count.
