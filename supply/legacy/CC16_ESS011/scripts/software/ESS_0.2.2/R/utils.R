subset_indexed_matrix <- function(data, index_list, error_message) {
  if (!is.list(index_list) || !length(index_list)) {
    stop("index_list must be a non-empty list.", call. = FALSE)
  }

  used_genes <- unlist(index_list, use.names = FALSE)
  used_genes <- sort(unique(used_genes[!is.na(used_genes) & used_genes > 0L]))
  if (!length(used_genes)) stop(error_message, call. = FALSE)

  old_to_new <- integer(nrow(data))
  old_to_new[used_genes] <- seq_along(used_genes)
  remap <- function(index) {
    positive <- !is.na(index) & index > 0L
    index[positive] <- old_to_new[index[positive]]
    index
  }

  list(
    data = as.matrix(data[used_genes, , drop = FALSE]),
    indices = lapply(index_list, remap),
    genes = used_genes
  )
}

make_index_matrix <- function(items, data_genes, complex_input = NULL,
                              cofactor_input = NULL, cofactor_cols = NULL,
                              keep_duplicates = FALSE) {
  pos <- match(data_genes, data_genes)
  names(pos) <- data_genes
  sets <- vector("list", length(items))

  for (i in seq_along(items)) {
    x <- items[[i]]
    genes <- character()
    if (!is.na(x) && nzchar(x)) {
      if (!is.null(complex_input) && x %in% rownames(complex_input)) {
        subunits <- unlist(complex_input[x, grep("^subunit", colnames(complex_input)), drop = FALSE],
                           use.names = FALSE)
        genes <- subunits[nzchar(subunits)]
      } else if (!is.null(cofactor_input) && x %in% rownames(cofactor_input)) {
        subunits <- unlist(cofactor_input[x, cofactor_cols, drop = FALSE], use.names = FALSE)
        genes <- subunits[nzchar(subunits)]
      } else {
        genes <- x
      }
    }
    idx <- if (keep_duplicates) {
      unname(pos[genes])
    } else {
      unname(pos[intersect(genes, data_genes)])
    }
    sets[[i]] <- idx[!is.na(idx)]
  }

  max_len <- max(1L, max(lengths(sets)))
  out <- matrix(0L, nrow = length(items), ncol = max_len)
  for (i in seq_along(sets)) {
    if (length(sets[[i]])) out[i, seq_along(sets[[i]])] <- as.integer(sets[[i]])
  }
  out
}

validate_nonnegative_finite_expression <- function(data) {
  values <- if (inherits(data, "sparseMatrix")) {
    if ("x" %in% methods::slotNames(data)) methods::slot(data, "x") else numeric()
  } else {
    data
  }
  if (any(!is.finite(values))) {
    stop(
      "Signaling expression matrix must contain only finite values; NA, NaN, and Inf are unsupported.",
      call. = FALSE
    )
  }
  if (any(values < 0)) {
    stop(
      "Signaling expression matrix must be non-negative; support-set pruning assumes expression values >= 0.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

check_supported_object <- function(object, type, population.size, distance.use) {
  if (!identical(type, "triMean")) {
    stop("CellChatESS currently supports only type = 'triMean'.", call. = FALSE)
  }
  if (!identical(population.size, FALSE)) {
    stop("CellChatESS currently supports only population.size = FALSE.", call. = FALSE)
  }
  if (!is.null(distance.use) && !identical(distance.use, FALSE)) {
    stop("CellChatESS currently supports only non-spatial RNA workflows with distance.use = NULL/FALSE.", call. = FALSE)
  }
  if (!identical(object@options$datatype, "RNA")) {
    stop("CellChatESS currently supports only object@options$datatype == 'RNA'.", call. = FALSE)
  }
}
