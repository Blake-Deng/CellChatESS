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
  direct_gene <- rep(TRUE, length(items))

  for (i in seq_along(items)) {
    x <- items[[i]]
    genes <- character()
    if (!is.na(x) && nzchar(x)) {
      if (!is.null(complex_input) && x %in% data_genes) {
        # CellChat gives measured gene names precedence over complex names.
        genes <- x
      } else if (!is.null(complex_input) && x %in% rownames(complex_input)) {
        direct_gene[i] <- FALSE
        subunits <- unlist(complex_input[x, grep("^subunit", colnames(complex_input)), drop = FALSE],
                           use.names = FALSE)
        genes <- as.character(subunits[!is.na(subunits) & nzchar(subunits)])
        missing <- setdiff(genes, data_genes)
        if (!length(genes) || length(missing)) {
          stop("Complex '", x, "' has no subunits or has missing expression rows: ",
               paste(missing, collapse = ", "), ".", call. = FALSE)
        }
      } else if (!is.null(complex_input)) {
        stop("Ligand/receptor '", x, "' is neither a measured gene nor a defined complex.",
             call. = FALSE)
      } else if (!is.null(cofactor_input) && x %in% rownames(cofactor_input)) {
        subunits <- unlist(cofactor_input[x, cofactor_cols, drop = FALSE], use.names = FALSE)
        genes <- subunits[!is.na(subunits) & nzchar(subunits)]
      } else if (is.null(cofactor_input)) {
        genes <- x
      }
    } else if (!is.null(complex_input)) {
      stop("Ligand/receptor names must be non-empty and non-missing.", call. = FALSE)
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
  if (!is.null(complex_input)) attr(out, "direct_gene") <- direct_gene
  out
}

validate_inference_input <- function(object, data, raw.use, nboot, seed.use, Kh, n) {
  scalar_flag <- function(x) is.logical(x) && length(x) == 1L && !is.na(x)
  if (!scalar_flag(raw.use)) stop("raw.use must be TRUE or FALSE.", call. = FALSE)
  whole <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x) &&
    x == floor(x) && abs(x) <= .Machine$integer.max
  if (!whole(nboot) || nboot < 1) stop("nboot must be a positive integer.", call. = FALSE)
  if (!whole(seed.use)) stop("seed.use must be a finite integer seed.", call. = FALSE)
  if (!is.numeric(Kh) || length(Kh) != 1L || !is.finite(Kh) || Kh <= 0 ||
      !is.numeric(n) || length(n) != 1L || !is.finite(n) || n <= 0) {
    stop("Kh and n must be finite positive numbers.", call. = FALSE)
  }
  if (length(dim(data)) != 2L || any(dim(data) == 0L) ||
      !(is.matrix(data) && is.numeric(data) || methods::is(data, "dMatrix"))) {
    stop("Expression data must be a non-empty numeric matrix or numeric Matrix object.", call. = FALSE)
  }
  genes <- rownames(data)
  if (is.null(genes) || anyNA(genes) || any(!nzchar(genes)) || anyDuplicated(genes)) {
    stop("Expression gene names must be present, non-empty, and unique.", call. = FALSE)
  }
  group <- object@idents
  if (!is.factor(group) || length(group) != ncol(data) || anyNA(group)) {
    stop("object@idents must be a non-missing factor with one label per expression column.", call. = FALSE)
  }
  if (!is.null(names(group)) && !is.null(colnames(data)) &&
      !identical(names(group), colnames(data))) {
    stop("Cell identities and expression columns must have the same cell order.", call. = FALSE)
  }
  if (nlevels(group) == 0L || nlevels(group) != length(unique(group))) {
    stop("Unused factor levels found in object@idents. Drop unused levels before running.", call. = FALSE)
  }
  if (as.double(nlevels(group))^2 > .Machine$integer.max) {
    stop("The current kernel requires K squared to fit in a signed 32-bit integer.", call. = FALSE)
  }
  fields <- c("ligand", "receptor", "agonist", "antagonist", "co_A_receptor", "co_I_receptor")
  if (!all(fields %in% colnames(object@LR$LRsig))) {
    stop("object@LR$LRsig lacks required ligand/receptor or modifier columns.", call. = FALSE)
  }
  invisible(TRUE)
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
