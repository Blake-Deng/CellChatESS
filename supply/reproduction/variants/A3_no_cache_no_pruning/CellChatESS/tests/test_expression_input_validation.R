library(CellChatESS)

validate <- getFromNamespace("validate_nonnegative_finite_expression", "CellChatESS")

expect_rejection <- function(x, expected_text) {
  error <- tryCatch({
    validate(x)
    NULL
  }, error = function(e) e)
  if (!inherits(error, "error")) {
    stop("Expected input validation to reject this matrix.")
  }
  if (!grepl(expected_text, conditionMessage(error), fixed = TRUE)) {
    stop("Unexpected validation message: ", conditionMessage(error))
  }
  invisible(TRUE)
}

stopifnot(isTRUE(validate(matrix(c(0, 1, 2, 0), nrow = 2L))))
expect_rejection(matrix(c(1, -0.1), nrow = 1L), "must be non-negative")
expect_rejection(matrix(c(1, NA_real_), nrow = 1L), "only finite values")
expect_rejection(matrix(c(1, NaN), nrow = 1L), "only finite values")
expect_rejection(matrix(c(1, Inf), nrow = 1L), "only finite values")
expect_rejection(matrix(c(1, -Inf), nrow = 1L), "only finite values")

if (requireNamespace("Matrix", quietly = TRUE)) {
  sparse_ok <- Matrix::Matrix(matrix(c(0, 1, 2, 0), nrow = 2L), sparse = TRUE)
  sparse_negative <- Matrix::Matrix(matrix(c(0, -0.1, 2, 0), nrow = 2L), sparse = TRUE)
  sparse_nan <- Matrix::Matrix(matrix(c(0, NaN, 2, 0), nrow = 2L), sparse = TRUE)
  stopifnot(isTRUE(validate(sparse_ok)))
  expect_rejection(sparse_negative, "must be non-negative")
  expect_rejection(sparse_nan, "only finite values")
}

cat("Expression input validation tests passed.\n")
