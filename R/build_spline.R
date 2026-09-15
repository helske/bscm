#' Build the spline matrix for Stan
#'
#' combines the spline basis the difference matrix and sum-to-zero constraint
#' @noRd
build_spline <- function(
  T_total,
  T_pre,
  spline_df = NULL,
  knot_spacing = NULL,
  type = "rw2",
  noncentered = TRUE
) {
  degree <- 3L
  T_pre <- max(T_pre)
  if (!is.null(knot_spacing)) {
    h <- knot_spacing
  } else {
    h <- (T_pre - 1L) / (spline_df - degree)
  }

  n <- ceiling((T_total - 1) / h)
  knots <- 1 + h * seq.int(-degree, n + degree)
  B <- splines::splineDesign(knots, seq_len(T_total), ord = degree + 1)

  d <- diag(ncol(B))
  L <- lower.tri(d, diag = TRUE)
  M <- B %*% L
  if (type == "rw2") {
    M <- M %*% L
    tc <- seq_len(T_total) - 0.5 * (T_pre + 1)
    a <- colSums(M[seq_len(T_pre), , drop = FALSE])
    b <- as.vector(tc[seq_len(T_pre)] %*% M[seq_len(T_pre), , drop = FALSE])
    Q <- qr.Q(qr(cbind(a, b, d)))[, -(1:2), drop = FALSE]
    A <- M %*% Q
  } else {
    a <- colSums(M[seq_len(T_pre), , drop = FALSE])
    Q <- qr.Q(qr(cbind(a, d)))[, -1]
    A <- M %*% Q
  }
  scale <- sqrt(exp(mean(log(rowSums(A[seq_len(T_pre), , drop = FALSE]^2)))))
  A <- A / scale

  list(
    A = A,
    D = ncol(A),
    noncentered_xi = noncentered,
    type = type,
    scale = scale
  )
}
