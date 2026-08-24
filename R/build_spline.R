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
    noncentered = TRUE,
    scale = "gm"
) {
  T_pre <- max(T_pre)
  if (!is.null(knot_spacing)) {
    h <- knot_spacing
  } else {
    h <- (T_pre - 1) / (spline_df - 3)
  }
  lb  <- 1 - 5 * h
  rb <- 1 + (ceiling((T_total - 1) / h) + 5) * h
  knots <- seq(lb + h, rb - h, by = h)
  
  B <- splines::bs(
    seq_len(T_total), knots = knots, intercept = TRUE, 
    Boundary.knots = c(lb, rb)
  )
  B <- B[, 5:(ncol(B) - 4), drop = FALSE]
  D <- ncol(B)
  d <- diag(D)
  L <- lower.tri(d, diag = TRUE)
  M <- B %*% L
  if (type == "rw2") {
    M <- M %*% L
  }
  a <- colSums(M[seq_len(T_pre), , drop = FALSE])
  Q <- qr.Q(qr(cbind(a, d)))[, -1]
  A <- M %*% Q
  if (scale == "gm") {
    A <- A / sqrt(exp(mean(log(rowSums(A[seq_len(T_pre), , drop = FALSE]^2)))))
  } else {
    A <- A / sqrt(mean(rowSums(A[seq_len(T_pre), , drop = FALSE]^2)))
  }
  list(
    A = A,
    D = ncol(A),
    noncentered_xi = noncentered,
    type = type,
    scale = scale
  )
}