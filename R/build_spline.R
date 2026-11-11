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
  if (h < 1) {
    msg <- if (is.null(knot_spacing)) "smaller df" else "larger knot_spacing"
    warning_(
      c(paste0(
        "Knot spacing of {signif(h, 3)} is less than one time step,",
        " so the spline has more coefficients than time points."
        ),
        i = "Consider using a {msg}."
      )
    )
  }
  n <- ceiling((T_total - 1) / h - 1e-8)
  knots <- 1 + h * seq.int(-degree, n + degree)
  knots[n + degree + 1L] <- max(knots[n + degree + 1L], T_total)
  B <- splines::splineDesign(knots, seq_len(T_total), ord = degree + 1)
  
  L <- lower.tri(diag(ncol(B)), diag = TRUE)
  M <- B %*% L
  if (type == "rw2") {
    M <- M %*% L
    tc <- seq_len(T_total) - 0.5 * (T_pre + 1)
    X <- cbind(1, tc)
  } else {
    X <- matrix(1, T_total, 1L)
  }
  pre <- seq_len(T_pre)
  X_pre <- X[pre, , drop = FALSE]
  A <- M[, -seq_len(ncol(X)), drop = FALSE]
  A <- A - X %*% solve(crossprod(X_pre), crossprod(X_pre, A[pre, , drop = FALSE]))
  
  scale <- sqrt(exp(mean(log(rowSums(A[pre, , drop = FALSE]^2)))))
  A <- A / scale
  
  list(
    A = A,
    D = ncol(A),
    noncentered_xi = noncentered,
    type = type,
    scale = scale
  )
}
