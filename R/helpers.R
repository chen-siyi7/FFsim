# Core primitives for the Flattened Franchise (FF) construction.
# Base R only.

#' Sublinear occupancy transform g_eta
#'
#' The transform that defines the Flattened Franchise. The group-component
#' occupancy \code{n} is mapped to \code{floor(max(1, n^eta))} before it is fed
#' to the Chinese Restaurant Table (CRT) draw. \code{eta = 1} is the identity
#' (the standard franchise); \code{eta < 1} compresses large occupancies.
#'
#' @param n Non-negative integer occupancy (scalar).
#' @param eta Flattening exponent in (0, 1].
#' @return An integer-valued occupancy to pass to \code{\link{rcrt}}.
#' @examples
#' g_eta(1000, 1)    # 1000
#' g_eta(1000, 0.5)  # 31
#' g_eta(1, 0.5)     # 1
#' @export
g_eta <- function(n, eta) {
  if (eta >= 1 || n <= 1) return(n)
  max(1, floor(n^eta))
}

#' Draw a Chinese Restaurant Table (CRT) variate
#'
#' Samples the number of occupied tables for \code{n} customers at concentration
#' \code{theta}, i.e. \eqn{m = \sum_{i=1}^{n} \mathrm{Bernoulli}(\theta/(\theta+i-1))}.
#' The Bernoulli draws use a uniform comparison rather than vectorized
#' \code{\link[stats]{rbinom}}, which can intermittently return \code{NA} on long
#' probability vectors.
#'
#' @param n Non-negative integer number of customers (scalar).
#' @param theta Positive concentration parameter (scalar).
#' @return A non-negative integer, the table count.
#' @examples
#' set.seed(1)
#' mean(replicate(2000, rcrt(500, 2)))  # ~ 2 * (log(500) + const)
#' @export
rcrt <- function(n, theta) {
  if (!is.finite(n) || n <= 0) return(0L)
  if (n == 1 || !is.finite(theta) || theta <= 0) return(1L)
  i <- seq_len(n)
  sum(stats::runif(n) < theta / (theta + i - 1))
}

#' Adjusted Rand index
#'
#' Permutation-invariant agreement between two label vectors.
#'
#' @param a,b Integer or factor label vectors of equal length.
#' @return The adjusted Rand index (numeric scalar); 1 for identical partitions.
#' @examples
#' adj_rand(c(1, 1, 2, 2), c(2, 2, 1, 1))  # 1
#' @export
adj_rand <- function(a, b) {
  tab <- table(a, b); n <- sum(tab)
  if (n < 2) return(NA_real_)
  c2 <- function(x) sum(x * (x - 1) / 2)
  idx <- c2(as.vector(tab)); ai <- c2(rowSums(tab)); bj <- c2(colSums(tab))
  expd <- ai * bj / c2(n); mx <- 0.5 * (ai + bj)
  if (mx - expd == 0) return(1)
  (idx - expd) / (mx - expd)
}

# internal: single Dirichlet draw from a (possibly unnormalized) shape vector
.rdirichlet1 <- function(shape) {
  shape <- pmax(shape, 1e-8)
  g <- stats::rgamma(length(shape), shape = shape, rate = 1)
  g <- pmax(g, 1e-300)
  g / sum(g)
}

# internal: vectorized categorical sampling via the Gumbel-max trick.
# logp is an n-by-K matrix of unnormalized log-probabilities.
.gumbel_argmax <- function(logp) {
  g <- -log(-log(matrix(stats::runif(length(logp)), nrow(logp), ncol(logp))))
  max.col(logp + g, ties.method = "first")
}

# internal: standard error across replicates (0 when only one replicate)
.se_of <- function(v, nr) {
  s <- stats::sd(v)
  if (is.na(s)) 0 else s / sqrt(nr)
}
