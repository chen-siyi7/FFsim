# Experiment B: the full finite weak-limit HDP Gibbs sampler on grouped
# Gaussian data, reporting clustering recovery (adjusted Rand index) vs eta.

#' Simulate grouped Gaussian data with one dominant group and satellites
#'
#' Generates data from a shared set of \code{K_true} univariate Gaussian
#' components. One dominant (anchor) group carries the high-prevalence
#' components; the satellite groups carry the minority component.
#'
#' @param regime Either \code{"separated"} or \code{"overlapping"}; sets the
#'   spacing between component means.
#' @param config A configuration list, see \code{\link{ff_default_config}}.
#' @param seed Integer RNG seed.
#' @return A list with \code{x} (observations), \code{group} (group index),
#'   \code{label} (true component), \code{J} (number of groups), and
#'   \code{mu_true} (true component means).
#' @examples
#' cfg <- ff_default_config(quick = TRUE)
#' d <- simulate_hdp_data("separated", cfg, seed = 1)
#' table(d$group, d$label)
#' @export
simulate_hdp_data <- function(regime, config, seed) {
  set.seed(seed); K <- config$K_true
  sp <- if (regime == "overlapping") config$sep_overlap else config$sep
  mu_true <- seq(-sp, sp, length.out = K)
  sizes <- c(config$n_dom, rep(config$n_sat, config$n_satellite)); J <- length(sizes)
  Pi <- rbind(config$pi_dom, matrix(config$pi_sat, config$n_satellite, K, byrow = TRUE))
  x <- numeric(0); group <- integer(0); label <- integer(0)
  for (j in seq_len(J)) {
    zj <- sample.int(K, sizes[j], replace = TRUE, prob = Pi[j, ])
    x <- c(x, stats::rnorm(sizes[j], mu_true[zj], config$sigma))
    group <- c(group, rep(j, sizes[j])); label <- c(label, zj)
  }
  list(x = x, group = group, label = label, J = J, mu_true = mu_true)
}

#' One run of the finite weak-limit HDP Gibbs sampler
#'
#' Runs the full sampler (component locations, instantiated group weights,
#' allocations, FF table counts, global weights) at a single \code{eta} and
#' returns the adjusted Rand index against the true allocations. The FF step is
#' the only difference from the standard sampler: occupancies are passed through
#' \code{\link{g_eta}} before the CRT draw.
#'
#' @param data Output of \code{\link{simulate_hdp_data}}.
#' @param eta Flattening exponent in (0, 1].
#' @param config A configuration list, see \code{\link{ff_default_config}}.
#' @param seed Integer RNG seed.
#' @return A list with \code{ari_mean} (mean ARI over retained draws) and
#'   \code{ari_draws} (the thinned ARI values).
#' @examples
#' \donttest{
#' cfg <- ff_default_config(quick = TRUE)
#' d <- simulate_hdp_data("separated", cfg, seed = 1)
#' ff_gibbs(d, eta = 0.6, cfg, seed = 1)$ari_mean
#' }
#' @export
ff_gibbs <- function(data, eta, config, seed) {
  set.seed(seed)
  x <- data$x; group <- data$group; N <- length(x); J <- data$J
  K <- config$K_max; sig2 <- config$sigma^2
  km <- suppressWarnings(stats::kmeans(x, centers = K, nstart = 5, iter.max = 50))
  theta <- as.numeric(km$centers); z <- km$cluster
  beta <- rep(1 / K, K)
  glist <- split(seq_len(N), group); ari <- numeric(0); keep <- 0L
  for (it in seq_len(config$n_iter)) {
    for (k in seq_len(K)) {
      idx <- which(z == k); nk <- length(idx)
      if (nk > 0) {
        prec <- 1 / config$tau0_sq + nk / sig2
        mn <- (config$mu0 / config$tau0_sq + sum(x[idx]) / sig2) / prec
        theta[k] <- stats::rnorm(1, mn, sqrt(1 / prec))
      } else theta[k] <- stats::rnorm(1, config$mu0, sqrt(config$tau0_sq))
    }
    nmat <- matrix(0L, J, K)
    tj <- table(group, factor(z, levels = seq_len(K)))
    nmat[as.integer(rownames(tj)), ] <- as.matrix(tj)
    loglik <- outer(x, theta, function(xx, th) stats::dnorm(xx, th, config$sigma, log = TRUE))
    for (j in seq_len(J)) {
      pij <- .rdirichlet1(config$alpha0_c * beta + nmat[j, ])
      rows <- glist[[j]]
      lp <- sweep(loglik[rows, , drop = FALSE], 2, log(pij + 1e-300), "+")
      z[rows] <- .gumbel_argmax(lp)
    }
    nmat <- matrix(0L, J, K)
    tj <- table(group, factor(z, levels = seq_len(K)))
    nmat[as.integer(rownames(tj)), ] <- as.matrix(tj)
    m_dot <- numeric(K)
    for (j in seq_len(J)) for (k in seq_len(K)) {
      njk <- nmat[j, k]
      if (njk > 0) m_dot[k] <- m_dot[k] + rcrt(g_eta(njk, eta), config$alpha0_c * beta[k])
    }
    beta <- .rdirichlet1(config$gamma_c / K + m_dot)
    if (it > config$burnin) {
      keep <- keep + 1L
      if (keep %% config$thin_ari == 0) ari <- c(ari, adj_rand(z, data$label))
    }
  }
  list(ari_mean = mean(ari, na.rm = TRUE), ari_draws = ari)
}

#' Clustering experiment (ARI versus eta)
#'
#' Sweeps \code{eta} for one or both separation regimes, averaging the adjusted
#' Rand index over replicate datasets.
#'
#' @param config A configuration list, see \code{\link{ff_default_config}}.
#' @param verbose Logical; print per-step progress messages.
#' @return A data frame with columns \code{regime}, \code{eta}, \code{ARI_mean},
#'   and \code{ARI_se}.
#' @examples
#' \donttest{
#' cfg <- ff_default_config(quick = TRUE)
#' ff_cluster(cfg, verbose = FALSE)
#' }
#' @export
ff_cluster <- function(config = ff_default_config(), verbose = TRUE) {
  regimes <- if (isTRUE(config$run_both)) c("separated", "overlapping") else "separated"
  raw <- list(); r <- 0L
  for (rep in seq_len(config$n_rep)) for (regime in regimes) {
    dat <- simulate_hdp_data(regime, config, seed = config$seed + rep)
    for (ei in seq_along(config$eta_grid)) {
      eta <- config$eta_grid[ei]
      if (verbose) message(sprintf("  [clustering: %s] rep %d/%d  eta = %.2f",
                                   regime, rep, config$n_rep, eta))
      a <- ff_gibbs(dat, eta, config, seed = config$seed + 1000 * rep + ei)$ari_mean
      r <- r + 1L
      raw[[r]] <- data.frame(rep = rep, regime = regime, eta = eta, ARI = a)
    }
  }
  raw <- do.call(rbind, raw); nr <- config$n_rep
  keys <- unique(raw[, c("regime", "eta")])
  keys <- keys[order(keys$regime, -keys$eta), ]
  keys$ARI_mean <- mapply(function(rg, e) mean(raw$ARI[raw$regime == rg & raw$eta == e]),
                          keys$regime, keys$eta)
  keys$ARI_se <- mapply(function(rg, e) .se_of(raw$ARI[raw$regime == rg & raw$eta == e], nr),
                        keys$regime, keys$eta)
  rownames(keys) <- NULL
  keys
}
