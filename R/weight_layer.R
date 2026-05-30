# Experiment A: the conditional weight layer.
# Holds the allocation table fixed and runs the (m | n, beta) -> (beta | m)
# Gibbs recursion, isolating the layer the FF theory describes.

# internal: fixed occupancy table plus the group-balanced MSE target
.build_occupancy <- function(config) {
  N <- rbind(config$dom_counts,
             matrix(rep(config$sat_counts, config$n_satellite),
                    config$n_satellite, config$K, byrow = TRUE))
  rownames(N) <- NULL
  props <- N / rowSums(N)
  list(N = N, beta0_target = colMeans(props), kstar = config$kstar)
}

# internal: one Gibbs chain at fixed occupancy N and given eta; returns the
# matrix of retained beta draws (rows = iterations kept, cols = components).
.weight_layer_chain <- function(N, eta, config, seed) {
  set.seed(seed)
  K <- config$K; J <- nrow(N)
  beta <- rep(1 / K, K)
  keep <- matrix(0, config$iters_w - config$burn_w, K); ki <- 0L
  for (it in seq_len(config$iters_w)) {
    m_dot <- numeric(K)
    for (j in seq_len(J)) for (k in seq_len(K)) {
      njk <- N[j, k]
      if (njk > 0)
        m_dot[k] <- m_dot[k] + rcrt(g_eta(njk, eta), config$alpha0_w * beta[k])
    }
    beta <- .rdirichlet1(config$gamma_w / K + m_dot)   # FF acts via m_dot above
    if (it > config$burn_w) { ki <- ki + 1L; keep[ki, ] <- beta }
  }
  keep
}

#' Weight-layer experiment (conditional on a fixed allocation table)
#'
#' Runs the table-count to global-weight Gibbs recursion on a fixed,
#' size-imbalanced occupancy table over the eta grid, averaged over replicate
#' chains. This isolates the layer the FF theory conditions on and exhibits the
#' eta^{-1/2} inflation of the posterior standard deviation of the
#' satellite-borne weight.
#'
#' @param config A configuration list, see \code{\link{ff_default_config}}.
#' @param verbose Logical; print per-step progress messages.
#' @return A list with components \code{summary} (a data frame with one row per
#'   eta: posterior mean, posterior SD and its SE, the group-balanced target,
#'   weight MSE and its SE, the FF/BGS SD ratio and its SE, and the
#'   \eqn{\eta^{-1/2}} reference) and \code{draws} (a named list of the
#'   replicate-1 posterior draws of the satellite weight, for density plots).
#' @examples
#' \donttest{
#' cfg <- ff_default_config(quick = TRUE)
#' res <- ff_weight_layer(cfg, verbose = FALSE)
#' res$summary
#' }
#' @export
ff_weight_layer <- function(config = ff_default_config(), verbose = TRUE) {
  occ <- .build_occupancy(config); ks <- occ$kstar
  E <- length(config$eta_grid); base_i <- which(abs(config$eta_grid - 1) < 1e-9)
  raw <- list(); ri <- 0L; draws_store <- list()
  for (rep in seq_len(config$n_rep)) {
    sdv <- numeric(E); block <- vector("list", E)
    for (ei in seq_len(E)) {
      eta <- config$eta_grid[ei]
      if (verbose) message(sprintf("  [weight layer] rep %d/%d  eta = %.2f",
                                   rep, config$n_rep, eta))
      draws <- .weight_layer_chain(occ$N, eta, config,
                                  seed = config$seed + 1000 * rep + ei)
      wk <- draws[, ks]; pm <- colMeans(draws)
      sdv[ei] <- stats::sd(wk)
      block[[ei]] <- data.frame(rep = rep, eta = eta, mean = mean(wk),
                                sd = stats::sd(wk),
                                mse = mean((pm - occ$beta0_target)^2))
      if (rep == 1) draws_store[[sprintf("%.2f", eta)]] <- wk
    }
    d <- do.call(rbind, block)
    d$sd_ratio <- d$sd / sdv[base_i]   # within-replicate matched baseline
    ri <- ri + 1L; raw[[ri]] <- d
  }
  raw <- do.call(rbind, raw); nr <- config$n_rep
  etas <- sort(unique(raw$eta), decreasing = TRUE)
  agg1  <- function(col) vapply(etas, function(e) mean(raw[[col]][raw$eta == e]), numeric(1))
  aggse <- function(col) vapply(etas, function(e) .se_of(raw[[col]][raw$eta == e], nr), numeric(1))
  summary <- data.frame(
    eta                   = etas,
    beta_kstar_post_mean  = agg1("mean"),
    beta_kstar_post_SD    = agg1("sd"),
    beta_kstar_post_SD_se = aggse("sd"),
    beta_kstar_target     = occ$beta0_target[ks],
    weight_MSE            = agg1("mse"),
    weight_MSE_se         = aggse("mse"),
    SD_ratio_vs_eta1      = agg1("sd_ratio"),
    SD_ratio_se           = aggse("sd_ratio"),
    eta_minus_half        = etas^(-0.5))
  list(summary = summary, draws = draws_store)
}
