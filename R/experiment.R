# Configuration and the top-level driver that runs both experiments and
# (optionally) writes tables and figures.

#' Default configuration for the FF simulation
#'
#' Returns the list of settings consumed by all driver functions. Edit the
#' returned list to change the data-generating process, the estimand, the
#' sampler, or the MCMC budget. Component 3 is the satellite-borne minority
#' component k*.
#'
#' @param quick Logical; if \code{TRUE}, return small settings for a fast
#'   end-to-end check (short chains, few replicates, one regime).
#' @return A named list of configuration values.
#' @examples
#' str(ff_default_config())
#' str(ff_default_config(quick = TRUE))
#' @export
ff_default_config <- function(quick = FALSE) {
  cfg <- list(
    eta_grid = c(1.0, 0.9, 0.8, 0.7, 0.6, 0.5),  # eta = 1 is the BGS baseline
    n_rep    = 10,
    seed     = 1,

    # Experiment A: weight layer (fixed allocation table)
    K           = 3,
    kstar       = 3,
    dom_counts  = c(1000, 1000, 0),   # dominant-group occupancy by component
    sat_counts  = c(2, 2, 17),        # each satellite-group occupancy
    n_satellite = 7,
    alpha0_w    = 5,                  # group concentration (weight layer)
    gamma_w     = 1,
    iters_w     = 4000,
    burn_w      = 1000,

    # Experiment B: full HDP Gibbs on grouped Gaussian data
    run_both    = TRUE,
    K_true      = 3,
    sep         = 5,
    sep_overlap = 2,
    sigma       = 1,
    n_dom       = 2000,
    n_sat       = 25,
    pi_dom      = c(0.50, 0.50, 0.00),
    pi_sat      = c(0.15, 0.15, 0.70),
    K_max       = 4,
    alpha0_c    = 2,
    gamma_c     = 1,
    mu0         = 0, tau0_sq = 100,
    n_iter      = 1500, burnin = 500, thin_ari = 20
  )
  if (quick) {
    cfg$eta_grid   <- c(1.0, 0.7, 0.5)
    cfg$n_rep      <- 2
    cfg$dom_counts <- c(200, 200, 0)
    cfg$iters_w    <- 600;  cfg$burn_w <- 150
    cfg$n_dom      <- 300
    cfg$n_iter     <- 250;  cfg$burnin <- 80
    cfg$run_both   <- FALSE
  }
  cfg
}

# internal: write the two figures (mean +/- SE) to outdir
.ff_plots <- function(A, B, outdir, config) {
  # weight layer
  grDevices::pdf(file.path(outdir, "ff_weight_layer.pdf"), width = 9, height = 8)
  op <- graphics::par(mfrow = c(2, 2), mar = c(4.3, 4.3, 2.6, 1))
  r <- A$summary[order(A$summary$eta), ]
  yl <- range(c(r$SD_ratio_vs_eta1 - r$SD_ratio_se,
                r$SD_ratio_vs_eta1 + r$SD_ratio_se, r$eta_minus_half))
  graphics::plot(r$eta, r$SD_ratio_vs_eta1, type = "b", pch = 19, ylim = yl,
                 xlab = expression(eta), ylab = "SD ratio  FF / BGS", main = "SD inflation")
  .ebars(r$eta, r$SD_ratio_vs_eta1, r$SD_ratio_se)
  graphics::lines(r$eta, r$eta_minus_half, lty = 2, col = "red")
  graphics::legend("topright", bty = "n", lty = c(1, 2), pch = c(19, NA),
                   col = c("black", "red"),
                   legend = c("observed", expression(eta^{-1/2})))
  yl <- range(c(r$beta_kstar_post_SD - r$beta_kstar_post_SD_se,
                r$beta_kstar_post_SD + r$beta_kstar_post_SD_se))
  graphics::plot(r$eta, r$beta_kstar_post_SD, type = "b", pch = 19, ylim = yl,
                 xlab = expression(eta),
                 ylab = expression("posterior SD of " * beta[k^"*"]), main = "Posterior SD")
  .ebars(r$eta, r$beta_kstar_post_SD, r$beta_kstar_post_SD_se)
  yl <- range(c(r$weight_MSE - r$weight_MSE_se, r$weight_MSE + r$weight_MSE_se))
  graphics::plot(r$eta, r$weight_MSE, type = "b", pch = 19, ylim = yl,
                 xlab = expression(eta), ylab = "weight MSE vs target", main = "Weight MSE")
  .ebars(r$eta, r$weight_MSE, r$weight_MSE_se)
  pick <- intersect(c("1.00", "0.80", sprintf("%.2f", min(config$eta_grid))), names(A$draws))
  cols <- c("black", "blue", "red")
  dens <- lapply(pick, function(k) stats::density(A$draws[[k]]))
  graphics::plot(NA, xlim = range(vapply(dens, function(d) range(d$x), numeric(2))),
                 ylim = range(vapply(dens, function(d) range(d$y), numeric(2))),
                 xlab = expression(beta[k^"*"]), ylab = "posterior density",
                 main = "Posterior of satellite weight (rep 1)")
  for (i in seq_along(dens)) graphics::lines(dens[[i]], col = cols[i], lwd = 2)
  graphics::abline(v = r$beta_kstar_target[1], lty = 3)
  graphics::legend("topright", bty = "n", lwd = 2, col = cols[seq_along(pick)],
                   legend = paste0("eta = ", pick))
  graphics::par(op); grDevices::dev.off()

  # clustering
  grDevices::pdf(file.path(outdir, "ff_clustering.pdf"), width = 6, height = 5)
  regs <- unique(B$regime); cols <- c("black", "blue")
  yl <- range(c(B$ARI_mean - B$ARI_se, B$ARI_mean + B$ARI_se))
  graphics::plot(NA, xlim = rev(range(config$eta_grid)), ylim = yl,
                 xlab = expression(eta), ylab = "adjusted Rand index",
                 main = "Clustering recovery vs eta")
  for (i in seq_along(regs)) {
    s <- B[B$regime == regs[i], ]; s <- s[order(-s$eta), ]
    graphics::lines(s$eta, s$ARI_mean, type = "b", pch = 19, col = cols[i])
    .ebars(s$eta, s$ARI_mean, s$ARI_se, col = cols[i])
  }
  graphics::legend("bottomleft", bty = "n", lwd = 1, pch = 19,
                   col = cols[seq_along(regs)], legend = regs)
  grDevices::dev.off()
  invisible(NULL)
}

# internal: vertical error bars
.ebars <- function(x, y, se, col = "black") {
  ok <- is.finite(se) & se > 0
  if (any(ok)) graphics::arrows(x[ok], (y - se)[ok], x[ok], (y + se)[ok],
                                angle = 90, code = 3, length = 0.03, col = col)
}

#' Run the full FF simulation study
#'
#' Runs both experiments over the eta grid and the configured replicates and,
#' optionally, writes the result tables (CSV) and figures (PDF) to a directory.
#' The returned numbers are produced by this code; the function does not
#' reproduce any specific published table.
#'
#' @param config A configuration list, see \code{\link{ff_default_config}}.
#' @param outdir Optional directory to write \code{ff_weight_layer.csv},
#'   \code{ff_clustering.csv} and, if \code{make_plots}, the matching PDFs. The
#'   directory is created if needed. If \code{NULL}, nothing is written.
#' @param make_plots Logical; write the figures when \code{outdir} is given.
#' @param verbose Logical; print progress messages.
#' @return Invisibly, a list with \code{weight_layer} (the Experiment A summary),
#'   \code{clustering} (the Experiment B summary), and \code{draws}.
#' @examples
#' \donttest{
#' res <- ff_run(ff_default_config(quick = TRUE), verbose = FALSE)
#' res$weight_layer
#' }
#' @export
ff_run <- function(config = ff_default_config(), outdir = NULL,
                   make_plots = TRUE, verbose = TRUE) {
  if (verbose) message(sprintf("Experiment A: weight layer (%d replicates)", config$n_rep))
  A <- ff_weight_layer(config, verbose = verbose)
  if (verbose) message(sprintf("Experiment B: full HDP Gibbs (%d replicates)", config$n_rep))
  B <- ff_cluster(config, verbose = verbose)
  if (!is.null(outdir)) {
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
    utils::write.csv(A$summary, file.path(outdir, "ff_weight_layer.csv"), row.names = FALSE)
    utils::write.csv(B, file.path(outdir, "ff_clustering.csv"), row.names = FALSE)
    if (make_plots) .ff_plots(A, B, outdir, config)
    if (verbose) message("Outputs written to: ", normalizePath(outdir))
  }
  invisible(list(weight_layer = A$summary, clustering = B, draws = A$draws))
}
