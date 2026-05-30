# Base-R tests (run directly by R CMD check; no testthat dependency required).
library(ffsim)

## g_eta
stopifnot(g_eta(1000, 1)   == 1000)
stopifnot(g_eta(1000, 0.5) == 31)      # floor(sqrt(1000))
stopifnot(g_eta(1, 0.5)    == 1)
stopifnot(g_eta(0, 0.5)    == 0)

## rcrt: degenerate guards and no NA on long vectors
stopifnot(rcrt(0, 2) == 0)
stopifnot(rcrt(1, 2) == 1)
stopifnot(rcrt(5, -1) == 1)            # non-positive concentration -> single table
set.seed(1)
v <- replicate(500, rcrt(800, 1.5))
stopifnot(!any(is.na(v)), all(v >= 1), all(v <= 800))

## adj_rand
stopifnot(abs(adj_rand(c(1, 1, 2, 2), c(2, 2, 1, 1)) - 1) < 1e-12)
stopifnot(is.finite(adj_rand(c(1, 2, 3, 4), c(1, 1, 2, 2))))

## internal Dirichlet draw is a valid probability vector
s <- ffsim:::.rdirichlet1(c(1, 2, 3))
stopifnot(abs(sum(s) - 1) < 1e-8, all(s >= 0))

## configuration
cfg <- ff_default_config(quick = TRUE)
stopifnot(is.list(cfg), cfg$K == 3, length(cfg$eta_grid) == 3)

## data generation
d <- simulate_hdp_data("separated", cfg, seed = 1)
stopifnot(length(d$x) == cfg$n_dom + cfg$n_satellite * cfg$n_sat)
stopifnot(all(d$label %in% seq_len(cfg$K_true)))

## weight layer: shape, matched baseline, finite SD ratios
wl <- ff_weight_layer(cfg, verbose = FALSE)
stopifnot(nrow(wl$summary) == length(cfg$eta_grid))
base <- wl$summary$SD_ratio_vs_eta1[wl$summary$eta == 1]
stopifnot(abs(base - 1) < 1e-12)
stopifnot(all(is.finite(wl$summary$SD_ratio_vs_eta1)),
          all(wl$summary$SD_ratio_vs_eta1 > 0))
stopifnot(abs(wl$summary$eta_minus_half[wl$summary$eta == 1] - 1) < 1e-12)

## one Gibbs run returns an ARI no greater than 1
a <- ff_gibbs(d, eta = 0.6, cfg, seed = 1)$ari_mean
stopifnot(is.finite(a), a <= 1 + 1e-9)

## full driver writes the expected files
out <- file.path(tempdir(), "ffsim_test_out")
res <- ff_run(cfg, outdir = out, make_plots = TRUE, verbose = FALSE)
stopifnot(file.exists(file.path(out, "ff_weight_layer.csv")),
          file.exists(file.path(out, "ff_clustering.csv")),
          file.exists(file.path(out, "ff_weight_layer.pdf")),
          file.exists(file.path(out, "ff_clustering.pdf")))
stopifnot(is.data.frame(res$weight_layer), is.data.frame(res$clustering))

cat("all ffsim tests passed\n")
