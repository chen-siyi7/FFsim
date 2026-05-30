# Flattened Franchise (FF) sim

Reference implementation of the Flattened Franchise (FF), a sublinear
occupancy transform applied at the Chinese Restaurant Table (CRT) step of a
hierarchical Dirichlet process (HDP) sampler. The single modification that
defines FF is

```
g_eta(n) = floor(max(1, n^eta)),   eta in (0, 1],
```

applied to the group-component occupancy `n` before the CRT draw. `eta = 1`
recovers the standard franchise (the BGS baseline); `eta < 1` compresses the
influence of large, dominant groups on the shared global weights.

## What it does

Two experiments, each averaged over replicates with standard-error bars:

* **Weight layer** (`ff_weight_layer`). Holds the allocation table fixed and
  runs the `(m | n, beta) -> (beta | m)` Gibbs recursion. This isolates the
  layer the FF theory conditions on and exhibits the `eta^(-1/2)` inflation of
  the minority-weight posterior standard deviation.
* **Clustering** (`ff_cluster`). A full finite weak-limit HDP Gibbs sampler on
  grouped Gaussian data, reporting the adjusted Rand index versus `eta`.

`ff_run()` runs both and can write the result tables (CSV) and figures (PDF).

> The numbers are produced by the code. This package is a runnable reference
> implementation of the FF mechanism; it does **not** reproduce any specific
> published table, and its default data-generating process, estimand, and
> sampler variant are meant to be edited to match a given study.

## Quick start

```r

# the defining transform and the CRT primitive
g_eta(1000, 0.6)
set.seed(1); rcrt(500, 2)

# a fast end-to-end run (small settings)
res <- ff_run(ff_default_config(quick = TRUE), outdir = "FF_sim_output")
res$weight_layer
res$clustering

# the full study (several minutes); writes CSVs and PDFs
res <- ff_run(ff_default_config(), outdir = "FF_sim_output")
```

## Aligning with a manuscript

* `ff_default_config()` returns every setting; edit it to change the design,
  the MCMC budget, or the estimand (`build_occupancy`'s group-balanced target).
* The `eta^(-1/2)` curve is the leading-order asymptote; at finite, simulable
  occupancies the observed SD ratio rises monotonically but stays below it.
* The clustering experiment uses an instantiated-weight sampler, so the FF
  table counts enter only the global-weight update. For a tighter coupling
  between `beta` and clustering, replace the assignment step with the collapsed
  direct-assignment predictive.

## License

MIT (c) 2026.
