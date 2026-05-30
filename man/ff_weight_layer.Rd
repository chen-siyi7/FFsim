\name{ff_weight_layer}
\alias{ff_weight_layer}
\title{Weight-layer experiment (conditional on a fixed allocation table)}
\usage{
ff_weight_layer(config = ff_default_config(), verbose = TRUE)
}
\arguments{
\item{config}{A configuration list, see \code{\link{ff_default_config}}.}
\item{verbose}{Logical; print per-step progress messages.}
}
\value{
A list with \code{summary} (one row per eta: posterior mean, posterior SD and
its SE, the group-balanced target, weight MSE and its SE, the FF/BGS SD ratio
and its SE, and the eta^(-1/2) reference) and \code{draws} (replicate-1
posterior draws of the satellite weight, for density plots).
}
\description{
Runs the table-count to global-weight Gibbs recursion on a fixed,
size-imbalanced occupancy table over the eta grid, averaged over replicate
chains. This isolates the layer the FF theory conditions on and exhibits the
eta^(-1/2) inflation of the minority-weight posterior standard deviation.
}
\examples{
\donttest{
res <- ff_weight_layer(ff_default_config(quick = TRUE), verbose = FALSE)
res$summary
}
}
