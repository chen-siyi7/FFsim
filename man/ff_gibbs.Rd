\name{ff_gibbs}
\alias{ff_gibbs}
\title{One run of the finite weak-limit HDP Gibbs sampler}
\usage{
ff_gibbs(data, eta, config, seed)
}
\arguments{
\item{data}{Output of \code{\link{simulate_hdp_data}}.}
\item{eta}{Flattening exponent in (0, 1].}
\item{config}{A configuration list, see \code{\link{ff_default_config}}.}
\item{seed}{Integer RNG seed.}
}
\value{
A list with \code{ari_mean} (mean ARI over retained draws) and \code{ari_draws}
(the thinned ARI values).
}
\description{
Runs the full sampler (component locations, instantiated group weights,
allocations, FF table counts, global weights) at a single \code{eta} and
returns the adjusted Rand index against the truth. The FF step is the only
difference from the standard sampler: occupancies pass through \code{\link{g_eta}}
before the CRT draw.
}
\examples{
\donttest{
cfg <- ff_default_config(quick = TRUE)
d <- simulate_hdp_data("separated", cfg, seed = 1)
ff_gibbs(d, eta = 0.6, cfg, seed = 1)$ari_mean
}
}
