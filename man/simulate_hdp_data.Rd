\name{simulate_hdp_data}
\alias{simulate_hdp_data}
\title{Simulate grouped Gaussian data with one dominant group and satellites}
\usage{
simulate_hdp_data(regime, config, seed)
}
\arguments{
\item{regime}{Either \code{"separated"} or \code{"overlapping"}; sets the
spacing between component means.}
\item{config}{A configuration list, see \code{\link{ff_default_config}}.}
\item{seed}{Integer RNG seed.}
}
\value{
A list with \code{x} (observations), \code{group} (group index), \code{label}
(true component), \code{J} (number of groups), and \code{mu_true} (true means).
}
\description{
Generates data from a shared set of \code{K_true} univariate Gaussian
components: one dominant (anchor) group carrying the high-prevalence components
and satellite groups carrying the minority component.
}
\examples{
cfg <- ff_default_config(quick = TRUE)
d <- simulate_hdp_data("separated", cfg, seed = 1)
table(d$group, d$label)
}
