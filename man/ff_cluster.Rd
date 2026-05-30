\name{ff_cluster}
\alias{ff_cluster}
\title{Clustering experiment (ARI versus eta)}
\usage{
ff_cluster(config = ff_default_config(), verbose = TRUE)
}
\arguments{
\item{config}{A configuration list, see \code{\link{ff_default_config}}.}
\item{verbose}{Logical; print per-step progress messages.}
}
\value{
A data frame with columns \code{regime}, \code{eta}, \code{ARI_mean}, and
\code{ARI_se}.
}
\description{
Sweeps \code{eta} for one or both separation regimes, averaging the adjusted
Rand index over replicate datasets.
}
\examples{
\donttest{
ff_cluster(ff_default_config(quick = TRUE), verbose = FALSE)
}
}
