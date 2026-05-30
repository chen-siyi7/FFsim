\name{ff_default_config}
\alias{ff_default_config}
\title{Default configuration for the FF simulation}
\usage{
ff_default_config(quick = FALSE)
}
\arguments{
\item{quick}{Logical; if \code{TRUE}, return small settings for a fast
end-to-end check (short chains, few replicates, one regime).}
}
\value{A named list of configuration values.}
\description{
Returns the list of settings consumed by all driver functions. Edit the
returned list to change the data-generating process, the estimand, the sampler,
or the MCMC budget. Component 3 is the satellite-borne minority component k*.
}
\examples{
str(ff_default_config())
str(ff_default_config(quick = TRUE))
}
