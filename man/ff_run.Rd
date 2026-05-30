\name{ff_run}
\alias{ff_run}
\title{Run the full FF simulation study}
\usage{
ff_run(config = ff_default_config(), outdir = NULL, make_plots = TRUE,
  verbose = TRUE)
}
\arguments{
\item{config}{A configuration list, see \code{\link{ff_default_config}}.}
\item{outdir}{Optional directory for \code{ff_weight_layer.csv},
\code{ff_clustering.csv} and, if \code{make_plots}, the matching PDFs. Created
if needed. If \code{NULL}, nothing is written.}
\item{make_plots}{Logical; write the figures when \code{outdir} is given.}
\item{verbose}{Logical; print progress messages.}
}
\value{
Invisibly, a list with \code{weight_layer} (Experiment A summary),
\code{clustering} (Experiment B summary), and \code{draws}.
}
\description{
Runs both experiments over the eta grid and the configured replicates and,
optionally, writes the result tables and figures to a directory. The returned
numbers are produced by this code; the function does not reproduce any specific
published table.
}
\examples{
\donttest{
res <- ff_run(ff_default_config(quick = TRUE), verbose = FALSE)
res$weight_layer
}
}
