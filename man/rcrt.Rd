\name{rcrt}
\alias{rcrt}
\title{Draw a Chinese Restaurant Table (CRT) variate}
\usage{
rcrt(n, theta)
}
\arguments{
\item{n}{Non-negative integer number of customers (scalar).}
\item{theta}{Positive concentration parameter (scalar).}
}
\value{A non-negative integer table count.}
\description{
Samples the number of occupied tables for \code{n} customers at concentration
\code{theta}, the sum over i of Bernoulli(theta / (theta + i - 1)). The
Bernoulli draws use a uniform comparison rather than vectorized
\code{\link[stats]{rbinom}}, which can intermittently return \code{NA} on long
probability vectors.
}
\examples{
set.seed(1)
mean(replicate(2000, rcrt(500, 2)))
}
