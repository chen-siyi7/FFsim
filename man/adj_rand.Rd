\name{adj_rand}
\alias{adj_rand}
\title{Adjusted Rand index}
\usage{
adj_rand(a, b)
}
\arguments{
\item{a, b}{Integer or factor label vectors of equal length.}
}
\value{The adjusted Rand index (numeric scalar); 1 for identical partitions.}
\description{
Permutation-invariant agreement between two label vectors.
}
\examples{
adj_rand(c(1, 1, 2, 2), c(2, 2, 1, 1))
}
