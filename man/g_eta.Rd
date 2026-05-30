\name{g_eta}
\alias{g_eta}
\title{Sublinear occupancy transform g_eta}
\usage{
g_eta(n, eta)
}
\arguments{
\item{n}{Non-negative integer occupancy (scalar).}
\item{eta}{Flattening exponent in (0, 1].}
}
\value{An integer-valued occupancy to pass to \code{\link{rcrt}}.}
\description{
The transform that defines the Flattened Franchise. The occupancy \code{n} is
mapped to \code{floor(max(1, n^eta))} before the Chinese Restaurant Table draw.
\code{eta = 1} is the identity (standard franchise); \code{eta < 1} compresses
large occupancies.
}
\examples{
g_eta(1000, 1)
g_eta(1000, 0.5)
g_eta(1, 0.5)
}
