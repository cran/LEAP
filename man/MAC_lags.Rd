\name{MAC_lags}
\alias{MAC_lags}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Internal function used by MAC_counter and MAC_perm
}
\description{
Performs the lag-based correlation analysis
}
\usage{
MAC_lags(data, max_lag_prop = 1/3, symmetric = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
Data being passed to the function
}
  \item{max_lag_prop}{
The largest proportion of your experiments that you want the size of the lag to be. Recommended not to go beyond 1/3. Default value is 1/3.
}
  \item{symmetric}{
Logical. TRUE indicates that a symmetric matrix is required for output MAC .csv file. Results in the absolute maximum value for each pair (i,j) and (j,i). Default value is FALSE.
}
}
\details{
See vignette for more details.
}
\value{
Returns MAC and associated lag matrices.
}
\examples{

x <- matrix(rnorm(6),2,3)
y <- MAC_lags(x)

\dontrun{MAC_results = MAC_lags(data=example_data, max_lag_prop=1/3, symmetric=F)}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
