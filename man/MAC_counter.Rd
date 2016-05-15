\name{MAC_counter}
\alias{MAC_counter}
\alias{max_lag_prop}
\alias{MAC_cutoff}
\alias{file_name}
\title{
Function to perform lag-based correlation anaylsis of single-cell sequencing data, sorted by pseudotime.
}
\description{
See vignette for more details.
}
\usage{
MAC_counter(data, max_lag_prop = 1/3, MAC_cutoff = 0.2, 
file_name = F, lag_matrix = T, symmetric = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A data matrix for which the rows are genes and the columns are experiments, sorted by their pseudotime.
}
  \item{max_lag_prop}{
The largest proportion of your experiments that you want the size of the lag to be. Recommended not to go beyond 1/3. Default value is 1/3.
}
  \item{MAC_cutoff}{
The lowest MAC desired to be shown in results. Default value is 0.2.
}
  \item{file_name}{
The name to be used for resulting .csv files. I.e., using file_name="mine" would create the file "MAC_mine.csv" and "lag_mine.csv"
}
  \item{lag_matrix}{
Logical. TRUE indicates that the resulting matrix of lags should be saved to a csv. Default value is TRUE.
}
  \item{symmetric}{
Logical. TRUE indicates that a symmetric matrix is required for output MAC .csv file. Results in the absolute maximum value for each pair (i,j) and (j,i). Default value is FALSE.
}
}
\details{
See vignette for more details.
}
\value{
Returns a dataset with four columns: Row gene index and Column gene index correspond to the indeces for the gene pair (i,j), Correlation is the maximum absolute correlation (MAC) achieved for the pair, and Lag is the lag at which the MAC occurred.
}
\examples{

x <- matrix(rnorm(6),2,3)
y <- MAC_counter(x)

\dontrun{MAC_results = MAC_counter(data=example_data)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
