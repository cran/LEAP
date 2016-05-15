\name{MAC_perm}
\alias{MAC_perm}
\alias{num_perms}
\alias{perm_file_name}
\alias{FDR_cutoffs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to perform a permutation analysis to determine a cutoff for significant MAC values.
}
\description{
See vignette for more details.
}
\usage{
 MAC_perm(data, MACs_observ, num_perms = 100, max_lag_prop = 1/3, 
 FDR_cutoffs = 101, perm_file_name = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A data matrix for which the rows are genes and the columns are experiments, sorted by their pseudotime.
}
  \item{MACs_observ}{
The resulting matrix of MACs from running MAC_counter on the dataset
}
  \item{num_perms}{
The number of permutations to be performed. Default is 100.  
}
  \item{max_lag_prop}{
The largest proportion of your experiments that you want the size of the lag to be. Recommended not to go beyond 1/3. Default value is 1/3.
}
  \item{FDR_cutoffs}{
The number of cutoffs between 0 and 1 to use for FDR analysis. Default value is 101, resulting in 0,0.01,0.02,...,0.98,0.99,1.
}
  \item{perm_file_name}{
The name to be used for resulting .csv file. I.e., using perm_file_name="mine" would create the file "perm_mine.csv"
}
}
\details{
See vignette for more details.
}
\value{
Returns a dataset with four columns: cors are the correlation cutoffs, MACs_observed are the number of observed correlations at that cutoff, MACs_ave_perm are the average number observed in the permuted datasets at that cutoff, and fdr is the false discovery rate (FDR) observed at that cutoff.
}
\examples{

x <- matrix(rnorm(6),2,3)
cor <- cor(x)
y <- MAC_perm(x, cor)

\dontrun{MAC_perm(data=example_data, MACs_observ=MAC_example)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
