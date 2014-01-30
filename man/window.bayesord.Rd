\name{window.bayesord}
\alias{window.bayesord}
\title{Subsetting function for objects of class "bayesord"}
\description{Can be used to extract and thin the MCMC outputs for objects of class "bayesord".}
\usage{
\method{window}{bayesord}(x, start = NA, end = NA, thin = NA, chains = NA, ...)
}
\arguments{
  \item{x}{an object of class "bayesord".}
  \item{start}{an integer specifying the length of the burn-in period if only a subset of the samples are to be returned.}
  \item{end}{an integer specifying the final iteration to be returned.}
  \item{thin}{an integer specifying the thinning rate of the chain.}
  \item{chains}{an integer specifying which chain to extract.}
  \item{\dots}{not used here.}
}
\details{Analogous to the \code{window} function in the \code{coda} package, but for "bayesord" objects.}
\value{
Returns a "bayesord" object according to subsetting criteria.
}
\references{McKinley et al., (2014), in review.}
\author{TJ McKinley}
\seealso{\code{\link{bayesord}}}
\examples{
# load in US data on religious fundamentalism against region of residence
# (data can be found in Agresti, "Analysis of Ordinal Categorical Data", (2010))
data(religion)

# fit non-proportional odds model
religion.NPO <- bayesord(Beliefs ~ Region, data = religion, niter = 2e+05, noutputsum = 2e+04, 
    nchains = 2, model.type = "NPO")
    
## check convergence plots (NOT RUN)
#plot(religion.NPO)
#summarise chains
print(religion.NPO)

# remove burnin and thin chain
religion.NPO <- window(religion.NPO, start = 50000, thin = 150)
#summarise thinned chains
print(religion.NPO)
summary(religion.NPO)
}
