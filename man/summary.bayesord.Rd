\name{summary.bayesord}
\alias{summary.bayesord}
\alias{print.summary.bayesord}
\title{
Produces posterior summaries resulting from a call to \code{\link{bayesord}}
}
\description{Produces posterior summaries resulting from a call to \code{\link{bayesord}}.}
\usage{
\method{summary}{bayesord}(object, type = c("reg", "rand", "all"), scale = c("logOR", "OR"), threshold = 0, topmodels = 5, digits = 2, ...)
\method{print}{summary.bayesord}(x, ...)
}
\arguments{
  \item{object}{an object of class "bayesord".}
  \item{x}{an object of class "summary.bayesord".}
  \item{type}{a character specifying whether to plot the regression parameters, the clustering parameters, or all of them.}
  \item{scale}{a character specifying whether to plot regression parameters as log(odds ratios) or odds-ratios.}
  \item{threshold}{a numeric between 0 and 1. Corresponds to a threshold applied to the posterior probabilities of association (PPA) such that only those variables with a PPA of >\code{threshold} are summarised. (Only necessary if variable selection is was used.)}
  \item{topmodels}{a positive integer specifying how many models to summarise. Lists them in decreasing PPA order.}
  \item{digits}{a positive integer controlling the number of digits used in rounding the output.}
  \item{\dots}{not used here.}
}
\details{Produces various posterior summaries. Also includes a measure of the effective sample size, as well as the ratio of the Monte Carlo error to the posterior standard deviation. If model selection was used, then this also produces summaries of the posterior probabilities of association (PPAs), in terms of both the posterior weight according to each model, as well as for each variable after averaging over all the models.}
\value{
Produces a list with elements:
\item{output}{a list, where each element of the list is a matrix of summary values for the various parameters of the model.}
\item{models}{a list of length \code{nchains}, where each element of the list is a matrix of summary values for the top models from each chain.}
\item{variables}{a list of length \code{nchains}, where each element of the list is a matrix of summary values for each of the variables from each chain.}
\item{threshold}{a numeric scalar recording the PPA threshold for printing.}
\item{topmodels}{an integer recording the number of top models to print.}
\item{digits}{a positive integer controlling the number of digits used in rounding the output.}
}
\references{McKinley et al., (2014), in review.}
\author{TJ McKinley}
\seealso{\code{\link{bayesord}}, \code{\link{plot.bayesord}}}
\examples{
# load in US data on religious fundamentalism against region of residence
# (data can be found in Agresti, "Analysis of Ordinal Categorical Data", (2010))
data(religion)

# fit non-proportional odds model
religion.NPO <- bayesord(Beliefs ~ Region, data = religion, niter = 2e+05, noutputsum = 2e+04, 
    nchains = 2, model.type = "NPO", start = 50000, thin = 150)
    
## check convergence plots (NOT RUN)
#plot(religion.NPO)

#summarise output
print(religion.NPO)
summary(religion.NPO)
}
