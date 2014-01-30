\name{plot.bayesord}
\alias{plot.bayesord}
\title{Plots summaries of MCMC runs resulting from a call to \code{\link{bayesord}}}
\description{
This function plots summaries of the MCMC runs resulting from a call to \code{\link{bayesord}}. It uses the \code{coda} package to produce trace and autocorrelation plots to help visually assess convergence and mixing of the chain(s).
}
\usage{
\method{plot}{bayesord}(x, which = c("trace", "AC"), type = c("reg", "cluster", "all"), scale = c("logOR", "OR"), ask = F, bystatus = F, ...)
}
\arguments{
  \item{x}{an object of class "bayesord"}
  \item{which}{a character specifying whether to produce trace plots or autocorrelation plots.}
  \item{type}{a character specifying whether to plot the regression parameters, the clustering parameters, or all of them.}
  \item{scale}{a character specifying whether to plot regression parameters as log(odds ratios) or odds-ratios.}
  \item{ask}{a logical specifying whether to ask before producing next set of plots.}
  \item{bystatus}{a logical specifying whether we want to produce trace plots stratified by status. See 'Details' section for more information.}
  \item{\dots}{not used here.}
}
\details{
With RJ-MCMC it is sometimes difficult to assess convergence due to dimension jumping. The \code{bystatus} option allows the trace plots to be stratified by status (i.e. PO or NPO). If variable selection is being used, then this will also extract the subset of samples that are 'included' in the model.
}
\author{TJ McKinley}
\seealso{\code{\link{bayesord}}}
\examples{
# load in US data on religious fundamentalism against region of residence
# (data can be found in Agresti, "Analysis of Ordinal Categorical Data", (2010))
data(religion)

# fit non-proportional odds model
religion.NPO <- bayesord(Beliefs ~ Region, data = religion, niter = 2e+05, noutputsum = 2e+04, 
    nchains = 2, model.type = "NPO", start = 50000, thin = 150)
    
# check convergence plots (set 'ask = TRUE' to pause plots)
plot(religion.NPO)
}
