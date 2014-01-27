# plot function for 'fitted.bayesord' objects
plot.fitted.bayesord <- function(object, variable = "all", ...) {
    if (class(object) != "fitted.bayesord") 
        stop("'object' is not a 'fitted.bayesord' object")
    plot.marginal.bayesord(object, variable = variable, ...)
}

