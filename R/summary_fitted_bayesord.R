# summary function for 'fitted.bayesord' object
summary.fitted.bayesord <- function(object, ...) {
    if (class(object) != "fitted.bayesord") 
        stop("'object' is not a 'fitted.bayesord' object")
    # calculate posterior summaries for proportions correctly predicted
    ind.summary <- sapply(as.list(1:length(object$fits)), function(i, x, y) {
        x <- x[[i]]
        y <- y[i]
        x <- x[y + 1, ]
        x
    }, x = object$fits, y = object$model.dat[, (ncol(object$model.dat) - 1)])
    ind.summary <- apply(ind.summary, 1, sum)/sum(object$var.info$orig.data$counts)
    ind.summary <- c(mean(ind.summary), quantile(ind.summary, probs = c(0.025, 0.975)))
    class(ind.summary) <- "summary.fitted.bayesord"
    ind.summary
}

