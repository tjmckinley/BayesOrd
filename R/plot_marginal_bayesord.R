# function to plot fitted distributions
plot.marginal.bayesord <- function(object, variable = "all", ...) {
    if (class(object) != "fitted.bayesord" & class(object) != "predict.bayesord") 
        stop("'object' is not a 'fitted.bayesord' or a 'predict.bayesord' object")
    if (!is.character(variable)) 
        stop("'variable' not a character")
    
    if (variable == "all") {
        # extract variables to plot from object
        variable <- object$var.info$variables[-1]
    } else {
        x <- match(variable, object$var.info$variables[-1])
        if (length(x[is.na(x)]) > 1) 
            stop(paste("One or more of 'variable' not found in data set", sep = ""))
    }
    x <- numeric(length(variable))
    for (i in 1:length(variable)) x[i] <- ifelse(!is.factor(object$var.info$orig.data[, 
        match(variable[i], colnames(object$var.info$orig.data))]), 0, 1)
    nonfactor <- variable[x == 0]
    if (length(nonfactor) > 0) 
        cat(paste("The following variables are non-categorical and have been removed:\n", 
            paste(nonfactor, collapse = "\n"), "\n"))
    variable <- variable[x == 1]
    if (length(variable) == 0) 
        stop("No non-categorical variables in list")
    
    # set up plotting parameters
    if (length(variable) > 1) {
        if (length(variable) == 2) 
            mfrow1 <- c(2, 1) else mfrow1 <- c(ceiling(length(variable)/2), 2)
    } else mfrow1 <- c(1, 1)
    par(mfrow = mfrow1)
    
    # plot each variable in turn
    data1 <- object$var.info$orig.data
    response <- data1[, 1]
    data1 <- data1[, -1]
    for (i in 1:length(variable)) {
        var <- data1[, match(c(variable[i], "counts"), colnames(data1))]
        var$response <- response
        
        #'aggregate' removes rows with zeros so as a fudge bind some dummy values to the
        # bottom before aggregating
        tempvar <- expand.grid(levels(var[, 1]), levels(response))
        tempvar <- data.frame(tempvar[, 1], rep(1, nrow(tempvar)), tempvar[, 2])
        colnames(tempvar) <- colnames(var)
        tempvar <- rbind(var, tempvar)
        temp1 <- aggregate(tempvar$counts, list(tempvar$response, tempvar[, 1]), 
            sum)
        temp1[, 3] <- temp1[, 3] - 1
        temp1 <- matrix(as.numeric(temp1$x), length(levels(var$response)))
        rownames(temp1) <- levels(var$response)
        colnames(temp1) <- levels(var[, 1])
        temp1 <- t(temp1)
        temp1 <- t(apply(temp1, 1, function(x) x/sum(x)))
        # plot output
        temp1 <- barplot(temp1, beside = T, legend = T, ylim = c(0, 1), main = variable[i])
        # add posterior model fits
        for (j in 1:length(levels(var[, 1]))) {
            temp <- object$fits[var[, 1] == levels(var[, 1])[j]]
            temp <- lapply(as.list(1:nrow(temp[[1]])), function(i, x) sapply(x, function(x, 
                i) x[i, ], i = i), x = temp)
            temp <- sapply(temp, function(x) apply(x, 1, sum))
            temp <- apply(temp, 1, function(x) x/sum(x))
            temp <- apply(temp, 1, function(x) c(mean(x), quantile(x, probs = c(0.025, 
                0.975))))
            
            # now add predict posterior summaries to plot
            for (k in 1:ncol(temp1)) {
                points(temp1[j, k], temp[1, k], pch = 20)
                lines(rep(temp1[j, k], 2), temp[2:3, k])
                lines(c(temp1[j, k] - 0.05, temp1[j, k] + 0.05), rep(temp[2, k], 
                  2))
                lines(c(temp1[j, k] - 0.05, temp1[j, k] + 0.05), rep(temp[3, k], 
                  2))
            }
        }
    }
    # reset graphical parameters
    par(mfrow = c(1, 1))
}

