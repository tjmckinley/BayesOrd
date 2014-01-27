# print function for 'bayesord' object
print.bayesord <- function(x, ...) {
    if (class(x) != "bayesord") 
        stop("'x' is not a 'bayesord' object")
    
    if (x$info$model.type == 0) 
        cat("\nPROPORTIONAL ODDS MODEL\n") else {
        if (x$info$model.type == 1) 
            cat("\nNON-PROPORTIONAL ODDS MODEL\n") else cat("\nPROPORTIONAL ODDS OR NON-PROPORTIONAL ODDS MODEL\n")
    }
    if (x$info$var.select == T) 
        cat("\nVariable selection specified\n")
    
    cat("\nCall:\n")
    print(x$formula)
    
    cat("\nData:\n")
    cat(paste(sum(x$var.info$orig.data$counts), "observations with", x$info$nbetagroup, 
        ifelse(x$info$nbetagroup > 1, "explanatory variables\n", "explanatory variable\n")))
    cat(paste("Response variable has", x$info$ntheta + 1, "ordered groups:", paste(levels(x$var.info$orig.data[, 
        1]), collapse = ", "), "\n"))
    
    cat("\nMCMC details:\n")
    cat(paste(x$info$nchains, "chains of", x$info$niter, "iterations\n"))
    
    cat("\nPrior information:\n")
    cat(paste("Regression variables: N(", x$info$mnb, ",", x$info$varb, ")\n", sep = ""))
    cat(paste("Cut-points: N(0,", x$info$vart, ")\n", sep = ""))
    if (x$info$RE == 1) {
        cat(paste(x$info$npsi, " random intercepts, each with: N(", x$info$mnpsi, 
            ",varp)\n", sep = ""))
        cat(paste("varp: G(", x$info$shvarp, ",", x$info$rtvarp, ") [i.e. mean=", 
            x$info$shvarp/x$info$rtvarp, ", variance=", x$info$shvarp/(x$info$rtvarp^2), 
            "]\n", sep = ""))
    }
}

