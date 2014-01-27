# print function for 'summary.fitted.bayesord' object
print.summary.fitted.bayesord <- function(x, digits = 2, ...) {
    if (class(x) != "summary.fitted.bayesord") 
        stop("'x' is not a 'summary.fitted.bayesord' object")
    
    x <- round(x, digits = digits)
    
    cat(paste("Posterior mean probability of correct classification (95% CI): ", 
        x[1], " (", paste(x[2:3], collapse = ","), ")\n", sep = ""))
}
