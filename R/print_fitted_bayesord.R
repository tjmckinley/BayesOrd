# print function for 'fitted.bayesord' object
print.fitted.bayesord <- function(x, ...) {
    if (class(x) != "fitted.bayesord") 
        stop("'x' is not a 'fitted.bayesord' object")
    
    x <- summary(x)
    print(x)
}

