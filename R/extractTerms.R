# function for extracting terms on either side of vertical bars in formula
extractTerms <- function(term) {
    if (!is.call(term)) 
        stop("Must have a valid formula")
    t <- as.character(term)
    exp <- t[[3]]
    exp <- strsplit(exp, " \\| ")[[1]]
    t[3] <- exp[1]
    t <- paste(t[2], t[1], t[3])
    RE <- NULL
    if (length(exp) > 1) {
        RE <- exp[2]
        RE <- gsub(" ", "", RE)
        RE <- strsplit(RE, "\\+")[[1]]
    }
    t <- list(form = as.formula(t), RE = RE)
    t
}

