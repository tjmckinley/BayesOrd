# window function for 'bayesord' objects
window.bayesord <- function(x, start = NA, end = NA, thin = NA, chains = NA, ...) {
    require(coda)
    if (class(x) != "bayesord") 
        stop("'x' is not a 'bayesord' object")
    if (!is.na(chains)) {
        x$beta <- x$beta[chains]
        x$theta <- x$theta[chains]
        x$status <- x$status[chains]
        x$status <- x$sdb[chains]
        if (!is.na(x$psi[1])) {
            x$psi <- x$psi[chains]
            x$varp <- x$varp[chains]
        }
        x$loglikelihood <- x$loglikelihood[chains]
        x$info$nchains <- length(chains)
    }
    if (is.na(start)) {
        if (is.na(end)) {
            if (is.na(thin)) {
                x <- x
            } else {
                for (i in 1:length(x)) {
                  if (is.mcmc.list(x[[i]])) 
                    x[[i]] <- window(x[[i]], thin = thin)
                }
            }
        } else {
            if (is.na(thin)) {
                for (i in 1:length(x)) {
                  if (is.mcmc.list(x[[i]])) 
                    x[[i]] <- window(x[[i]], end = end)
                }
            } else {
                for (i in 1:length(x)) {
                  if (is.mcmc.list(x[[i]])) 
                    x[[i]] <- window(x[[i]], end = end, thin = thin)
                }
            }
        }
    } else {
        if (is.na(end)) {
            if (is.na(thin)) {
                for (i in 1:length(x)) {
                  if (is.mcmc.list(x[[i]])) 
                    x[[i]] <- window(x[[i]], start = start)
                }
            } else {
                for (i in 1:length(x)) {
                  if (is.mcmc.list(x[[i]])) 
                    x[[i]] <- window(x[[i]], start = start, thin = thin)
                }
            }
        } else {
            if (is.na(thin)) {
                for (i in 1:length(x)) {
                  if (is.mcmc.list(x[[i]])) 
                    x[[i]] <- window(x[[i]], start = start, end = end)
                }
            } else {
                for (i in 1:length(x)) {
                  if (is.mcmc.list(x[[i]])) 
                    x[[i]] <- window(x[[i]], start = start, end = end, thin = thin)
                }
            }
        }
    }
    x
}

