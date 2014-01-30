# extract posterior 'fitted' distributions for 'bayesord' object
fitted.bayesord <- function(object, multi = F, mc.cores = NA, cluster = TRUE, 
    ...) {
    if (class(object) != "bayesord") 
        stop("'object' is not a 'bayesord' object")
    
    # extract posteriors
    beta <- as.matrix(object$beta)
    theta <- as.matrix(object$theta)
    
    # augment beta posterior if PO model
    if (object$info$model.type == 0) {
        temp.lev <- colnames(beta)
        beta <- matrix(rep(as.vector(beta), ncol(theta)), nrow = nrow(beta))
        colnames(beta) <- as.vector(apply(cbind(matrix(rep(temp.lev, ncol(theta)), 
            nrow = ncol(theta), byrow = T), 1:ncol(theta)), 1, function(x) {
            y <- x[-length(x)]
            y <- paste(y, "_", x[length(x)], sep = "")
            y
        }))
    }
    if (!is.na(object$psi[1]) & cluster == TRUE) {
        psi <- as.matrix(object$psi)
        mcmc <- cbind(beta, theta, psi)
    } else {
        psi <- NA
        mcmc <- cbind(beta, theta)
    }
    nbeta <- ncol(beta)
    ntheta <- ncol(theta)
    
    if (multi == T) {
        if(require(multicore))
        {
            x.list <- as.list(as.data.frame(t(object$model.dat)))
            
            #check number of cores
            if(is.na(mc.cores)) mc.cores <- multicore:::detectCores()
            
            # produce predicted posteriors for individuals in data set
            pred <- mclapply(x.list, function(dat, mcmc, nbeta, ntheta) {
                id <- dat[length(dat)]
                dat <- dat[-length(dat)]
                response <- dat[length(dat)]
                dat <- dat[-length(dat)]
                counts <- dat[length(dat)]
                dat <- dat[-length(dat)]
                if (!is.na(psi[1])) {
                    psi <- mcmc[, (nbeta + ntheta + 1):ncol(mcmc)]
                    psi <- psi[, id]
                } else psi <- numeric(nrow(mcmc))
                
                mcmc <- cbind(mcmc[, 1:(nbeta + ntheta)], psi)
                
                pred <- apply(mcmc, 1, function(mcmc, dat, counts, nbeta, ntheta) {
                    psi <- mcmc[length(mcmc)]
                    mcmc <- mcmc[-length(mcmc)]
                    beta <- mcmc[1:nbeta]
                    theta <- mcmc[-(1:nbeta)]
                    
                    beta <- matrix(beta, nrow = ntheta, byrow = T)
                    beta <- apply(beta, 1, function(x, dat) x * dat, dat = dat)
                    if (!is.null(nrow(beta))) 
                      beta <- apply(beta, 2, sum)
                    beta <- beta + psi
                    beta <- theta - beta
                    gamma <- exp(beta)/(1 + exp(beta))
                    p <- c(gamma[1], gamma[-1] - gamma[-length(gamma)], 1 - gamma[length(gamma)])
                    # now simulate a response score
                    g <- rmultinom(1, size = counts, prob = p)
                    # g<-apply(g,2,function(x) which(x==1)-1)
                    # g<-summary(factor(g,levels=0:(length(p)-1)))
                    g
                }, dat = dat, counts = counts, nbeta = nbeta, ntheta = ntheta)
                list(pred)
            }, mcmc = mcmc, nbeta = nbeta, ntheta = ntheta, mc.cores = mc.cores)
        } else {
            stop("Parallel processing requires multicore library. Please install or set multi = FALSE")
        }
    } else {
        # produce predicted posteriors for individuals in data set
        pred <- apply(object$model.dat, 1, function(dat, mcmc, nbeta, ntheta) {
            id <- dat[length(dat)]
            dat <- dat[-length(dat)]
            response <- dat[length(dat)]
            dat <- dat[-length(dat)]
            counts <- dat[length(dat)]
            dat <- dat[-length(dat)]
            if (!is.na(psi[1])) {
                psi <- mcmc[, (nbeta + ntheta + 1):ncol(mcmc)]
                psi <- psi[, id]
            } else psi <- numeric(nrow(mcmc))
            
            mcmc <- cbind(mcmc[, 1:(nbeta + ntheta)], psi)
            
            pred <- apply(mcmc, 1, function(mcmc, dat, counts, nbeta, ntheta) {
                psi <- mcmc[length(mcmc)]
                mcmc <- mcmc[-length(mcmc)]
                beta <- mcmc[1:nbeta]
                theta <- mcmc[-(1:nbeta)]
                
                beta <- matrix(beta, nrow = ntheta, byrow = T)
                beta <- apply(beta, 1, function(x, dat) x * dat, dat = dat)
                if (!is.null(nrow(beta))) 
                  beta <- apply(beta, 2, sum)
                beta <- beta + psi
                beta <- theta - beta
                gamma <- exp(beta)/(1 + exp(beta))
                p <- c(gamma[1], gamma[-1] - gamma[-length(gamma)], 1 - gamma[length(gamma)])
                # now simulate a response score
                g <- rmultinom(1, size = counts, prob = p)
                # g<-apply(g,2,function(x) which(x==1)-1)
                # g<-summary(factor(g,levels=0:(length(p)-1)))
                g
            }, dat = dat, counts = counts, nbeta = nbeta, ntheta = ntheta)
            list(pred)
        }, mcmc = mcmc, nbeta = nbeta, ntheta = ntheta)
    }
    pred <- lapply(pred, function(x) x[[1]])
    # output predictions
    pred <- list(fits = pred, model.dat = object$model.dat, var.info = object$var.info)
    class(pred) <- "fitted.bayesord"
    pred
}
