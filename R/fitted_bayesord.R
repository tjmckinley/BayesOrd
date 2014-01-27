# extract posterior 'fitted' distributions for 'bayesord' object
fitted.bayesord <- function(x, multi = F, mc.cores = getOption("cores"), RE = TRUE, 
    ...) {
    require(coda)
    if (class(x) != "bayesord") 
        stop("'x' is not a 'bayesord' object")
    
    # extract posteriors
    beta <- as.matrix(x$beta)
    theta <- as.matrix(x$theta)
    
    # augment beta posterior if PO model
    if (x$info$model.type == 0) {
        temp.lev <- colnames(beta)
        beta <- matrix(rep(as.vector(beta), ncol(theta)), nrow = nrow(beta))
        colnames(beta) <- as.vector(apply(cbind(matrix(rep(temp.lev, ncol(theta)), 
            nrow = ncol(theta), byrow = T), 1:ncol(theta)), 1, function(x) {
            y <- x[-length(x)]
            y <- paste(y, "_", x[length(x)], sep = "")
            y
        }))
    }
    if (!is.na(x$psi[1]) & RE == TRUE) {
        psi <- as.matrix(x$psi)
        mcmc <- cbind(beta, theta, psi)
    } else {
        psi <- NA
        mcmc <- cbind(beta, theta)
    }
    nbeta <- ncol(beta)
    ntheta <- ncol(theta)
    
    if (multi == T) {
        require(multicore)
        x.list <- as.list(as.data.frame(t(x$model.dat)))
        
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
                if (length(p[p < 0]) > 0) 
                  browser()
                # now simulate a response score
                g <- rmultinom(1, size = counts, prob = p)
                # g<-apply(g,2,function(x) which(x==1)-1)
                # g<-summary(factor(g,levels=0:(length(p)-1)))
                g
            }, dat = dat, counts = counts, nbeta = nbeta, ntheta = ntheta)
            list(pred)
        }, mcmc = mcmc, nbeta = nbeta, ntheta = ntheta, mc.cores = mc.cores)
    } else {
        # produce predicted posteriors for individuals in data set
        pred <- apply(x$model.dat, 1, function(dat, mcmc, nbeta, ntheta) {
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
                if (length(p[p < 0]) > 0) 
                  browser()
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
    pred <- list(fits = pred, model.dat = x$model.dat, var.info = x$var.info)
    class(pred) <- "fitted.bayesord"
    pred
}
