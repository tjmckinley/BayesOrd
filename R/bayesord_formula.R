# set formula object for 'bayesord' method than incorporates a random intercept
bayesord.formula <- function(formula, data = list(), nchains = 1, multi = F, model.type = c("PO", 
    "NPO", "both"), var.select = TRUE, niter = 50000, nsavecoda = 1000, mnb = 0, 
    varb = 1000, maxsdb = 20, fixed = FALSE, vart = 1, mnpsi = 0, shvarp = 0.01, 
    rtvarp = 0.01, propsdb = 1, propsdt = 1, propsdp = 1, propsdvarp = 1, runtraining = FALSE, 
    nitertrain = 1000, start = NA, end = NA, thin = NA, ...) {
    if (!is.call(formula)) 
        stop("Must have a valid formula")
    if (length(nchains) > 1 | !is.numeric(nchains)) 
        stop("'nchains' is not a numeric scalar")
    if (length(multi) > 1 | !is.logical(multi)) 
        stop("'multi' is not a logical value")
    if (model.type[1] != "PO" & model.type[1] != "NPO" & model.type[1] != "both") 
        stop("'model.type' is incorrectly specified")
    if (length(niter) > 1 | !is.numeric(niter)) 
        stop("'niter' is not a numeric scalar")
    if (length(nsavecoda) > 1 | !is.numeric(nsavecoda)) 
        stop("'nsavecoda' is not a numeric scalar")
    if (length(mnb) > 1 | !is.numeric(mnb)) 
        stop("'mnb' is not a numeric scalar")
    if (length(varb) > 1 | !is.numeric(varb)) 
        stop("'varb' is not a numeric scalar")
    if (length(maxsdb) > 1 | !is.numeric(maxsdb)) 
        stop("'maxsdb' is not a numeric scalar")
    if (length(fixed) > 1 | !is.logical(fixed)) 
        stop("'fixed' is not a logical value")
    if (length(vart) > 1 | !is.numeric(vart)) 
        stop("'vart' is not a numeric scalar")
    if (length(mnpsi) > 1 | !is.numeric(mnpsi)) 
        stop("'mnpsi' is not a numeric scalar")
    if (length(shvarp) > 1 | !is.numeric(shvarp)) 
        stop("'shvarp' is not a numeric scalar")
    if (length(rtvarp) > 1 | !is.numeric(rtvarp)) 
        stop("'rtvarp' is not a numeric scalar")
    if (length(propsdb) > 1 | !is.numeric(propsdb)) 
        stop("'propsdb' is not a numeric scalar")
    if (length(propsdt) > 1 | !is.numeric(propsdt)) 
        stop("'propsdt' is not a numeric scalar")
    if (length(propsdp) > 1 | !is.numeric(propsdp)) 
        stop("'propsdp' is not a numeric scalar")
    if (length(propsdvarp) > 1 | !is.numeric(propsdvarp)) 
        stop("'propsdvarp' is not a numeric scalar")
    if (length(runtraining) > 1 | !is.logical(runtraining)) 
        stop("'runtraining' is not a logical value")
    if ((length(start) > 1 | !is.numeric(start)) & !is.na(start)) 
        stop("'start' is not a numeric scalar or NA")
    if (length(nitertrain) > 1 | !is.numeric(nitertrain)) 
        stop("'nitertrain' is not a numeric scalar")
    if ((length(end) > 1 | !is.numeric(end)) & !is.na(end)) 
        stop("'end' is not a numeric scalar or NA")
    if ((length(thin) > 1 | !is.numeric(thin)) & !is.na(thin)) 
        stop("'thin' is not a numeric scalar or NA")
    
    # extract formula
    form <- extractTerms(formula)
    # extract random intercepts term if required
    RE <- form[[2]]
    form <- form[[1]]
    mf <- model.frame(formula = form, data = data, na.action = na.fail)
    # check there are no columns called 'RE' or 'counts'
    temp <- match(c("RE", "counts"), attr(mf, "names"))
    if (length(temp[!is.na(temp)]) > 0) 
        stop("Can't name variables 'RE' or 'counts'")
    
    # create vector for random intercepts if required
    if (is.null(RE)) 
        rand.int.vec <- NA else {
        rand.int.vec <- data[, match(RE, colnames(data)), drop = F]
        for (j in 1:ncol(rand.int.vec)) if (!is.factor(rand.int.vec[, j])) 
            stop("Random intercepts term is not a factor")
        for (j in 1:ncol(rand.int.vec)) rand.int.vec[, j] <- as.numeric(rand.int.vec[, 
            j])
    }
    
    # condense data into succinct form for model
    data <- data[, match(attr(mf, "names"), colnames(data))]
    if (!is.null(RE)) 
        data$RE <- factor(rand.int.vec)
    # now aggregate data
    data <- aggregate(rep(1, nrow(data)), data, table)
    data[, ncol(data)] <- as.numeric(data[, ncol(data)])
    colnames(data)[ncol(data)] <- "counts"
    # create design matrix from standard formula
    mf <- model.frame(formula = form, data = data, na.action = na.fail)
    # remove intercept if not already removed
    if (attr(attr(mf, "terms"), "intercept") == 0) 
        stop("Must have intercept in formula")
    x <- model.matrix(form, data = mf)
    
    # set up response variable
    y <- model.response(mf)
    if (!is.factor(y)) 
        stop("Response variable is not a factor")
    
    # now create matrices for variable selection
    variables <- attr(attr(mf, "terms"), "term.labels")
    variables <- variables[variables != "RE"]
    variables <- variables[variables != "counts"]
    nvariables <- length(variables)
    
    xassign <- attr(x, "assign")[-1]
    xassign <- rev((length(xassign):1)[!duplicated(rev(xassign))])
    xassign <- c(0, xassign)  #add nbetagroup to end of this later
    
    # remove intercept from design matrix
    x <- x[, -1]
    
    # extract variable information for use in marginal posterior plotting
    var.info <- as.character(attr(attr(mf, "terms"), "variables"))[-1]
    interactions <- attr(attr(mf, "terms"), "order")
    intpresent <- ifelse(length(interactions[interactions > 1]) > 0, 1, 0)
    var.info <- list(variables = var.info, orig.data = data, interactions = intpresent, 
        assign = attributes(x)$assign[-1])
    
    # dealing with interaction effects
    if (intpresent == 1) {
        # create intfactor matrix
        temp <- variables
        temp <- lapply(as.list(temp), function(x) strsplit(x, ":")[[1]])
        temp1 <- matrix(0, length(temp), length(temp))
        for (j in 1:length(temp)) {
            if (length(temp[[j]]) > 1) {
                for (i in 1:(j - 1)) {
                  if (length(temp[[i]]) < length(temp[[j]])) {
                    temp2 <- match(temp[[j]], temp[[i]])
                    if (length(temp2[!is.na(temp2)]) > 0) 
                      temp1[i, j] <- 1
                  }
                }
            }
        }
        write.table(temp1, "intfactor.txt", row.names = F, col.names = F, sep = "\t", 
            quote = F)
        temp <- sapply(temp, length) - 1
        write.table(temp, "interaction.txt", row.names = F, col.names = F, sep = "\t", 
            quote = F)
        temp1 <- c((1:length(temp))[!duplicated(temp)] - 1, length(temp))
        write.table(temp1, "intstart.txt", row.names = F, col.names = F, sep = "\t", 
            quote = F)
        maxinteraction <- max(temp)
        rm(temp, temp1, temp2)
    } else maxinteraction <- 0
    
    # now run ordinal regression model
    
    # load 'coda' library
    require(coda)
    
    # set options for writing out data
    options(scipen = 999)
    
    # extract data in correct format
    N <- length(y)
    ntheta <- length(levels(y)) - 1
    
    # extract response variable
    response <- as.numeric(y) - 1
    # set random intercepts terms
    if (!is.null(RE)) {
        id <- as.numeric(data$RE)
        if (length(id) == 0) 
            stop("ID variable is not present")
        RE <- 1
    } else {
        id <- rep(1, nrow(x))
        RE <- 0
    }
    id <- as.numeric(id)
    npsi <- length(unique(id))
    
    # extract explanatory variables and convert to correct format
    dat <- x
    temp.lev <- colnames(dat)
    
    # produce string of linear terms
    nbeta <- length(temp.lev) * ntheta
    nbetagroup <- length(temp.lev)
    
    # produce argument list for likelihood function
    dat1 <- data.frame(dat, counts = data$counts, response, id)
    
    # sort dat1 according to RE terms
    if (RE == 1) {
        var.info$orig.data <- var.info$orig.data[sort.list(dat1[, ncol(dat1)]), ]
        dat1 <- dat1[sort.list(dat1[, ncol(dat1)]), ]
        psicount <- cumsum(tapply(dat1[, ncol(dat1)], factor(dat1[, ncol(dat1)], 
            levels = unique(dat1[, ncol(dat1)])), length))
        write.table(psicount, "psicount.txt", row.names = F, col.names = F, sep = "\t", 
            quote = F)
    }
    # browser() compile C code
    system("cp ../ccode/* .")
    system("gcc -Wall main.c index.c loglikelihood.c loglikelihoodpsi.c moveNPO.c moveNPOtoPO.c movePO.c movePOtoNPO.c moveinctoexc.c moveexctoinc.c validitycheck.c movetheta.c movesdb.c -lgsl -lgslcblas -lm -O3 -o MCMCmulti")
    system("rm *.c functions.h")
    
    # now output data file
    write.table(dat1, "data.dat", row.names = F, col.names = F, sep = "\t", quote = F)
    
    # assign data file for use in model checking
    if (RE == 0) 
        id <- NA
    
    # assign model type
    if (model.type[1] == "PO") 
        model.type <- 0 else {
        if (model.type[1] == "NPO") 
            model.type <- 1 else model.type <- 2
    }
    
    # object for output
    info <- as.data.frame(t(as.matrix(c(n = N, niter = niter, nchains = nchains, 
        nbetagroup = nbetagroup, ntheta = ntheta, mnb = mnb, varb = varb, maxsdb = maxsdb, 
        fixed = ifelse(fixed == TRUE, 1, 0), vart = vart, mnpsi = mnpsi, shvarp = shvarp, 
        rtvarp = rtvarp, npsi = npsi, RE = RE, model.type = model.type, var.select = var.select, 
        runtraining = runtraining, nitertrain = ifelse(runtraining == T, nitertrain, 
            NA)))))
    
    # now write priors file
    priors <- data.frame(n = N, niter = niter, nsavecoda = nsavecoda, nvariables = nvariables, 
        nbeta = nbeta, nbetagroup = nbetagroup, ntheta = ntheta, mnb = mnb, varb = varb, 
        maxsdb = maxsdb, fixed = ifelse(fixed == TRUE, 1, 0), vart = vart, propsdb = propsdb, 
        propsdt = propsdt, mnpsi = mnpsi, shvarp = shvarp, rtvarp = rtvarp, propsdp = propsdp, 
        propsdvarp = propsdvarp, npsi = npsi, RE = RE, model.type = model.type, varselect = ifelse(var.select == 
            TRUE, 1, 0), maxinteraction = maxinteraction, runtraining = ifelse(runtraining == 
            TRUE, 1, 0), nitertrain = nitertrain)
    priors <- data.frame(t(as.matrix(priors)), colnames(priors))
    write.table(priors, "priors.txt", row.names = F, col.names = F, quote = F, sep = "\t")
    
    # write assignment file
    write.table(xassign, "xassign.txt", row.names = F, col.names = F, quote = F, 
        sep = "\t")
    
    if (multi == T) {
        library(multicore)
        # run model
        mclapply(as.list(1:nchains), function(x) {
            seed <- floor(2147483647 * runif(1, 0, 1))
            temp <- paste("./MCMCmulti ", x, " ", seed, sep = "")  #,' > output',x,'.txt',sep='')
            system(temp)
        }, mc.cores = nchains, mc.set.seed = T)
    } else {
        for (i in 1:nchains) {
            seed <- floor(2147483647 * runif(1, 0, 1))
            temp <- paste("./MCMCmulti ", i, " ", seed, sep = "")  #' > output',i,'.txt',sep='')
            system(temp)
        }
    }
    print("Finished MCMC run")
    
    options(scipen = 5)
    
    # read in and process data
    beta <- list(NULL)
    theta <- list(NULL)
    psi <- list(NULL)
    status <- list(NULL)
    sdb <- list(NULL)
    varp <- list(NULL)
    loglikelihood <- list(NULL)
    for (i in 1:nchains) {
        # read in table of posterior samples
        beta[[i]] <- read.table(paste("codaMCMC", i, ".txt", sep = ""), header = T)
        theta[[i]] <- as.mcmc(beta[[i]][, (nbeta + 1):(nbeta + ntheta)])
        if (model.type == 2 | var.select == TRUE) {
            status[[i]] <- beta[[i]][, (nbeta + ntheta + 1):(nbeta + ntheta + nbetagroup)]
            if (is.null(ncol(status[[i]]))) 
                status[[i]] <- matrix(status[[i]], ncol = 1)
            colnames(status[[i]]) <- paste("status_", temp.lev, sep = "")
            status[[i]] <- as.mcmc(status[[i]])
        } else status[[i]] <- NA
        psi[[i]] <- as.mcmc(beta[[i]][, (nbeta + ntheta + nbetagroup + 1):(nbeta + 
            ntheta + nbetagroup + npsi)])
        if (fixed == FALSE) {
            sdb[[i]] <- as.mcmc(beta[[i]][, (nbeta + ntheta + nbetagroup + npsi + 
                1):(2 * nbeta + ntheta + nbetagroup + npsi)])
            varp[[i]] <- as.mcmc(beta[[i]][, (2 * nbeta + ntheta + nbetagroup + npsi + 
                1)])
            loglikelihood[[i]] <- as.mcmc(beta[[i]][, (2 * nbeta + ntheta + nbetagroup + 
                npsi + 2)])
        } else {
            sdb[[i]] <- NA
            varp[[i]] <- as.mcmc(beta[[i]][, (nbeta + ntheta + nbetagroup + npsi + 
                1)])
            loglikelihood[[i]] <- as.mcmc(beta[[i]][, (nbeta + ntheta + nbetagroup + 
                npsi + 2)])
        }
        
        if (model.type == 0) {
            # remove extraneous columns if proportional odds model
            beta[[i]] <- beta[[i]][, 1:nbetagroup]
            if (is.null(ncol(beta[[i]]))) 
                beta[[i]] <- matrix(beta[[i]], ncol = 1)
            colnames(beta[[i]]) <- temp.lev
            
            if (fixed == FALSE) {
                sdb[[i]] <- as.mcmc(as.matrix(sdb[[i]])[, 1:nbetagroup])
            }
        } else {
            beta[[i]] <- beta[[i]][, 1:nbeta]
            if (is.null(ncol(beta[[i]]))) 
                beta[[i]] <- matrix(beta[[i]], ncol = 1)
            colnames(beta[[i]]) <- as.vector(apply(cbind(matrix(rep(temp.lev, ntheta), 
                nrow = ntheta, byrow = T), 1:ntheta), 1, function(x) {
                y <- x[-length(x)]
                y <- paste(y, "_", x[length(x)], sep = "")
                y
            }))
        }
        beta[[i]] <- as.mcmc(beta[[i]])
    }
    if (model.type == 2 | var.select == TRUE) 
        status <- as.mcmc.list(status) else status <- NA
    if (fixed == FALSE) 
        sdb <- as.mcmc.list(sdb) else sdb <- NA
    if (RE == 0) {
        varp <- NA
        psi <- NA
    } else {
        varp <- mcmc.list(varp)
        psi <- mcmc.list(psi)
    }
    ans <- list(beta = mcmc.list(beta), theta = mcmc.list(theta), status = status, 
        psi = psi, sdb = sdb, varp = varp, loglikelihood = mcmc.list(loglikelihood), 
        info = info, model.dat = dat1, var.info = var.info)
    # remove text files
    system("rm data.dat priors.txt xassign.txt codaMCMC* MCMCmulti")
    if (RE == 1) 
        system("rm psicount.txt")
    if (intpresent == 1) 
        system("rm interaction.txt intfactor.txt intstart.txt")
    # set class
    class(ans) <- "bayesord"
    # thin if necessary
    if (missing(start)) 
        start <- NA
    if (missing(end)) 
        end <- NA
    if (missing(thin)) 
        thin <- NA
    ans <- window.bayesord(ans, start = start, end = end, thin = thin)
    # output call and formula for printing
    ans$call <- match.call()
    ans$formula <- formula
    ans
}

