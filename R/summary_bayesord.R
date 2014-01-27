# summary function for 'bayesord' object
summary.bayesord <- function(object, type = c("reg", "rand", "all"), scale = c("logOR", 
    "OR"), threshold = 0, topmodels = 5, digits = 2, ...) {
    require(coda)
    if (class(object) != "bayesord") 
        stop("'object' is not a 'bayesord' object")
    if (type[1] != "reg" & type[1] != "rand" & type[1] != "all") 
        stop("'type' is incorrect")
    if (scale[1] != "logOR" & scale[1] != "OR") 
        stop("'scale' is incorrect")
    if (!is.numeric(topmodels) | length(topmodels) > 1) 
        stop("'topmodels' is incorrect")
    if (!is.numeric(digits) | length(digits) > 1) 
        stop("'digits' is incorrect")
    if (!is.numeric(threshold) | length(threshold) > 1) 
        stop("'threshold' is incorrect")
    
    reg <- ifelse(type[1] == "reg" | type[1] == "all", 1, 0)
    cutp <- ifelse(type[1] == "reg" | type[1] == "all", 1, 0)
    randint <- ifelse(type[1] == "rand" | type[1] == "all", 1, 0)
    varp <- ifelse(type[1] == "rand" | type[1] == "all", 1, 0)
    sdb <- ifelse(type[1] == "reg" | type[1] == "all", 1, 0)
    loglikelihood <- ifelse(type[1] == "reg" | type[1] == "all", 1, 0)
    
    status <- ifelse(type[1] == "reg" | type[1] == "all", 1, 0)
    
    # adjust possibilities by model.type
    if (object$info$model.type != 2 & object$info$var.select != 1) 
        status <- 0
    if (object$info$fixed == 1) 
        sdb <- 0
    
    output <- list(NULL)
    ind <- 1
    if (reg == 1) {
        reg <- as.matrix(object$beta)
        temp <- reg
        if (scale[1] == "OR") 
            temp <- exp(temp)
        temp <- as.mcmc(temp)
        ess <- effectiveSize(temp)
        temp <- summary(temp)
        temp1 <- cbind(temp$statistics[, 1:2], temp$quantiles[, c(3, 1, 5)], ess, 
            temp$statistics[, 4]/temp$statistics[, 2])
        rownames(temp1) <- sapply(as.list(colnames(reg)), function(x, y) paste(y, 
            "(", x, ")", sep = ""), y = scale[1])
        colnames(temp1) <- c("Mean", "SD", "Median", "2.5%", "97.5%", "ESS", "MC error/SD")
        if (status == 1) {
            # #if variable selection has been used then average across non-zero models
            # #produce PPAs for each variable status1<-as.matrix(object$status) #produce
            # model averaged results reg<-as.matrix(object$beta)
            # temp1<-matrix(NA,ncol(reg),7) k<-1 for(j in 1:ncol(reg)) {
            # temp<-reg[status1[,k]<2,j] if(scale[1]=='OR') temp<-exp(temp)
            # if(length(temp)>1) { temp<-as.mcmc(temp) ess<-effectiveSize(temp)
            # temp<-summary(temp)
            # temp1[j,]<-c(temp$statistics[1:2],temp$quantiles[c(3,1,5)],ess,temp$statistics[4]/temp$statistics[2])
            # } if(k%%ncol(status1)==0) k<-1 else k<-k+1 }
            # rownames(temp1)<-sapply(as.list(colnames(reg)),function(x,y)
            # paste(y,'(',x,')',sep=''),y=scale[1])
            # colnames(temp1)<-c('Mean','SD','Median','2.5%','97.5%','ESS','MC error/SD')
            
            status <- list(NULL)
            for (i in 1:length(object$status)) {
                # produce models via PPAs
                status[[i]] <- as.matrix(object$status[[i]])
                lev <- colnames(status[[i]])
                lev <- sapply(as.list(lev), function(x) {
                  x <- strsplit(x, "status_")[[1]]
                  x[x != ""]
                })
                # extract unique models
                status[[i]] <- table(factor(apply(status[[i]], 1, function(x) paste(x, 
                  collapse = ""))))
                status[[i]] <- status[[i]]/sum(status[[i]])
                status1 <- matrix(NA, length(status[[i]]), length(lev) + 1)
                for (j in 1:length(status[[i]])) status1[j, ] <- c(strsplit(names(status[[i]])[j], 
                  "")[[1]], status[[i]][j])
                status1[, -ncol(status1)] <- t(apply(status1[, -ncol(status1), drop = F], 
                  1, function(x) {
                    x[x == "0"] <- "PO"
                    x[x == "1"] <- "NPO"
                    x[x == "2"] <- "0"
                    x
                  }))
                status1 <- as.data.frame(status1)
                colnames(status1) <- c(lev, "PPAs")
                status[[i]] <- status1
                status[[i]] <- status[[i]][sort.list(status[[i]]$PPAs, decreasing = T), 
                  ]
            }
            
            status1 <- list(NULL)
            for (i in 1:length(object$status)) {
                status1[[i]] <- as.matrix(object$status[[i]])
                status1[[i]] <- apply(status1[[i]], 2, function(x) table(factor(x, 
                  levels = c("0", "1", "2"))))
                status1[[i]] <- apply(status1[[i]], 2, function(x) x/sum(x))
                rownames(status1[[i]])[rownames(status1[[i]]) == "0"] <- "PO"
                rownames(status1[[i]])[rownames(status1[[i]]) == "1"] <- "NPO"
                rownames(status1[[i]])[rownames(status1[[i]]) == "2"] <- "Excluded"
                colnames(status1[[i]]) <- lev
            }
        } else {
            # reg<-as.matrix(object$beta) temp<-reg if(scale[1]=='OR') temp<-exp(temp)
            # temp<-as.mcmc(temp) ess<-effectiveSize(temp) temp<-summary(temp)
            # temp1<-cbind(temp$statistics[,1:2],temp$quantiles[,c(3,1,5)],ess,temp$statistics[,4]/temp$statistics[,2])
            # rownames(temp1)<-sapply(as.list(colnames(reg)),function(x,y)
            # paste(y,'(',x,')',sep=''),y=scale[1])
            # colnames(temp1)<-c('Mean','SD','Median','2.5%','97.5%','ESS','MC error/SD')
            
            status <- NA
            status1 <- NA
        }
        if (object$info$model.type != 0) {
            regno <- as.numeric(matrix(1:(object$info$nbetagroup * object$info$ntheta), 
                object$info$ntheta, byrow = T))
            temp1 <- temp1[regno, ]
        }
        output[[ind]] <- temp1
        ind <- ind + 1
    } else {
        output[[ind]] <- NA
        ind <- ind + 1
    }
    if (cutp == 1) {
        cutp <- as.matrix(object$theta)
        temp <- cutp
        temp <- as.mcmc(temp)
        ess <- effectiveSize(temp)
        temp <- summary(temp)
        temp1 <- cbind(temp$statistics[, 1:2], temp$quantiles[, c(3, 1, 5)], ess, 
            temp$statistics[, 4]/temp$statistics[, 2])
        colnames(temp1) <- c("Mean", "SD", "Median", "2.5%", "97.5%", "ESS", "MC error/SD")
        output[[ind]] <- temp1
        ind <- ind + 1
    } else {
        output[[ind]] <- NA
        ind <- ind + 1
    }
    if (sdb == 1) {
        sdb <- as.matrix(object$sdb)
        temp <- sdb
        temp <- as.mcmc(temp)
        ess <- effectiveSize(temp)
        temp <- summary(temp)
        temp1 <- cbind(temp$statistics[, 1:2], temp$quantiles[, c(3, 1, 5)], ess, 
            temp$statistics[, 4]/temp$statistics[, 2])
        colnames(temp1) <- c("Mean", "SD", "Median", "2.5%", "97.5%", "ESS", "MC error/SD")
        if (object$info$model.type != 0) {
            regno <- as.numeric(matrix(1:(object$info$nbetagroup * object$info$ntheta), 
                object$info$ntheta, byrow = T))
            temp1 <- temp1[regno, ]
        }
        output[[ind]] <- temp1
        ind <- ind + 1
    } else {
        output[[ind]] <- NA
        ind <- ind + 1
    }
    if (varp == 1) {
        if (!is.na(object$varp[1])) {
            varp <- as.matrix(object$varp)
            temp <- varp
            temp <- as.mcmc(temp)
            ess <- effectiveSize(temp)
            temp <- summary(temp)
            temp1 <- cbind(temp$statistics[, 1:2], temp$quantiles[, c(3, 1, 5)], 
                ess, temp$statistics[, 4]/temp$statistics[, 2])
            colnames(temp1) <- c("Mean", "SD", "Median", "2.5%", "97.5%", "ESS", 
                "MC error/SD")
            output[[ind]] <- temp1
            ind <- ind + 1
        }
    } else {
        output[[ind]] <- NA
        ind <- ind + 1
    }
    if (randint == 1) {
        if (!is.na(object$psi[1])) {
            randint <- as.matrix(object$psi[[i]])
            randint <- cbind(randint, as.matrix(object$varp))
            temp <- randint
            temp <- as.mcmc(temp)
            ess <- effectiveSize(temp)
            temp <- summary(temp)
            temp1 <- cbind(temp$statistics[, 1:2], temp$quantiles[, c(3, 1, 5)], 
                ess, temp$statistics[, 4]/temp$statistics[, 2])
            colnames(temp1) <- c("Mean", "SD", "Median", "2.5%", "97.5%", "ESS", 
                "MC error/SD")
            output[[ind]] <- temp1
            ind <- ind + 1
        } else {
            if (type != "all") 
                stop("Can't plot as no random intercepts in this model")
        }
    } else {
        output[[ind]] <- NA
        ind <- ind + 1
    }
    if (loglikelihood == 1) {
        loglikelihood <- as.matrix(object$loglikelihood)
        colnames(loglikelihood) <- "loglikelihood"
        temp <- loglikelihood
        temp <- as.mcmc(temp)
        ess <- effectiveSize(temp)
        temp <- summary(temp)
        temp1 <- matrix(c(temp$statistics[1:2], temp$quantiles[c(3, 1, 5)], ess, 
            temp$statistics[4]/temp$statistics[2]), nrow = 1)
        colnames(temp1) <- c("Mean", "SD", "Median", "2.5%", "97.5%", "ESS", "MC error/SD")
        output[[ind]] <- temp1
        ind <- ind + 1
    } else {
        output[[ind]] <- NA
        ind <- ind + 1
    }
    
    # output summary objects for printing
    output <- list(output = output, models = status, variables = status1, threshold = threshold, 
        topmodels = topmodels, digits = digits)
    class(output) <- "summary.bayesord"
    output
}
