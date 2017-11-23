# plot function for 'bayesord' x
plot.bayesord <- function(x, which = c("conv", "AC"), type = c("reg", "cluster", 
    "all"), scale = c("logOR", "OR"), trace = T, density = T, ask = T, bystatus = F, ...) {
    if (class(x) != "bayesord") 
        stop("'x' is not a 'bayesord' x")
    if (which[1] != "conv" & which[1] != "AC") 
        stop("'which' is incorrect")
    if (type[1] != "reg" & type[1] != "cluster" & type[1] != "all") 
        stop("'type' is incorrect")
    if (scale[1] != "logOR" & scale[1] != "OR") 
        stop("'scale' is incorrect")
    if (!is.logical(trace[1]) | length(trace) > 1) 
        stop("'trace' is not a logical value")
    if (!is.logical(density[1]) | length(density) > 1) 
        stop("'density' is not a logical value")
    if (!is.logical(bystatus[1]) | length(bystatus) > 1) 
        stop("'bystatus' is not a logical value")
        
    if (which[1] == "conv" & (trace[1] == F & density[1] == F))
        stop("Either 'trace' or 'density' must be TRUE if 'which = conv'")
    
    reg <- ifelse(type[1] == "reg" | type[1] == "all", 1, 0)
    cutp <- ifelse(type[1] == "reg" | type[1] == "all", 1, 0)
    status <- ifelse(type[1] == "reg" | type[1] == "all", 1, 0)
    clusterint <- ifelse(type[1] == "cluster" | type[1] == "all", 1, 0)
    varp <- ifelse(type[1] == "cluster" | type[1] == "all", 1, 0)
    sdb <- ifelse(type[1] == "reg" | type[1] == "all", 1, 0)
    loglikelihood <- ifelse(type[1] == "reg" | type[1] == "all", 1, 0)
    
    conv <- ifelse(which[1] == "conv", 1, 0)
    CI <- ifelse(which[1] == "CI", 1, 0)
    AC <- ifelse(which[1] == "AC", 1, 0)
    
    if (conv == 1 | AC == 1) {
        output <- list(NULL)
        ind <- 1
        if (reg == 1) {
            reg <- list(NULL)
            for (i in 1:length(x$beta)) {
                reg[[i]] <- as.matrix(x$beta[[i]])
                if (scale[1] == "OR") 
                  reg[[i]] <- exp(reg[[i]])
            }
            output[[ind]] <- reg
            ind <- ind + 1
            
            # if we wish to plot traces by status
            if (bystatus == TRUE) {
                if (is.na(x$status[1])) 
                  stop("Can't plot by status")
                # extract names
                variables <- colnames(as.matrix(x$beta[[1]]))
                # extract statuses
                status <- list(NULL)
                for (i in 1:length(x$beta)) status[[i]] <- as.matrix(x$status[[i]])
                # bind multiple chains together for each variable
                reg <- output[[ind - 1]]
                n <- nrow(reg[[1]])
                reg <- do.call("rbind", reg)
                # bind statuses together for each variable
                status <- do.call("rbind", status)
                # now convert to correct format for plotting
                output <- list(NULL)
                l <- 1
                for (i in 1:ncol(reg)) {
                  # print(paste('i=',i,'l=',l))
                  temp <- matrix(reg[, i], nrow = n)
                  tempstatus <- matrix(status[, l], nrow = n)
                  l <- l + 1
                  if (l%%(ncol(status) + 1) == 0) 
                    l <- 1
                  temp3 <- list(NULL)
                  for (j in 1:ncol(temp)) {
                    temp2 <- list(NULL)
                    for (k in 0:2) temp2[[k + 1]] <- temp[tempstatus[, j] == k, j]
                    temp3[[j]] <- temp2
                  }
                  # reorder
                  temp2 <- list(NULL)
                  for (k in 1:3) temp2[[k]] <- lapply(temp3, function(y, k) y[[k]], 
                    k = k)
                  for (k in 1:3) {
                    temp2[[k]] <- lapply(as.list(1:length(temp2[[k]])), function(i, 
                      y) {
                      y <- y[[i]]
                      y <- cbind(y, rep(i, length(y)))
                      y
                    }, y = temp2[[k]])
                    temp2[[k]] <- do.call("rbind", temp2[[k]])
                  }
                  output[[i]] <- temp2
                }
                # set up plotting parameters
                if(trace == T & density == T)
                {
                    maxp <- length(output) * 2
                    if (maxp >= 4) 
                      mfrow1 <- c(4, 2) else mfrow1 <- c(maxp, 2)
                }
                else
                {
                    maxp <- length(output)
                    if(maxp >= 4)
                        mfrow1 <- c(4, 1) else mfrow1 <- c(maxp, 1)
                }
                cols1 <- c("black", "red", "blue", "green", "yellow", "purple")
                par(mfrow = mfrow1)
                # produce plots
                l <- 1
                for (i in 1:length(output)) {
                  for (k in 1:2) {
                    # plot trace
                    temp <- output[[i]][[k]]
                    if (nrow(temp) > 1) {
                      temp <- cbind(1:nrow(temp), temp)
                      if(trace == T)
                      {
                          plot(NULL, xlim = range(temp[, 1]), ylim = range(temp[, 2]), 
                            type = "l", main = paste(variables[i], ": ", ifelse(k == 
                              1, "PO", "NPO"), sep = ""), xlab = "Index", ylab = "Value")
                          for (j in unique(temp[, 3])) lines(temp[temp[, 3] == j, 1], 
                            temp[temp[, 3] == j, 2], col = cols1[j])
                      }
                      if(density == T)
                      {
                          # plot density
                          plot(density(temp[, 2]), main = paste(variables[i], ": ", ifelse(k == 
                            1, "PO", ifelse(k == 2, "NPO", "Excluded")), sep = ""), xlab = "Value", 
                            ylab = "Density")
                      }
                    } else {
                      if(trace == T)
                      {
                          plot(0, 0, xaxt = "n", yaxt = "n", xlab = "", ylab = "", col = "white", 
                            main = paste(variables[i], ": ", ifelse(k == 1, "PO", "NPO"), 
                              sep = ""))
                          text(0, 0, "<= 1 Sample")
                      }
                      if(density == T)
                      {
                          plot(0, 0, xaxt = "n", yaxt = "n", xlab = "", ylab = "", col = "white", 
                            main = paste(variables[i], ": ", ifelse(k == 1, "PO", "NPO"), 
                              sep = ""))
                          text(0, 0, "<= 1 Sample")
                      }
                    }
                    l <- l + 1
                    if (ask == T) {
                      if ((l - 1)%%4 == 0 && (l - 1) != maxp) 
                        readline("Press any key to continue:")
                    }
                  }
                }
                par(mfrow = c(1, 1))
                return(cat(""))
            }
        }
        if (cutp == 1) {
            cutp <- list(NULL)
            for (i in 1:length(x$theta)) cutp[[i]] <- as.matrix(x$theta[[i]])
            output[[ind]] <- cutp
            ind <- ind + 1
        }
        if (status == 1) {
            if (!is.na(x$status[1])) {
                status <- list(NULL)
                for (i in 1:length(x$status)) status[[i]] <- as.matrix(x$status[[i]])
                output[[ind]] <- status
                ind <- ind + 1
            }
        }
        if (sdb == 1) {
            if (!is.na(x$sdb[1])) {
                sdb <- list(NULL)
                for (i in 1:length(x$sdb)) sdb[[i]] <- as.matrix(x$sdb[[i]])
                output[[ind]] <- sdb
                ind <- ind + 1
            }
        }
        if (varp == 1) {
            if (!is.na(x$varp[1])) {
                varp <- list(NULL)
                for (i in 1:length(x$varp)) varp[[i]] <- as.matrix(x$varp[[i]])
                output[[ind]] <- varp
                ind <- ind + 1
            }
        }
        if (clusterint == 1) {
            if (!is.na(x$psi[1])) {
                clusterint <- list(NULL)
                for (i in 1:length(x$psi)) {
                  clusterint[[i]] <- as.matrix(x$psi[[i]])
                  clusterint[[i]] <- cbind(clusterint[[i]], as.matrix(x$varp[[i]]))
                }
                output[[ind]] <- clusterint
                ind <- ind + 1
            } else {
                if (type != "all") 
                  stop("Can't plot as no cluster intercepts in this model")
            }
        }
        if (loglikelihood == 1) {
            loglikelihood <- list(NULL)
            for (i in 1:length(x$loglikelihood)) loglikelihood[[i]] <- as.matrix(x$loglikelihood[[i]])
            loglikelihood <- lapply(loglikelihood, function(y) {
                colnames(y) <- "loglikelihood"
                y
            })
            output[[ind]] <- loglikelihood
            ind <- ind + 1
        }
        
        # now convert list of outputs into a combined 'mcmc.list'
        output1 <- list(NULL)
        for (i in 1:length(output[[1]])) output1[[i]] <- as.mcmc(do.call("cbind", 
            lapply(output, function(y, i) y[[i]], i = i)))
        output1 <- as.mcmc.list(output1)
        if (conv == 1) 
            plot(output1, trace = trace, density = density, ask = ask)
        if (AC == 1) 
            autocorr.plot(output1, ask = ask)
    }
}

