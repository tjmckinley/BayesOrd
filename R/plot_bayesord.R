# plot function for 'bayesord' object
plot.bayesord <- function(object, which = c("trace", "AC"), type = c("reg", "rand", 
    "all"), scale = c("logOR", "OR"), ask = F, bystatus = F, threshold = 0.5, ...) {
    require(coda)
    if (class(object) != "bayesord") 
        stop("'object' is not a 'bayesord' object")
    if (which[1] != "trace" & which[1] != "AC") 
        stop("'which' is incorrect")
    if (type[1] != "reg" & type[1] != "rand" & type[1] != "all") 
        stop("'type' is incorrect")
    if (scale[1] != "logOR" & scale[1] != "OR") 
        stop("'scale' is incorrect")
    if (!is.logical(bystatus[1]) | length(bystatus) > 1) 
        stop("'bystatus' is not a logical value")
    
    reg <- ifelse(type[1] == "reg" | type[1] == "all", 1, 0)
    cutp <- ifelse(type[1] == "reg" | type[1] == "all", 1, 0)
    status <- ifelse(type[1] == "reg" | type[1] == "all", 1, 0)
    randint <- ifelse(type[1] == "rand" | type[1] == "all", 1, 0)
    varp <- ifelse(type[1] == "rand" | type[1] == "all", 1, 0)
    sdb <- ifelse(type[1] == "reg" | type[1] == "all", 1, 0)
    loglikelihood <- ifelse(type[1] == "reg" | type[1] == "all", 1, 0)
    
    trace <- ifelse(which[1] == "trace", 1, 0)
    CI <- ifelse(which[1] == "CI", 1, 0)
    AC <- ifelse(which[1] == "AC", 1, 0)
    
    if (trace == 1 | AC == 1) {
        output <- list(NULL)
        ind <- 1
        if (reg == 1) {
            reg <- list(NULL)
            for (i in 1:length(object$beta)) {
                reg[[i]] <- as.matrix(object$beta[[i]])
                if (scale[1] == "OR") 
                  reg[[i]] <- exp(reg[[i]])
            }
            output[[ind]] <- reg
            ind <- ind + 1
            
            # if we wish to plot traces by status
            if (bystatus == TRUE) {
                if (is.na(object$status[1])) 
                  stop("Can't plot by status")
                # extract names
                variables <- colnames(as.matrix(object$beta[[1]]))
                # extract statuses
                status <- list(NULL)
                for (i in 1:length(object$beta)) status[[i]] <- as.matrix(object$status[[i]])
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
                  for (k in 1:3) temp2[[k]] <- lapply(temp3, function(x, k) x[[k]], 
                    k = k)
                  for (k in 1:3) {
                    temp2[[k]] <- lapply(as.list(1:length(temp2[[k]])), function(i, 
                      x) {
                      x <- x[[i]]
                      x <- cbind(x, rep(i, length(x)))
                      x
                    }, x = temp2[[k]])
                    temp2[[k]] <- do.call("rbind", temp2[[k]])
                  }
                  output[[i]] <- temp2
                }
                # set up plotting parameters
                maxp <- length(output) * 2
                if (maxp >= 4) 
                  mfrow1 <- c(4, 2) else mfrow1 <- c(maxp, 2)
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
                      plot(NULL, xlim = range(temp[, 1]), ylim = range(temp[, 2]), 
                        type = "l", main = paste(variables[i], ": ", ifelse(k == 
                          1, "PO", "NPO"), sep = ""), xlab = "Index", ylab = "Value")
                      for (j in unique(temp[, 3])) lines(temp[temp[, 3] == j, 1], 
                        temp[temp[, 3] == j, 2], col = cols1[j])
                      # plot density
                      plot(density(temp[, 2]), main = paste(variables[i], ": ", ifelse(k == 
                        1, "PO", ifelse(k == 2, "NPO", "Excluded")), sep = ""), xlab = "Value", 
                        ylab = "Density")
                    } else {
                      plot(0, 0, xaxt = "n", yaxt = "n", xlab = "", ylab = "", col = "white", 
                        main = paste(variables[i], ": ", ifelse(k == 1, "PO", "NPO"), 
                          sep = ""))
                      text(0, 0, "None")
                      plot(0, 0, xaxt = "n", yaxt = "n", xlab = "", ylab = "", col = "white", 
                        main = paste(variables[i], ": ", ifelse(k == 1, "PO", "NPO"), 
                          sep = ""))
                      text(0, 0, "None")
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
            for (i in 1:length(object$theta)) cutp[[i]] <- as.matrix(object$theta[[i]])
            output[[ind]] <- cutp
            ind <- ind + 1
        }
        if (status == 1) {
            if (!is.na(object$status[1])) {
                status <- list(NULL)
                for (i in 1:length(object$status)) status[[i]] <- as.matrix(object$status[[i]])
                output[[ind]] <- status
                ind <- ind + 1
            }
        }
        if (sdb == 1) {
            if (!is.na(object$sdb[1])) {
                sdb <- list(NULL)
                for (i in 1:length(object$sdb)) sdb[[i]] <- as.matrix(object$sdb[[i]])
                output[[ind]] <- sdb
                ind <- ind + 1
            }
        }
        if (varp == 1) {
            if (!is.na(object$varp[1])) {
                varp <- list(NULL)
                for (i in 1:length(object$varp)) varp[[i]] <- as.matrix(object$varp[[i]])
                output[[ind]] <- varp
                ind <- ind + 1
            }
        }
        if (randint == 1) {
            if (!is.na(object$psi[1])) {
                randint <- list(NULL)
                for (i in 1:length(object$psi)) {
                  randint[[i]] <- as.matrix(object$psi[[i]])
                  randint[[i]] <- cbind(randint[[i]], as.matrix(object$varp[[i]]))
                }
                output[[ind]] <- randint
                ind <- ind + 1
            } else {
                if (type != "all") 
                  stop("Can't plot as no random intercepts in this model")
            }
        }
        if (loglikelihood == 1) {
            loglikelihood <- list(NULL)
            for (i in 1:length(object$loglikelihood)) loglikelihood[[i]] <- as.matrix(object$loglikelihood[[i]])
            loglikelihood <- lapply(loglikelihood, function(x) {
                colnames(x) <- "loglikelihood"
                x
            })
            output[[ind]] <- loglikelihood
            ind <- ind + 1
        }
        
        # now convert list of outputs into a combined 'mcmc.list'
        output1 <- list(NULL)
        for (i in 1:length(output[[1]])) output1[[i]] <- as.mcmc(do.call("cbind", 
            lapply(output, function(x, i) x[[i]], i = i)))
        output1 <- as.mcmc.list(output1)
        if (trace == 1) 
            plot(output1, ask = ask)
        if (AC == 1) 
            autocorr.plot(output1, ask = ask)
    }
    
    # reset indicators reg<-ifelse(type[1]=='reg' | type[1]=='all',1,0)
    # randint<-ifelse(type[1]=='rand' | type[1]=='all',1,0) if(CI==1) { if(reg==1) {
    # if(object$var.info$interactions==1) cat(paste('Be careful with interpretation
    # of log(ORs) and ORs due to interaction effect\n')) if(scale[1]=='logOR') {
    # if(CI==1) { if(trace==1) dev.new() #reorder for plotting
    # regno<-as.numeric(matrix(1:(object$info$nbetagroup*object$info$ntheta),object$info$ntheta,byrow=T))
    # #produce PPAs for each variable status1<-as.matrix(object$status)
    # reg<-as.matrix(object$beta) if(!is.na(status1[1])) { #produce model averaged
    # results reg<-as.matrix(object$beta) temp1<-matrix(NA,ncol(reg),3) k<-1 for(j in
    # 1:ncol(reg)) { temp<-reg[,j]#[status1[,k]<2,j] temp<-temp[!is.na(temp)]
    # if(length(temp)>0) temp1[j,]<-c(mean(temp),quantile(temp,probs=c(0.025,0.975)))
    # if(k%%ncol(status1)==0) k<-1 else k<-k+1 } } else
    # temp1<-t(apply(reg,2,function(x) c(mean(x),quantile(x,probs=c(0.025,0.975)))))
    # rownames(temp1)<-sapply(as.list(colnames(reg)),function(x,y)
    # paste(y,'(',x,')',sep=''),y=scale[1]) colnames(temp1)<-c('Mean','2.5%','97.5%')
    # if(length(regno)==nrow(temp1)) reg<-temp1[regno,] else reg<-temp1
    # if(!is.na(status1[1])) { status1<-apply(status1,2,function(x)
    # length(x[x!=2])/length(x)) plot.CI(reg,'Log(odds
    # ratios)',0,status1,threshold=threshold) } else plot.CI(reg,'Log(odds
    # ratios)',0,threshold=threshold) } } else { if(CI==1) { if(trace==1) dev.new()
    # #reorder for plotting
    # regno<-as.numeric(matrix(1:(object$info$nbetagroup*object$info$ntheta),object$info$ntheta,byrow=T))
    # #produce PPAs for each variable status1<-as.matrix(object$status)
    # reg<-as.matrix(object$beta) if(!is.na(status1[1])) { #produce model averaged
    # results reg<-as.matrix(object$beta) temp1<-matrix(NA,ncol(reg),3) k<-1 for(j in
    # 1:ncol(reg)) { temp<-reg[status1[,k]<2,j] temp<-temp[!is.na(temp)]
    # temp<-exp(temp) if(length(temp)>0)
    # temp1[j,]<-c(mean(temp),quantile(temp,probs=c(0.025,0.975)))
    # if(k%%ncol(status1)==0) k<-1 else k<-k+1 } } else
    # temp1<-t(apply(reg,2,function(x) c(mean(x),quantile(x,probs=c(0.025,0.975)))))
    # rownames(temp1)<-sapply(as.list(colnames(reg)),function(x,y)
    # paste(y,'(',x,')',sep=''),y=scale[1]) colnames(temp1)<-c('Mean','2.5%','97.5%')
    # reg<-temp1[regno,] if(!is.na(status1[1])) {
    # status1<-apply(status1,2,function(x) length(x[x!=2])/length(x))
    # plot.CI(reg,'Odds ratios',1,status1,threshold=threshold) } else
    # plot.CI(reg,'Odds ratios',1,threshold=threshold) } } } if(randint==1) {
    # if(!is.na(object$psi[1])) { if(trace==1 | reg==1) dev.new()
    # temp<-as.matrix(object$psi) temp1<-apply(temp,2,function(x)
    # c(mean(x),quantile(x,probs=c(0.025,0.975)))) rownames(temp1)<-colnames(temp)
    # colnames(temp1)<-c('Mean','2.5%','97.5%') plot.CI(temp1,'Random intercepts',0)
    # } else { if(type!='all') stop('Can't plot as no random intercepts in this
    # model') } } }
}

