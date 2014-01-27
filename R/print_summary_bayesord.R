# print function for 'summary.bayesord' object
print.summary.bayesord <- function(x, ...) {
    if (class(x) != "summary.bayesord") 
        stop("'x' is not a 'summary.bayesord' object")
    # set up summary parameters
    threshold <- x$threshold
    digits <- x$digits
    topmodels <- x$topmodels
    
    cat("\n")
    if (!is.na(x$output[[1]][1])) {
        x1 <- signif(x$output[[1]], digits = digits)
        if (is.list(x$variables)) {
            x$variables <- lapply(as.list(1:length(x$variables)), function(i, x) {
                x[[i]] <- data.frame(chain = rep(i, nrow(x[[i]])), structure = rownames(x[[i]]), 
                  x[[i]])
            }, x = x$variables)
            x$variables <- do.call("rbind", x$variables)
            ymean <- matrix(NA, 3, ncol(x$variables) - 2)
            ysd <- matrix(NA, 3, ncol(x$variables) - 2)
            z <- c("PO", "NPO", "Excluded")
            for (i in 1:3) {
                ymean[i, ] <- apply(x$variables[x$variables$structure == z[i], -(1:2)], 
                  2, mean)
                ysd[i, ] <- apply(x$variables[x$variables$structure == z[i], -(1:2)], 
                  2, sd)
            }
            rownames(ymean) <- z
            rownames(ysd) <- z
            ymean <- t(ymean)
            ysd <- t(ysd)
            ymean <- signif(ymean, digits = digits)
            ysd <- signif(ysd, digits = digits)
            y <- matrix(NA, ncol(x$variables) - 2, 3)
            temp <- rep(NA, 3)
            for (i in 1:(ncol(x$variables) - 2)) {
                for (j in 1:3) temp[j] <- paste(ymean[i, j], " (", ysd[i, j], ")", 
                  sep = "")
                y[i, ] <- temp
            }
            rownames(y) <- colnames(x$variables)[-(1:2)]
            colnames(y) <- z
            variables <- t(y)
            
            # now print only those variables with PPA>threshold
            if (threshold > 0) {
                status <- 1 - ymean[, 3]
                ntheta <- nrow(x1)/length(status)
                # extract only variables with PPA above threshold
                if (length(status) < nrow(x1)) {
                  status1 <- rep(1:length(status), each = ntheta)
                  status1[rep(ifelse(status > threshold, 1, 0), each = ntheta) == 
                    0] <- NA
                  x1 <- x1[!is.na(status1), ]
                  # status<-status[status>threshold]
                } else {
                  x1 <- x1[status > threshold, ]
                  # status<-status[status>threshold]
                }
            }
            cat("########### COEFFICIENTS ###########\n")
            if (threshold > 0) 
                cat(paste("########### PPAs > ", threshold, " ###########\n", sep = ""))
            print(x1)
            cat("\n")
            
            # now produce top models and variable PPAs
            topmodels <- min(topmodels, nrow(x$models))
            cat(paste("Top", topmodels, "models according to PPA:\n"))
            y <- do.call("rbind", x$models)
            y$group <- factor(apply(y[-ncol(y)], 1, function(x) paste(x, collapse = "+")))
            y$PPAs <- as.numeric(as.character(y$PPAs))
            y <- aggregate(y$PPAs, list(y$group), function(x, n) {
                if (length(x) != n) 
                  x <- c(x, rep(0, n - length(x)))
                c(mean(x), sd(x))
            }, n = length(x$models))
            y <- data.frame(model = y[, 1], mean = y[, 2][, 1], sd = y[, 2][, 2])
            y <- y[sort.list(y[, 2], decreasing = T), ]
            y <- y[1:topmodels, ]
            y[, -1] <- signif(y[, -1], digits = digits)
            y <- cbind(do.call("rbind", lapply(as.list(as.character(y[, 1])), function(x) strsplit(x, 
                "\\+")[[1]])), y[, -1])
            y <- t(y)
            y[y == "0"] <- ""
            rownames(y) <- c(colnames(variables), "Mean PPA", "SD PPA")
            colnames(y) <- paste("M", 1:topmodels, sep = "")
            y <- as.data.frame(y)
            cat("########## MODELS ##########\n")
            print(y, quote = F)
            cat("\n")
            cat("########## VARIABLES ##########\n")
            print(t(variables), quote = F)
            cat("\n")
        } else {
            cat("########### COEFFICIENTS ###########\n")
            print(x1)
            cat("\n")
        }
    }
    if (!is.na(x$output[[2]][1])) {
        cat("############ CUT-POINTS ############\n")
        print(signif(x$output[[2]], digits = digits))
        cat("\n")
    }
    if (!is.na(x$output[[3]][1])) {
        # now print only those variables with PPA>threshold
        if (is.list(x$variables) & threshold > 0) {
            x2 <- rownames(x$output[[1]])
            x2 <- sapply(as.list(x2), function(x) strsplit(x, "OR\\(")[[1]][2])
            x2 <- sapply(as.list(x2), function(x) strsplit(x, "\\)")[[1]][1])
            x2 <- paste("SD_", x2, sep = "")
            x1 <- signif(x$output[[3]], digits = digits)
            rownames(x1) <- x2
            status <- 1 - ymean[, 3]
            ntheta <- nrow(x1)/length(status)
            # extract only variables with PPA above threshold
            if (length(status) < nrow(x1)) {
                status1 <- rep(1:length(status), each = ntheta)
                status1[rep(ifelse(status > threshold, 1, 0), each = ntheta) == 0] <- NA
                x1 <- x1[!is.na(status1), ]
                # status<-status[status>threshold]
            } else {
                x1 <- x1[status > threshold, ]
                # status<-status[status>threshold]
            }
            cat("########### STD. DEV. ###########\n")
            cat(paste("########### PPAs > ", threshold, " ###########\n", sep = ""))
            print(x1)
            cat("\n")
        } else {
            x2 <- rownames(x1)
            x2 <- sapply(as.list(x2), function(x) strsplit(x, "OR\\(")[[1]][2])
            x2 <- sapply(as.list(x2), function(x) strsplit(x, "\\)")[[1]][1])
            x2 <- paste("SD_", x2, sep = "")
            x1 <- signif(x$output[[3]], digits = digits)
            rownames(x1) <- x2
            cat("########### STD. DEV. ###########\n")
            print(x1)
            cat("\n")
        }
    }
    if (!is.na(x$output[[4]][1])) {
        cat("############# ST. DEV. #############\n")
        print(signif(x$output[[4]], digits = digits))
        cat("\n")
    }
    if (!is.na(x$output[[5]][1])) {
        cat("########## RANDOM EFFECTS ##########\n")
        print(signif(x$output[[5]], digits = digits))
        cat("\n")
    }
    if (!is.na(x$output[[6]][1])) {
        cat("########## LOG-LIKELIHOOD ##########\n")
        print(signif(x$output[[6]], digits = digits))
        cat("\n")
    }
}
