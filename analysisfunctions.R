#functions to run Bayesian ordinal regression model

#generic function
bayesord<-function(formula, ...) UseMethod("bayesord")

#function for extracting terms on either side of vertical bars in formula
extract_terms<-function(term)
{
	if(!is.call(term)) stop("Must have a valid formula")
	t<-as.character(term)
	exp<-t[[3]]
	exp<-strsplit(exp," \\| ")[[1]]
	t[3]<-exp[1]
	t<-paste(t[2],t[1],t[3])
	RE<-NULL
	if(length(exp)>1)
	{
		RE<-exp[2]
		RE <- gsub(" ", "", RE)
		RE <- strsplit(RE, "\\+")[[1]]
	}
	t<-list(form=as.formula(t),RE=RE)
	t
}

#set formula object for "bayesord" method than incorporates a random intercept
bayesord.formula<-function(formula,data=list(),nchains=1,multi=F,model.type=c("PO","NPO","both"),var.select=TRUE,niter=50000,nsavecoda=1000,mnb=0,varb=1000,maxsdb=20,fixed=FALSE,vart=1,mnpsi=0,shvarp=0.01,rtvarp=0.01,propsdb=1,propsdt=1,propsdp=1,propsdvarp=1,runtraining=FALSE,nitertrain=1000,start=NA,end=NA,thin=NA, ...)
{
	if(!is.call(formula)) stop("Must have a valid formula")
	if(length(nchains)>1 | !is.numeric(nchains)) stop("'nchains' is not a numeric scalar")
	if(length(multi)>1 |!is.logical(multi)) stop("'multi' is not a logical value")
	if(model.type[1]!="PO" & model.type[1]!="NPO" & model.type[1]!="both") stop("'model.type' is incorrectly specified")
	if(length(niter)>1 | !is.numeric(niter)) stop("'niter' is not a numeric scalar")
	if(length(nsavecoda)>1 | !is.numeric(nsavecoda)) stop("'nsavecoda' is not a numeric scalar")
	if(length(mnb)>1 | !is.numeric(mnb)) stop("'mnb' is not a numeric scalar")
	if(length(varb)>1 | !is.numeric(varb)) stop("'varb' is not a numeric scalar")
	if(length(maxsdb)>1 | !is.numeric(maxsdb)) stop("'maxsdb' is not a numeric scalar")
	if(length(fixed)>1 |!is.logical(fixed)) stop("'fixed' is not a logical value")
	if(length(vart)>1 | !is.numeric(vart)) stop("'vart' is not a numeric scalar")
	if(length(mnpsi)>1 | !is.numeric(mnpsi)) stop("'mnpsi' is not a numeric scalar")
	if(length(shvarp)>1 | !is.numeric(shvarp)) stop("'shvarp' is not a numeric scalar")
	if(length(rtvarp)>1 | !is.numeric(rtvarp)) stop("'rtvarp' is not a numeric scalar")
	if(length(propsdb)>1 | !is.numeric(propsdb)) stop("'propsdb' is not a numeric scalar")
	if(length(propsdt)>1 | !is.numeric(propsdt)) stop("'propsdt' is not a numeric scalar")
	if(length(propsdp)>1 | !is.numeric(propsdp)) stop("'propsdp' is not a numeric scalar")
	if(length(propsdvarp)>1 | !is.numeric(propsdvarp)) stop("'propsdvarp' is not a numeric scalar")
	if(length(runtraining)>1 |!is.logical(runtraining)) stop("'runtraining' is not a logical value")
	if((length(start)>1 | !is.numeric(start)) & !is.na(start)) stop("'start' is not a numeric scalar or NA")
	if(length(nitertrain)>1 | !is.numeric(nitertrain)) stop("'nitertrain' is not a numeric scalar")
	if((length(end)>1 | !is.numeric(end)) & !is.na(end)) stop("'end' is not a numeric scalar or NA")
	if((length(thin)>1 | !is.numeric(thin)) & !is.na(thin)) stop("'thin' is not a numeric scalar or NA")
	
	#extract formula	
	form<-extract_terms(formula)
#	browser()
	#extract random intercepts term if required
	RE<-form[[2]]
	form<-form[[1]]
	mf<-model.frame(formula=form,data=data,na.action=na.fail)
	#check there are no columns called "RE" or "counts"
	temp<-match(c("RE","counts"),attr(mf,"names"))
	if(length(temp[!is.na(temp)])>0) stop("Can't name variables 'RE' or 'counts'")
	
	#create vector for random intercepts if required
	if(is.null(RE)) rand.int.vec<-NA
	else
	{
		rand.int.vec<-data[,match(RE,colnames(data)), drop = F]
		for(j in 1:ncol(rand.int.vec)) if(!is.factor(rand.int.vec[, j])) stop("Random intercepts term is not a factor")
		for(j in 1:ncol(rand.int.vec)) rand.int.vec[, j] <- as.numeric(rand.int.vec[, j])
	}
	
	#condense data into succinct form for model
	data<-data[,match(attr(mf,"names"),colnames(data))]
	if(!is.null(RE)) data$RE<-factor(rand.int.vec)
	#now aggregate data
	data<-aggregate(rep(1,nrow(data)),data,table)
	data[,ncol(data)]<-as.numeric(data[,ncol(data)])
	colnames(data)[ncol(data)]<-"counts"
	#create design matrix from standard formula
	mf<-model.frame(formula=form,data=data,na.action=na.fail)
	#remove intercept if not already removed
	if(attr(attr(mf,"terms"),"intercept")==0) stop("Must have intercept in formula")
	x<-model.matrix(form,data=mf)
	
	#set up response variable
	y<-model.response(mf)
	if(!is.factor(y)) stop("Response variable is not a factor")
	
	#now create matrices for variable selection
	variables<-attr(attr(mf,"terms"),"term.labels")
	variables<-variables[variables!="RE"]
	variables<-variables[variables!="counts"]
	nvariables<-length(variables)
	
	xassign<-attr(x,"assign")[-1]
	xassign<-rev((length(xassign):1)[!duplicated(rev(xassign))])
	xassign<-c(0,xassign) #add nbetagroup to end of this later
	
	#remove intercept from design matrix
	x<-x[,-1]
	
	#extract variable information for use in marginal posterior plotting
	var.info<-as.character(attr(attr(mf,"terms"),"variables"))[-1]
	interactions<-attr(attr(mf,"terms"),"order")
	intpresent<-ifelse(length(interactions[interactions>1])>0,1,0)
	var.info<-list(variables=var.info,orig.data=data,interactions=intpresent,assign=attributes(x)$assign[-1])
		
	#dealing with interaction effects
	if(intpresent==1)
	{
		#create intfactor matrix
		temp<-variables
		temp<-lapply(as.list(temp),function(x) strsplit(x,":")[[1]])
		temp1<-matrix(0,length(temp),length(temp))
		for(j in 1:length(temp))
		{
			if(length(temp[[j]])>1)
			{
				for(i in 1:(j-1))
				{
					if(length(temp[[i]])<length(temp[[j]]))
					{
						temp2<-match(temp[[j]],temp[[i]])
						if(length(temp2[!is.na(temp2)])>0) temp1[i,j]<-1
					}
				}
			}
		}
		write.table(temp1,"intfactor.txt",row.names=F,col.names=F,sep="\t",quote=F)
		temp<-sapply(temp,length)-1
		write.table(temp,"interaction.txt",row.names=F,col.names=F,sep="\t",quote=F)
		temp1<-c((1:length(temp))[!duplicated(temp)]-1,length(temp))
		write.table(temp1,"intstart.txt",row.names=F,col.names=F,sep="\t",quote=F)
		maxinteraction<-max(temp)
		rm(temp,temp1,temp2)
	}
	else maxinteraction<-0
			
	#now run ordinal regression model
	
	#load 'coda' library
	require(coda)
	
	#set options for writing out data	
	options(scipen=999)
			
	#extract data in correct format
	N<-length(y)
	ntheta<-length(levels(y))-1
		
	#extract response variable
	response<-as.numeric(y)-1
	#set random intercepts terms
	if(!is.null(RE))
	{
		id<-as.numeric(data$RE)
		if(length(id)==0) stop("ID variable is not present")
		RE<-1
	}
	else
	{
		id<-rep(1,nrow(x))
		RE<-0
	}
	id<-as.numeric(id)
	npsi<-length(unique(id))
	
	#extract explanatory variables and convert to correct format
	dat<-x
	temp.lev<-colnames(dat)
		
	#produce string of linear terms
	nbeta<-length(temp.lev)*ntheta
	nbetagroup<-length(temp.lev)
	
	#produce argument list for likelihood function
	dat1<-data.frame(dat,counts=data$counts,response,id)
		
	#sort dat1 according to RE terms
	if(RE==1)
	{
		var.info$orig.data<-var.info$orig.data[sort.list(dat1[,ncol(dat1)]),]
		dat1<-dat1[sort.list(dat1[,ncol(dat1)]),]
		psicount<-cumsum(tapply(dat1[,ncol(dat1)],factor(dat1[,ncol(dat1)],levels=unique(dat1[,ncol(dat1)])),length))
		write.table(psicount,"psicount.txt",row.names=F,col.names=F,sep="\t",quote=F)
	}		
#	browser()
	#compile C code
	system("cp ../ccode/* .")
	system("gcc -Wall main.c index.c loglikelihood.c loglikelihoodpsi.c moveNPO.c moveNPOtoPO.c movePO.c movePOtoNPO.c moveinctoexc.c moveexctoinc.c validitycheck.c movetheta.c movesdb.c -lgsl -lgslcblas -lm -O3 -o MCMCmulti")
	system("rm *.c functions.h")
			
	#now output data file
	write.table(dat1,"data.dat",row.names=F,col.names=F,sep="\t",quote=F)
	
	#assign data file for use in model checking
	if(RE==0) id<-NA
	
	#assign model type
	if(model.type[1]=="PO") model.type<-0
	else
	{
		if(model.type[1]=="NPO") model.type<-1
		else model.type<-2
	}
		
	#object for output
	info<-as.data.frame(t(as.matrix(c(n=N,niter=niter,nchains=nchains,nbetagroup=nbetagroup,ntheta=ntheta,mnb=mnb,varb=varb,maxsdb=maxsdb,fixed=ifelse(fixed==TRUE,1,0),vart=vart,mnpsi=mnpsi,shvarp=shvarp,rtvarp=rtvarp,npsi=npsi,RE=RE,model.type=model.type,var.select=var.select,runtraining=runtraining,nitertrain=ifelse(runtraining==T,nitertrain,NA)))))
	
	#now write priors file
	priors<-data.frame(n=N,niter=niter,nsavecoda=nsavecoda,nvariables=nvariables,nbeta=nbeta,nbetagroup=nbetagroup,ntheta=ntheta,mnb=mnb,varb=varb,maxsdb=maxsdb,fixed=ifelse(fixed==TRUE,1,0),vart=vart,propsdb=propsdb,propsdt=propsdt,mnpsi=mnpsi,shvarp=shvarp,rtvarp=rtvarp,propsdp=propsdp,propsdvarp=propsdvarp,npsi=npsi,RE=RE,model.type=model.type,varselect=ifelse(var.select==TRUE,1,0),maxinteraction=maxinteraction,runtraining=ifelse(runtraining==TRUE,1,0),nitertrain=nitertrain)
	priors<-data.frame(t(as.matrix(priors)),colnames(priors))
	write.table(priors,"priors.txt",row.names=F,col.names=F,quote=F,sep="\t")
	
	#write assignment file
	write.table(xassign,"xassign.txt",row.names=F,col.names=F,quote=F,sep="\t")
			
	if(multi==T)
	{
		library(multicore)
		#run model
		mclapply(as.list(1:nchains),function(x)
		{
			seed<-floor(2147483647*runif(1,0,1))
			temp<-paste("./MCMCmulti ",x," ",seed,sep="")#," > output",x,".txt",sep="")
			system(temp)
		},mc.cores=nchains,mc.set.seed=T)
	}
	else
	{
		for(i in 1:nchains)
		{
			seed<-floor(2147483647*runif(1,0,1))
			temp<-paste("./MCMCmulti ",i," ",seed,sep="")#" > output",i,".txt",sep="")
			system(temp)
		}
	}
	print("Finished MCMC run")
	
	options(scipen=5)
	
	#read in and process data
	beta<-list(NULL)
	theta<-list(NULL)
	psi<-list(NULL)
	status<-list(NULL)
	sdb<-list(NULL)
	varp<-list(NULL)
	loglikelihood<-list(NULL)
	for(i in 1:nchains) 
	{
		#read in table of posterior samples
		beta[[i]]<-read.table(paste("codaMCMC",i,".txt",sep=""),header=T)			
		theta[[i]]<-as.mcmc(beta[[i]][,(nbeta+1):(nbeta+ntheta)])
		if(model.type==2 | var.select==TRUE)
		{
			status[[i]]<-beta[[i]][,(nbeta+ntheta+1):(nbeta+ntheta+nbetagroup)]
			if(is.null(ncol(status[[i]]))) status[[i]]<-matrix(status[[i]],ncol=1)
			colnames(status[[i]])<-paste("status_",temp.lev,sep="")
			status[[i]]<-as.mcmc(status[[i]])
		}
		else status[[i]]<-NA
		psi[[i]]<-as.mcmc(beta[[i]][,(nbeta+ntheta+nbetagroup+1):(nbeta+ntheta+nbetagroup+npsi)])
		if(fixed==FALSE)
		{
			sdb[[i]]<-as.mcmc(beta[[i]][,(nbeta+ntheta+nbetagroup+npsi+1):(2*nbeta+ntheta+nbetagroup+npsi)])
			varp[[i]]<-as.mcmc(beta[[i]][,(2*nbeta+ntheta+nbetagroup+npsi+1)])
			loglikelihood[[i]]<-as.mcmc(beta[[i]][,(2*nbeta+ntheta+nbetagroup+npsi+2)])
		}
		else
		{
			sdb[[i]]<-NA
			varp[[i]]<-as.mcmc(beta[[i]][,(nbeta+ntheta+nbetagroup+npsi+1)])
			loglikelihood[[i]]<-as.mcmc(beta[[i]][,(nbeta+ntheta+nbetagroup+npsi+2)])
		}
		
		if(model.type==0)
		{
			#remove extraneous columns if proportional odds model
			beta[[i]]<-beta[[i]][,1:nbetagroup]
			if(is.null(ncol(beta[[i]]))) beta[[i]]<-matrix(beta[[i]],ncol=1)
			colnames(beta[[i]])<-temp.lev
			
			if(fixed==FALSE)
			{
				sdb[[i]]<-as.mcmc(as.matrix(sdb[[i]])[,1:nbetagroup])
			}
		}
		else 
		{
			beta[[i]]<-beta[[i]][,1:nbeta]
			if(is.null(ncol(beta[[i]]))) beta[[i]]<-matrix(beta[[i]],ncol=1)
			colnames(beta[[i]])<-as.vector(apply(cbind(matrix(rep(temp.lev,ntheta),nrow=ntheta,byrow=T),1:ntheta),1,function(x){y<-x[-length(x)];y<-paste(y,"_",x[length(x)],sep="");y}))
		}
		beta[[i]]<-as.mcmc(beta[[i]])
	}
	if(model.type==2 | var.select==TRUE) status<-as.mcmc.list(status)
	else status<-NA
	if(fixed==FALSE) sdb<-as.mcmc.list(sdb)
	else sdb<-NA
	if(RE==0)
	{
		varp<-NA
		psi<-NA
	}
	else
	{
		varp<-mcmc.list(varp)
		psi<-mcmc.list(psi)
	}
	ans<-list(beta=mcmc.list(beta),theta=mcmc.list(theta),status=status,psi=psi,sdb=sdb,varp=varp,loglikelihood=mcmc.list(loglikelihood),info=info,model.dat=dat1,var.info=var.info)
	#remove text files
	system("rm data.dat priors.txt xassign.txt codaMCMC* MCMCmulti")
	if(RE==1) system("rm psicount.txt")
	if(intpresent==1) system("rm interaction.txt intfactor.txt intstart.txt")
	#set class
	class(ans)<-"bayesord"
	#thin if necessary
	if(missing(start)) start<-NA
	if(missing(end)) end<-NA
	if(missing(thin)) thin<-NA
	ans<-window.bayesord(ans,start=start,end=end,thin=thin)
	#output call and formula for printing
	ans$call<-match.call()
	ans$formula<-formula
	ans
}

#print function for "bayesord" object
print.bayesord<-function(x, ...)
{
	if(class(x)!="bayesord") stop("'x' is not a 'bayesord' object")
	
	if(x$info$model.type==0) cat("\nPROPORTIONAL ODDS MODEL\n")
	else
	{
		if(x$info$model.type==1) cat("\nNON-PROPORTIONAL ODDS MODEL\n")
		else cat("\nPROPORTIONAL ODDS OR NON-PROPORTIONAL ODDS MODEL\n")
	}
	if(x$info$var.select==T) cat("\nVariable selection specified\n")
	
	cat("\nCall:\n")
	print(x$formula)
	
	cat("\nData:\n")
	cat(paste(sum(x$var.info$orig.data$counts),"observations with",x$info$nbetagroup,ifelse(x$info$nbetagroup>1,"explanatory variables\n","explanatory variable\n")))
	cat(paste("Response variable has",x$info$ntheta+1,"ordered groups:",paste(levels(x$var.info$orig.data[,1]),collapse=", "),"\n"))
	
	cat("\nMCMC details:\n")
	cat(paste(x$info$nchains,"chains of",x$info$niter,"iterations\n"))
		
	cat("\nPrior information:\n")
	cat(paste("Regression variables: N(",x$info$mnb,",",x$info$varb,")\n",sep=""))
	cat(paste("Cut-points: N(0,",x$info$vart,")\n",sep=""))
	if(x$info$RE==1)
	{
		cat(paste(x$info$npsi," random intercepts, each with: N(",x$info$mnpsi,",varp)\n",sep=""))
		cat(paste("varp: G(",x$info$shvarp,",",x$info$rtvarp,") [i.e. mean=",x$info$shvarp/x$info$rtvarp,", variance=",x$info$shvarp/(x$info$rtvarp^2),"]\n",sep=""))
	}	
}

#summary function for "bayesord" object
summary.bayesord<-function(object,type=c("reg","rand","all"),scale=c("logOR","OR"),threshold=0,topmodels=5,digits=2, ...)
{
	require(coda)
	if(class(object)!="bayesord") stop("'object' is not a 'bayesord' object")
	if(type[1]!="reg" & type[1]!="rand" & type[1]!="all") stop("'type' is incorrect")
	if(scale[1]!="logOR" & scale[1]!="OR") stop("'scale' is incorrect")
	if(!is.numeric(topmodels) | length(topmodels)>1) stop("'topmodels' is incorrect")
	if(!is.numeric(digits) | length(digits)>1) stop("'digits' is incorrect")
	if(!is.numeric(threshold) | length(threshold)>1) stop("'threshold' is incorrect")
	
	reg<-ifelse(type[1]=="reg" | type[1]=="all",1,0)
	cutp<-ifelse(type[1]=="reg" | type[1]=="all",1,0)
	randint<-ifelse(type[1]=="rand" | type[1]=="all",1,0)
	varp<-ifelse(type[1]=="rand" | type[1]=="all",1,0)
	sdb<-ifelse(type[1]=="reg" | type[1]=="all",1,0)
	loglikelihood<-ifelse(type[1]=="reg" | type[1]=="all",1,0)
		
	status<-ifelse(type[1]=="reg" | type[1]=="all",1,0)
	
	#adjust possibilities by model.type
	if(object$info$model.type!=2 & object$info$var.select!=1) status<-0
	if(object$info$fixed==1) sdb<-0
			
	output<-list(NULL)
	ind<-1
	if(reg==1)
	{
		reg<-as.matrix(object$beta)
		temp<-reg
		if(scale[1]=="OR") temp<-exp(temp)
		temp<-as.mcmc(temp)
		ess<-effectiveSize(temp)
		temp<-summary(temp)
		temp1<-cbind(temp$statistics[,1:2],temp$quantiles[,c(3,1,5)],ess,temp$statistics[,4]/temp$statistics[,2])
		rownames(temp1)<-sapply(as.list(colnames(reg)),function(x,y) paste(y,"(",x,")",sep=""),y=scale[1])
		colnames(temp1)<-c("Mean","SD","Median","2.5%","97.5%","ESS","MC error/SD")
		if(status==1)
		{
#			#if variable selection has been used then average across non-zero models
#			
#			#produce PPAs for each variable
#			status1<-as.matrix(object$status)
#			#produce model averaged results
#			reg<-as.matrix(object$beta)
#			temp1<-matrix(NA,ncol(reg),7)
#			k<-1
#			for(j in 1:ncol(reg))
#			{
#				temp<-reg[status1[,k]<2,j]
#				if(scale[1]=="OR") temp<-exp(temp)
#				if(length(temp)>1)
#				{
#					temp<-as.mcmc(temp)
#					ess<-effectiveSize(temp)
#					temp<-summary(temp)
#					temp1[j,]<-c(temp$statistics[1:2],temp$quantiles[c(3,1,5)],ess,temp$statistics[4]/temp$statistics[2])
#				}
#				if(k%%ncol(status1)==0) k<-1
#				else k<-k+1
#			}
#			rownames(temp1)<-sapply(as.list(colnames(reg)),function(x,y) paste(y,"(",x,")",sep=""),y=scale[1])
#			colnames(temp1)<-c("Mean","SD","Median","2.5%","97.5%","ESS","MC error/SD")
			
			status<-list(NULL)
			for(i in 1:length(object$status))
			{
				#produce models via PPAs
				status[[i]]<-as.matrix(object$status[[i]])
				lev<-colnames(status[[i]])
				lev<-sapply(as.list(lev),function(x){x<-strsplit(x,"status_")[[1]];x[x!=""]})
				#extract unique models
				status[[i]]<-table(factor(apply(status[[i]],1,function(x) paste(x,collapse=""))))
				status[[i]]<-status[[i]]/sum(status[[i]])
				status1<-matrix(NA,length(status[[i]]),length(lev)+1)
				for(j in 1:length(status[[i]])) status1[j,]<-c(strsplit(names(status[[i]])[j],"")[[1]],status[[i]][j])
				status1[,-ncol(status1)]<-t(apply(status1[,-ncol(status1), drop = F],1,function(x)
				{
					x[x=="0"]<-"PO"
					x[x=="1"]<-"NPO"
					x[x=="2"]<-"0"
					x
				}))
				status1<-as.data.frame(status1)
				colnames(status1)<-c(lev,"PPAs")
				status[[i]]<-status1
				status[[i]]<-status[[i]][sort.list(status[[i]]$PPAs,decreasing=T),]
			}
			
			status1<-list(NULL)
			for(i in 1:length(object$status))
			{
				status1[[i]]<-as.matrix(object$status[[i]])
				status1[[i]]<-apply(status1[[i]],2,function(x) table(factor(x,levels=c("0","1","2"))))
				status1[[i]]<-apply(status1[[i]],2,function(x) x/sum(x))
				rownames(status1[[i]])[rownames(status1[[i]])=="0"]<-"PO"
				rownames(status1[[i]])[rownames(status1[[i]])=="1"]<-"NPO"
				rownames(status1[[i]])[rownames(status1[[i]])=="2"]<-"Excluded"
				colnames(status1[[i]])<-lev
			}
		}
		else
		{
#			reg<-as.matrix(object$beta)
#			temp<-reg
#			if(scale[1]=="OR") temp<-exp(temp)
#			temp<-as.mcmc(temp)
#			ess<-effectiveSize(temp)
#			temp<-summary(temp)
#			temp1<-cbind(temp$statistics[,1:2],temp$quantiles[,c(3,1,5)],ess,temp$statistics[,4]/temp$statistics[,2])
#			rownames(temp1)<-sapply(as.list(colnames(reg)),function(x,y) paste(y,"(",x,")",sep=""),y=scale[1])
#			colnames(temp1)<-c("Mean","SD","Median","2.5%","97.5%","ESS","MC error/SD")
							
			status<-NA
			status1<-NA
		}
		if(object$info$model.type!=0)
		{
			regno<-as.numeric(matrix(1:(object$info$nbetagroup*object$info$ntheta),object$info$ntheta,byrow=T))
			temp1<-temp1[regno,]
		}
		output[[ind]]<-temp1
		ind<-ind+1
	}
	else
	{
		output[[ind]]<-NA
		ind<-ind+1
	}
	if(cutp==1)
	{
		cutp<-as.matrix(object$theta)
		temp<-cutp
		temp<-as.mcmc(temp)
		ess<-effectiveSize(temp)
		temp<-summary(temp)
		temp1<-cbind(temp$statistics[,1:2],temp$quantiles[,c(3,1,5)],ess,temp$statistics[,4]/temp$statistics[,2])
		colnames(temp1)<-c("Mean","SD","Median","2.5%","97.5%","ESS","MC error/SD")
		output[[ind]]<-temp1
		ind<-ind+1
	}
	else
	{
		output[[ind]]<-NA
		ind<-ind+1
	}
	if(sdb==1)
	{
		sdb<-as.matrix(object$sdb)
		temp<-sdb
		temp<-as.mcmc(temp)
		ess<-effectiveSize(temp)
		temp<-summary(temp)
		temp1<-cbind(temp$statistics[,1:2],temp$quantiles[,c(3,1,5)],ess,temp$statistics[,4]/temp$statistics[,2])
		colnames(temp1)<-c("Mean","SD","Median","2.5%","97.5%","ESS","MC error/SD")
		if(object$info$model.type!=0)
		{
			regno<-as.numeric(matrix(1:(object$info$nbetagroup*object$info$ntheta),object$info$ntheta,byrow=T))
			temp1<-temp1[regno,]
		}
		output[[ind]]<-temp1
		ind<-ind+1
	}
	else
	{
		output[[ind]]<-NA
		ind<-ind+1
	}
	if(varp==1)
	{
		if(!is.na(object$varp[1]))
		{
			varp<-as.matrix(object$varp)
			temp<-varp
			temp<-as.mcmc(temp)
			ess<-effectiveSize(temp)
			temp<-summary(temp)
			temp1<-cbind(temp$statistics[,1:2],temp$quantiles[,c(3,1,5)],ess,temp$statistics[,4]/temp$statistics[,2])
			colnames(temp1)<-c("Mean","SD","Median","2.5%","97.5%","ESS","MC error/SD")
			output[[ind]]<-temp1
			ind<-ind+1
		}
	}
	else
	{
		output[[ind]]<-NA
		ind<-ind+1
	}
	if(randint==1)
	{
		if(!is.na(object$psi[1]))
		{
			randint<-as.matrix(object$psi[[i]])
			randint<-cbind(randint,as.matrix(object$varp))
			temp<-randint
			temp<-as.mcmc(temp)
			ess<-effectiveSize(temp)
			temp<-summary(temp)
			temp1<-cbind(temp$statistics[,1:2],temp$quantiles[,c(3,1,5)],ess,temp$statistics[,4]/temp$statistics[,2])
			colnames(temp1)<-c("Mean","SD","Median","2.5%","97.5%","ESS","MC error/SD")
			output[[ind]]<-temp1
			ind<-ind+1
		}
		else
		{
			if(type!="all") stop("Can't plot as no random intercepts in this model")
		}
	}
	else
	{
		output[[ind]]<-NA
		ind<-ind+1
	}
	if(loglikelihood==1)
	{
		loglikelihood<-as.matrix(object$loglikelihood)
		colnames(loglikelihood)<-"loglikelihood"
		temp<-loglikelihood
		temp<-as.mcmc(temp)
		ess<-effectiveSize(temp)
		temp<-summary(temp)
		temp1<-matrix(c(temp$statistics[1:2],temp$quantiles[c(3,1,5)],ess,temp$statistics[4]/temp$statistics[2]),nrow=1)
		colnames(temp1)<-c("Mean","SD","Median","2.5%","97.5%","ESS","MC error/SD")
		output[[ind]]<-temp1
		ind<-ind+1
	}
	else
	{
		output[[ind]]<-NA
		ind<-ind+1
	}
	
	#output summary objects for printing
	output<-list(output=output,models=status,variables=status1,threshold=threshold,topmodels=topmodels,digits=digits)
	class(output)<-"summary.bayesord"
	output
}

#print function for "summary.bayesord" object
print.summary.bayesord<-function(x, ...)
{
	if(class(x)!="summary.bayesord") stop("'x' is not a 'summary.bayesord' object")
	#set up summary parameters
	threshold<-x$threshold
	digits<-x$digits
	topmodels<-x$topmodels
	
	cat("\n")
	if(!is.na(x$output[[1]][1]))
	{
		x1<-signif(x$output[[1]],digits=digits)
		if(is.list(x$variables))
		{	
			x$variables<-lapply(as.list(1:length(x$variables)),function(i,x)
			{
				x[[i]]<-data.frame(chain=rep(i,nrow(x[[i]])),structure=rownames(x[[i]]),x[[i]])
			},x=x$variables)
			x$variables<-do.call("rbind",x$variables)
			ymean<-matrix(NA,3,ncol(x$variables)-2)
			ysd<-matrix(NA,3,ncol(x$variables)-2)
			z<-c("PO","NPO","Excluded")
			for(i in 1:3)
			{
				ymean[i,]<-apply(x$variables[x$variables$structure==z[i],-(1:2)],2,mean)
				ysd[i,]<-apply(x$variables[x$variables$structure==z[i],-(1:2)],2,sd)
			}
			rownames(ymean)<-z
			rownames(ysd)<-z
			ymean<-t(ymean)
			ysd<-t(ysd)
			ymean<-signif(ymean,digits=digits)
			ysd<-signif(ysd,digits=digits)
			y<-matrix(NA,ncol(x$variables)-2,3)
			temp<-rep(NA,3)
			for(i in 1:(ncol(x$variables)-2))
			{	
				for(j in 1:3) temp[j]<-paste(ymean[i,j]," (",ysd[i,j],")",sep="")
				y[i,]<-temp
			}
			rownames(y)<-colnames(x$variables)[-(1:2)]
			colnames(y)<-z
			variables<-t(y)
			
			#now print only those variables with PPA>threshold
			if(threshold>0)
			{
				status<-1-ymean[,3]
				ntheta<-nrow(x1)/length(status)
				#extract only variables with PPA above threshold
				if(length(status)<nrow(x1))
				{
					status1<-rep(1:length(status),each=ntheta)
					status1[rep(ifelse(status>threshold,1,0),each=ntheta)==0]<-NA
					x1<-x1[!is.na(status1),]
#					status<-status[status>threshold]
				}
				else
				{
					x1<-x1[status>threshold,]
#					status<-status[status>threshold]
				}
			}
			cat("########### COEFFICIENTS ###########\n")
			if(threshold>0) cat(paste("########### PPAs > ",threshold," ###########\n",sep=""))
			print(x1)
			cat("\n")
			
			#now produce top models and variable PPAs
			topmodels<-min(topmodels,nrow(x$models))
			cat(paste("Top",topmodels,"models according to PPA:\n"))
			y<-do.call("rbind",x$models)
			y$group<-factor(apply(y[-ncol(y)],1,function(x) paste(x,collapse="+")))
			y$PPAs<-as.numeric(as.character(y$PPAs))
			y<-aggregate(y$PPAs,list(y$group),function(x,n)
			{
				if(length(x)!=n) x<-c(x,rep(0,n-length(x)))
				c(mean(x),sd(x))
			},n=length(x$models))
			y<-data.frame(model=y[,1],mean=y[,2][,1],sd=y[,2][,2])
			y<-y[sort.list(y[,2],decreasing=T),]
			y<-y[1:topmodels,]
			y[,-1]<-signif(y[,-1],digits=digits)
			y<-cbind(do.call("rbind",lapply(as.list(as.character(y[,1])),function(x) strsplit(x,"\\+")[[1]])),y[,-1])
			y<-t(y)
			y[y=="0"]<-""
			rownames(y)<-c(colnames(variables),"Mean PPA","SD PPA")
			colnames(y)<-paste("M",1:topmodels,sep="")
			y<-as.data.frame(y)
			cat("########## MODELS ##########\n")
			print(y,quote=F)
			cat("\n")
			cat("########## VARIABLES ##########\n")
			print(t(variables),quote=F)
			cat("\n")
		}
		else
		{
			cat("########### COEFFICIENTS ###########\n")
			print(x1)
			cat("\n")
		}
	}
	if(!is.na(x$output[[2]][1]))
	{
		cat("############ CUT-POINTS ############\n")
		print(signif(x$output[[2]],digits=digits))
		cat("\n")
	}
	if(!is.na(x$output[[3]][1]))
	{
		#now print only those variables with PPA>threshold
		if(is.list(x$variables) & threshold>0)
		{
			x2<-rownames(x$output[[1]])
			x2<-sapply(as.list(x2),function(x) strsplit(x,"OR\\(")[[1]][2])
			x2<-sapply(as.list(x2),function(x) strsplit(x,"\\)")[[1]][1])
			x2<-paste("SD_",x2,sep="")
			x1<-signif(x$output[[3]],digits=digits)
			rownames(x1)<-x2
			status<-1-ymean[,3]
			ntheta<-nrow(x1)/length(status)
			#extract only variables with PPA above threshold
			if(length(status)<nrow(x1))
			{
				status1<-rep(1:length(status),each=ntheta)
				status1[rep(ifelse(status>threshold,1,0),each=ntheta)==0]<-NA
				x1<-x1[!is.na(status1),]
#					status<-status[status>threshold]
			}
			else
			{
				x1<-x1[status>threshold,]
#					status<-status[status>threshold]
			}
			cat("########### STD. DEV. ###########\n")
			cat(paste("########### PPAs > ",threshold," ###########\n",sep=""))
			print(x1)
			cat("\n")
		}
		else
		{
			x2<-rownames(x1)
			x2<-sapply(as.list(x2),function(x) strsplit(x,"OR\\(")[[1]][2])
			x2<-sapply(as.list(x2),function(x) strsplit(x,"\\)")[[1]][1])
			x2<-paste("SD_",x2,sep="")
			x1<-signif(x$output[[3]],digits=digits)
			rownames(x1)<-x2
			cat("########### STD. DEV. ###########\n")
			print(x1)
			cat("\n")
		}
	}
	if(!is.na(x$output[[4]][1]))
	{
		cat("############# ST. DEV. #############\n")
		print(signif(x$output[[4]],digits=digits))
		cat("\n")
	}
	if(!is.na(x$output[[5]][1]))
	{
		cat("########## RANDOM EFFECTS ##########\n")
		print(signif(x$output[[5]],digits=digits))
		cat("\n")
	}
	if(!is.na(x$output[[6]][1]))
	{
		cat("########## LOG-LIKELIHOOD ##########\n")
		print(signif(x$output[[6]],digits=digits))
		cat("\n")
	}
}

#function for plotting credible intervals
#plot.CI<-function(x, ylab1, comp.line, status, threshold=0.5, ...)
#{
##	if(!is.mcmc(x) & !is.mcmc.list(x)) stop("'x' function in 'plot.CI' is not an 'mcmc' or 'mcmc.list' object")
#	if(missing(ylab1)) stop("'ylab1' argument is missing in 'plot.CI'")
#	
##	x<-as.matrix(x)
##	x<-t(apply(x,2,function(x) c(mean(x),quantile(x,probs=c(0.025,0.975)))))
#	if(!missing(status))
#	{
#		ntheta<-nrow(x)/length(status)
#		#extract only variables with PPA above threshold
#		if(length(status)<nrow(x))
#		{
#			status1<-rep(1:length(status),each=ntheta)
#			status1[rep(ifelse(status>threshold,1,0),each=ntheta)==0]<-NA
#			x<-x[!is.na(status1),]
#			status<-status[status>threshold]
#		}
#		else
#		{
#			x<-x[status>threshold,]
#			status<-status[status>threshold]
#		}
#	}

#	#extend margins
#	par(mar=c(8, 4, 4, 2)+0.1)
#	ylim1<-range(x[!is.na(x)])
#	if(!missing(status)) ylim1[2]<-ylim1[2]*1.1
#	#plot CIs
#	plot(1:nrow(x),x[,1],ylim=ylim1,pch=20,xlab="",ylab=ylab1,xaxt="n",main=ifelse(!missing(status),paste("Posterior means and 95% CIs (PPAs > ",threshold,")",sep=""),"Posterior means and 95% CIs"))
#	for(i in 1:nrow(x)) 
#	{
#		lines(c(i,i),x[i,2:3])
#		lines(c(i-0.05,i+0.05),rep(x[i,2],2))
#		lines(c(i-0.05,i+0.05),rep(x[i,3],2))
#	}
#	if(!missing(comp.line)) abline(h=comp.line,lty=2)
#	#add dividing lines if necessary
#	lev<-rownames(x)
#	lev<-sapply(as.list(lev),function(x)
#	{
#		x<-strsplit(x,"_")[[1]]
#		x<-x[-length(x)]
#		if(length(x)>1) x<-paste(x,collapse="_")
#		x
#	})
#	lev<-sapply(as.list(lev),function(x)
#	{
#		x<-strsplit(x,"\\(")[[1]][-1]
#		if(length(x)>1) x<-paste(x,collapse="(")
#		x
#	})
##	lev<-sapply(as.list(lev),function(x)
##	{
##		x<-strsplit(x,"\\)")[[1]]
##		x<-x[-length(x)]
##		if(length(x)>1) x<-paste(x,collapse=")")
##		x
##	})
#	lev<-unique(lev)
#	if(!missing(status))
#	{
#		if(length(status)<nrow(x))
#		{
#			for(i in 1:(length(status)-1)) abline(v=(i*ntheta+0.5),lty=2)
#			axis(1,at=1:nrow(x),labels=rep(1:ntheta,length(status)),las=2)
##			lev<-names(status)
##			lev<-sapply(as.list(lev),function(x){x<-strsplit(x,"status_")[[1]];x[x!=""]})
#			for(i in 1:length(lev))
#			{
#				mtext(lev[i],side=1,line=2,at=(i-1+0.5)*ntheta+0.5,las=2)
#				text((i-1+0.5)*ntheta+0.5,ylim1[2],labels=signif(status[i],digits=2),font=2)
#			}
#		}
#		else
#		{
#			axis(1,at=1:nrow(x),labels=lev,las=2)
#			text(1:nrow(x),rep(ylim1[2],nrow(x)),labels=signif(status,digits=2),font=2)
#		}
#	}
#	else
#	{
#		if(length(lev)<nrow(x))
#		{
#			ntheta<-nrow(x)/length(lev)
#			for(i in 1:(length(lev)-1)) abline(v=(i*ntheta+0.5),lty=2)
#			axis(1,at=1:nrow(x),labels=rep(1:ntheta,length(lev)),las=2)
#			for(i in 1:length(lev)) mtext(lev[i],side=1,line=2,at=(i-1+0.5)*ntheta+0.5,las=2)
#		}
#		else axis(1,at=1:nrow(x),labels=lev,las=2)
#	}
#	#reduce margins
#	par(mar=c(5, 4, 4, 2)+0.1)
#}

#plot function for "bayesord" object
plot.bayesord<-function(object,which=c("trace","AC"),type=c("reg","rand","all"),scale=c("logOR","OR"),ask=F,bystatus=F,threshold=0.5, ...)
{
	require(coda)
	if(class(object)!="bayesord") stop("'object' is not a 'bayesord' object")
	if(which[1]!="trace" & which[1]!="AC") stop("'which' is incorrect")
	if(type[1]!="reg" & type[1]!="rand" & type[1]!="all") stop("'type' is incorrect")
	if(scale[1]!="logOR" & scale[1]!="OR") stop("'scale' is incorrect")
	if(!is.logical(bystatus[1]) | length(bystatus)>1) stop("'bystatus' is not a logical value")
	
	reg<-ifelse(type[1]=="reg" | type[1]=="all",1,0)
	cutp<-ifelse(type[1]=="reg" | type[1]=="all",1,0)
	status<-ifelse(type[1]=="reg" | type[1]=="all",1,0)
	randint<-ifelse(type[1]=="rand" | type[1]=="all",1,0)
	varp<-ifelse(type[1]=="rand" | type[1]=="all",1,0)
	sdb<-ifelse(type[1]=="reg" | type[1]=="all",1,0)
	loglikelihood<-ifelse(type[1]=="reg" | type[1]=="all",1,0)
	
	trace<-ifelse(which[1]=="trace",1,0)
	CI<-ifelse(which[1]=="CI",1,0)
	AC<-ifelse(which[1]=="AC",1,0)
	
	if(trace==1 | AC==1)
	{
		output<-list(NULL)
		ind<-1
		if(reg==1)
		{
			reg<-list(NULL)
			for(i in 1:length(object$beta))
			{
				reg[[i]]<-as.matrix(object$beta[[i]])
				if(scale[1]=="OR") reg[[i]]<-exp(reg[[i]])		
			}
			output[[ind]]<-reg
			ind<-ind+1
			
			#if we wish to plot traces by status
			if(bystatus==TRUE)
			{
				if(is.na(object$status[1])) stop("Can't plot by status")
				#extract names
				variables<-colnames(as.matrix(object$beta[[1]]))
				#extract statuses
				status<-list(NULL)
				for(i in 1:length(object$beta)) status[[i]]<-as.matrix(object$status[[i]])
				#bind multiple chains together for each variable
				reg<-output[[ind-1]]
				n<-nrow(reg[[1]])
				reg<-do.call("rbind",reg)
				#bind statuses together for each variable
				status<-do.call("rbind",status)
				#now convert to correct format for plotting
				output<-list(NULL)
				l<-1
				for(i in 1:ncol(reg))
				{
#					print(paste("i=",i,"l=",l))
					temp<-matrix(reg[,i],nrow=n)
					tempstatus<-matrix(status[,l],nrow=n)
					l<-l+1
					if(l%%(ncol(status)+1)==0) l<-1
					temp3<-list(NULL)
					for(j in 1:ncol(temp))
					{
						temp2<-list(NULL)
						for(k in 0:2) temp2[[k+1]]<-temp[tempstatus[,j]==k,j]
						temp3[[j]]<-temp2
					}
					#reorder
					temp2<-list(NULL)
					for(k in 1:3) temp2[[k]]<-lapply(temp3,function(x,k) x[[k]],k=k)
					for(k in 1:3)
					{
						temp2[[k]]<-lapply(as.list(1:length(temp2[[k]])),function(i,x)
						{
							x<-x[[i]]
							x<-cbind(x,rep(i,length(x)))
							x
						},x=temp2[[k]])
						temp2[[k]]<-do.call("rbind",temp2[[k]])
					}
					output[[i]]<-temp2
				}
				#set up plotting parameters
				maxp<-length(output)*2
				if(maxp>=4) mfrow1<-c(4,2)
				else mfrow1<-c(maxp,2)
				cols1<-c("black","red","blue","green","yellow","purple")
				par(mfrow=mfrow1)
				#produce plots
				l<-1
				for(i in 1:length(output))
				{
					for(k in 1:2)
					{
						#plot trace
						temp<-output[[i]][[k]]
						if(nrow(temp)>1)
						{
							temp<-cbind(1:nrow(temp),temp)
							plot(NULL,xlim=range(temp[,1]),ylim=range(temp[,2]),type="l",main=paste(variables[i],": ",ifelse(k==1,"PO","NPO"),sep=""),xlab="Index",ylab="Value")
							for(j in unique(temp[,3])) lines(temp[temp[,3]==j,1],temp[temp[,3]==j,2],col=cols1[j])
							#plot density
							plot(density(temp[,2]),main=paste(variables[i],": ",ifelse(k==1,"PO",ifelse(k==2,"NPO","Excluded")),sep=""),xlab="Value",ylab="Density")
						}
						else
						{
							plot(0,0,xaxt="n",yaxt="n",xlab="",ylab="",col="white",main=paste(variables[i],": ",ifelse(k==1,"PO","NPO"),sep=""))
							text(0,0,"None")
							plot(0,0,xaxt="n",yaxt="n",xlab="",ylab="",col="white",main=paste(variables[i],": ",ifelse(k==1,"PO","NPO"),sep=""))
							text(0,0,"None")
						}
						l<-l+1
						if(ask==T)
						{
							if((l-1)%%4==0 && (l-1)!=maxp) readline("Press any key to continue:")
						}
					}	
				}
				par(mfrow=c(1,1))
				return(cat(""))
			}
		}
		if(cutp==1)
		{
			cutp<-list(NULL)
			for(i in 1:length(object$theta)) cutp[[i]]<-as.matrix(object$theta[[i]])
			output[[ind]]<-cutp
			ind<-ind+1
		}
		if(status==1)
		{
			if(!is.na(object$status[1]))
			{
				status<-list(NULL)
				for(i in 1:length(object$status)) status[[i]]<-as.matrix(object$status[[i]])
				output[[ind]]<-status
				ind<-ind+1
			}
		}
		if(sdb==1)
		{
			if(!is.na(object$sdb[1]))
			{
				sdb<-list(NULL)
				for(i in 1:length(object$sdb)) sdb[[i]]<-as.matrix(object$sdb[[i]])
				output[[ind]]<-sdb
				ind<-ind+1
			}
		}
		if(varp==1)
		{
			if(!is.na(object$varp[1]))
			{
				varp<-list(NULL)
				for(i in 1:length(object$varp)) varp[[i]]<-as.matrix(object$varp[[i]])
				output[[ind]]<-varp
				ind<-ind+1
			}
		}
		if(randint==1)
		{
			if(!is.na(object$psi[1]))
			{
				randint<-list(NULL)
				for(i in 1:length(object$psi))
				{
					randint[[i]]<-as.matrix(object$psi[[i]])
					randint[[i]]<-cbind(randint[[i]],as.matrix(object$varp[[i]]))
				}
				output[[ind]]<-randint
				ind<-ind+1
			}
			else
			{
				if(type!="all") stop("Can't plot as no random intercepts in this model")
			}
		}
		if(loglikelihood==1)
		{
			loglikelihood<-list(NULL)
			for(i in 1:length(object$loglikelihood)) loglikelihood[[i]]<-as.matrix(object$loglikelihood[[i]])
			loglikelihood<-lapply(loglikelihood,function(x){colnames(x)<-"loglikelihood";x})
			output[[ind]]<-loglikelihood
			ind<-ind+1
		}		
		
		#now convert list of outputs into a combined "mcmc.list"
		output1<-list(NULL)
		for(i in 1:length(output[[1]]))	output1[[i]]<-as.mcmc(do.call("cbind",lapply(output,function(x,i) x[[i]],i=i)))
		output1<-as.mcmc.list(output1)
		if(trace==1) plot(output1,ask=ask)
		if(AC==1) autocorr.plot(output1,ask=ask)
	}
	
	#reset indicators
#	reg<-ifelse(type[1]=="reg" | type[1]=="all",1,0)
#	randint<-ifelse(type[1]=="rand" | type[1]=="all",1,0)
#	if(CI==1)
#	{
#		if(reg==1)
#		{
#			if(object$var.info$interactions==1) cat(paste("Be careful with interpretation of log(ORs) and ORs due to interaction effect\n"))
#			if(scale[1]=="logOR")
#			{
#				if(CI==1)
#				{
#					if(trace==1) dev.new()
#					#reorder for plotting
#					regno<-as.numeric(matrix(1:(object$info$nbetagroup*object$info$ntheta),object$info$ntheta,byrow=T))					
#					#produce PPAs for each variable
##					status1<-as.matrix(object$status)
#					reg<-as.matrix(object$beta)
##					if(!is.na(status1[1]))
##					{
##						#produce model averaged results
###						reg<-as.matrix(object$beta)
##						temp1<-matrix(NA,ncol(reg),3)
##						k<-1
##						for(j in 1:ncol(reg))
##						{
##							temp<-reg[,j]#[status1[,k]<2,j]
##							temp<-temp[!is.na(temp)]
##							if(length(temp)>0) temp1[j,]<-c(mean(temp),quantile(temp,probs=c(0.025,0.975)))
##							if(k%%ncol(status1)==0) k<-1
##							else k<-k+1
##						}
##					}
##					else 
#					temp1<-t(apply(reg,2,function(x) c(mean(x),quantile(x,probs=c(0.025,0.975)))))
#					rownames(temp1)<-sapply(as.list(colnames(reg)),function(x,y) paste(y,"(",x,")",sep=""),y=scale[1])
#					colnames(temp1)<-c("Mean","2.5%","97.5%")
#					if(length(regno)==nrow(temp1)) reg<-temp1[regno,]
#					else reg<-temp1
#					if(!is.na(status1[1]))
#					{
#						status1<-apply(status1,2,function(x) length(x[x!=2])/length(x))
#						plot.CI(reg,"Log(odds ratios)",0,status1,threshold=threshold)
#					}
#					else plot.CI(reg,"Log(odds ratios)",0,threshold=threshold)
#				}
#			}	
#			else
#			{	
#				if(CI==1)
#				{
#					if(trace==1) dev.new()
#					#reorder for plotting
#					regno<-as.numeric(matrix(1:(object$info$nbetagroup*object$info$ntheta),object$info$ntheta,byrow=T))					
#					#produce PPAs for each variable
#					status1<-as.matrix(object$status)
#					reg<-as.matrix(object$beta)
#					if(!is.na(status1[1]))
#					{
#						#produce model averaged results
#						reg<-as.matrix(object$beta)
#						temp1<-matrix(NA,ncol(reg),3)
#						k<-1
#						for(j in 1:ncol(reg))
#						{
#							temp<-reg[status1[,k]<2,j]
#							temp<-temp[!is.na(temp)]
#							temp<-exp(temp)
#							if(length(temp)>0) temp1[j,]<-c(mean(temp),quantile(temp,probs=c(0.025,0.975)))
#							if(k%%ncol(status1)==0) k<-1
#							else k<-k+1
#						}
#					}
#					else temp1<-t(apply(reg,2,function(x) c(mean(x),quantile(x,probs=c(0.025,0.975)))))
#					rownames(temp1)<-sapply(as.list(colnames(reg)),function(x,y) paste(y,"(",x,")",sep=""),y=scale[1])
#					colnames(temp1)<-c("Mean","2.5%","97.5%")
#					reg<-temp1[regno,]
#					if(!is.na(status1[1]))
#					{
#						status1<-apply(status1,2,function(x) length(x[x!=2])/length(x))
#						plot.CI(reg,"Odds ratios",1,status1,threshold=threshold)
#					}
#					else plot.CI(reg,"Odds ratios",1,threshold=threshold)
#				}
#			}
#		}
#		if(randint==1)
#		{
#			if(!is.na(object$psi[1]))
#			{
#				if(trace==1 | reg==1) dev.new()
#				temp<-as.matrix(object$psi)
#				temp1<-apply(temp,2,function(x) c(mean(x),quantile(x,probs=c(0.025,0.975))))
#				rownames(temp1)<-colnames(temp)
#				colnames(temp1)<-c("Mean","2.5%","97.5%")
#				plot.CI(temp1,"Random intercepts",0)
#			}
#			else
#			{
#				if(type!="all") stop("Can't plot as no random intercepts in this model")
#			}
#		}
#	}
}

#extract posterior "fitted" distributions for "bayesord" object
fitted.bayesord<-function(x,multi=F,mc.cores=getOption("cores"),RE=TRUE, ...)
{
	require(coda)
	if(class(x)!="bayesord") stop("'x' is not a 'bayesord' object")
		
	#extract posteriors
	beta<-as.matrix(x$beta)
	theta<-as.matrix(x$theta)
	
	#augment beta posterior if PO model
	if(x$info$model.type==0)
	{
		temp.lev<-colnames(beta)
		beta<-matrix(rep(as.vector(beta),ncol(theta)),nrow=nrow(beta))
		colnames(beta)<-as.vector(apply(cbind(matrix(rep(temp.lev,ncol(theta)),nrow=ncol(theta),byrow=T),1:ncol(theta)),1,function(x){y<-x[-length(x)];y<-paste(y,"_",x[length(x)],sep="");y}))
	}
	if(!is.na(x$psi[1]) & RE==TRUE)
	{
		psi<-as.matrix(x$psi)
		mcmc<-cbind(beta,theta,psi)
	}
	else
	{
		psi<-NA
		mcmc<-cbind(beta,theta)
	}
	nbeta<-ncol(beta)
	ntheta<-ncol(theta)
		
	if(multi==T)
	{
		require(multicore)
		x.list<-as.list(as.data.frame(t(x$model.dat)))
	
		#produce predicted posteriors for individuals in data set
		pred<-mclapply(x.list,function(dat,mcmc,nbeta,ntheta)
		{
			id<-dat[length(dat)]
			dat<-dat[-length(dat)]
			response<-dat[length(dat)]
			dat<-dat[-length(dat)]
			counts<-dat[length(dat)]
			dat<-dat[-length(dat)]
			if(!is.na(psi[1])) 
			{
				psi<-mcmc[,(nbeta+ntheta+1):ncol(mcmc)]
				psi<-psi[,id]
			}
			else psi<-numeric(nrow(mcmc))
		
			mcmc<-cbind(mcmc[,1:(nbeta+ntheta)],psi)
		
			pred<-apply(mcmc,1,function(mcmc,dat,counts,nbeta,ntheta)
			{	
				psi<-mcmc[length(mcmc)]
				mcmc<-mcmc[-length(mcmc)]
				beta<-mcmc[1:nbeta]
				theta<-mcmc[-(1:nbeta)]
			
				beta<-matrix(beta,nrow=ntheta,byrow=T)
				beta<-apply(beta,1,function(x,dat) x*dat,dat=dat)
				if(!is.null(nrow(beta))) beta<-apply(beta,2,sum)
				beta<-beta+psi
				beta<-theta-beta
				gamma<-exp(beta)/(1+exp(beta))
				p<-c(gamma[1],gamma[-1]-gamma[-length(gamma)],1-gamma[length(gamma)])
				if(length(p[p<0])>0) browser()
				#now simulate a response score
				g<-rmultinom(1,size=counts,prob=p)
#				g<-apply(g,2,function(x) which(x==1)-1)
#				g<-summary(factor(g,levels=0:(length(p)-1)))
				g
			},dat=dat,counts=counts,nbeta=nbeta,ntheta=ntheta)
			list(pred)
		},mcmc=mcmc,nbeta=nbeta,ntheta=ntheta,mc.cores=mc.cores)
	}
	else
	{
		#produce predicted posteriors for individuals in data set
		pred<-apply(x$model.dat,1,function(dat,mcmc,nbeta,ntheta)
		{
			id<-dat[length(dat)]
			dat<-dat[-length(dat)]
			response<-dat[length(dat)]
			dat<-dat[-length(dat)]
			counts<-dat[length(dat)]
			dat<-dat[-length(dat)]
			if(!is.na(psi[1])) 
			{
				psi<-mcmc[,(nbeta+ntheta+1):ncol(mcmc)]
				psi<-psi[,id]
			}
			else psi<-numeric(nrow(mcmc))
		
			mcmc<-cbind(mcmc[,1:(nbeta+ntheta)],psi)
		
			pred<-apply(mcmc,1,function(mcmc,dat,counts,nbeta,ntheta)
			{
				psi<-mcmc[length(mcmc)]
				mcmc<-mcmc[-length(mcmc)]
				beta<-mcmc[1:nbeta]
				theta<-mcmc[-(1:nbeta)]
			
				beta<-matrix(beta,nrow=ntheta,byrow=T)
				beta<-apply(beta,1,function(x,dat) x*dat,dat=dat)
				if(!is.null(nrow(beta))) beta<-apply(beta,2,sum)
				beta<-beta+psi
				beta<-theta-beta
				gamma<-exp(beta)/(1+exp(beta))
				p<-c(gamma[1],gamma[-1]-gamma[-length(gamma)],1-gamma[length(gamma)])
				if(length(p[p<0])>0) browser()
				#now simulate a response score
				g<-rmultinom(1,size=counts,prob=p)
#				g<-apply(g,2,function(x) which(x==1)-1)
#				g<-summary(factor(g,levels=0:(length(p)-1)))
				g
			},dat=dat,counts=counts,nbeta=nbeta,ntheta=ntheta)
			list(pred)
		},mcmc=mcmc,nbeta=nbeta,ntheta=ntheta)
	}
	pred<-lapply(pred,function(x) x[[1]])
	#output predictions
	pred<-list(fits=pred,model.dat=x$model.dat,var.info=x$var.info)
	class(pred)<-"fitted.bayesord"
	pred
}

#print function for "fitted.bayesord" object
print.fitted.bayesord<-function(x, ...)
{
	if(class(x)!="fitted.bayesord") stop("'x' is not a 'fitted.bayesord' object")
	
	x<-summary(x)
	print(x)
}

#summary function for "fitted.bayesord" object
summary.fitted.bayesord<-function(object, ...)
{
	if(class(object)!="fitted.bayesord") stop("'object' is not a 'fitted.bayesord' object")
	#calculate posterior summaries for proportions correctly predicted
	ind.summary<-sapply(as.list(1:length(object$fits)),function(i,x,y)
	{
		x<-x[[i]]
		y<-y[i]
		x<-x[y+1,]
		x
	},x=object$fits,y=object$model.dat[,(ncol(object$model.dat)-1)])
	ind.summary<-apply(ind.summary,1,sum)/sum(object$var.info$orig.data$counts)
	ind.summary<-c(mean(ind.summary),quantile(ind.summary,probs=c(0.025,0.975)))
	class(ind.summary)<-"summary.fitted.bayesord"
	ind.summary
}

#print function for "summary.fitted.bayesord" object
print.summary.fitted.bayesord<-function(x,digits=2, ...)
{
	if(class(x)!="summary.fitted.bayesord") stop("'x' is not a 'summary.fitted.bayesord' object")
	
	x<-round(x,digits=digits)
	
	cat(paste("Posterior mean probability of correct classification (95% CI): ",x[1]," (",paste(x[2:3],collapse=","),")\n",sep=""))
}

##predict function for "bayesord" object
#predict.bayesord<-function(object,newdata=NULL,multi=FALSE, ...)
#{
#	if(class(object)!="bayesord") stop("'object' is not a 'bayesord' object")
#	if(is.null(newdata)) stop("'newdata' not specified")
#	
#	#since model has been fitted using formula 
#	form<-extract_terms(object$formula)
#	form<-form[[1]]
#	x<-model.matrix(form,newdata)
#	#swap model.dat over to newdata and remove RE terms
#	object$model.dat<-x
#	object$psi<-NA
#	object$varp<-NA
#	x<-fitted(object,predict=TRUE,multi=multi)
#	x
#}

#plot function for "fitted.bayesord" objects
plot.fitted.bayesord<-function(object, variable="all", ...)
{
	if(class(object)!="fitted.bayesord") stop("'object' is not a 'fitted.bayesord' object")
	plot.marginal.bayesord(object,variable=variable, ...)
}

##plot function for "predict.bayesord" objects
#plot.predict.bayesord<-function(object, variable="all", ...)
#{
#	if(class(object)!="predict.bayesord") stop("'object' is not a 'predict.bayesord' object")
#	plot.marginal.bayesord(object,variable=variable, ...)
#}
	
#function to plot fitted distributions
plot.marginal.bayesord<-function(object,variable="all", ...)
{
	if(class(object)!="fitted.bayesord" & class(object)!="predict.bayesord") stop("'object' is not a 'fitted.bayesord' or a 'predict.bayesord' object")
	if(!is.character(variable)) stop("'variable' not a character")
	
	if(variable=="all")
	{
		#extract variables to plot from object
		variable<-object$var.info$variables[-1]
	}
	else
	{
		x<-match(variable,object$var.info$variables[-1])
		if(length(x[is.na(x)])>1) stop(paste("One or more of 'variable' not found in data set",sep=""))
	}
	x<-numeric(length(variable))
	for(i in 1:length(variable)) x[i]<-ifelse(!is.factor(object$var.info$orig.data[,match(variable[i],colnames(object$var.info$orig.data))]),0,1)
	nonfactor<-variable[x==0]
	if(length(nonfactor)>0) cat(paste("The following variables are non-categorical and have been removed:\n",paste(nonfactor,collapse="\n"),"\n"))
	variable<-variable[x==1]
	if(length(variable)==0) stop("No non-categorical variables in list")
	
	#set up plotting parameters
	if(length(variable)>1)
	{
		if(length(variable)==2) mfrow1<-c(2,1)
		else mfrow1<-c(ceiling(length(variable)/2),2)
	}
	else mfrow1<-c(1,1)
	par(mfrow=mfrow1)
	
	#plot each variable in turn
	data1<-object$var.info$orig.data
	response<-data1[,1]
	data1<-data1[,-1]
	for(i in 1:length(variable))
	{	
		var<-data1[,match(c(variable[i],"counts"),colnames(data1))]
		var$response<-response
		
		#'aggregate' removes rows with zeros so as a fudge bind some dummy values to the
		#bottom before aggregating
		tempvar<-expand.grid(levels(var[,1]),levels(response))
		tempvar<-data.frame(tempvar[,1],rep(1,nrow(tempvar)),tempvar[,2])
		colnames(tempvar)<-colnames(var)
		tempvar<-rbind(var,tempvar)
		temp1<-aggregate(tempvar$counts,list(tempvar$response,tempvar[,1]),sum)
		temp1[,3]<-temp1[,3]-1
		temp1<-matrix(as.numeric(temp1$x),length(levels(var$response)))
		rownames(temp1)<-levels(var$response)
		colnames(temp1)<-levels(var[,1])
		temp1<-t(temp1)
		temp1<-t(apply(temp1,1,function(x) x/sum(x)))
		#plot output
		temp1<-barplot(temp1,beside=T,legend=T,ylim=c(0,1),main=variable[i])		
		#add posterior model fits
		for(j in 1:length(levels(var[,1])))
		{
			temp<-object$fits[var[,1]==levels(var[,1])[j]]
			temp<-lapply(as.list(1:nrow(temp[[1]])),function(i,x) sapply(x,function(x,i) x[i,],i=i),x=temp)
			temp<-sapply(temp,function(x) apply(x,1,sum))
			temp<-apply(temp,1,function(x) x/sum(x))
			temp<-apply(temp,1,function(x) c(mean(x),quantile(x,probs=c(0.025,0.975))))
			
			#now add predict posterior summaries to plot
			for(k in 1:ncol(temp1))
			{
				points(temp1[j,k],temp[1,k],pch=20)
				lines(rep(temp1[j,k],2),temp[2:3,k])
				lines(c(temp1[j,k]-0.05,temp1[j,k]+0.05),rep(temp[2,k],2))
				lines(c(temp1[j,k]-0.05,temp1[j,k]+0.05),rep(temp[3,k],2))
			}
		}
	}
	#reset graphical parameters
	par(mfrow=c(1,1))
}

#window function for "bayesord" objects
window.bayesord<-function(x,start=NA,end=NA,thin=NA,chains=NA, ...)
{
	require(coda)
	if(class(x)!="bayesord") stop("'x' is not a 'bayesord' object")
	if(!is.na(chains))
	{
		x$beta<-x$beta[chains]
		x$theta<-x$theta[chains]
		x$status<-x$status[chains]
		x$status<-x$sdb[chains]
		if(!is.na(x$psi[1]))
		{	
			x$psi<-x$psi[chains]
			x$varp<-x$varp[chains]
		}
		x$loglikelihood<-x$loglikelihood[chains]
		x$info$nchains<-length(chains)
	}
	if(is.na(start))
	{
		if(is.na(end))
		{
			if(is.na(thin))
			{	
				x<-x
			}
			else
			{
				for(i in 1:length(x))
				{
					if(is.mcmc.list(x[[i]])) x[[i]]<-window(x[[i]],thin=thin)
				}
			}
		}
		else
		{
			if(is.na(thin))
			{
				for(i in 1:length(x))
				{
					if(is.mcmc.list(x[[i]])) x[[i]]<-window(x[[i]],end=end)
				}
			}
			else
			{
				for(i in 1:length(x))
				{
					if(is.mcmc.list(x[[i]])) x[[i]]<-window(x[[i]],end=end,thin=thin)
				}
			}
		}
	}
	else
	{
		if(is.na(end))
		{
			if(is.na(thin))
			{
				for(i in 1:length(x))
				{
					if(is.mcmc.list(x[[i]])) x[[i]]<-window(x[[i]],start=start)
				}
			}
			else
			{
				for(i in 1:length(x))
				{
					if(is.mcmc.list(x[[i]])) x[[i]]<-window(x[[i]],start=start,thin=thin)
				}
			}
		}
		else
		{
			if(is.na(thin))
			{
				for(i in 1:length(x))
				{
					if(is.mcmc.list(x[[i]])) x[[i]]<-window(x[[i]],start=start,end=end)
				}
			}
			else
			{
				for(i in 1:length(x))
				{
					if(is.mcmc.list(x[[i]])) x[[i]]<-window(x[[i]],start=start,end=end,thin=thin)
				}
			}
		}
	}
	x
}




