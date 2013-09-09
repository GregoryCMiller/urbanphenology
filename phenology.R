#' Phenology calculation functions and summaries

PRISM_PPT <- function(info,nc){
  
  #PRISM Precip
  #======================
  
  apri <- apply(get.var.ncdf(nc, 'APRI'), 3, rbind)
  uyear <- unique(info$dates$year)
  APRI <- array(NA, dim=c(nrow(apri), length(uyear)))
  for (i in seq(along=uyear)){
    idx <- info$dates$year == uyear[i]
    APRI[,i] <- apply(apri[,idx], 1, sum, na.rm=T)
  }
  
  df <- data.frame(year=as.vector(outer( rep(1,nrow(apri)), uyear)))
  
  df$city <- info$abbrev
  subs <- get.var.ncdf(nc, info$landcover$name)[,,1]
  df$cat <- as.vector(outer(subs, rep(1,length(uyear))))
  df$cat <- factor(df$cat,levels=info$coverkey$CODE,labels=info$coverkey$NAME)
  df$catcode <- as.integer(df$cat)
  df$ANNPRISMPPT <- as.vector(APRI)
  sdf <- subset(df, catcode %in% c(1,5))
  
  #write.table(df, file='PRISM_pix.csv', sep=',', col.names=NA)
  write.table(sdf, file='PRISM_pix_15.csv', sep=',', col.names=NA)
  
}
#' Calculate Phenology metrics with growing season precipitation
#'  
#' SOS,EOS : start / end of season (highest and lowest derivative)
#' PEAKDAY : day of year of highest value
#' AMP     : annual amplitude (max-min)
#' SUM     : integrated evi
#' 
PhenMetrics <- function(info, nc, varname, n_plot=200){	
	cat('\nGrowing Season Phenology Metrics...')
	#browser()
	
	evi <- apply(get.var.ncdf(nc, 'evi'), 3, rbind)
	
	pdf('evihist.pdf')
	hist(evi,breaks=seq(-33000,33000,100))
	abline(v=0,col='red')
	dev.off()
	
	#arbitrarily remove values not ideal implementation
	evi[evi < 0] <- NA
	
	smooth <- apply(get.var.ncdf(nc, varname), 3, rbind)
	deriv <- cbind(rep(NA,nrow(smooth)), t(apply(smooth,1,diff)))
	
	uyear <- unique(info$dates$rainyear)
	SOS <- array(NA,dim=c(nrow(smooth),length(uyear)))
	EOS <- array(NA,dim=dim(SOS))
	SUM <- array(NA,dim=dim(SOS))
	AMP <- array(NA,dim=dim(SOS))
	PEAKDAY <- array(NA,dim=dim(SOS))
	for (i in seq(along=uyear)){
		idx <- info$dates$rainyear == uyear[i]
		dates <- info$dates$integer[idx]
		SOS[,i] <- dates[as.numeric(apply(deriv[,idx], 1, which.max))]
		EOS[,i] <- dates[as.numeric(apply(deriv[,idx], 1, which.min))]
		SUM[,i] <- apply(evi[,idx], 1, sum, na.rm=T)
		AMP[,i] <- apply(evi[,idx], 1, max, na.rm=T) - apply(evi[,idx], 1, min, na.rm=T)
		PEAKDAY[, i] <- dates[as.numeric(apply(evi[,idx], 1, which.max))]		
	}
	
	if (n_plot > 0){
		PlotPhen(n_plot=n_plot,x=info$dates$date, smoy=smooth,
				yearidx=as.numeric(info$dates$year_date),yearidx2=as.numeric(info$dates$rainyear_origin), SOS=SOS, EOS=EOS)
	}
	#browser()
	df <- data.frame(
			year    = as.vector(outer( rep(1,nrow(smooth)), uyear)), 
			SOSdate = as.vector(SOS), 
			EOSdate = as.vector(EOS),
			SOS     = as.numeric(format(as.Date(as.vector(SOS),origin='1970-01-01'),'%j')), 
			EOS     = as.numeric(format(as.Date(as.vector(EOS),origin='1970-01-01'),'%j')),
			PEAKDAY = as.numeric(format(as.Date(as.vector(PEAKDAY),origin='1970-01-01'),'%j')),
			SUM     = as.vector(SUM),
			AMP     = as.vector(AMP)
	)
	
	df$LOS <- df$EOS - df$SOS
	
	# create landcover, prcp, income variables

	subs <- get.var.ncdf(nc, info$landcover$name)[,,1]
	df$cat <- as.vector(outer(subs, rep(1,length(uyear))))
	df$cat <- factor(df$cat,levels=info$coverkey$CODE,labels=info$coverkey$NAME)
	
	#precip
	prcp <- tapply(get.var.ncdf(nc, 'PRCP')[1,1,],info$dates$rainyear,sum)
	df$PRCP = as.vector(outer(rep(1, nrow(smooth)), prcp))
	
	#income variable
	inco <- as.vector(get.var.ncdf(nc, 'FAVINC0'))
	inco[inco == -2147483647] <- NA
	inco[inco == 0] <- NA
	df$inco = as.vector(outer(inco, rep(1,length(uyear)))) / 1000
	incobreaks <- c(0,quantile(df$inco[df$cat == 'DevLow'],c(0.33,0.66),na.rm=T),10^6)
	df$inco <- cut(df$inco,breaks=incobreaks)
	
	# remove where LOS <= 0
	df <- subset(df, LOS > 0)
	df <- subset(df, EOS > SOS)
	df <- subset(df, EOSdate > SOSdate)
	
	#Get Means
	#=========
	bylist 			    <- list(year=df$year, cat=df$cat)
	phen 				<- aggregate(df$LOS, by=bylist, FUN=length)
	phen$catcode        <- as.integer(phen$cat)
	phen$annualprcp 	<- aggregate(df$PRCP, by=bylist, FUN=mean, na.rm=T)$x
	phen$los_mean 		<- aggregate(df$LOS, by=bylist, FUN=mean, na.rm=T)$x
	phen$los_sd 		<- aggregate(df$LOS, by=bylist, FUN=sd, na.rm=T)$x
	phen$sos_mean 		<- aggregate(df$SOSdate, by=bylist, FUN=mean, na.rm=T)$x
	phen$sos_sd 		<- aggregate(df$SOSdate, by=bylist, FUN=sd, na.rm=T)$x
	phen$eos_mean 		<- aggregate(df$EOSdate, by=bylist, FUN=mean, na.rm=T)$x
	phen$eos_sd 		<- aggregate(df$EOSdate, by=bylist, FUN=sd, na.rm=T)$x
	phen$peakday_mean 	<- aggregate(df$PEAKDAY, by=bylist, FUN=mean, na.rm=T)$x
	phen$peakday_sd 	<- aggregate(df$PEAKDAY, by=bylist, FUN=sd, na.rm=T)$x
	phen$sum_mean 	    <- aggregate(df$SUM, by=bylist, FUN=mean, na.rm=T)$x
	phen$sum_sd 	    <- aggregate(df$SUM, by=bylist, FUN=sd, na.rm=T)$x
	phen$amp_mean 	    <- aggregate(df$AMP, by=bylist, FUN=mean, na.rm=T)$x
	phen$amp_sd 	    <- aggregate(df$AMP, by=bylist, FUN=sd, na.rm=T)$x	

	#get growing season precip
	start <- findInterval(phen$sos_mean, info$dates$integer)
	end <- findInterval(phen$eos_mean, info$dates$integer)
	cumprcp <- cumsum(get.var.ncdf(nc, 'PRCP')[1,1,])
	phen$growingprcp <- cumprcp[end] - cumprcp[start]
	phen$seasonstartdate <- info$dates$date[start]
	phen$seasonenddate <- info$dates$date[end]
	
	#write 
	write.table(phen, file='PhenologySummary_by_cat.csv', sep=',', col.names=NA)
	
	# scatter plot of prcp vs mean los with sd error bars
	#browser()
	yearcol <- rainbow(length(unique(phen$year)))[as.integer(factor(phen$year))]
	pdf('phen_growseason_prcp_v_los.pdf')
	plot(phen$growingprcp, phen$los_mean,
			xlim=c(0, 1.1*max(phen$growingprcp)),
			pch=as.integer(phen$cat),
			col=yearcol,
			xlab='Growing Season Precip',
			ylab='Length of Season'
	)
	legend('bottomright', as.character(unique(phen$year)), fill=rainbow(length(unique(phen$year))))
	legend('topright', as.character(info$coverkey$NAME), pch=1:length(levels(phen$cat)))
	dev.off()
	
	# by cat and income
	bylist 			  <- list(year=df$year, cat=df$cat, inco=df$inco)
	phen 			  <- aggregate(df$LOS, by=bylist, FUN=length)
	phen$catcode      <- as.integer(phen$cat)
	phen$annualprcp   <- aggregate(df$PRCP, by=bylist, FUN=mean, na.rm=T)$x
	phen$los_mean     <- aggregate(df$LOS, by=bylist, FUN=mean, na.rm=T)$x
	phen$los_sd       <- aggregate(df$LOS, by=bylist, FUN=sd, na.rm=T)$x
	phen$sos_mean     <- aggregate(df$SOSdate, by=bylist, FUN=mean, na.rm=T)$x
	phen$sos_sd       <- aggregate(df$SOSdate, by=bylist, FUN=sd, na.rm=T)$x
	phen$eos_mean     <- aggregate(df$EOSdate, by=bylist, FUN=mean, na.rm=T)$x
	phen$eos_sd       <- aggregate(df$EOSdate, by=bylist, FUN=sd, na.rm=T)$x
	phen$peakday_mean <- aggregate(df$PEAKDAY, by=bylist, FUN=mean, na.rm=T)$x
	phen$peakday_sd   <- aggregate(df$PEAKDAY, by=bylist, FUN=sd, na.rm=T)$x
	phen$sum_mean 	  <- aggregate(df$SUM, by=bylist, FUN=mean, na.rm=T)$x
	phen$sum_sd 	  <- aggregate(df$SUM, by=bylist, FUN=sd, na.rm=T)$x
	phen$amp_mean 	  <- aggregate(df$AMP, by=bylist, FUN=mean, na.rm=T)$x
	phen$amp_sd 	  <- aggregate(df$AMP, by=bylist, FUN=sd, na.rm=T)$x
	
	start <- findInterval(phen$sos_mean, info$dates$integer)
	end <- findInterval(phen$eos_mean, info$dates$integer)
	cumprcp <- cumsum(get.var.ncdf(nc, 'PRCP')[1,1,])
	phen$growingprcp <- cumprcp[end] - cumprcp[start]
	
	write.table(phen, file='PhenologySummary_by_catinco.csv', sep=',', col.names=NA)
	
	#phenology by neighborhood age
	#bef80 <- get.var.ncdf(nc, 'BLTBEF80')[,,1]
	#bef80 <- get.var.ncdf(nc, 'BLTBEF80')[,,1]
	
}

#' apply kmeans clustering to the smoothed time series
PhenClusters <- function(info, nc, varname, nclust, samples=100, sampsize=1000){
  require(cluster)
  #browser()
	start <- Sys.time()
	cat('\nK-Medoids Clustering...')
	
	data <- apply(get.var.ncdf(nc, varname), 3, rbind)	
	clust <- rep(NA,nrow(data))
	keeprows <- apply(!is.na(data), 1, sum) > 0
	data <- data[keeprows,]
	
	# cluster fit
	fit <- clara(x=data, k=nclust, samples=samples,
			sampsize=min(sampsize, nrow(data)),
			stand=FALSE
	)
	#reorder / relabel clusters
	clust[keeprows] <- fit$clustering
	means <- aggregate(data, by=list(clust[keeprows]), FUN=mean, na.rm=T)[,-1]
	sdev <- aggregate(data, by=list(clust[keeprows]), FUN=sd, na.rm=T)[,-1]
	cluorder <- order(rowSums(means))
	medoids <- fit$medoids[cluorder,]
	means <- means[cluorder,]
	sdev <- sdev[cluorder,]
	clust <- cluorder[clust]
	n <- table(clust[keeprows])
	
	#plot
	pdf('phenCluster.pdf',height=1.5*(nclust+2),width=12)
	par(mfrow=c(nclust+2, 1), mai=c(0.5, 0.5, 0.1, 0.1))
	
	subs <- as.vector(get.var.ncdf(nc, info$landcover$name)[,,1])
	subs <- factor(subs, levels=info$coverkey$CODE, labels=info$coverkey$NAME)
	classtable <- as.matrix(prop.table(table(clust,subs),margin=2))
	classtable <- classtable[,!is.na(colnames(classtable))]
	idx <- as.logical(round(seq(0,1,length.out=ncol(classtable))))
	barplot(classtable[,idx], beside=T, col=rainbow(nclust))
	barplot(classtable[,!idx], beside=T, col=rainbow(nclust))
	
	xvals <- info$dates$date
	for (i in seq(nclust)){
		plot(x=xvals, y=means[i,], 
				col=rainbow(nclust)[i], 
				ylim=c(0,max(means,na.rm=T)),
				ylab='smoothed evi',
				xlab='date',
				main=paste('n =', n[i])
		)
		
		points(xvals, medoids[i,], cex=.5)
		points(xvals, means[i,] + sdev[i,], pch='+', cex=.5)
		points(xvals, means[i,] - sdev[i,], pch='+', cex=.5)
		abline(v=info$dates$year_date,col=grey(.5))
	
	}
	dev.off()
	
	im <- matrix(clust, nrow=nc$dim$x$len, ncol=nc$dim$y$len)
	im <- im[, seq(ncol(im), 1, -1)]
	
	pdf('phen_img.pdf',
			height=2+(ncol(im)*info$scale), 
			width=2+(nrow(im)*info$scale)
	)
	par(mfrow=c(1,1), mai=c(1,1,1,1),omi=c(0,0,0,0))
	image(z=im, 
			x=get.var.ncdf(nc,'xcoord'), 
			y=get.var.ncdf(nc,'ycoord'), 
			xlab=att.get.ncdf(nc,'xcoord','units')$value, 
			ylab=att.get.ncdf(nc,'ycoord','units')$value,
			col=rainbow(nclust), 
			breaks=seq(0,nclust),
			useRaster=T
	)
	
	dev.off()
	cat(format(round(Sys.time() - start,1)))
	
	#wavelet correlation? 
    #====================
	#prcp <- get.var.ncdf(nc, 'PRCP')[1,1,]
	#wlet3 <- t(apply(means,1,WaveletCorr,y=prcp,lag=3))
	#wlet6 <- t(apply(means,1,WaveletCorr,y=prcp,lag=6))
	
	#pdf('phen_prcp_corr.pdf',height=1.5*(nclust+1),width=12)
	#par(mfrow=c(nclust, 1), mai=c(0.5, 0.5, 0.1, 0.1))
	
	#for (i in seq(nclust)){
	#	plot(x=xvals, y=wlet6[i,], 
	#			col=rainbow(nclust)[i], 
	#			ylim=c(min(wlet,na.rm=T),max(wlet,na.rm=T)),
	#			ylab='spearman wavelet correlation',
	#			xlab='date',
	#			main=paste('n =', n[i])
	#	)
	#	
	#	#points(x=xvals,y=wlet6[i,],col=rainbow(nclust)[i],pch='+')
	#	abline(h=0)
	#	abline(v=info$dates$year_date,col=grey(.5))
	#	
	#}
	#dev.off()
}

WaveletCorr <- function(x,y,lag){
	cc <- array(NA,dim=length(x))
	for (i in seq(along=x)){
		if ((i-lag) >= 1){
			c <- cor.test(x[(i-lag):i],y[(i-lag):i], method='spearman')
			cc[i] <- c$estimate 
		}
	}
	return(cc)
}
#' Plot individual time series and indicate phenology metrics
PlotPhen <- function(x, smoy, yearidx, SOS, EOS, n_plot, main='', yearidx2=NULL){	
	cat('\n\tPlot Phenology Time Series...')
	#browser()
	pdf('phenplots.pdf')
	par(mfrow=c(4,1))
	cols <- c(rgb(0,1,0,alpha=.4), rgb(1,0,0,alpha=.4))
	maxy <- max(smoy,na.rm=T)
	rawy <- smoy
	for (i in sample(seq(nrow(smoy)),n_plot)){
		plot(NA, xlim=c(min(x), max(x)), ylim=c(0,maxy))
		
		col <- rep(rgb(0,1,0,alpha=.4),length(SOS[i,]))
		col[SOS[i,] > EOS[i,]] <- rgb(1, 0, 0, alpha=.4)
		rect(SOS[i,], 0, EOS[i,], maxy, col=col)
		abline(v=SOS[i,],col='green',lw=2)
		abline(v=EOS[i,],col='red',lw=2)
	
		abline(v=yearidx, col='black',lty=3)
		abline(v=yearidx2, col='blue',lty=1)
		points(x, rawy[i,], col='black', pch='+', cex=.5)
		lines(x, smoy[i,], col='blue')
		
	}
	dev.off()
	cat('ok')
}

#' Create a new loess smoothed variable from an existing variable. 
#' Loess applied to each slice in dim 3 [i,j,:]  
LOESS_Var <- function(info, nc, srcvar, newvar, units, missval, span=5, maxna=0.2){
	cat('\nloess smooth', srcvar,'...')
	start <- Sys.time()

	#create new derived variable if does not already exist
	if (!(newvar %in% names(nc$var))){ 
		cat('\n\tcreate var')	
		close.ncdf( nc )
		nc <- open.ncdf(info$netcdf, write=TRUE) # open for writing and create variable
		var <- var.def.ncdf(newvar, units, list(nc$dim[['x']], nc$dim[['y']], nc$dim[['t']]), missval)
		nc <- var.add.ncdf( nc, var )
	}
	close.ncdf(nc) # close the current file
	
	# get parent variable and apply smooth function
	nc <- open.ncdf(info$netcdf,write=TRUE)
	srcdata <- get.var.ncdf(nc, srcvar)
	srcdata[srcdata == missval] <- NA
	rowdata <- apply(srcdata,3,rbind)
	newdata <- t(apply(rowdata, 1, LOESS_smooth, x=info$dates$integer, span=span/ncol(srcdata) ))
	newdata <- array(newdata,dim=dim(srcdata))
	
	put.var.ncdf( nc, 'sevi', newdata)
	close.ncdf(nc)
	
	# reopen updated file for reading 
	nc <- open.ncdf(info$netcdf)
	cat(format(Sys.time() - start))
	return(nc)
}

#' gaussian local regression fit (LOESS)
LOESS_smooth <- function(y, x, span){
	nacount <- sum(is.na(y))
	if (nacount > 40){
		return(y*NA)
	}
	fit <- loess(y ~ x, 
			na.action=na.exclude,
			family='gaussian',
			degree=1,
			span=span
	)
	
	pred <- as.integer(predict(fit, newvals=x))
	return(pred)
}
