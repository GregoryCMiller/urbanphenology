#Copyright (C) 2013 Greg Miller <gmill002@gmail.com>
#' Netcdf utilities 
library(yaml)
library(ncdf)

missval <- c(-32767,-2147483647)

#' Create date information data frame
#'   date       : date object (1989-10-31) 
#'   integer    : days since 1970-01-01 (5123)
#'   doy        : day of year (156)
#'   year       : year only (1989)
#'   year_date  : jan1 of year (1989-01-01)
#'   idx        : ordinal index (17)
#' 
#'  Rainyear goes from MM-DD in prev year to same date in current year
#'   rainyear_start: rainyear start of that year
#' 
DateInfo <- function(nc, rainyearmonth=7, rainyearday=1){
	dates <- as.Date(as.character(get.var.ncdf(nc,'tcoord')),format='%Y%j')
	df <- data.frame(date = dates)
	
	df$integer = as.numeric(dates)
	df$doy = as.numeric(format(dates,'%j'))
	df$year = as.numeric(format(dates,'%Y'))
	df$year_date <- as.Date(ISOdate(df$year,1,1))
	df$full_idx = seq(length(dates))
	df$year_idx <- floor(df$doy / 16)
	
	uyear <- c(min(df$year - 1), unique(df$year), max(df$year)+1)
	rainyears <- as.Date(ISOdate(uyear, rainyearmonth, rainyearday))
	df$rainyear_origin <- as.Date(cut(df$date, rainyears))
	df$rainyear <- as.numeric(format(df$rainyear_origin,'%Y'))
	df$rainyear_doy <- as.numeric(df$date - df$rainyear_origin)
	
	return(df)
}

#' load a var, replace NAs, mask subset
GetNCVar <- function(nc,vname,sub=NULL,start=NA,count=NA){
  x <- get.var.ncdf(nc,vname,start,count)
  x[x %in% missval] <- NA
  if (!is.null(sub)){
    x[!sub] <- NA
  }
  return(x)
}

#' Landcover Map
CoverMap <- function(info,nc,x0,name){
  #browser()
  ck <- info$coverkey[order(info$coverkey$CODE), ]
  
  x <- matrix(x0, nrow=dim(x0)[2], ncol=dim(x0)[1])
  x <- x[,seq(ncol(x),1,-1)]
  
  pdf(paste(name,'map.pdf',sep=''),
		  height=2+(ncol(x)*info$scale), 
		  width=2+(nrow(x)*info$scale)
  )
  par(mai=c(1,1,1,1),omi=c(0,0,0,0))
  
  image(x=get.var.ncdf(nc,'xcoord'),
		y=get.var.ncdf(nc,'ycoord'),
		z=x,
        main=paste(info$name,name),
        col=as.character(ck$COL), 
        breaks=c(-1,ck$CODE),
        xlab=att.get.ncdf(nc,'xcoord','units')$value, 
		ylab=att.get.ncdf(nc,'ycoord','units')$value,
		useRaster=T
  )
  dev.off()
  
  pdf(paste(name,'map_legend.pdf',sep=''))
  plot(NA,xlim=c(0,1),ylim=c(0,1),axes=FALSE)
  legend('bottomleft', as.character(ck$NAME), fill=as.character(ck$COL), bg=NA,ncol=2)
  dev.off()
  }

#' Calculate correlation coefs
CorrCoefs <- function(info,nc,sub,maxn=2000){
  cat('\nCor Coefs')
  #print(info$excludevars)
  vars <- names(nc$var)[!(names(nc$var) %in% info$excludevars)]
  #print(vars)
  ccdim <- c(length(vars),length(info$coverkey$NAME))
  cc <- array(NA,dim=ccdim)
  n <- array(NA,dim=ccdim)
  p <- array(NA,dim=ccdim)
  
  rownames(cc) <- vars
  rownames(n) <- vars
  rownames(p) <- vars
  
  colnames(cc) <-info$coverkey$NAME
  colnames(n) <- info$coverkey$NAME
  colnames(p) <- info$coverkey$NAME
  
  x <- GetNCVar(nc,'evi')
  
  for (r in seq(along=vars)) {
    vname <- vars[r]
	  y <- GetNCVar(nc, vname, NULL)
	
    for (d in seq(length(info$coverkey$NAME))){
      dval <- info$landcover$uvals[d]
      y0 <- y[sub == dval]
      x0 <- x[sub == dval]
	  idx <- !is.na(y0) & !is.na(x0)
	  x0 <- x0[idx]
	  y0 <- y0[idx]
	  if (length(x0) > maxn){
	  	sampidx <- sample(1:length(x0), maxn, replace=F)
		x0 <- x0[sampidx]
	  	y0 <- y0[sampidx]
		}
		
      grpname <- as.character(info$coverkey$NAME[which(dval==info$coverkey$CODE)])
      n[r,d] <- length(x0)
	  if ( n[r,d] > 10){
		c <- cor.test(x0, y0, method='spearman')
		cc[r,d] <- c$estimate 
		p[r,d] <- format(round(c$p.value,8),scientific=FALSE)
	  }
    }
  }
  info$corcoef <- cc
  write.table(cc,'corrcoef.csv',sep=',',col.names=NA)
  write.table(p,'corr_p.csv',sep=',',col.names=NA)
  write.table(n,'corr_n.csv',sep=',',col.names=NA)
  write.table(info$varsummary,'var_summary.csv',sep=',',col.names=NA)
  return(info)
}

#' Point wise correlation maps
PointWiseCor <- function(info,nc,xname,yname){
  #browser()
  mapname <- paste('pointcorr',xname,yname,sep='_')
  cat('\n', mapname,sep='')
  x <- GetNCVar(nc,xname)
  y <- GetNCVar(nc,yname)
  
  cc <- array(NA,dim=c(dim(x)[1],dim(x)[2]))
  for (i in seq(dim(x)[1])){
    for (j in seq(dim(x)[2])){
      x0 <- x[i,j,]
      y0 <- y[i,j,]
	  if (sum(!is.na(y0) & !is.na(x0)) > 4){  
        cc[i,j] <- cor(x0, y0, use="complete.obs", method='spearman')
      }
    }
  }
  
  cc <- matrix(cc,nrow=dim(cc)[2],ncol=dim(cc)[1])
  cc <- cc[,seq(ncol(cc),1,-1)]
  write.table(cc, paste(mapname, '.csv', sep=''), col.names=NA)
  
  pdf(paste(mapname,'map.pdf',sep='_'),
		  height=2+(nrow(x)*info$scale), 
		  width=2+(ncol(x)*info$scale)
  
  )
  par(mai=c(1,1,1,1),omi=c(0,0,0,0))
  breaks <- seq(-1, 1, .05)
  col <- rainbow(length(breaks)-1)
  image(z=cc,
		  x=get.var.ncdf(nc,'xcoord'),
		  y=get.var.ncdf(nc,'ycoord'), 
		  main=paste('Pointwise Spearman R (',xname,':',yname,')'),
		  xlab=att.get.ncdf(nc,'xcoord','units')$value, 
		  ylab=att.get.ncdf(nc,'ycoord','units')$value,
		  col=col,
		  breaks=breaks,
		  useRaster=T
  )
  
  dev.off()
  
  pdf(paste(mapname,'map_legend.pdf',sep='_'))
  plot(NA,xlim=c(0,1),ylim=c(0,1),axes=FALSE)
  legend('bottomleft', as.character(cut(breaks,breaks)),
		  fill=col, ncol=2, bg=NULL)
  dev.off()
}

#' 
GetVarsAsDataFrame <- function(vnames,nc){
  d = vector('list',length=length(vnames))
  for (vname in vnames){
    d[[which(vname == vnames)]] <- as.vector(get.var.ncdf(nc, vname))  
  }
  names(d) <- vnames
  D <- as.data.frame(d)
  return(D)
}

#' 
RunCorMx <- function(nc,info,subs) {
	outdir <- file.path(info$outpath,'CorrMx')
	dir.create(outdir,showWarnings=FALSE)
	setwd(outdir)
	
	retval <- GenericExpandDoCall(
		FUN='OneCorrelationMatrix',
		constants=list(
				nc=nc,
				info=info,
				sub=subs
		),
		varlist=list(
				subval = c(20,40,50,70,81,82), #info$landcover$uvals	
				xname	= 'evi',
				yname = c('PRCP','PRCPA3','PRCPA6'),
				downsample = c(10,12),#,3,2,1),
				maxn 	= 10**5
		)
	)
	setwd(info$outpath)
}

#' copy to target, clear before if necessary               
Copy.dir <- function(old,new) {
	unlink(new,recursive=T)
	dir.create(new)
	file.copy(file.path(old,dir(old)), file.path(new,dir(old)))
}

#' generate a 2d meshgrid
meshgrid <- function(a,b) {
	list(
			x=outer(b*0,a,FUN="+"),
			y=outer(b,a*0,FUN="+")
	)
} 

#' Execute a function for all combinations of parameters in varlist
#' varlist should be a list of named arguments possibly vectors
#' constants is a list of named args that are the same for each iteration
GenericExpandDoCall <- function(FUN, constants, varlist){
  varlist$stringsAsFactors <- FALSE
  varframe <- do.call('expand.grid', varlist)
  for (i in seq(nrow(varframe))){
    do.call(FUN, c(constants, as.list(varframe[i,])))
  }
}

#' Fisher's z
#' Convert Pearson's r's to the normally distributed variable z
fisherRtoZ <- function(r){
  0.5 * (log(1+r) - log(1-r)) 
}

#' copy everything including or except the *.nc file
CopyWithoutData <- function(studydir){
	outdir <- paste(studydir,'_nodata',sep='')
	unlink(outdir,recursive=T)
	dir.create(outdir)
	files <- dir(studydir,'(csv|yaml|jpg|pdf)$')
	newfiles <- file.copy(file.path(studydir,files),file.path(outdir,files))  
}
