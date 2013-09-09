homedir <- 'I:\\urbphen\\Gregs_Scripts_NLCDpatch'
source(file.path(homedir,'Util.R'))
source(file.path(homedir,'phenology.R'))
source(file.path(homedir,'CorrelationMatrix.R'))

# Manual control of analysis functions
COVERMAP  <- !TRUE
PHENOLOGY <- !TRUE
PRISM     <- TRUE
CORCOEFS  <- !TRUE
CLUSTER   <- FALSE
POINTCOR  <- FALSE

#' Main urban phenology function
MainUrbPhen <- function(studydir){
	cat('\nMain:',studydir)
	setwd(studydir)
  
  # Load saved parameters
  info <- yaml.load_file(file.path(studydir, 'info.yaml'))
  
  # open netcdf data file
  nc <- open.ncdf(info$netcdf)
  
  # dates
	info$dates <- DateInfo(nc,info$rainyearmonth, info$rainyearday)
  
  # Landcover  
	info$coverkey <- read.table(info$nlcdcol,header=T)
	subs <- GetNCVar(nc, info$landcover$name)
	info$landcover$uvals <- sort(unique(as.vector(subs)))
	info$landcover$uvals <- info$landcover$uvals[info$landcover$uvals != -1]
	
	#phenology 
	if (PHENOLOGY) {
		if ( (!('sevi' %in% names(nc$var))) | FALSE ) {
			nc <- LOESS_Var(info, nc, srcvar='evi', newvar='sevi', units='evi X 10000', missval=-32767, span=5)
		}	
    
    PhenMetrics(info, nc, varname='sevi')
		
    if (CLUSTER){
			PhenClusters(info, nc, varname='sevi', nclust=8)
		}
	}
	
  #Land Cover Map
	if (COVERMAP) {
		CoverMap(info, nc, subs[,,1], info$landcover$name)
	}
	
	#pointwise Correlation
	if (POINTCOR) {
		PointWiseCor(info, nc, 'evi','PRCP')
		PointWiseCor(info, nc, 'evi','PRCPA3')
		PointWiseCor(info, nc, 'evi','PRCPA6')
	}
	
	#correlation coefs 
	if (CORCOEFS){
		info <- CorrCoefs(info,nc,subs)
	}
	
	#Wavelet Correlation mx
	#RunCorMx(nc,info,subs)
	
  # Prism annual data
  if (PRISM){
	  PRISM_PPT(info, nc)
  }
  
  #write(as.yaml(info),info$yaml)
	CopyWithoutData(studydir)
	cat('\n\nR UrbanPhen Completed')
	return(info)
}

if (length(commandArgs(TRUE))==1){
  MainUrbPhen(commandArgs(TRUE)[1])
} else {
  MainUrbPhen(file.path("I:","urbphen","netcdf","abq"))
  MainUrbPhen(file.path("I:","urbphen","netcdf","ie"))
  MainUrbPhen(file.path("I:","urbphen","netcdf","lax"))
  MainUrbPhen(file.path("I:","urbphen","netcdf","phx"))
  MainUrbPhen(file.path("I:","urbphen","netcdf","las"))
  MainUrbPhen(file.path("I:","urbphen","netcdf","lbb"))
  MainUrbPhen(file.path("I:","urbphen","netcdf","aus"))
}