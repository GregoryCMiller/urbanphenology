library(yaml)
library(abind)
library(reshape2)

CombineTable <- function(info,value='coef'){	
	for (i in seq(along=info)){
		ccfile <- file.path(info[[i]]['outpath'],'corrcoef.csv')
		if (file.exists(ccfile)){
			df <- read.table(ccfile, sep=',', header=T, row.names=1)
			df$area <- info[[i]]$name
			df$var <- rownames(df)
			ldf <- melt(df,id.vars=c('area','var'), value.name=value,variable.name='type')
			if (i == 1){
				ccdf <- ldf
			} else {
				ccdf <- rbind(ccdf, ldf)
			}
		}
	}
	wide <- dcast(ccdf, type+var ~ area, value.var="coef")
	return(wide)
}
OrganizeLagVars <- function(info,wide){
	
	lagtab <- list()
	types <- unique(wide$type)
	for (l in seq(along=info[[1]]$lag)) {
		for (t in seq(along=types)){
			type <- types[t]
			swide <- wide[wide$type == type,]
			lag <- info[[1]]$lag[[l]]
			#sdf <- swide[swide$var %in% lag$names,]
			sdf <- swide[grep(paste(lag$var,lag$type,'[0-9]',sep=''),swide$var),]
			sdf <- rbind(swide[swide$var==lag$var,],sdf)
			sdf[1,'var'] <- paste(sdf[1,'var'],lag$type,'0',sep='') 
			sdf$times <- as.numeric(gsub('[A-Z]','',sdf$var))
			sdf <- sdf[order(sdf$times),]
			lag$tab <- sdf
			
			lagtab[[length(lagtab)+1]] <- lag
			
		}
	}
	
	write('','reorg.csv')
	for (i in seq(along=lagtab)){
	
		write.table(lagtab[[i]]$tab,'reorg.csv',append=TRUE,sep=',',col.names=NA)
		write('','reorg.csv',append=TRUE)
		
	}
	
	return(lagtab)
}
PlotLagVars <- function(lagtab,df){
	pdf('lagplots.pdf')
	#browser()
	for (i in seq(along=lagtab)){
		lag <- lagtab[[i]]
		tab <- lag$tab
		stab <- tab[,!names(tab) %in% c('times','var','type')]
		plot(NA,xlim=c(0,max(tab$times)),
				ylim=c(-1,1), 
				main=paste(tab$type[1],lag$var,sep='-'),
				xlab=lag$xlab,
				ylab='Spearman R')
		
		abline(h=0)
		for (a in seq(along=names(stab))){
			points(tab$times,stab[,a],col=df$col[a])
			lines(tab$times,stab[,a],col=df$col[a])	
		}	
		
		legend('bottomright',df$name,fill=df$col)
	}
	dev.off()
}
GetLagInfo <- function(info){	
	for (i in seq(along=info)){
		for (j in seq(along=info[[i]]$lag)){
			info[[i]]$lag[[j]]$type <- 'L'	
			info[[i]]$lag[[j]]$xlab <- 'Lag Time (x16 days)'
		}
		for (j in seq(along=info[[i]]$accumulate)){
			info[[i]]$accumulate[[j]]$type <- 'A'
			info[[i]]$lag <- c(info[[i]]$accumulate[j], info[[i]]$lag)
			info[[i]]$lag[[j]]$xlab <- 'Accumulate Time (x16 days)'
		}		
		for (j in seq(along=info[[i]]$lag)){
			ij <- info[[i]]$lag[[j]]
			info[[i]]$lag[[j]]$times <- c(0,ij$times)
			info[[i]]$lag[[j]]$names <- c(ij$var, paste(ij$var, 'L', ij$times, sep=''))
		}
	}
	return(info)
}

GetInfo <- function(folders){
	info <- lapply(file.path(folders, 'info.yaml'), yaml.load_file)
	return(info)
}

GetStudyDF <- function(info, folders){
	df <- data.frame(folder=folders)
	params <- c('abbrev','name','climate_type','climate_rank')
	for (n in params) {
		df[n] <- rep(NA,length(info))
		for (i in seq(along=info)){
			df[i,n] <- info[[i]][n]		
		}
	}
	df <- df[order(df$climate_rank),]
	return(df)	
}
  
topcompare <- function(folders, outdir){
	dir.create(outdir,show=FALSE)
	setwd(outdir)
	info <- GetInfo(folders)
	df <- GetStudyDF(info, folders)	
	df$col <- rainbow(nrow(df))
	
	wide <- CombineTable(info)
	info <- GetLagInfo(info)
	lagtab <- OrganizeLagVars(info,wide)
	PlotLagVars(lagtab,df)
	return(info)
}

abbrev <- c("phx", "lbb", "las","abq","aus","lax",'ie') #,, 
#abbrev <- paste(abbrev,'_nodata',sep='')
folders <- file.path("I:\\urbphen\\netcdf",abbrev)
outdir <- "I:\\urbphen\\netcdf\\compare"

info <- topcompare(folders,outdir)
