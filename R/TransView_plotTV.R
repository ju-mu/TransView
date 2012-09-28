

#' Plots overview of peak locations
#' @param ... 
#' @param peaks 
#' @param gtf 
#' @param scale 
#' @param cluster 
#' @param control 
#' @param interpolate 
#' @param show_names 
#' @param label_size 
#' @param zero_alpha 
#' @param colr 
#' @param colr_df 
#' @param colour_spread 
#' @param key_limit 
#' @param set_zero 
#' @param rowv 
#' @param ex_windows 
#' @param gclust 
#' @param norm_readc 
#' @param no_key 
#' @param stranded_peak 
#' @param ck_size 
#' @param remove_lowex 
#' @param verbose 
#' @returnType 
#' @return 
#' @author Julius Muller
#' @export
plotTV<-function ( ..., regions, gtf=NA, scale="global", cluster="none", control = F, interpolate = 1, 
		show_names=T, label_size=1, zero_alpha=0.5, colr=c("white","blue", "red"), colr_df="redgreen",
		colour_spread=c(0.05,0.05), key_limit="auto", key_limit_rna="auto", set_zero="center", rowv=NA,
		ex_windows=5, gclust="peaks", norm_readc=T, no_key=F, stranded_peak=T, ck_size=c(2,1), remove_lowex=0, verbose=1 ) 
{
	stopifnot(is.numeric(set_zero) || set_zero=="center")
	argList<-list(...)

	
	ttlRNA<-c();ttl<-c()
	
	if(class(regions)[1]!="GRanges"){
		if(class(regions)[1]=="character"){
			trefs<-length(unique(regions))
			regions<-as.data.frame(gtf[which(values(gtf)$transcript_id %in% regions),])
			tpeaks<-length(unique(regions$transcript_id))
			if(tpeaks!=tpeaks){
				if(tpeaks==0){stop("No identifier in column 'transcript_id' of the gtf is matching the regions!")
				}else{warning(paste(tpeaks-length(unique(regions$transcript_id )),"transcript_id's not found in GTF"))}
			}
			rg<-split(regions,f=regions$transcript_id)
			rg<-lapply(rg,function(x){c(as.character(x$seqnames)[1],min(x$start),max(x$end),sum(x$width),as.character(x$strand)[1],as.character(x$transcript_id)[1])})
			regions<-as.data.frame(do.call("rbind",rg),stringsAsFactors=F)
			colnames(regions)<-c("seqnames","start","end","width","strand","transcript_id")
		}else if(class(regions)[1]!="logical" || !is.na(regions))stop("regions must be of class 'GRanges' or 'character'")
	}else{
		regions<-as.data.frame(regions)
		regions$seqnames<-as.character(regions$seqnames)
		regions$strand<-as.character(regions$strand)
		tpeaks<-nrow(regions)
	}
	if(class(gtf)[1]!="GRanges"){
		if(class(gtf)[1]!="logical" || !is.na(gtf))stop("gtf must be of class 'GRanges'")
	}
	
	
	if(tpeaks<2)stop("At least 2 rows have to be present in regions")
	
	if(is.numeric(rowv)){
		if(cluster!="none"){rowv<-NA
		}else cluster<-1
	}else if(!is.na(rowv))stop("rowv has to be a numeric vector")
	if(!is.character(key_limit) & (!is.numeric(key_limit) | length(key_limit)!=2))stop("key_limit must be a numeric vector of length 2")
	if(!is.character(key_limit_rna) & (!is.numeric(key_limit_rna) | length(key_limit_rna)!=2))stop("key_limit_rna must be a numeric vector of length 2")
	if(!(gclust %in% c("expression","peaks","both")))stop("Argument gclust must be either 'expression','peaks' or 'both'")
	if(!(scale %in% c("global","individual")))stop("Argument scale must be either 'global' or 'individual'")
	tcvg<-c();rcvg<-c();hmap<-NA
	for (arg in argList) {
		
		if (!.is.dc(arg)){
			if (length( dim(arg)) == 2 && is.numeric(arg)){
				if(!is.na(hmap))stop("Only one matrix allowed!")
				hmap<-arg
				if(nrow(hmap)!=tpeaks)stop("A matrix has to have the same amount of rows like regions")
				next
			}else stop("Data sets must be of any number of class 'DensityContainer' and maximally one matrix")
		} 
		if(spliced(arg)){
			if(class(gtf)[1]!="GRanges")stop("Expression data detected but no GTF found! Please re-run with gtf2gr")
			rcvg <- c(rcvg, ifelse(norm_readc,fmapmass(arg),1))
			if(!any(colnames(regions) %in% "transcript_id"))stop("'transcript_id' column is missing in regions. RNA-Seq can not be associated to regions.")
			regions$transcript_id<-as.character(regions$transcript_id)
			ttlRNA<-c(ttlRNA,ex_name(arg))
		}else{
			tcvg <- c(tcvg, ifelse(norm_readc,fmapmass(arg),1))
			ttl<-c(ttl,ex_name(arg))
		}
	}
	
	hmapc<-ifelse(length( dim(hmap)) == 2 && is.numeric(hmap),1,0)
	if (hmapc && length(rcvg)>0)stop("Expression data can be passed as one matrix or one or many DensityContainer but not both!")
	
	if(!is.logical(control) && length(control)!=(length(argList)-hmapc))stop("If control is provided, it must match the amount of experiments.")
	
	if(!is.null(tcvg))nvec <- tcvg/min(tcvg)
	if(!is.null(rcvg))rvec <- rcvg/min(rcvg)
	argc <- 0; argcRNA <- 0
	plotmat<-list();scalevec<-list();key_limits<-list()
	plotmatRNA<-list()
	if(hmapc)plotmatRNA<-hmap
	scalevecRNA<-list();key_limitsRNA<-list()
	
	if(verbose>0)message("Fetching densities...")
	
	for (arg in argList) {
		if(.is.dc(arg)){
			if (!is.logical(control)){
				if(!.is.dc(control[[argc+argcRNA+1]]))stop("Input must be of class 'DensityContainer'")
				ctrl<-control[[argc+argcRNA+1]]
			}else{ctrl=F}
			if(spliced(arg)){
				argcRNA <- argcRNA+1
				
				dsts <- sliceNT(arg, gtf=gtf, tnames=regions$transcript_id, control = ctrl,treads_norm = norm_readc)
				if(remove_lowex)glens<-unlist(lapply(dsts,length))
				
				dsts<-.gene2window(dsts,ex_windows,window_fun="approx")
				
				dsts<- do.call(rbind, dsts)/rvec[argcRNA]#matrix for fast plotting normalized by total reads
				plotmatRNA[[argcRNA]] <- dsts
				
			}else{
				argc <- argc + 1
				dsts <- sliceN(arg, regions, control = ctrl,treads_norm = norm_readc)
				if(stranded_peak && "strand" %in% colnames(regions)){#flip negative strand at the tss
					oord<-which(regions$strand == "-")
					dsts[oord]<-lapply(dsts[oord],rev)
				}
				usize <- unique(sapply(dsts, length))/interpolate
				if (length(usize) > 1)stop("The query regions must have equal length")
				if(interpolate>1)dsts<-.gene2window(dsts,usize,window_fun="approx")
				
				dsts<- do.call(rbind, dsts)/nvec[argc]
				plotmat[[argc]] <- dsts
				
				scalevec[[argc]]<-as.vector(plotmat[[argc]])#for keys and scale finding
				kh<-quantile( scalevec[[argc]][scalevec[[argc]]>mean(scalevec[[argc]])[1]],prob=c(colour_spread[1],1-colour_spread[1]))
				key_limits[[argc]]<-c(floor(min(kh)),ceiling(max(kh)))
			}
		}
	}
	if(argcRNA>0 && (cluster!="none"||remove_lowex))rnaj<-do.call("cbind",plotmatRNA)
	
	### ###
	
	if(verbose>0)message("Plotting...")
	
	#### Remove not expressed ####
	if(remove_lowex>0 && argcRNA>0){
		lowex<-which((rowSums(rnaj)/glens)<remove_lowex)
		if(length(lowex)>0){
			regions<-regions[-lowex,]
			message(sprintf("%d genes did not pass the expression threshold",length(lowex)))
			rnaj<-rnaj[-lowex,]
			for(x in 1:argcRNA)plotmatRNA[[x]] <- plotmatRNA[[x]][-lowex,]
			if(argc>0){
				for(x in 1:argc){
					plotmat[[x]] <- plotmat[[x]][-lowex,]
					kh<-quantile( scalevec[[x]][scalevec[[x]]>mean(scalevec[[x]])[1]],prob=c(colour_spread[1],1-colour_spread[1]))
					key_limits[[argc]]<-c(floor(min(kh)),ceiling(max(kh)))
				}
			}
		}
		tpeaks<-nrow(regions)
	}
	
	### ###
	
	
	
	#### CLUSTER ####
	
	if(cluster!="none"){
		if(argc && !hmapc && !argcRNA){cob<-do.call("cbind",plotmat)
		}else if(argcRNA && !argc && !hmapc){cob<-do.call("cbind",plotmatRNA)                   
		}else if(hmapc && !argcRNA && !argc){cob<-plotmatRNA
		}else if(gclust=="peaks"){cob<-do.call("cbind",plotmat)
		}else if(gclust=="expression" && !hmapc){cob<-do.call("cbind",plotmatRNA)
		}else if(gclust=="expression" && hmapc){ cob<-plotmatRNA                                         
		}else if(gclust=="both" && !hmapc){ cob<-cbind(do.call("cbind",lapply(plotmat,.row_z_score)),do.call("cbind",lapply(plotmatRNA,.row_z_score)))   
		}else if(gclust=="both"){ cob<-cbind(do.call("cbind",lapply(plotmat,.row_z_score)),.row_z_score(plotmatRNA))
		}
		
		cob<-.row_z_score(cob)#do kmeans clustering only on z scores
		if(is.numeric(cluster)){
			#do kmeans clustering only on z scores
			kclust<-kmeans(.row_z_score(cob),cluster)$cluster
		}else if(substr(cluster,1,3)=="hc_"){
			if(substr(cluster,4,5)=="sp"){dend<-as.dendrogram(hclust(as.dist(1-cor(t(cob), method="spearman"))))
			}else if(substr(cluster,4,5)=="pe"){dend<-as.dendrogram(hclust(as.dist(1-cor(t(cob), method="pearson"))))
			}else if(substr(cluster,4,5)=="rm"){dend<-as.dendrogram(hclust(dist(rowMeans(cob))))
			}else stop("Clustering not implemented:",cluster) 
		}else stop("Clustering not implemented:",cluster) 
	}
	
	### ###
	
	#### Re scale expression ####
	
	if(argcRNA>0){
		cob<-do.call("cbind",plotmatRNA)
		cob<-.row_z_score(cob)
		vcob<-as.vector(cob)#vector for keys and scale finding
		xa<-quantile( vcob[abs(vcob)>mean(vcob)],prob=c(colour_spread[2],1-colour_spread[2]))
		maa<-abs(xa[which(abs(xa)==max(abs(xa)))[1]])
		for(x in 1:argcRNA){
			key_limitsRNA[[x]]<-c(floor(-maa),ceiling(maa))
			plotmatRNA[[x]] <- cob[,(1+(x-1)*ex_windows):(ex_windows+((x-1)*ex_windows))]
			scalevecRNA[[x]]<-as.vector(plotmatRNA[[x]])#for keys and scale finding
		}
	} else if (hmapc){
		plotmatRNA<-.row_z_score(plotmatRNA)
		scalevecRNA<-as.vector(plotmatRNA)#vector for keys and scale finding
		xa<-quantile( scalevecRNA[abs(scalevecRNA)>mean(scalevecRNA)],prob=c(colour_spread[2],1-colour_spread[2]))
		maa<-abs(xa[which(abs(xa)==max(abs(xa)))[1]])
		key_limitsRNA<-c(floor(-maa),ceiling(maa))
	}
	### ###
	
	
	#### LAYOUT ####
	
	op <- par(no.readonly = TRUE)
	on.exit(par(op))
	
	lhei <- c(ck_size[1], 8)
	
	lwid<-rep(10,argc+argcRNA+hmapc)
	lmax<-(2*(argc+argcRNA+hmapc))
	uorder<-1:(argc+argcRNA+hmapc);lorder<-(argc+argcRNA+hmapc+1):lmax
	if(show_names){
		if(argc>0){
			lwid<-c(1,lwid)#reserve some space for peak names in first column
			uorder<-c(0,uorder)#dont plot upper panel
			lmax<-lmax+1
			lorder<-(argc+argcRNA+1):lmax
		}
		if(argcRNA>0){
			lwid<-c(lwid,1)#reserve some space for gene names in last column
			uorder<-c(uorder,0)#dont plot upper panel
			lmax<-lmax+1
			lorder<-(argc+argcRNA+1):lmax
		} else if(hmapc){
			lwid<-c(lwid,1)#reserve some space for gene names in last column
			uorder<-c(uorder,0)#dont plot upper panel
			lmax<-lmax+1
			lorder<-(argc+hmapc+1):lmax
		}
	}
	
	if(cluster!="none"){
		uorder<-c(0,uorder)
		lmax<-lmax+1
		lorder<-(argc+argcRNA+hmapc+1):lmax
		if(substr(cluster,1,3)=="hc_"){
			lwid<-c(3,lwid)
		} else {#names are plotted last so positions have to be flipped for kmeans
			lwid<-c(1,lwid)
			if(argc>1 && show_names){
				lorder[1:2]<-rev(lorder[1:2])
				lwid[1]<-1
			}
		}
	}  
	lmat <- rbind(uorder, lorder)
	below<-rep(0,length(uorder))
	below[which(uorder>0)]<-(max(lorder)+1):(max(lorder)+length(which(uorder>0)))
	lmat<-rbind(lmat,below)
	lhei<-c(lhei,0.1)
	
	if(no_key){
		lhei<-lhei[2:3]
		mu<-max(uorder)
		lorder[lorder>0]<-lorder[lorder>0]-mu
		below[below>0]<-below[below>0]-mu
		lmat<-rbind(lorder,below)
	}
	
	layout(lmat, widths = lwid, heights = lhei)
	
	
#layout.show(nf) 
	### ###
	
	
	
	#### PLOT KEYS AND SET SCALES ####
	
	if(scale=="global"){
		if(argc>0){
			gmin<-min(unlist(key_limits))
			gmax<-max(unlist(key_limits))
			for (x in 1:argc)key_limits[[x]]<-c(floor(gmin),ceiling(gmax))
		}
		if(argcRNA>0){
			gmin<-min(unlist(key_limitsRNA))
			gmax<-max(unlist(key_limitsRNA))
			for (x in 1:argcRNA)key_limitsRNA[[x]]<-c(floor(gmin),ceiling(gmax))
		}
	}
	
	shrink<-function(x,ext){x[x < ext[1]] <- ext[1];x[x > ext[2]] <- ext[2];x}
	rowMax<-function(x){apply(x,1,max)}
	if(!no_key)par(mar = c(1.5, 3/ck_size[2], 1.5,3/ck_size[2]), cex = 0.45)#c(bottom, left, top, right)
	col<-list();breaks<-list()
	
	for (argn in 1:argc) {
		if(argc==0)break
		rmax<-min(rowMax(plotmat[[argn]]))
		
		if(!is.character(key_limit)){key_limits[[argn]]<-key_limit
		}else if(key_limits[[argn]][1]>rmax && rmax>=1 && scale!="global")key_limits[[argn]][1]<-rmax-1
		
		plotmat[[argn]]<-shrink(plotmat[[argn]],key_limits[[argn]])
		col[[argn]] <- colorpanel(key_limits[[argn]][2]-key_limits[[argn]][1], colr[1],colr[2],colr[3])
		if(!no_key){
			scalevec[[argn]]<-shrink(scalevec[[argn]],key_limits[[argn]])
			breaks <- seq(key_limits[[argn]][1], key_limits[[argn]][2], length = length(col[[argn]])+1)
			image( matrix(breaks, ncol = 1), col = col[[argn]],breaks=breaks, xaxt = "n", yaxt = "n",ylim=c(0,1))
			
			hc <- hist(scalevec[[argn]], plot = F, breaks = breaks)$counts
			ktitle<-sprintf("Reads %s > %d","",key_limits[[argn]][1])
			kpeak<-2*max(hc[2:(length(hc)-1)])
			if(hc[1]>kpeak){#correct for very skewed distributions
				ktitle<-sprintf("Reads %s > %d","",key_limits[[argn]][1]+1)
				hc[1]<-0
			}
			if(hc[length(hc)]>kpeak){#correct for very skewed distributions
				ktitle<-sprintf("Reads %s > %d",paste(key_limits[[argn]][2]-1,"> x"),key_limits[[argn]][1]+1)
				hc[length(hc)]<-0
			}
			if(ttl[argn]!="NA")title(ttl[argn])
			mtext(side = 1, ktitle, line = 1.2,cex=label_size*0.5)
			hy <- c(hc, hc[length(hc)])
			lines(seq(0,1,length.out=length(breaks)), hy/max(hy) * 0.95, lwd = 2, type = "s", col = "cyan")
			lv <- round(seq(key_limits[[argn]][1], key_limits[[argn]][2], length = 5))
			axis(1, at = seq(0,1,length.out=length(lv)), labels = lv, line = -1, tick = 0,font=2,cex.axis=label_size*0.7)
			axis(2, at = pretty(hy)/max(hy) * 0.8, pretty(hy),cex.axis=label_size*0.5)
			mtext(side = 2, "Count", line = 2,cex=label_size*0.5)#font=1
		}
	}
	colres<-100
	colRNA <- colorpanel(colres, "blue","white", "red")
	breaksRNA<-list()
	for (argn in 1:argcRNA) {
		if(argcRNA==0)break
		
		rmax<-min(rowMax(plotmatRNA[[argn]]))
		
		if(!is.character(key_limit_rna)){key_limitsRNA[[argn]]<-key_limit_rna
		}else if(key_limitsRNA[[argn]][1]>rmax && rmax>=1 && scale!="global")key_limitsRNA[[argn]][1]<-rmax-1
		plotmatRNA[[argn]]<-shrink(plotmatRNA[[argn]],key_limitsRNA[[argn]])
		if(!no_key){
			scalevecRNA[[argn]]<-shrink(scalevecRNA[[argn]],key_limitsRNA[[argn]])
			breaks <- seq(key_limitsRNA[[argn]][1], key_limitsRNA[[argn]][2], length = colres+1)
			image( matrix(breaks, ncol = 1), col = colRNA,breaks=breaks, xaxt = "n", yaxt = "n",ylim=c(0,1))
			
			hc <- hist(scalevecRNA[[argn]], plot = F, breaks = breaks)$counts
			kpeak<-3*max(hc[2:(length(hc)-1)])
			if(ttlRNA[argn]!="NA")title(ttlRNA[argn])
			mtext(side = 1, "Z-Score", line = 1.2,cex=label_size*0.5)
			hy <- c(hc, hc[length(hc)])
			lines(seq(0,1,length.out=length(breaks)), hy/max(hy) * 0.95, lwd = 2, type = "s", col = "cyan")
			lv <- round(seq(key_limitsRNA[[argn]][1], key_limitsRNA[[argn]][2], length = 5))
			axis(1, at = seq(0,1,length.out=length(lv)), labels = lv, line = -1, tick = 0,font=2,cex.axis=label_size*0.7)
			axis(2, at = pretty(hy)/max(hy) * 0.8, pretty(hy),cex.axis=label_size*0.7)
			mtext(side = 2, "Count", line = 2,cex=label_size*0.5)#font=1
		}
	}
	
	if(hmapc){# heatmap
		if(!is.character(key_limit_rna))key_limitsRNA<-key_limit_rna
		plotmatRNA<-shrink(plotmatRNA,key_limitsRNA)
		colres<-100
		if(colr_df=="redgreen")colRNA=greenred(100)#looks better for microarrays?
		if(!no_key){
			scalevecRNA<-shrink(scalevecRNA,key_limitsRNA)
			breaks <- seq(key_limitsRNA[1], key_limitsRNA[2], length = colres+1)
			image( matrix(breaks, ncol = 1), col = colRNA,breaks=breaks, xaxt = "n", yaxt = "n",ylim=c(0,1))
			
			hc <- hist(scalevecRNA, plot = F, breaks = breaks)$counts
			kpeak<-3*max(hc[2:(length(hc)-1)])
			
			mtext(side = 1, "Z-Score", line = 1.2,cex=label_size*0.5)
			hy <- c(hc, hc[length(hc)])
			lines(seq(0,1,length.out=length(breaks)), hy/max(hy) * 0.95, lwd = 2, type = "s", col = "cyan")
			lv <- round(seq(key_limitsRNA[1], key_limitsRNA[2], length = 5))
			axis(1, at = seq(0,1,length.out=length(lv)), labels = lv, line = -1, tick = 0,font=2,cex.axis=label_size*0.7)
			axis(2, at = pretty(hy)/max(hy) * 0.8, pretty(hy),cex.axis=label_size*0.7)
			mtext(side = 2, "Count", line = 2,cex=label_size*0.5)#font=1
		}
	}
	
	### ###
	
	
	#### PLOT CLUSTER ####
	par(mar=c(0,0,1,0))#c(bottom, left, top, right)
	if(is.numeric(rowv)){
		if(length(rowv)!=tpeaks)stop("rowv has to have the same row number as the data")
		kclust<-rowv
		cluster<-length(unique(rowv))
	}
	if(cluster!="none"){
		if(is.numeric(cluster)){
			morder<-order(kclust)
			if(is.numeric(rowv) & "Order" %in% colnames(regions))morder<-regions$Order
			kcout<-kclust
			cpos<-1
			kcol<-rainbow(cluster)
			for(k in 1:cluster)kclust[which(kclust==k)]<-kcol[k]
			image(rbind(seq(1,0,-1/(tpeaks-1))),col = kclust[morder], axes = F)
		}else if(substr(cluster,1,3)=="hc_"){
			plot(dend, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none",edgePar = list(lwd=1))
			morder<-order.dendrogram(dend)
		}
	}
	### ###
	
	
	#### PEAK NAMES ####
	if(argc>0 && show_names){
		par(mar=c(0,0,1,0))#c(bottom, left, top, right)
		plot(x=c(0.5,0.5),y=c(0,tpeaks),type="n",axes=F,ylim=c(-0.5,tpeaks-0.5),yaxs="i")
		text(.5,seq((tpeaks-1),0), labels =  rownames(plotmat[[argc]]), xpd = F,cex = label_size)
	}
	### ###
	
	#### PLOT MAIN FIGURE ####
	par(mar = c(0, 1, 1, 1))#c(bottom, left, top, right)
	rotate = function(mat) t(mat[nrow(mat):1,,drop=FALSE])
	for (argn in 1:argc) {
		if(argc==0)break
		if(cluster!="none")plotmat[[argn]]<-plotmat[[argn]][morder,]
		image(rotate(plotmat[[argn]]), col = col[[argn]], useRaster = T, axes = F)
		lplot<-dim(plotmat[[argn]])[2]*interpolate
		if(set_zero=="center")set_zero<-lplot/2
		#set_zero<-set_zero/interpolate
		#axis(1, at = seq(0,1,length.out=5), labels = round(seq(-set_zero,lplot-set_zero,lplot/4)), cex.axis=label_size, line = -1, tick = 0,font=2)
		lines(c(set_zero/lplot,set_zero/lplot),c(0,1),col=rgb(0,0,0,alpha=zero_alpha),lwd=3,lty=1)
	}
	for (argn in 1:argcRNA) {
		if(argcRNA==0)break
		if(cluster!="none")plotmatRNA[[argn]]<-plotmatRNA[[argn]][morder,]
		image(rotate(plotmatRNA[[argn]]), col = colRNA, useRaster = T, axes = F)
	}
	if(hmapc){
		if(cluster!="none")plotmatRNA<-plotmatRNA[morder,]
		image(rotate(plotmatRNA), col = colRNA, useRaster = T, axes = F)
		#axis(1, at= seq(0,1,length.out=length(colnames(hmap))), labels = colnames(hmap), cex.axis=label_size, line = -1, tick = 0,font=1)
	}
	### ###
	
	#### GENE NAMES ####
	if(show_names){
		if(argcRNA>0){
			par(mar=c(0,0,1,1))#c(bottom, left, top, right)
			plot(x=c(0.5,0.5),y=c(0,tpeaks),type="n",axes=F,ylim=c(-0.5,tpeaks-0.5),yaxs="i")
			text(y=seq((tpeaks-1),0), .5, labels =  gsub(".y","",rownames(plotmatRNA[[argcRNA]])), xpd = F,cex = label_size)
		}else if(hmapc){
			par(mar=c(0,0,1,1))#c(bottom, left, top, right)
			plot(x=c(0.5,0.5),y=c(0,tpeaks),type="n",axes=F,ylim=c(-0.5,tpeaks-0.5),yaxs="i")
			text(y=seq((tpeaks-1),0), .5, labels = rownames(plotmatRNA), xpd = F,cex = label_size)
		}
	}
	### ###
	
	#### PLOT AXIS LABELS ####
	par(mar = c(0, 0, 0, 0))#c(bottom, left, top, right)
	
	for (argn in 1:argc) {
		if(argc==0)break
		plot.new()
		text(seq(0,1,length.out=5),0.5, labels = round(seq(-set_zero,lplot-set_zero,lplot/4)), cex=label_size, font=2)
	}
	for (argn in 1:argcRNA) {
		if(argcRNA==0)break
		plot.new()
		text(c(0,1),0.5, labels = c("5'","3'"), cex=label_size, font=2)
	}
	if(hmapc){
		plot.new()
		corec<-(1/ncol(hmap))/2.2
		text(seq(0+corec,1-corec,length.out=ncol(hmap)),.5, labels = colnames(hmap), cex=label_size, font=2)
	}
	### ###
	
	if(is.numeric(cluster) & !is.numeric(rowv)){invisible(cbind(regions,"Cluster"=kcout,"Order"=morder))
	}else if(cluster!="none" & !is.numeric(rowv)){invisible(regions[morder,])
	} else {invisible(regions)}
}


