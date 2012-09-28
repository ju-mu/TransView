
library("pasillaBamSubset")

exbam<-dir(system.file("extdata", package="TransView"),full=T,patt="bam$")
filename<-system.file("extdata", "test.sam", package="TransView")

test_parseReads <- function(){
	test_data<-parseReads(filename, spliced=TRUE,verbose=0,hwindow=20,compression=20)
	test_stat<-tvStats(test_data)
	checkEqualsNumeric(c(1, NA, NA), histogram(test_data)[1:3])
	checkEqualsNumeric(1.34783, round(test_stat$gcoverage,digits=5))
}


exls<-dir(system.file("extdata", package="TransView"),full=T,patt="xls$")
peaks<-macs2gr(exls,psize=500)
peak<-as.data.frame(peaks[18,])
peak$seqnames<-as.character(peak$seqnames)

test_parseVariations <- function(){
	
	exden.ctrl<-parseReads(exbam[1],verbose=0,hwindow=20,compression=20)
	exden.chip<-parseReads(exbam[2],verbose=0,hwindow=20,compression=20)
	test_stat<-tvStats(exden.chip)
	checkIdentical(c("1"=49, "2"=40, "3"=21, "4"=28, "5"=23, "6"=21, "7"=20, "8"=14, "9"= 7,"10"=4), histogram(exden.chip)[1:10])
	x1<-slice1(exden.chip,peak$seqnames,peak$start,peak$end)[300:310]
	x2<-slice1(exden.ctrl,peak$seqnames,peak$start,peak$end)[300:310]
	x3<-slice1(exden.chip,peak$seqnames,peak$start,peak$end,control=exden.ctrl,treads_norm=F)[300:310]
	checkIdentical(0.00023, round(test_stat$gcoverage,digits=5))
	checkEqualsNumeric(c(77 ,74 ,74 ,73 ,71 ,71 ,69 ,72 ,76 ,74 ,75),x1)
	checkEqualsNumeric(c( 7, 7, 7 ,7, 7, 7, 5, 5, 5, 5, 5),x2)
	checkEqualsNumeric(c(70, 67, 67, 66, 64, 64, 64, 67, 71, 69, 70),x3)
}


test_parseVariations2 <- function(){
	
	exden.ctrl<-parseReads(exbam[1],verbose=0,hwindow=20,compression=20)
	exden.chip<-parseReads(exbam[2],min_quality=10,extendreads=10,read_stranded=-1,max_dups=4,verbose=0,hwindow=20,compression=20)
	test_stat<-tvStats(exden.chip)
	checkEqualsNumeric(c(  72, 40, 37, 28, 18 ,14, 13,  7 , 5,10), histogram(exden.chip)[1:10])
	x1<-slice1(exden.chip,peak$seqnames,peak$start,peak$end)[300:310]
	x2<-slice1(exden.ctrl,peak$seqnames,peak$start,peak$end)[300:310]
	x3<-slice1(exden.chip,peak$seqnames,peak$start,peak$end,control=exden.ctrl,treads_norm=T,input_method="/")[300:310]
	checkIdentical(0.0001037, round(test_stat$gcoverage,digits=7))
	checkEqualsNumeric(c(81,78,77,77,73,73,70,71,71,70,73),x1)
	checkEqualsNumeric(c( 7, 7, 7 ,7, 7, 7, 5, 5, 5, 5, 5),x2)
	checkIdentical(c(3.35755, 3.30378, 3.28540, 3.28540, 3.20945, 3.20945, 3.56478, 3.58496, 3.58496, 3.56478, 3.62449),round(x3,digits=5))
}

exgtf<-dir(system.file("extdata", package="TransView"),full=T,patt="gtf.gz$")
GTF.mm9<-gtf2gr(exgtf[2])
GTF.dm3<-gtf2gr(exgtf[1])
fn.pas_paired<-untreated1_chr4()

test_parse_T <- function(){
	exden.exprs<-parseReads(fn.pas_paired,spliced=T,min_quality=30,read_stranded=1,max_dups=5,verbose=0)
	x1<-slice1T(exden.exprs,values(GTF.dm3)$transcript_id[50],gtf=GTF.dm3,concatenate=T,stranded=T)[2735:2745]
	checkEqualsNumeric(c(33, 33, 33, 34, 34, 34, 36, 36, 36, 36, 35),x1)
	x2<-sliceNT(exden.exprs,unique(values(GTF.dm3)$transcript_id[101:150]),gtf=GTF.dm3,concatenate=F,stranded=F)$"NM_001014696.2"  [31:40]
	checkEqualsNumeric(c( 8,  8, 11, 12, 12, 13, 13, 13, 12, 13),x2)
}


test_annotate<- function(){
	apeaks<-annotatePeaks(peaks=peaks,gtf=GTF.mm9)
	checkIdentical(values(apeaks)$transcript_id[5],"NM_011031")
}
