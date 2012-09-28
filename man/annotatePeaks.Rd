\name{annotatePeaks}
\alias{annotatePeaks}

\title{
Associates peaks to TSS
}

\description{
A convenience function to associate a genomic position to a TSS provided by a gtf file. 
}

\usage{
annotatePeaks(peaks, gtf, limit=10e3, remove_unmatched=T, unifyBy=F, unify_fun="sum_by_length", min_genelength=300)
}

\arguments{
  \item{peaks}{
A \linkS4class{GRanges} object. 
}
  \item{gtf}{
A \linkS4class{GRanges} object with a meta data column \sQuote{transcript_id} and \sQuote{exon_id} like e.g. from \code{gtf2gr}. 
}
  \item{limit}{
Maximal distance in positive and negative direction for a peak - TSS association in base pairs.  
}
  \item{remove_unmatched}{
If TRUE, only TSS associated peaks will be returned.
}
  \item{unifyBy}{
If a transcript has multiple isoforms, the peak will be associated to all of them. In order to reduce the mapping to one unique isoform, a \code{DensityContainer} can be provided. The choice of the returned isoform will be made based on unify_fun. If the character string \sQuote{first} is provided, the resulting column transcript_id will contain arbitrarily the first ID found.  
}
  \item{unify_fun}{
A function which will choose the isoform in case of non unique peak - TSS associations. Defaults to the isoform with the highest sum of the reads divided by the length in base pairs \code{function(x){sum(x)/length(x)}}.
}
  \item{min_genelength}{
Genes with a total sum of all exons smaller than this value will be excluded from the output.
}
}

\details{
Convenience function to annotate a \linkS4class{GRanges} object having one row per peak from e.g. \code{macs2gr}. The resulting peak - TSS associations can be customized by the restricting the distance and resolving multiple matches using unify_fun.
}

\value{
\linkS4class{GRanges} object with row names according to the peak names provided and an added or updated meta data column \sQuote{transcript_id} with the associated transcript IDs.
}
\author{Julius Muller \email{ju-mu@alumni.ethz.ch}}
\examples{

exgtf<-dir(system.file("extdata", package="TransView"),full=TRUE,patt="gtf.gz$")[2]
exls<-dir(system.file("extdata", package="TransView"),full=TRUE,patt="xls$")

GTF<-gtf2gr(exgtf)
peaks<-macs2gr(exls,psize=500)
apeaks<-annotatePeaks(peaks=peaks,gtf=GTF)

}


