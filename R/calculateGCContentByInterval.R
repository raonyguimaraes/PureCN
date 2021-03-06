#' Calculates GC content by interval
#' 
#' Uses \code{scanFa} from the Rsamtools package to retrieve GC content of
#' intervals in a reference FASTA file.
#' 
#' 
#' @param interval.file File specifying the intervals. Interval is expected in
#' first column in format CHR:START-END.  Instead of a file, a \code{GRanges}
#' object can be provided. This allows the use of BED files for example. Note
#' that GATK interval files are 1-based (first position of the genome is 1).
#' Other formats like BED files are often 0-based. The \code{import} function
#' will automatically convert to 1-based \code{GRanges}.
#' @param reference.file Reference FASTA file.
#' @param output.file Optionally, write GC content file.
#' @param \dots Additional parameters passed to the \code{read.delim} function
#' that reads the \code{interval.file}.
#' @return Returns GC content by interval.
#' @author Markus Riester
#' @examples
#' 
#' reference.file <- system.file("extdata", "ex2_reference.fa", 
#'     package="PureCN", mustWork = TRUE)
#' interval.file <- system.file("extdata", "ex2_intervals.txt", 
#'     package="PureCN", mustWork = TRUE)
#' bed.file <- system.file("extdata", "ex2_intervals.bed", 
#'     package="PureCN", mustWork = TRUE)
#' calculateGCContentByInterval(interval.file, reference.file, 
#'     output.file="gc_file.txt")
#' 
#' intervals <- import(bed.file)
#' calculateGCContentByInterval(intervals, reference.file, 
#'     output.file="gc_file.txt")
#' 
#' @export calculateGCContentByInterval
#' @importFrom rtracklayer import
#' @importFrom Biostrings letterFrequency
calculateGCContentByInterval <- function(interval.file, reference.file,
output.file = NULL, ... ) {
    if (class(interval.file)=="GRanges") {
        interval.gr <- interval.file
        interval <- data.frame(Target=paste0(seqnames(interval.gr),":", 
            start(interval.gr),"-", end(interval.gr)))
    } else {    
        interval <- read.delim(interval.file, as.is=TRUE, ...)
        colnames(interval)[1] <- "Target"
        pos <- as.data.frame(do.call(rbind, strsplit(interval$Target, ":|-")),
            stringsAsFactors = FALSE)
        interval.gr <- GRanges(seqnames = pos[,1], 
            IRanges(start = as.numeric(pos[,2]), end = as.numeric(pos[,3])))
    }    
    if (min(start(interval.gr)) < 1) {
        .stopUserError("Interval coordinates should start at 1, not at 0.")
    }    
    x <- scanFa(reference.file, interval.gr)
    GC.count <- letterFrequency(x,"GC")
    all.count <- letterFrequency(x,"ATGC")
    gc <- data.frame(
        Target=interval$Target,
        gc_bias=as.vector(ifelse(all.count==0,NA,GC.count/all.count))
    )
    if (!is.null(output.file)) {
        write.table(gc, file=output.file, row.names=FALSE, quote=FALSE, 
            sep="\t")
    }    
    invisible(gc)
}
