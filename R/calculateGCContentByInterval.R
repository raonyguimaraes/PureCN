calculateGCContentByInterval <- structure(
function(# Calculates GC content by interval
### Uses scanFa from the Rsamtools package to retrieve GC 
### content of intervals in a reference FASTA file.
interval.file,
### File specifying the intervals. Interval is expected in 
### first column in format CHR:START-END.
reference.file,
### Reference Fasta file.
output.file=NULL
### Optionally, write GC content file. 
) {
    interval <- read.delim(interval.file, as.is=TRUE)
    colnames(interval)[1] <- "Target"
    pos <- as.data.frame(do.call(rbind, strsplit(interval$Target, ":|-")), 
        stringsAsFactors = FALSE)
    interval.gr <- GRanges(seqnames = pos[,1], 
        IRanges(start = as.numeric(pos[,2]), end = as.numeric(pos[,3])))
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
### Returns GC content by interval.
}, ex=function() {
reference.file <- system.file("extdata", "ex2_reference.fa", 
    package="PureCN", mustWork = TRUE)
interval.file <- system.file("extdata", "ex2_intervals.txt", 
    package="PureCN", mustWork = TRUE)
calculateGCContentByInterval(interval.file, reference.file, 
    output.file="gc_file.txt")
}) 