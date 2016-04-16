correctCoverageBias <- structure(function(# Correct for GC bias
### Takes as input coverage data in GATK format (or data read by readCoverageGatk)
### and a mapping file for GC content, and then uses a loess normalization for
### bias correction. Largely follows the GC correction of the TitanCNA package.
gatk.coverage.file, 
### Exon coverage file as produced by GATK.
gc.gene.file,
### GC gene file, which provides GC content for each exon in the coverage files. 
### First column in format CHR:START-END. Second column GC content (0 to 1). 
### Third column provides gene symbols, which are optional, but used in
### runAbsoluteCN to generate gene level calls.
output.file=NULL
### Optionally, write minimal GATK coverage file with GC corrected coverage.
) {
    if (is.character(gatk.coverage.file)) {
        tumor  <- readCoverageGatk(gatk.coverage.file)
    } else {
        tumor <- gatk.coverage.file
    }    
    
    gc <- read.delim(gc.gene.file)
    if (!identical(as.character(gc[,1]), as.character(tumor[,1]))) {
        stop("Interval files in gatk.coverage.file and gc.gene.file different.")
    }    

    # taken from TitanCNA
    gc$valid <- TRUE
    gc$valid[tumor$average.coverage <= 0 | gc$gc_bias < 0] <- FALSE
    gc$ideal <- TRUE
    routlier <- 0.01
    range <- quantile(tumor$average.coverage[gc$valid], prob = c(0, 1 - routlier), 
        na.rm = TRUE)
    doutlier <- 0.001
    domain <- quantile(gc$gc_bias[gc$valid], prob = c(doutlier, 1 - doutlier), 
        na.rm = TRUE)
    gc$ideal[!gc$valid | 
        tumor$average.coverage <= range[1] |
        tumor$average.coverage > range[2] | 
        gc$gc_bias < domain[1] | 
        gc$gc_bias > domain[2]] <- FALSE

    rough <- loess(tumor$average.coverage[gc$ideal] ~ gc$gc_bias[gc$ideal], span = 0.03)
    i <- seq(0, 1, by = 0.001)
    final <- loess(predict(rough, i) ~ i, span = 0.3)
    cor.gc <- predict(final, gc$gc_bias)
    cor.gc.factor <- cor.gc/mean(tumor$average.coverage, na.rm=TRUE)

    tumor$average.coverage <- tumor$average.coverage / cor.gc.factor
    tumor$coverage <- tumor$coverage / cor.gc.factor
    if (!is.null(output.file)) {
        tmp <- tumor[, c("probe", "coverage", "average.coverage")]
        colnames(tmp) <- c("Target", "total_coverage", "average_coverage")
        write.table(tmp, file=output.file, row.names=FALSE, quote=FALSE)
    }    
    tumor
### GC normalized coverage
}, ex=function() {
gatk.normal.file <- system.file("extdata", "example_normal.txt", 
    package="PureCN")
gc.gene.file <- system.file("extdata", "example_gc.gene.file.txt", 
    package="PureCN")
coverage <- correctCoverageBias(gatk.normal.file, gc.gene.file)
})