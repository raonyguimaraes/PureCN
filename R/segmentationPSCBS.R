segmentationPSCBS <-
structure(function(# PSCBS segmentation
### Segmentation function. Uses the PSCBS package.
### This function is called via the fun.segmentation argument of runAbsoluteCN.
### The arguments are passed via args.segmentation.
normal, 
### GATK coverage file for normal sample.
tumor,  
### GATK coverage file for tumor sample.
log.ratio, 
### Copy number log-ratios, one for each exon in coverage file.
plot.cnv, 
### Segmentation plots.
coverage.cutoff, 
### Minimum coverage in both normal and tumor.
sampleid=sampleid,
### Sample id, used in output files.
exon.weight.file=NULL,
### Can be used to assign weights to exons. NOT SUPPORTED YET.
flavor="tcn&dh",
### Flavor value for PSBCS. See segmentByNonPairedPSCBS.
tauA=0.03,
### tauA argument for PSCBS. See segmentByNonPairedPSCBS.
vcf=NULL,
### Optional VCF object with germline allelic ratios.
tumor.id.in.vcf=1,
### Id of tumor in case multiple samples are stored in VCF.
verbose=TRUE,
### Verbose output.
...
### Additional parameters passed to the segmentByNonPairedPSCBS function.
) {
    if (!requireNamespace("PSCBS", quietly = TRUE)) {
        stop("segmentationPSCBS requires the PSCBS package.")
    }

    exon.weights <- NULL
    if (!is.null(exon.weight.file)) {
        exon.weights <- read.delim(exon.weight.file, as.is=TRUE)
        exon.weights <- exon.weights[match(as.character(tumor[,1]), exon.weights[,1]),2]
        if (verbose) message("Exon weights found, but currently not supported by PSCBS.")
    }
    well.covered.exon.idx = ((normal$average.coverage > coverage.cutoff) & (tumor$average.coverage > coverage.cutoff)) | ((normal$average.coverage > 1.5 * coverage.cutoff) &  (tumor$average.coverage > 0.5 * coverage.cutoff))
    #MR: fix for missing chrX/Y 
    well.covered.exon.idx[is.na(well.covered.exon.idx)] <- FALSE
    tumor <- tumor[well.covered.exon.idx,]
    log.ratio <- log.ratio[well.covered.exon.idx]
    exon.gr <- GRanges(seqnames=tumor$chr, IRanges(start=tumor$probe_start, end=tumor$probe_end))
    ov <- findOverlaps(vcf, exon.gr)
    d.f <- cbind(tumor[subjectHits(ov),], CT=2^(log.ratio+1)[subjectHits(ov)], betaT=unlist(geno(vcf[queryHits(ov)])$FA[,tumor.id.in.vcf]), x=start(vcf[queryHits(ov)]) )
    d.f.2 <- cbind(tumor[-subjectHits(ov),], CT=2^(log.ratio+1)[-subjectHits(ov)], betaT=NA, x=tumor$probe_start[-subjectHits(ov)] )
    d.f.3 <- rbind(d.f, d.f.2)
    d.f.3 <- d.f.3[order(.strip.chr.name(d.f.3$chr), d.f.3$x),]
    d.f <- d.f.3
    seg <- PSCBS::segmentByNonPairedPSCBS(CT=d.f$CT, betaT=d.f$betaT, chromosome=.strip.chr.name(d.f$chr), x=d.f$x, tauA=tauA, flavor=flavor,...)
    .PSCBSoutput2DNAcopy(seg, sampleid)
### A list with elements seg and size. "seg" contains the segmentation, 
### "size" the size of all segments in base pairs.    
},ex=function() {
gatk.normal.file <- system.file("extdata", "example_normal.txt", package="PureCN")
gatk.tumor.file <- system.file("extdata", "example_tumor.txt", package="PureCN")
vcf.file <- system.file("extdata", "example_vcf.vcf", package="PureCN")
gc.gene.file <- system.file("extdata", "example_gc.gene.file.txt", package="PureCN")

ret <-runAbsoluteCN(gatk.normal.file=gatk.normal.file, gatk.tumor.file=gatk.tumor.file, 
   vcf.file=vcf.file, sampleid='Sample1', gc.gene.file=gc.gene.file, fun.segmentation=segmentationPSCBS)
})    

    
.PSCBSoutput2DNAcopy <- function(seg, sampleid) {
    sx <- cbind(ID=sampleid, seg$output[complete.cases(seg$output),])
    sx <- sx[,c("ID", "chromosome", "tcnStart", "tcnEnd", "tcnNbrOfLoci", "tcnMean")]
    colnames(sx) <- c("ID", "chrom", "loc.start",  "loc.end", "num.mark", "seg.mean")
    sx$seg.mean <- log2(sx$seg.mean/2)
    list(seg=sx, size=sx$loc.end-sx$loc.start+1)
}
    