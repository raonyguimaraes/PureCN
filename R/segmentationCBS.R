segmentationCBS <-
structure(function(# CBS segmentation
### The default segmentation function. This function is called via the 
### fun.segmentation argument of runAbsoluteCN. The arguments are passed
### via args.segmentation.
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
### Can be used to assign weights to exons.
alpha=0.005,
### Alpha value for CBS, see documentation for the segment function.
vcf=NULL,
### Optional VCF object with germline allelic ratios.
tumor.id.in.vcf=1,
### Id of tumor in case multiple samples are stored in VCF.
verbose=TRUE
### Verbose output.
) {
    exon.weights <- NULL
    if (!is.null(exon.weight.file)) {
        exon.weights <- read.delim(exon.weight.file, as.is=TRUE)
        exon.weights <- exon.weights[match(as.character(tumor[,1]), exon.weights[,1]),2]
        if (verbose) message("Exon weights found, will use weighted CBS.")
    }
    x <- .CNV.analyze2(normal, tumor, logR=log.ratio, plot.cnv=plot.cnv, coverage.cutoff=coverage.cutoff, sampleid=sampleid, alpha=alpha, weights=exon.weights, verbose=verbose) 
    if (!is.null(vcf)) {
        x <- .pruneByVCF(x, vcf, tumor.id.in.vcf)
    }
    idx.enough.markers <- x$cna$output$num.mark > 1
    ##value<< A list with elements
    xx <- list(
        seg=x$cna$output[idx.enough.markers,], ##<< The segmentation.
        size=x$cnv$size[idx.enough.markers]  ##<< The size of all segments (in base pairs).
    )
##end<<
},ex=function() {
gatk.normal.file <- system.file("extdata", "example_normal.txt", package="PureCN")
gatk.tumor.file <- system.file("extdata", "example_tumor.txt", package="PureCN")
vcf.file <- system.file("extdata", "example_vcf.vcf", package="PureCN")
gc.gene.file <- system.file("extdata", "example_gc.gene.file.txt", package="PureCN")

# speed-up the runAbsoluteCN by using the stored grid-search.
# (purecn.example.output$candidates).
data(purecn.example.output)

# The max.candidate.solutions is set to a very low value only to speed-up this example.
# This is not a good idea for real samples.
ret <-runAbsoluteCN(gatk.normal.file=gatk.normal.file, gatk.tumor.file=gatk.tumor.file, 
   vcf.file=vcf.file, sampleid='Sample1', gc.gene.file=gc.gene.file, 
   candidates=purecn.example.output$candidates, max.candidate.solutions=2,
   fun.segmentation=segmentationCBS, args.segmentation=list(alpha=0.001))
})    

# looks at breakpoints, and if p-value is higher than max.pval, merge unless there is evidence based
# on germline SNPs
.pruneByVCF <- function(x, vcf, tumor.id.in.vcf, min.size=5, max.pval=0.00001, iterations=3, debug=FALSE) {
    seg <- segments.p(x$cna)
    for (iter in 1:iterations) {
        seg.gr <- GRanges(seqnames=.add.chr.name(seg$chrom), IRanges(start=seg$loc.start, end=seg$loc.end))
        ov <- findOverlaps(seg.gr, vcf)
        ar <- sapply(geno(vcf)$FA[,tumor.id.in.vcf], function(x) x[1])
        ar.r <- ifelse(ar>0.5, 1-ar, ar)
        merged <- rep(FALSE, nrow(seg))
        for (i in 2:nrow(seg)) {
            # don't try to merge chromosomes or very significant breakpoints
            if (is.na(seg$pval[i-1]) || seg$pval[i-1]<max.pval) next
            # don't merge when we have no germline data for segments    
            if (!(i %in% queryHits(ov) && (i-1) %in% queryHits(ov))) next
            ar.i <- list(ar.r[subjectHits(ov)][queryHits(ov)==i],ar.r[subjectHits(ov)][queryHits(ov)==i-1])
            if (length(ar.i[[1]]) < min.size || length(ar.i[[2]]) < min.size) next
            if (merged[i-1]) next
            
            p.t <- t.test(ar.i[[1]], ar.i[[2]])$p.value
            if (p.t>0.2) {
                merged[i] <- TRUE
                x$cna$output$seg.mean[i-1] <- weighted.mean(c(seg$seg.mean[i],seg$seg.mean[i-1]),w=c(seg$num.mark[i],seg$num.mark[i-1]))
                x$cnv$size[i-1] <- x$cnv$size[i]+x$cnv$size[i-1]
                x$cna$output$num.mark[i-1] <- seg$num.mark[i]+seg$num.mark[i-1]
                x$cna$output$loc.end[i-1] <- seg$loc.end[i]
                seg$pval[i-1] <- seg$pval[i]
            }
            if (debug) message(paste(i, "LR diff:", abs(seg$seg.mean[i]-seg$seg.mean[i-1]), "Size: ", seg$num.mark[i-1], "PV:", p.t, "PV bp:",seg$pval[i-1] , "Merged:", merged[i],"\n", sep=" "))
        }
        x$cna$output <- x$cna$output[!merged,]
        x$cnv <- x$cnv[!merged,]
        seg <- seg[!merged,]
    }
    x
}
    

# ExomeCNV version without the x11() calls 
.CNV.analyze2 <-
function(normal, tumor, logR=NULL, coverage.cutoff=15, normal.chrs=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"), normal.chr=normal.chrs, c=0.5, write.file=FALSE, file=NULL, weights=NULL, doDNAcopy=TRUE, sdundo=0.5, undo.splits="sdundo", smooth=TRUE, alpha=0.01, sampleid=NULL, plot.cnv=TRUE, verbose=TRUE) {
    `%+%` <- function(x,y) paste(x,y,sep="")
    normal.chrs = intersect(levels(normal$chr), normal.chrs)

    # first, do it for exons with enough coverage. MR: added less stringent cutoff in case normal looks great. these could be homozygous deletions in high purity samples
    well.covered.exon.idx = ((normal$average.coverage > coverage.cutoff) & (tumor$average.coverage > coverage.cutoff)) | ((normal$average.coverage > 1.5 * coverage.cutoff) &  (tumor$average.coverage > 0.5 * coverage.cutoff))
    #MR: fix for missing chrX/Y 
    well.covered.exon.idx[is.na(well.covered.exon.idx)] <- FALSE

    if (verbose) message(paste("Removing", sum(!well.covered.exon.idx), "low coverage exons."))
    if (is.null(logR)) norm.log.ratio = .calcLogRatio(normal, tumor, verbose)
    else norm.log.ratio = logR

    if (doDNAcopy) {

        CNA.obj = CNA(norm.log.ratio[well.covered.exon.idx], .strip.chr.name(normal$chr[well.covered.exon.idx]), (normal$probe_start[well.covered.exon.idx] + normal$probe_end[well.covered.exon.idx])/2, data.type="logratio", sampleid=sampleid)
        smoothed.CNA.obj = if (smooth) smooth.CNA(CNA.obj) else CNA.obj
        if (!is.null(weights)) { 
            weights <- weights[well.covered.exon.idx]
            # MR: this shouldn't happen. In doubt, count them as maximum (assuming that poorly performing exons are down-weighted)
            weights[is.na(weights)] <- max(weights, na.rm=TRUE)
            segment.smoothed.CNA.obj = segment(smoothed.CNA.obj, undo.splits=undo.splits, undo.SD=sdundo, verbose=ifelse(verbose, 1, 0), alpha=alpha,weights=weights)
        } else {
            segment.smoothed.CNA.obj = segment(smoothed.CNA.obj, undo.splits=undo.splits, undo.SD=sdundo, verbose=ifelse(verbose, 1, 0), alpha=alpha)
        }        

        if (plot.cnv) {
            if (write.file && !is.null(file)) png(file, width=2000, height=1000, units="px")
            else if (write.file) png("CNV detection for exons with > " %+% coverage.cutoff %+% " coverage.png", width=2000, height=1000, units="px")
            plot(segment.smoothed.CNA.obj, plot.type="s")
            if (write.file) dev.off()

            if (write.file && !is.null(file)) png("allchr." %+% file, width=2000, height=1000, units="px")
            else if (write.file) png("all chromosome CNV detection for exons with > " %+% coverage.cutoff %+% " coverage.png", width=2000, height=1000, units="px") 
            plot(segment.smoothed.CNA.obj, plot.type="w")
            abline(h=log2(c + (1-c)*c(1,3,4,5)/2), col="purple")
            if (write.file) dev.off()
        }

        cnv = .get.proper.cnv.positions(normal[well.covered.exon.idx,], print(segment.smoothed.CNA.obj))

        return(list(cnv=cnv, cna=segment.smoothed.CNA.obj, logR=norm.log.ratio))

    } else {

        logR.mean = mean(norm.log.ratio[well.covered.exon.idx])
        logR.sd = sd(norm.log.ratio[well.covered.exon.idx])
        logR.min = min(norm.log.ratio[well.covered.exon.idx])
        logR.max = max(norm.log.ratio[well.covered.exon.idx])

        if (plot.cnv) {
            if (write.file && !is.null(file)) png(file, width=2000, height=1000, units="px")
            else if (write.file) png("CNV detection for exons with > " %+% coverage.cutoff %+% " coverage.noDNAcopy.png", width=2000, height=1000, units="px")
            par(mfrow=c(4,6))
            for (chr in levels(normal$chr)) {
                plot((normal$probe_start[well.covered.exon.idx & normal$chr==chr] + normal$probe_end[well.covered.exon.idx & normal$chr==chr])/2, norm.log.ratio[well.covered.exon.idx & normal$chr==chr], pch="*", pc=20, ylim=c(logR.min, logR.max), main=chr, xlab="position", ylab="log ratio")
                abline(h=logR.mean + logR.sd, col="red")
                abline(h=logR.mean - logR.sd, col="red")
                abline(h=0, col="gray")
            }
            if (write.file) dev.off()
        }

        return(list(logR=norm.log.ratio))

    }
}

.get.proper.cnv.positions <- function (exons, cnv) 
{
    chr.hash <- NULL    
    data(chr.hash, envir=environment())
    `%+%` <- function(x, y) paste(x, y, sep = "")
    order.by.chr = order(.strip.chr.name(exons$chr))
    exons = exons[order.by.chr, ]
    cnv$chr = as.character(chr.hash[cnv$chrom, "chr"])
    cnv$probe = "cnv" %+% as.character(1:nrow(cnv))
    end.idx = cumsum(cnv$num.mark)
    start.idx = c(1, 1 + end.idx[-length(end.idx)])
    cnv$probe_start = exons$probe_start[start.idx]
    cnv$probe_end = exons$probe_end[end.idx]
    cnv$size = cnv$probe_end - cnv$probe_start + 1
    sum.chunk = function(i, colName) {
        sum(exons[start.idx[i]:end.idx[i], colName])
    }
    cnv$targeted.base = sapply(1:nrow(cnv), sum.chunk, colName = "targeted.base")
    cnv$sequenced.base = sapply(1:nrow(cnv), sum.chunk, colName = "sequenced.base")
    cnv$coverage = sapply(1:nrow(cnv), sum.chunk, colName = "coverage")
    cnv$average.coverage = cnv$coverage/cnv$targeted.base
    cnv$base.with..10.coverage = sapply(1:nrow(cnv), sum.chunk, 
        colName = "base.with..10.coverage")
    return(cnv)
}
