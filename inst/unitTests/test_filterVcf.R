test_filterVcf <- function() {
    vcf.file <- system.file("extdata", "example_vcf.vcf", package = "PureCN")
    vcf <- readVcf(vcf.file, "hg19")
    x <- filterVcfBasic(vcf)

    # TODO: remove in 1.8
    d.f <- data.frame(id=head(names(vcf)), Count=1)
    write.csv(d.f, file="snpbl.csv", row.names=FALSE)

    y <- filterVcfBasic(vcf, snp.blacklist="snpbl.csv")
    rtracklayer::export(head(rowRanges(vcf)), con="snpbl2.bed", format="bed")
    # TODO: remove in 1.8
    rtracklayer::export(head(rowRanges(vcf)), con="snpbl3.csv", format="bed")
    z <- filterVcfBasic(vcf, snp.blacklist="snpbl2.bed")
    zz <- filterVcfBasic(vcf, snp.blacklist="snpbl3.csv")

    checkEqualsNumeric(6, nrow(x$vcf)-nrow(y$vcf))
    checkEqualsNumeric(6, nrow(x$vcf)-nrow(z$vcf))
    checkEqualsNumeric(6, nrow(x$vcf)-nrow(zz$vcf))
    checkEqualsNumeric(6, sum(d.f$id %in% names(x$vcf)))
    checkEqualsNumeric(0, sum(d.f$id %in% names(y$vcf)))
    checkEqualsNumeric(0, sum(d.f$id %in% names(z$vcf)))
    checkEqualsNumeric(0, sum(d.f$id %in% names(zz$vcf)))
    gc.gene.file <- system.file("extdata", "example_gc.gene.file.txt", 
            package = "PureCN")
    
    # check wrong stats.file
    vcfMutectFilter <- filterVcfMuTect(vcf, stats.file=gc.gene.file)
    checkEqualsNumeric(nrow(x$vcf), nrow(vcfMutectFilter$vcf))

    # check stats file missing the failure_reasons column
    filename <- "statstest.txt"
    cat("#testfile\n", file=filename)
    d.f <- data.frame(head(rowRanges(vcf)))[,1:2]
    colnames(d.f) <- c("contig", "position")
    suppressWarnings(write.table(d.f, file=filename, append=TRUE, 
        quote=FALSE, row.names=FALSE, sep="\t"))
    vcfMutectFilter <- filterVcfMuTect(vcf, stats.file=filename)
    checkEqualsNumeric(6, nrow(vcfMutectFilter$vcf))
}    
