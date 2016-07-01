test_filterVcf <- function() {
    vcf.file <- system.file("extdata", "example_vcf.vcf", package = "PureCN")
    vcf <- readVcf(vcf.file, "hg19")
    x <- filterVcfBasic(vcf)

    d.f <- data.frame(id=head(names(vcf)), Count=1)
    write.csv(d.f, file="snpbl.csv", row.names=FALSE)
    y <- filterVcfBasic(vcf, snp.blacklist="snpbl.csv")
    d.f2 <- as.data.frame(head(rowRanges(vcf)))
    write.csv(d.f2, file="snpbl2.csv", row.names=FALSE)
    write.table(d.f2, file="snpbl3.csv", row.names=FALSE,sep="\t")
    z <- filterVcfBasic(vcf, snp.blacklist="snpbl2.csv")
    zz <- filterVcfBasic(vcf, snp.blacklist="snpbl3.csv")

    checkEqualsNumeric(6, nrow(x$vcf)-nrow(y$vcf))
    checkEqualsNumeric(6, nrow(x$vcf)-nrow(z$vcf))
    checkEqualsNumeric(6, nrow(x$vcf)-nrow(zz$vcf))
    checkEqualsNumeric(6, sum(d.f$id %in% names(x$vcf)))
    checkEqualsNumeric(0, sum(d.f$id %in% names(y$vcf)))
    checkEqualsNumeric(0, sum(d.f$id %in% names(z$vcf)))
    checkEqualsNumeric(0, sum(d.f$id %in% names(zz$vcf)))
}    