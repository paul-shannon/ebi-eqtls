library(RUnit)
library(ebi.eqtls)
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_ctor()
   test_catalogContents()
   test_get.eQTLs.smallRegion()
   test_get.eQTLs.bigRegion.with.chunks()

} # runTests
#----------------------------------------------------------------------------------------------------
test_ctor <- function()
{
    message(sprintf("--- test_ctor"))

    ee <- ebi.eqtls$new()
    checkTrue(all(c("R6", "ebi.eqtls") %in% class(ee)))

    tbl.cat <- ee$getCatalog()
    checkTrue(nrow(tbl.cat) >= 498)
    checkEquals(ncol(tbl.cat), 12)

} # test_ctor
#----------------------------------------------------------------------------------------------------
test_catalogContents <- function()
{
    message(sprintf("--- test_catalogContents"))
    ee <- ebi.eqtls$new()
    tbl.cat <- ee$getCatalog()
       # get all the unique ge (geneexpression) studies
    table(tbl.cat$unique_id)
    tbl.brain <- subset(tbl.cat, quant_method=="ge" & grepl("brain", unique_id, ignore.case=TRUE))
    checkTrue(nrow(tbl.brain) >= 30)
    checkEquals(length(grep("GTEx_V8", tbl.brain$unique_id)), 13)

} # test_catalogContents
#----------------------------------------------------------------------------------------------------
test_get.eQTLs.smallRegion <- function()
{
    message(sprintf("--- test_get.eQTLs.smallRegion"))

    ee <- ebi.eqtls$new()

       #--------------------------------------------------------------
       # get eQTLs in the proximal promoter of TFAM.  hg38 is assumed.
       # all genes, simplify
       #--------------------------------------------------------------

    tbl.eqtls <- ee$get.eQTLsByLocationAndStudyID(chrom="chr10",
                                                  start=58384760,
                                                  end=58385954,
                                                  studyIDs=c("GTEx_V8.Brain_Cortex", "GTEx_V8.Brain_Cerebellum"),
                                                  targetGene=NA,
                                                  simplify=TRUE)
    target.genes <- unique(tbl.eqtls$gene)
    checkTrue(length(target.genes) > 5)
    studies <- unique(tbl.eqtls$id)
    checkEquals(length(studies), 2)
    rsids <- unique(tbl.eqtls$rsid)
    checkTrue(length(rsids) > 5)
    checkEquals(colnames(tbl.eqtls), c("rsid","pvalue","gene","total.alleles","beta","id","chrom","hg38"))

       #-------------------------------------------
       # hand-check one rsid
       #-------------------------------------------

    rsid.oi <- "rs12247015"
    checkTrue(rsid.oi %in% tbl.eqtls$rsid)  # actually shows up 18 times, different target genes, both studies
    snp.info <- as.list(subset(tbl.eqtls, rsid==rsid.oi & gene=="TFAM" & id=="GTEx_V8.Brain_Cerebellum"))
    checkEquals(snp.info$rsid, "rs12247015")
    checkEqualsNumeric(snp.info$pvalue, 2.32911e-11, tol=1e-9)
    checkEquals(snp.info$gene, "TFAM")
    checkEquals(snp.info$total.alleles, 418)
    checkEqualsNumeric(snp.info$beta, -0.257902, tol=1e-5)
    checkEquals(snp.info$id, "GTEx_V8.Brain_Cerebellum")
    checkEquals(snp.info$chrom, "10")
    checkEquals(snp.info$hg38, 58385319)

       #--------------------------------------------------------------
       # get eQTLs in the proximal promoter of TFAM.
       # all genes, do not simplify
       #--------------------------------------------------------------

     tbl.eqtls <- ee$get.eQTLsByLocationAndStudyID(chrom="chr10",
                                                   start=58384760,
                                                   end=58385954,
                                                   studyIDs=c("GTEx_V8.Brain_Cortex", "GTEx_V8.Brain_Cerebellum"),
                                                   targetGene=NA,
                                                   simplify=FALSE)
     target.genes <- unique(tbl.eqtls$gene)
     checkTrue(length(target.genes) > 5)
     studies <- unique(tbl.eqtls$id)
     checkEquals(length(studies), 2)
     rsids <- unique(tbl.eqtls$rsid)
     checkTrue(length(rsids) > 5)
     checkEquals(colnames(tbl.eqtls),
                 c("qtl_id","study_id.QTL","qtl_group.QTL","rsid.QTL","chromosome.QTL","position.QTL",
                   "pvalue.QTL","condition_label.QTL","tissue_label.QTL","molecular_trait_id.QTL",
                   "gene","ac.QTL","ref.QTL","beta.QTL","variant.QTL","an.QTL","median_tpm.QTL",
                   "condition.QTL","r2.QTL","alt.QTL","type.QTL","maf.QTL","tissue.QTL","gene.QTL","id"))


       #-------------------------------------------
       # specify a targetGene, simplify
       #-------------------------------------------

    tbl.eqtls <- ee$get.eQTLsByLocationAndStudyID(chrom="chr10",
                                                  start=58384760,
                                                  end=58385954,
                                                  studyIDs=c("GTEx_V8.Brain_Cortex", "GTEx_V8.Brain_Cerebellum"),
                                                  targetGene="TFAM",
                                                  simplify=TRUE)

    checkTrue(nrow(tbl.eqtls) > 10)   # 12 on (7 jun 2022)
    checkTrue(all(tbl.eqtls$gene=="TFAM"))

       #-------------------------------------------
       # specify a targetGene, don't simplify
       #-------------------------------------------

    tbl.eqtls <- ee$get.eQTLsByLocationAndStudyID(chrom="chr10",
                                                  start=58384760,
                                                  end=58385954,
                                                  studyIDs="GTEx_V8.Brain_Cerebellum",
                                                  targetGene="TFAM",
                                                  simplify=FALSE)
    checkEquals(dim(tbl.eqtls), c(6,25))
    checkTrue(all(tbl.eqtls$gene == "TFAM"))
    checkEquals(colnames(tbl.eqtls),
                c("qtl_id","study_id.QTL","qtl_group.QTL","rsid.QTL","chromosome.QTL","position.QTL",
                  "pvalue.QTL","condition_label.QTL","tissue_label.QTL","molecular_trait_id.QTL",
                  "gene","ac.QTL","ref.QTL","beta.QTL","variant.QTL","an.QTL","median_tpm.QTL",
                  "condition.QTL","r2.QTL","alt.QTL","type.QTL","maf.QTL","tissue.QTL","gene.QTL","id"))

} # test_get.eQTLs.smallRegion
#----------------------------------------------------------------------------------------------------
test_get.eQTLs.bigRegion.with.chunks <- function()
{
    message(sprintf("--- test_get.eQTLs.bigRegion.with.chunks"))

    ee <- ebi.eqtls$new()

       # a silly small start
    tbl.small <- ee$fetch.eqtls.in.chunks(chrom="chr10",
                                          start=58384760,
                                            end=58385954,
                                          study="GTEx_V8.Brain_Cortex",
                                          simplify=TRUE,
                                          chunk.size=2000)
    checkTrue(nrow(tbl.small) > 45)


       # now a 40k region, with 10k chunks
    tbl.bigger <- ee$fetch.eqtls.in.chunks(chrom="chr10",
                                          start=58384760,
                                            end=58424760,
                                          study="GTEx_V8.Brain_Cortex",
                                          simplify=TRUE,
                                          chunk.size=10000)
    checkTrue(nrow(tbl.bigger) > 1000)


} # test_get.eQTLs.bigRegion.with.chunks
#----------------------------------------------------------------------------------------------------
test_create.snpLocs.table <- function()
{
    message(sprintf("--- test_create.snpLocs.table"))

    ee <- ebi.eqtls$new()
    tbl.snpLocs <- ee$create.snpLocs.table(chrom="chr10",
                                           start=58384760,
                                           end=583859540)

} # test_create.snpLocs.table
#----------------------------------------------------------------------------------------------------
if(!interactive())
    runTests()
