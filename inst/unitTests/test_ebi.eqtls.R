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
       # get eQTLs in the proximal promoter of TFAM.  hg38 is assumed.
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
if(!interactive())
    runTests()
