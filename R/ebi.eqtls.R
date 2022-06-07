# library(EnsDb.Hsapiens.v79)
#' @title ebi.eqtls
#' @description A template for building documented, tested R6 classes
#' @name ebi.eqtls
#'
#' @import catalogueR
#' @import EnsDb.Hsapiens.v79
#' @import plyr
#' @import SNPlocs.Hsapiens.dbSNP155.GRCh38
#' @import GenomicRanges
#'
#' @export

ebi.eqtls = R6Class("ebi.eqtls",

    #--------------------------------------------------------------------------------
    private = list(study.id=NULL,
                   catalog=NULL

                   ),

    #--------------------------------------------------------------------------------
    public = list(
         #' @description
         #' Creates a new instance of this [R6][R6::R6Class] class.
         #' @param study.id character, the ebi-archived eqtl study name
         #' @return a new instance of ebi.eqtls
        initialize = function(study.id=NA){
            private$study.id <- study.id
            private$catalog <-  get(data(meta, package="catalogueR"))
            },
        #------------------------------------------------------------
        #' @description accessor for the object's id field
        #' @return the current value of the id member
        getStudyId = function(){
            private$study.id
            },
        #------------------------------------------------------------
        #' @description retrieves the full ebi eqtl catalog
        #' @return data.frame
        getCatalog = function(){
            private$catalog
            },
        #------------------------------------------------------------
        #' @description gets all eQTLs in the region, or just those for the targetGene,
        #' @param chrom character
        #' @param start numeric
        #' @param end numeric
        #' @param studyIDs character from the catalog's unique_id field
        #' @param targetGene character default NA restricts eQTLs if supplied
        #' @param simplify logical return the full EBI table, or just a few crucial columsnn
        #' @return data.frame
        get.eQTLsByLocationAndStudyID = function(chrom, start, end, studyIDs, targetGene=NA, simplify){
            tbls <- list()
            for(study.id in studyIDs){
                message(sprintf("--- fetching %s (ge)", study.id))
                tryCatch({
                    suppressWarnings({tbl <- eQTL_Catalogue.fetch(unique_id=study.id,
                                                                  quant_method="ge",
                                                                  method="REST",
                                                                  chrom = sub("chr", "", chrom),
                                                                  bp_lower=start,
                                                                  bp_upper=end,
                                                                  verbose=TRUE)})
                    message(sprintf("%s %d-%d, %d", study.id, start, end, nrow(tbl)))
                    if(nrow(tbl) > 0){
                       tbl$id <- study.id
                       tbls[[study.id]] <- tbl
                       }
                   }, # tryCatch main block
                error = function(e){
                    message(sprintf("eQTL_Catalogue.fetch failed on study %s", study.id))
                    print(e)
                }) # tryCatch conclusion
               } # for study.id
               # create the table from the list of tables, map ensg to symbols, simplify
            tbl.out <- do.call(rbind.fill, tbls)
            rownames(tbl.out) <- NULL
            if(is.null(tbl.out))
                return(data.frame())
            new.order <- order(tbl.out$pvalue.QTL, decreasing=FALSE)
            tbl.out <- as.data.frame(tbl.out[new.order,])
            if(simplify){
                coi <- c("rsid.QTL", "pvalue.QTL", "gene_id.QTL", "an.QTL", "beta.QTL", "id",
                         "chromosome.QTL", "position.QTL")
                tbl.out <- tbl.out[, coi]
                colnames(tbl.out) <- c("rsid", "pvalue", "gene", "total.alleles", "beta", "id", "chrom", "hg38")
            } else {  # even though not simplifying, for uniformity we need "gene" as the column name
                colnames(tbl.out)[grep("gene_id.QTL", colnames(tbl.out))] <- "gene"
                }
            if(!is.na(targetGene)){  # find its ensg, subset on that, then replace
                ensg <- as.character(AnnotationDbi::mapIds(EnsDb.Hsapiens.v79, targetGene, "GENEID", "SYMBOL"))
                if(!ensg %in% tbl.out$gene){   # these eQTLs are not for the target gene
                    message(sprintf("ADv$geteQTLsByLocationAndStudyID, %d eqtls found, but none for %s",
                                    nrow(tbl.out), targetGene))
                    return(data.frame())
                    }
                tbl.out <- subset(tbl.out, gene==ensg)
                tbl.out$gene <- targetGene
                message(sprintf("--- %d variants for %s, corrected from %s", nrow(tbl.out), targetGene, ensg))
            } else {
                map <- AnnotationDbi::mapIds(EnsDb.Hsapiens.v79, tbl.out$gene, "SYMBOL", "GENEID")
                tbl.map <- data.frame(ensg=names(map), symbol=as.character(map), stringsAsFactors=FALSE)
                na.indices <- which(is.na(tbl.map$symbol))
                length(na.indices)
                tbl.map$symbol[na.indices] <- tbl.map$ensg[na.indices]
                tbl.out$gene <- tbl.map$symbol
                }
            rownames(tbl.out) <- NULL
            invisible(as.data.frame(tbl.out))
            }, # get.eQTLsByLocationAndStudyID

        #---------------------------------------------------------------------------------------------------
        #' @description gets all eQTLs in the region, chunking requests as specified
        #' @param chrom character
        #' @param start numeric
        #' @param end numeric
        #' @param study character from the catalog's unique_id field
        #' @param simplify logical return the full EBI table, or just a few crucial columsnn
        #' @param chunk.size integer often 10000, below the limit at which EBI REST fails
        #' @return data.frame
        fetch.eqtls.in.chunks = function(chrom, start, end, study, simplify, chunk.size){
            stopifnot(length(study) == 1)
            roi.width <- 1 + end - start
            if(roi.width <= chunk.size){
                message(sprintf("--- just one chunk"))
                tbl <- self$get.eQTLsByLocationAndStudyID(chrom, start, end, study, simplify=simplify)
            } else {
                boundaries.needed <- 1 + (roi.width %/% chunk.size)
                starts <- as.integer(seq(from=start, to=end, length.out=boundaries.needed))
                ends <- starts[-1]
                starts <- starts[-(length(starts))]
                tbls <- list()
                intervals <- length(starts)
                message(sprintf("==== pgc1a fetch.eqtls, %d chunks", intervals))
                for(i in seq_len(intervals)){
                    message(sprintf("--- fetching chunk %2d/%d for %s", i, intervals, study))
                    tbl.chunk <- self$get.eQTLsByLocationAndStudyID(chrom,
                                                                    as.integer(starts[i]),
                                                                    as.integer(ends[i]),
                                                                    study,
                                                                    targetGene=NA,
                                                                    simplify=simplify)
                    tbls[[i]] <- tbl.chunk
                    } # for i
                tbl <- do.call(rbind, tbls)
                } # else
            invisible(tbl)
            }, # fetch.eqtls.in.chunks
        #' @description create a table of hg19 and hg38 locs for all rsids in region
        #' @param chrom character
        #' @param start numeric
        #' @param end numeric
        create.snpLocs.table = function(chrom, start, end){
            require(GenomicRanges)
            require(SNPlocs.Hsapiens.dbSNP155.GRCh38)
            message(sprintf("initializing hg38 snpLocs, may take a minute"))
            t1 <- system.time(x <- snpsById(SNPlocs.Hsapiens.dbSNP155.GRCh38, "rs769450"))
            message(sprintf("dbSNP155 hg38: %5.2f", t1[["elapsed"]]))

            gr <- GRanges(seqnames=chrom, #sub("chr", "", chromosome),
                          IRanges(start=start, end=end))
            seqlevelsStyle(gr) <- "NCBI"
            gr.snps <- snpsByOverlaps(SNPlocs.Hsapiens.dbSNP155.GRCh38, gr)
            seqlevelsStyle(gr.snps) <- "UCSC"   # needed for rtracklayer liftover
            chain.file <- "hg38ToHg19.over.chain.gz"
            if(!file.exists(chain.file)){
                system(sprintf("curl -O http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/%s",
                               chain.file))
                system(sprintf("gunzip %s", chain.file))
                }
            chain <- import.chain(sub(".gz", "", chain.file, fixed=TRUE))
            x <- liftOver(gr.snps, chain)
            gr.hg19 <- unlist(x)
            tbl.snpLocs.hg19 <- as.data.frame(gr.hg19)[, c("seqnames", "start", "RefSNP_id")]
            colnames(tbl.snpLocs.hg19) <- c("chrom", "hg19", "rsid")

            tbl.snpLocs.hg38 <- as.data.frame(gr.snps)[, c(1,2,4)]
            colnames(tbl.snpLocs.hg38) <- c("chrom", "hg38", "rsid")
            tbl.snpLocs.hg38$chrom <- as.character(tbl.snpLocs.hg38$chrom)
            dim(tbl.snpLocs.hg38)

            tbl.snpLocs <- merge(tbl.snpLocs.hg38, tbl.snpLocs.hg19[, c("rsid", "hg19")], by=c("rsid"), all=TRUE)
            rownames(tbl.snpLocs) <- tbl.snpLocs$rsid
            coi <- c("chrom", "hg19", "hg38", "rsid")
            tbl.snpLocs <- tbl.snpLocs[, coi]
            rownames(tbl.snpLocs) <- NULL
            tbl.snpLocs
            } # create.snpLocs.table
        #----------------------------------------------------------------------------------------------------
       ) # public

    ) # class
#--------------------------------------------------------------------------------
