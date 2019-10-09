## ---- include = FALSE, warning = FALSE-----------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  warning = FALSE,
  message = FALSE
)

## ----setup---------------------------------------------------------------
require(intePareto)

## ------------------------------------------------------------------------
data("promoter", package = "intePareto")
file.bam <- system.file("extdata", "SRR925640.bam", 
                        package = "intePareto")
# note: this is just a subsampled test bam file to show how it works
b2c <- bam2counts(bamFile = file.bam,
                  region = promoter,
                  fragLength = 180)
tail(b2c)
b2r <- bam2rpm(bamFile = file.bam,
               region = promoter,
               fragLength = 180)
tail(b2r)


## ------------------------------------------------------------------------
data("test_rna_meta", package = "intePareto")
test_rna_meta

## ------------------------------------------------------------------------
data("test_chip_meta", package = "intePareto")
test_chip_meta

## ------------------------------------------------------------------------
# get the exact place of the tsv files for RNA-Seq and 
# bam files for ChIP-Seq data.
test_rna_meta$files <- system.file("extdata", 
                                   paste0(test_rna_meta$SRR,".tsv"), 
                                   package = "intePareto")


test_chip_meta$files <- system.file("extdata", 
                                    paste0(test_chip_meta$SRR,".bam"), 
                                    package = "intePareto")


# match of RNA-Seq and ChIP-Seq data on the gene level 
# through "weighted.mean" strategy.
res <- doMatch(rnaMeta = test_rna_meta,
               chipMeta = test_chip_meta,
               region = "promoter",
               method = "weighted.mean",
               ensemblDataset = "mmusculus_gene_ensembl")

## ------------------------------------------------------------------------
data("res", package = "intePareto")
df_final <- doIntegration(res = res,
                          ref = "wild.type",
                          type = "apeglm", 
                          apeAdapt = FALSE)

head(df_final)


## ------------------------------------------------------------------------
# chosse the first 3 fronts
objective <- data.frame(mark = c("z.H3K27ac","z.H3K4me3"), 
                        obj = c("max","max"), stringsAsFactors = FALSE)
nr.fronts <- 3
res_final <- doPareto(df_final = df_final, 
                      objective = objective, 
                      nr.fronts = nr.fronts)
head(res_final)
  

## ------------------------------------------------------------------------
sessionInfo()

