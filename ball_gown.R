library(ballgown)

args <- commandArgs(trailingOnly = TRUE)

sampleids = paste0(args[1:length(args)-1],collapse = ',' ) 
output = args[length(args)]

## generate the  expression table for each sample
sampleids = '06_ballgown/test_BC-sham-ctr-19-02-14/'
bg <-  ballgown(samples = sampleids , meas='FPKM')

print('process sample data')

print(sampleNames(bg))
gene_expression = ballgown::gexpr(bg)
gene_table = as.tibble(gene_expression)
gene_table$geneid = rownames(gene_expression)

geneID = tibble::tibble(genename = ballgown::geneNames(bg) ,
                        geneid =ballgown::geneIDs(bg) ) %>%
  dplyr::filter(! genename =='.' ) %>% dplyr::distinct()

gene_table <- dplyr::inner_join(gene_table ,geneID )
readr::write_csv(gene_table, path = output)
