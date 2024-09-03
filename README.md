# Generating SNP-Gene-Network using R Studio and Cytoscape Tool

## loading the library
library(biomaRt)
library(RCy3)
library(tidyr)
library(dplyr)
library(stringr)

## listing the databases
mymarts <- listMarts()

## selection of database
ensembl=useMart("ENSEMBL_MART_SNP")

## listing of datasets 
mydatasets <-listDatasets(ensembl)

## selection of dataset
ensembl=useMart("ENSEMBL_MART_SNP",
                dataset = "hsapiens_snp")

## listing of filters
myfilters <- listFilters(ensembl)

## making vector for filter
filter1 <- 'snp_filter'

## importing the SNPS list
variants <- read.table("variant_list.txt")

## listing of attributes
myattributes<- listAttributes(ensembl)

## making a variable for the attributes
att1 <- c('refsnp_id', 'ensembl_gene_name','snp','associated_gene')

searchResults <-getBM(att1, 
                      filters = filter1,
                      values = variants$V1, 
                      mart= ensembl)
## saving the results
write.csv(searchResults, file = "searchResults.csv", row.names = FALSE)

# Filtering data 
## Read the CSV file
results<- read.csv("searchResults.csv", stringsAsFactors = FALSE)

## Separate genes based on comma and semicolon, then filter unique genes
data_long <- results %>%
  separate_rows(associated_gene, sep = ",") %>%  # Separated by comma
  separate_rows(associated_gene, sep = ";") %>%  # Separated by semicolon
  mutate(genes = str_trim(associated_gene)) %>%  # Removing any leading or trailing whitespace
  distinct(refsnp_id, genes)  # Filter unique combinations of variant and gene

## save file
write.csv(data_long, "separated_unique_genes.csv", row.names = FALSE, quote = FALSE)

## Reading the File
unique_genes <- read.csv("separated_unique_genes.csv", stringsAsFactors = FALSE)

## Connecting to Cytoscape
cytoscapePing()

## Creating network and identifying nodes
createNetworkFromDataFrames(
  nodes = data.frame(id = unique(c(unique_genes$refsnp_id, unique_genes$genes))),
  edges = data.frame(source = unique_genes$refsnp_id, target = unique_genes$genes),
  title = "SNP-Gene Network"
)
