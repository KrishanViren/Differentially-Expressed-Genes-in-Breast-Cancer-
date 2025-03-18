# Install Bioconductor Manager if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install necessary bioinformatics packages
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "enrichplot", "ReactomePA", "AnnotationDbi"))


# Load the DEGs file
deg <- read.csv("C:/Users/danny/Downloads/Significant_DEGs.csv", header = TRUE)

# View the first few rows
head(deg)

# Convert Gene Symbols to Entrez IDs
entrez_ids <- select(org.Hs.eg.db, 
                     keys = deg$GeneSymbol, 
                     keytype = "SYMBOL", 
                     columns = "ENTREZID")

# Merge Entrez IDs with DEGs dataset
deg <- merge(deg, entrez_ids, by.x = "GeneSymbol", by.y = "SYMBOL")

# Check if conversion was successful
head(deg)

library(AnnotationDbi)
library(org.Hs.eg.db)


methods(select)

entrez_ids <- AnnotationDbi::select(
  x = org.Hs.eg.db,      # Explicitly specify the database
  keys = deg$GeneSymbol, # Your gene symbols
  keytype = "SYMBOL", 
  columns = "ENTREZID"
)


# Ensure your DEGs dataset is loaded
deg <- read.csv("C:/Users/danny/Downloads/Significant_DEGs.csv", header = TRUE)

# Convert Gene Symbols to Entrez IDs
entrez_ids <- AnnotationDbi::select(org.Hs.eg.db, 
                                    keys = deg$GeneSymbol, 
                                    keytype = "SYMBOL", 
                                    columns = "ENTREZID")

# Merge Entrez IDs with DEGs dataset
deg <- merge(deg, entrez_ids, by.x = "GeneSymbol", by.y = "SYMBOL")

# View results
head(deg)


# Perform KEGG pathway enrichment
kegg_results <- enrichKEGG(gene = deg$ENTREZID, 
                           organism = "hsa", 
                           keyType = "kegg", 
                           pvalueCutoff = 0.05, 
                           qvalueCutoff = 0.05)

# View top KEGG pathways
head(kegg_results@result)

# Visualize KEGG results using a dot plot
dotplot(kegg_results, showCategory = 15, title = "KEGG Pathway Enrichment")



# Perform GO enrichment (Biological Process)
go_results <- enrichGO(gene = deg$ENTREZID, 
                       OrgDb = org.Hs.eg.db, 
                       keyType = "ENTREZID", 
                       ont = "BP", 
                       pAdjustMethod = "BH", 
                       pvalueCutoff = 0.05, 
                       qvalueCutoff = 0.05)

# View top GO terms
head(go_results@result)

# Visualize GO results
dotplot(go_results, showCategory = 15, title = "GO Enrichment (Biological Process)")



# Perform Reactome Pathway Analysis
reactome_results <- enrichPathway(gene = deg$ENTREZID, 
                                  organism = "human", 
                                  pvalueCutoff = 0.05)

# View top Reactome pathways
head(reactome_results@result)

# Plot Reactome enrichment results
dotplot(reactome_results, showCategory = 15, title = "Reactome Pathway Enrichment")




write.csv(kegg_results@result, "KEGG_Enrichment_Results.csv", row.names = FALSE)
write.csv(go_results@result, "GO_Enrichment_Results.csv", row.names = FALSE)
write.csv(reactome_results@result, "Reactome_Enrichment_Results.csv", row.names = FALSE)

# KEGG Bar Plot
barplot(kegg_results, showCategory = 10, title = "Top 10 KEGG Pathways")

# GO Bar Plot
barplot(go_results, showCategory = 10, title = "Top 10 GO Terms (Biological Process)")

# Reactome Bar Plot
barplot(reactome_results, showCategory = 10, title = "Top 10 Reactome Pathways")



library(clusterProfiler)
ls("package:clusterProfiler")


library(ReactomePA)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ReactomePA", force = TRUE)
library(ReactomePA)
ls("package:ReactomePA")


