# Load required libraries
library(DESeq2)
library(data.table)

# Define the file path
file_path <- "C:/Users/danny/Downloads/TCGA_Count_Matrix_0.csv"

# Read the count matrix efficiently
count_matrix <- fread(file_path, data.table=FALSE)

# Convert first column (gene IDs) to row names
rownames(count_matrix) <- count_matrix[,1]  # Set gene IDs as row names
count_matrix <- count_matrix[,-1]  # Remove first column after setting row names

# Check dimensions
dim(count_matrix)  # Should return (genes, 1231)

# Ensure column names are unique
colnames(count_matrix) <- make.names(colnames(count_matrix), unique=TRUE)

# Create sample metadata (Modify condition assignment as needed)
num_samples <- ncol(count_matrix)  # Should be 1231

# Assign conditions dynamically (Modify this if you have a proper metadata file)
condition_labels <- rep(c("Tumor", "Normal"), length.out = num_samples)

# Create sample_info dataframe
sample_info <- data.frame(
  sample = colnames(count_matrix),
  condition = condition_labels
)

# Set row names
rownames(sample_info) <- sample_info$sample

# Verify metadata
dim(sample_info)  # Should return (1231, 2)
table(sample_info$condition)  # Check counts of each condition

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = sample_info, design = ~ condition)

# Run differential expression analysis
dds <- DESeq(dds)
res <- results(dds)

# Save results
write.csv(res, "DESeq2_Results.csv")

# View top differentially expressed genes (sorted by adjusted p-value)
head(res[order(res$padj), ])