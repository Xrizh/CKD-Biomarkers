options(scipen = 999)

library(GEOquery)
library(preprocessCore)
library(sva)
library(randomForest)
library(pheatmap)
library(ggplot2)
library(limma)
library(AnnotationDbi)
library(clusterProfiler)
library(org.Hs.eg.db)
library(glmnet)
library(caret)
library(pROC)
library(ComplexHeatmap)
library(circlize)
library(illuminaHumanv4.db)
library(hgu133plus2.db)
library(reshape2)
library(DOSE)
library(ggrepel)



#GSE180394 (external)
# Define the GEO accession numbers for the datasets
accession_numbers <- c("GSE35488", "GSE32591", "GSE66494", "GSE47184")

# Create an empty list to store the datasets, expression data, and labels
datasets <- list()

# Loop through the accession numbers and fetch the datasets
for (accession in accession_numbers) {
  # Fetch the dataset from GEO
  dataset <- getGEO(accession, GSEMatrix = TRUE)
  
  # Store the dataset in the list
  datasets[[accession]] <- dataset
  
  # Extract the expression data table
  expression_data <- exprs(dataset[[1]])
  
  # Store the expression data in the list as well
  datasets[[paste0(accession, "_expression")]] <- expression_data
}

# Now you can access the datasets using their accession numbers
gse35488 <- datasets[["GSE35488"]]
gse32591 <- datasets[["GSE32591"]]
gse66494 <- datasets[["GSE66494"]]
gse47184 <- datasets[["GSE47184"]]

gse35488_expression <- datasets[["GSE35488_expression"]]
gse32591_expression <- datasets[["GSE32591_expression"]]
gse66494_expression <- datasets[["GSE66494_expression"]]
gse47184_expression <- gse47184$`GSE47184-GPL14663_series_matrix.txt.gz`@assayData$exprs

#Labels for GSE47184
# Extract the source name column
source_name <- gse47184$`GSE47184-GPL14663_series_matrix.txt.gz`@phenoData@data$`disease status:ch1`

# Initialize an empty vector to store the labels
labelsGSE47184 <- numeric(length(source_name))

# Assign labels based on source name
for (i in 1:length(source_name)) {
  if (is.na(source_name[i])) {
    labelsGSE47184[i] <- 0  # Control
  } else if (source_name[i] == "Diabetic Nephropathy" || 
             source_name[i] == "Minimal Change Disease" || 
             source_name[i] == "Thin Membrane Disease" || 
             source_name[i] == "Tumor Nephrectomy" || 
             source_name[i] == "Focal and Segmental Glomerulosclerosis" || 
             source_name[i] == "Hypertensive nephropathy" || 
             source_name[i] == "IgA Nephropathy" ||
             source_name[i] == "Membranous Glomerulonephritis") {
    labelsGSE47184[i] <- 1  # Positive for CKD
  } else {
    labelsGSE47184[i] <- NA  # Unknown label
  }
}


# Map ENTREZIDs to gene symbols
gene_symbols <- mapIds(org.Hs.eg.db,
                       keys = rownames(gse47184_expression),
                       keytype = "ENTREZID",
                       column = "SYMBOL",
                       multiVals = "first") # Takes the first symbol if there are multiple

# Identify the genes with NA symbols (unmapped) and remove them
valid_indices <- which(!is.na(gene_symbols))
gene_symbols <- gene_symbols[valid_indices]
gse47184_expression <- gse47184_expression[valid_indices, ]

# Update the row names of the expression matrix to the corresponding gene symbols
rownames(gse47184_expression) <- gene_symbols


#Labels for GSE35488
# Extract the source name column
source_name <- gse35488$GSE35488_series_matrix.txt.gz@phenoData@data$source_name_ch1

# Initialize an empty vector to store the labels
labelsGSE35488 <- numeric(length(source_name))

# Assign labels based on source name
for (i in 1:length(source_name)) {
  if (source_name[i] == "Kidney biopsy from human IgAN") {
    labelsGSE35488[i] <- 1  # Positive for CKD
  } else if (source_name[i] == "Kidney biopsy from human pre-transplant living donor") {
    labelsGSE35488[i] <- 0  # Control
  } else {
    labelsGSE35488[i] <- NA  # Unknown label
  }
}

#Labels for GSE32591
# Extract the source name column
source_name <- gse32591$GSE32591_series_matrix.txt.gz@phenoData@data$`disease status:ch1`

# Initialize an empty vector to store the labels
labelsGSE32591 <- numeric(length(source_name))

# Assign labels based on source name
for (i in 1:length(source_name)) {
  if (source_name[i] == "LN patient") {
    labelsGSE32591[i] <- 1  # Positive for CKD
  } else if (source_name[i] == "control") {
    labelsGSE32591[i] <- 0  # Control
  } else {
    labelsGSE32591[i] <- NA  # Unknown label
  }
}

#Labels for GSE66494
# Extract the source name column
source_name <- gse66494$GSE66494_series_matrix.txt.gz@phenoData@data$`disease status:ch1`

# Initialize an empty vector to store the labels
labelsGSE66494 <- numeric(length(source_name))

# Assign labels based on source name
for (i in 1:length(source_name)) {
  if (source_name[i] == "chronic kidney disease (CKD)") {
    labelsGSE66494[i] <- 1  # Positive for CKD
  } else if (source_name[i] == "normal kidney") {
    labelsGSE66494[i] <- 0  # Control
  } else {
    labelsGSE66494[i] <- NA  # Unknown label
  }
}

# Map ENTREZIDs to gene symbols
gene_symbols <- mapIds(org.Hs.eg.db,
                       keys = rownames(gse35488_expression),
                       keytype = "ENTREZID",
                       column = "SYMBOL",
                       multiVals = "first") # Takes the first symbol if there are multiple

# Identify the genes with NA symbols (unmapped) and remove them
valid_indices <- which(!is.na(gene_symbols))
gene_symbols <- gene_symbols[valid_indices]
gse35488_expression <- gse35488_expression[valid_indices, ]
gse32591_expression <- gse32591_expression[valid_indices, ]

# Update the row names of the expression matrix to the corresponding gene symbols
rownames(gse35488_expression) <- gene_symbols
rownames(gse32591_expression) <- gene_symbols

# Fetch the GEO dataset by accession number
geo_data <- getGEO("GPL6480")

# Access the expression data (data table)
data_table <- geo_data@dataTable@table

# Create a data frame from the platform data
database <- data.frame(ID = geo_data@dataTable@table$ID, GENE = geo_data@dataTable@table$GENE, GENE_SYMBOL = geo_data@dataTable@table$GENE_SYMBOL)
rownames(gse66494_expression) <- database$GENE_SYMBOL[match(rownames(gse66494_expression), database$ID)]

# Combine expression datasets

# Identify common genes across all three datasets
common_genes <- Reduce(intersect, list(rownames(gse35488_expression), rownames(gse32591_expression), rownames(gse66494_expression), rownames(gse47184_expression)))

# Subset each dataset to only include the common genes
gse35488_expression_common <- gse35488_expression[common_genes, ]
gse32591_expression_common <- gse32591_expression[common_genes, ]
gse66494_expression_common <- gse66494_expression[common_genes, ]
gse47184_expression_common <- gse47184_expression[common_genes, ]

# Sort rows by row names for each dataset
gse35488_expression_common <- gse35488_expression_common[order(rownames(gse35488_expression_common)), ]
gse32591_expression_common <- gse32591_expression_common[order(rownames(gse32591_expression_common)), ]
gse66494_expression_common <- gse66494_expression_common[order(rownames(gse66494_expression_common)), ]
gse47184_expression_common <- gse47184_expression_common[order(rownames(gse47184_expression_common)), ]

# Now, rbind these datasets
all_expression_data <- cbind(gse35488_expression_common, gse32591_expression_common, gse66494_expression_common, gse47184_expression_common)

# Combine labels
all_labels <- c(labelsGSE35488, labelsGSE32591, labelsGSE66494, labelsGSE47184)

# Remove the last patient from the all_expression_data
all_expression_data <- all_expression_data[, -ncol(all_expression_data)]

# Remove the last label from all_labels
all_labels <- all_labels[-length(all_labels)]

# Normalize entire dataset
all_expression_data_normalized <- normalize.quantiles(all_expression_data)

# Set the row and column names for the normalized data
rownames(all_expression_data_normalized) <- rownames(all_expression_data)
colnames(all_expression_data_normalized) <- colnames(all_expression_data)

# Log2 transform the normalized data
all_expression_data_log2 <- log2(all_expression_data_normalized)

# Remove rows with any NA values
all_expression_data_log2 <- all_expression_data_log2[!apply(is.na(all_expression_data_log2), 1, any), ]

# Define the number of samples based on the number of columns in the expression data
num_samples <- ncol(all_expression_data_log2)

# Shuffle the expression data and labels in the same order
set.seed(123)  # For reproducibility
random_order <- sample(1:num_samples)
all_expression_data_log2_shuffled <- all_expression_data_log2[, random_order]
all_labels_shuffled <- all_labels[random_order]

# Update batch_indicator for the shuffled data
batch_indicator <- factor(rep(c("Batch1", "Batch2"), each = num_samples / 2))
batch_indicator_shuffled <- batch_indicator[random_order]

# Perform ComBat on the entire dataset
corrected_data <- ComBat(dat = all_expression_data_log2_shuffled, batch = batch_indicator_shuffled)

#----------------------------------------------------

# Melt the dataframes for ggplot
melted_data_before <- melt(all_expression_data_log2_shuffled)
melted_data_before$Correction <- 'After'

melted_data_after <- melt(corrected_data)
melted_data_after$Correction <- 'Before'

# Combine melted dataframes
all_melted_data <- rbind(melted_data_before, melted_data_after)

# Assign dataset identities based on column names of expression data
dataset_assignments <- c(rep("GSE35488", ncol(gse35488_expression_common)),
                         rep("GSE32591", ncol(gse32591_expression_common)),
                         rep("GSE66494", ncol(gse66494_expression_common)),
                         rep("GSE47184", ncol(gse47184_expression_common)))

names(dataset_assignments) <- colnames(all_expression_data_log2_shuffled)

# Add dataset information to the melted data
all_melted_data$Dataset <- dataset_assignments[all_melted_data$Var2]

# Create the box plot with the dataset as both the fill and color variable
ggplot(all_melted_data, aes(x = Var2, y = value, fill = Dataset, color = Dataset)) +
  geom_boxplot(outlier.shape = NA) + # hide outliers for cleaner look
  facet_wrap(~Correction, ncol = 1, scales = "free_x") +
  coord_cartesian(ylim = c(.5, 5.5)) +  # This line zooms in on the y-axis without removing data outside this range
  theme_minimal() +
  theme(axis.text.x = element_blank(), # remove x-axis labels
        axis.ticks.x = element_blank()) + # remove x-axis ticks
  labs(title = "Gene Distribution Before and After Batch Effect Correction",
       y = "Expression Value",
       x = "Sample") +
  scale_fill_manual(values = c("GSE35488" = "brown2", "GSE32591" = "olivedrab3", "GSE66494" = "dodgerblue", "GSE47184" = "magenta2")) +
  scale_color_manual(values = c("GSE35488" = "brown2", "GSE32591" = "olivedrab3", "GSE66494" = "dodgerblue", "GSE47184" = "magenta2"))

#----------------------------------------------------

# First, calculate variance for each gene in the corrected data
gene_variance <- apply(corrected_data, 1, var)

# Identify duplicated gene symbols
duplicated_genes <- which(duplicated(rownames(corrected_data)))

# If there are any duplicated genes, handle them
if (length(duplicated_genes) > 0) {
  # Create a data frame with gene symbols and their corresponding variance
  gene_variance_df <- data.frame(gene = rownames(corrected_data), variance = gene_variance, stringsAsFactors = FALSE)
  
  # Order the data frame by variance, so highest variance comes first
  gene_variance_df <- gene_variance_df[order(-gene_variance_df$variance), ]
  
  # Remove duplicated gene symbols, keeping only the one with the highest variance
  unique_genes <- !duplicated(gene_variance_df$gene)
  
  # Subset the corrected data to only include the unique genes
  corrected_data <- corrected_data[gene_variance_df$gene[unique_genes], ]
}

# Function to split data
split_data <- function(data, labels, num_train, num_train_pos, num_train_neg) {
  # Get indices of positive and negative cases
  pos_indices <- which(labels == 1)
  neg_indices <- which(labels == 0)
  
  # Randomly sample specified number of positive and negative cases for the training set
  train_pos_indices <- sample(pos_indices, num_train_pos)
  train_neg_indices <- sample(neg_indices, num_train_neg)
  train_indices <- c(train_pos_indices, train_neg_indices)
  
  # The rest go to the test set
  test_indices <- setdiff(1:length(labels), train_indices)
  
  # Split the data and labels
  train_data <- data[, train_indices]
  test_data <- data[, test_indices]
  train_labels <- labels[train_indices]
  test_labels <- labels[test_indices]
  
  return(list(train_data = train_data, test_data = test_data, 
              train_labels = train_labels, test_labels = test_labels))
}

# Set the number of patients you want for the training dataset
num_train <- 174 # You can adjust this
num_train_pos <- 141 # Half are positive
num_train_neg <- 33 # Half are negative

# Split the data
split_result <- split_data(corrected_data, all_labels_shuffled, num_train, num_train_pos, num_train_neg)

train_data <- split_result$train_data
test_data <- split_result$test_data
train_labels <- split_result$train_labels
test_labels <- split_result$test_labels

# Append labels to patient names (column names)
colnames(train_data) <- paste0(colnames(train_data), ifelse(train_labels == 1, "_1", "_0"))
colnames(test_data) <- paste0(colnames(test_data), ifelse(test_labels == 1, "_1", "_0"))

# Shuffle the datasets
set.seed(123) # for reproducibility
# Shuffle the training dataset and its labels
train_shuffled_indices <- sample(length(train_labels))
train_data_shuffled <- train_data[, train_shuffled_indices]
train_labels_shuffled <- train_labels[train_shuffled_indices]

# Shuffle the test dataset and its labels
test_shuffled_indices <- sample(length(test_labels))
test_data_shuffled <- test_data[, test_shuffled_indices]
test_labels_shuffled <- test_labels[test_shuffled_indices]

corrected_train_gene_symbols <- train_data_shuffled
corrected_test_gene_symbols <- test_data_shuffled
status_Trainshuffled <- train_labels_shuffled
status_Testshuffled <- test_labels_shuffled
corrected_train_gene_symbols <- as.matrix(corrected_train_gene_symbols)
corrected_test_gene_symbols <- as.matrix(corrected_test_gene_symbols)

#----------------------------------------------------

# Create a design matrix for the linear model
design_matrix <- model.matrix(~ status_Trainshuffled)

# Perform differential expression analysis using limma
fit <- lmFit(corrected_train_gene_symbols, design_matrix)
fit <- eBayes(fit)

# Criteria for important DEGs
DEG_p_value <- 0.05
DEG_coefficient <- .75  # Adjust this value as per your requirement

# Create a data frame for the volcano plot
volcano_data <- data.frame(
  Gene = rownames(fit$coefficients),
  Coefficient = fit$coefficients[, "status_Trainshuffled"],
  PValue = fit$p.value[, "status_Trainshuffled"]
)

# Color-coding based on significance criteria
volcano_data$Color <- ifelse(volcano_data$PValue < DEG_p_value & abs(volcano_data$Coefficient) > DEG_coefficient, 
                             ifelse(volcano_data$Coefficient > 0, "green", "red"), 
                             "grey")

# Create the volcano plot
volcano_plot <- ggplot(volcano_data, aes(x = Coefficient, y = -log10(PValue), color = Color)) +
  geom_point() +
  scale_color_identity() +
  labs(x = "Log2 Fold Change", y = "-log10(FDR)") +
  theme_minimal() +
  # Adding dotted lines for the thresholds
  geom_hline(yintercept = -log10(DEG_p_value), linetype = "dotted", color = "black") + # for the p-value threshold
  geom_vline(xintercept = DEG_coefficient, linetype = "dotted", color = "black") +     # for the positive log2 fold change threshold
  geom_vline(xintercept = -DEG_coefficient, linetype = "dotted", color = "black")      # for the negative log2 fold change threshold


# Labeling genes that are significant based on p-value and coefficient
significant_genes <- subset(volcano_data, PValue < DEG_p_value & (Coefficient > DEG_coefficient | Coefficient < -DEG_coefficient))
print(volcano_plot + geom_text(data = significant_genes, aes(label = Gene), vjust = 1, hjust = 1, size = 3))

# Extract genes that meet the criteria from the volcano_data
selected_genes <- volcano_data$Gene[volcano_data$PValue < DEG_p_value & abs(volcano_data$Coefficient) > DEG_coefficient]

# Use these selected genes to subset the corrected_train_gene_symbols data
regulated_genes <- corrected_train_gene_symbols[rownames(corrected_train_gene_symbols) %in% selected_genes, ]

# Create a heatmap
'pheatmap(regulated_gene_expression,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         clustering_distance_cols = "euclidean",
         clustering_distance_rows = "euclidean",
         main = "Heatmap of Up/Down-Regulated Genes",
         fontsize = 8)'
#----------------------------------------------------

# Extract the patient type from column names
patient_type <- status_Trainshuffled

# Create the HeatmapAnnotation for patient type
ha <- HeatmapAnnotation(
  df = data.frame(type = as.factor(patient_type)),
  annotation_height = unit(4, "mm"),
  annotation_name_gp = gpar(fontsize = 7)
)

# Draw the heatmap with annotation
Heatmap(
  regulated_genes,
  name = "expression",
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_names_gp = gpar(fontsize = 7),       # Set the font size for column names
  row_names_gp = gpar(fontsize = 7),       # Set the font size for row names
  top_annotation = ha,
  cluster_columns = FALSE, # We're not clustering columns; instead, we order them by patient_type
  column_order = order(patient_type) # Order columns based on patient_type
)

#----------------------------------------------------

# Extracting DEGs based on coefficient value
DEGs <- volcano_data$Gene[volcano_data$PValue < DEG_p_value & abs(volcano_data$Coefficient) > DEG_coefficient]

# Transpose the datasets
transposed_train_data <- t(corrected_train_gene_symbols)
transposed_test_data <- t(corrected_test_gene_symbols)

# Filter to DEGs
train_subset <- transposed_train_data[, DEGs]
test_subset <- transposed_test_data[, DEGs]

# Check for any missing genes in the test dataset and fill them with median values from the training data
missing_genes <- setdiff(DEGs, colnames(test_subset))
for (gene in missing_genes) {
  test_subset[, gene] <- median(train_subset[, gene], na.rm = TRUE)
}
#----------------------------------------------------

# Enrichment Analysis

# For GO and KEGG analysis:
ego <- enrichGO(gene = DEGs,
                OrgDb = org.Hs.eg.db,
                keyType = "SYMBOL",
                ont = "BP",
                pAdjustMethod = "BH",
                qvalueCutoff = 0.05)

# Convert gene symbols to Entrez IDs
geneList_entrez <- bitr(DEGs, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
ekk <- enrichKEGG(gene = geneList_entrez$ENTREZID,
                  organism = 'hsa',
                  keyType = 'kegg',
                  pAdjustMethod = 'BH',
                  qvalueCutoff = 0.05)

# For DO analysis:
edo <- enrichDO(gene = geneList_entrez$ENTREZID, pAdjustMethod = "BH", qvalueCutoff = 0.05)

# Function to generate a formatted plot for enrichment results
generate_formatted_plot <- function(enrichment_result, enrichment_type) {
  
  # Extract necessary data
  df <- as.data.frame(enrichment_result)
  df <- df[order(df$qvalue), ]  # Sorting based on adjusted p-values
  
  # Determine the title based on enrichment type
  title_text <- switch(enrichment_type,
                       GO = "Gene Ontology Enrichment",
                       KEGG = "KEGG Pathway Enrichment",
                       DO = "Disease Ontology Enrichment",
                       "Enrichment Analysis")  # default title
  
  # Generate the ggplot
  p <- ggplot(df, aes(x = GeneRatio, y = Description, size = Count, color = qvalue)) +
    geom_point(alpha = 0.8) + 
    scale_size_continuous(range = c(3, 15)) +  # Adjust the range for dot sizes as needed
    scale_color_gradient(low = "red", high = "blue", limits = c(0, 0.04), 
                         breaks = c(0.01, 0.02, 0.03, 0.04), name = "q value") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "right") +
    labs(title = title_text, x = "Gene ratio", y = NULL)
  
  print(p)
}

# Visualizing GO, KEGG, and DO enrichment results using the formatted plot:
generate_formatted_plot(ego, "GO")  # GO terms
#generate_formatted_plot(ekk, "KEGG")  # KEGG pathways
# Extract the data from the 'ekk' object
ekk_data <- ekk@result

# Filter the data to only include pathways with Count > 2
filtered_ekk_data <- subset(ekk_data, Count > 1)

ekk_plot <- ggplot(filtered_ekk_data, aes_string(x="GeneRatio", y="Description")) + 
  geom_point(aes_string(size="Count", colour="p.adjust")) + 
  scale_colour_gradient(low = "red", high = "blue") +  # Adjust color scheme to red-blue gradient
  labs(title = "KEGG Enrichment Analysis", x = "Gene ratio", y = "Pathway/Disease")

# Display the plot
print(ekk_plot)

# Visualizing DO enrichment results using circular plots:
generate_formatted_plot(edo, "DO")  # DO terms


#----------------------------------------------------

# Lasso Regression
x <- as.matrix(train_subset)  # Predictors 
y <- as.factor(status_Trainshuffled)  # Response variable

# Lasso with cross-validation
lasso_fit <- glmnet(x, y, family = "binomial", alpha = 1)
cv.lasso <- cv.glmnet(x, y, family = "binomial", alpha = 1)
best_lambda <- cv.lasso$lambda.min
coef_lasso <- predict(lasso_fit, type = "coefficients", s = best_lambda)

min_coef <- min(cv.lasso$glmnet.fit$beta)
max_coef <- max(cv.lasso$glmnet.fit$beta)

# Plotting Lasso Regression curve for DEGs
plot(cv.lasso$glmnet.fit, xvar="lambda", label=TRUE, ylim=c(min_coef, max_coef))
#title(paste("Lasso Regression Curve of", length(DEGs), "DEGs"))

#----------------------------------------------------

#Isolating top genes
plot(cv.lasso)

# Extracting the optimal number of genes from the cv.lasso object
optimal_genes_count <- cv.lasso$glmnet.fit$df[which.min(cv.lasso$cvm)]

# Convert the coefficients into a data frame
coef_df <- as.data.frame(as.matrix(coef_lasso))

# Remove the (Intercept) row
coef_df <- coef_df[-which(rownames(coef_df) == "(Intercept)"), , drop = FALSE]

# Ensure you only have one column now which represents the coefficients for each gene
gene_coef <- coef_df[[1]]

# Sort the genes by absolute coefficient value
sorted_genes <- coef_df %>%
  dplyr::arrange(desc(abs(gene_coef))) %>%
  dplyr::slice(1:optimal_genes_count)  # Select top genes based on cv.lasso

# Extracting genes
top_genes <- rownames(sorted_genes)

# Subset the original dataset for the top genes
top_genes_subset <- train_subset[, top_genes]

# Check the dimensions to ensure we have only optimal # of genes
#dim(top_genes_subset)

#----------------------------------------------------

# Using RF with RFE for feature selection

# Set up repeated 10-fold CV
train_control <- trainControl(method = "repeatedcv", number = 10, repeats = 5)

# Control for RFE
ctrl <- rfeControl(functions=rfFuncs, method=train_control$method)

# Perform RFE
results_rfe <- rfe(top_genes_subset, status_Trainshuffled, sizes=c(1:length(DEGs)), rfeControl=ctrl)
#results_rfe <- rfe(train_subset, status_Trainshuffled, sizes=c(1:length(DEGs)), rfeControl=ctrl)


# Extracting the number of genes and RMSE
num_genes <- results_rfe$results$Variables
rmse_values <- results_rfe$results$RMSE

# Assuming binary classification, calculating accuracy
accuracy <- 1 - rmse_values

# Identify the position of the highest accuracy
max_accuracy_idx <- which.max(accuracy)

# Define a cyan/medium blue color
cyan_blue <- rgb(0, 0.6, 1)

# Plotting the 10-fold CV of RF-RFE with hollow dots and blue connecting lines
plot(num_genes, accuracy, type="c", xlab="Number of Variables (Genes)", 
     ylab="Accuracy (Cross-Validation)", col=cyan_blue) #10-fold CV for Signature Gene Combination of RF-RFE
points(num_genes, accuracy, pch=21, col=cyan_blue, lwd=1)

# Highlight the point with the highest accuracy
points(num_genes[max_accuracy_idx], accuracy[max_accuracy_idx], col=cyan_blue, pch=19, lwd=1.5)

# Indicate the number of variables used at the point with the highest accuracy
text(num_genes[max_accuracy_idx], accuracy[max_accuracy_idx], 
     labels=paste(num_genes[max_accuracy_idx], "genes"), pos=1, offset=0.5, col=cyan_blue)

'pos = 1: Below
pos = 2: Left
pos = 3: Above
pos = 4: Right'

# Extract optimal number of genes
optimal_num_genes <- num_genes[max_accuracy_idx]

# Sort genes by importance
importance_data <- results_rfe$fit$importance
sorted_genes <- rownames(importance_data[order(importance_data[,"IncNodePurity"], decreasing = TRUE), ])

# Select top x genes
top_genes <- sorted_genes[1:optimal_num_genes]

# Subset data to include only top x genes
selected_train_subset <- train_subset[, top_genes]

# Selecting candidate genes from RF-RFE
selected_genes <- predictors(results_rfe)

# Train a RandomForest model using these genes
set.seed(123)
rf_model_selected <- randomForest(selected_train_subset, status_Trainshuffled)

# Plotting the importance scores for the selected subset in random forests
gene_importance_selected <- importance(rf_model_selected)[,1]
sorted_importance_selected <- gene_importance_selected[order(-gene_importance_selected)]

# Rescale importance scores to range [0, 100]
min_score <- min(gene_importance_selected)
max_score <- max(gene_importance_selected)
scaled_importance_selected <- (gene_importance_selected - min_score) * (100 - 0) / (max_score - min_score) + 0

# Sort the scaled importance scores
sorted_scaled_importance_selected <- scaled_importance_selected[order(-scaled_importance_selected)]

# Plotting the scaled importance scores for the selected subset in random forests
barplot(sorted_scaled_importance_selected, las=2, cex.names=0.7, ylim=c(0,100), ylab="Importance")#Scaled Gene Importance Scores for Selected Subset in Random Forests

importantgenes <- names(scaled_importance_selected)[scaled_importance_selected >25]
print(importantgenes)

#----------------------------------------------------

#Incorporate important genes
# Define the genes of interest
genes_of_interest <- importantgenes

# Subset the data to only include these genes
selected_genes_subset <- train_subset[, colnames(train_subset) %in% genes_of_interest]

#Grid search for the mtry parameter and calculate accuracy using repeated 5 times 10-fold cross-validation:
# Set up repeated 10-fold CV
train_control <- trainControl(method = "repeatedcv", number = 10, repeats = 5)

# Define the search grid
tuneGrid <- expand.grid(.mtry=1:length(genes_of_interest))

# Perform the grid search
set.seed(123)
rf_gridsearch <- train(selected_genes_subset, status_Trainshuffled, method="rf", tuneGrid=tuneGrid, trControl=train_control)

# Extract the best mtry
best_mtry <- rf_gridsearch$bestTune$mtry
print(best_mtry)
#----------------------------------------------------

# Cross-Validation on Training Data Only
num_folds <- 5
fold_size <- floor(nrow(selected_genes_subset) / num_folds)
all_fold_accuracies <- numeric(num_folds)
all_gene_importances <- list()
all_fold_roc <- list()  # Store ROC curves here

# Cross-Validation on Training Data Only
par(mfrow=c(1,5)) # Arrange plots in a 1x5 grid

calculate_roc <- function(predictions, true_labels) {
  thresholds <- predictions
  tpr <- numeric(length(thresholds))
  fpr <- numeric(length(thresholds))
  
  cat("\n----- Debugging ROC Calculation -----\n")
  cat("Number of unique thresholds:", length(thresholds), "\n")
  
  for (i in 1:length(thresholds)) {
    predicted_classes <- ifelse(predictions >= thresholds[i], 1, 0)
    true_positives <- sum(predicted_classes == 1 & true_labels == 1)
    false_positives <- sum(predicted_classes == 1 & true_labels == 0)
    true_negatives <- sum(predicted_classes == 0 & true_labels == 0)
    false_negatives <- sum(predicted_classes == 0 & true_labels == 1)
    
    tpr[i] <- true_positives / (true_positives + false_negatives)
    fpr[i] <- false_positives / (false_positives + true_negatives)
    
    # Print out the intermediate calculations for debugging
    cat("\nThreshold:", thresholds[i], "\n")
    cat("True Positives:", true_positives, "\n")
    cat("False Positives:", false_positives, "\n")
    cat("True Negatives:", true_negatives, "\n")
    cat("False Negatives:", false_negatives, "\n")
    cat("Sensitivity (TPR):", tpr[i], "\n")
    cat("1-Specificity (FPR):", fpr[i], "\n")
  }
  
  # Order the ROC curve by increasing values of FPR
  ordered_indices <- order(fpr, tpr)
  fpr <- fpr[ordered_indices]
  tpr <- tpr[ordered_indices]
  
  list(fpr=fpr, tpr=tpr)
}

# Function to compute AUC from the manually calculated ROC data
calculate_auc <- function(roc_data) {
  # Trapezoidal rule for integration
  auc = sum(diff(roc_data$fpr) * (roc_data$tpr[-1] + roc_data$tpr[-length(roc_data$tpr)]) / 2)
  return(auc)
}

for (k in 1:num_folds) {
  cv_indices <- ((k-1)*fold_size + 1):(k*fold_size)
  
  # Shuffle the test set
  shuffled_test_indices <- sample(cv_indices)
  
  cv_test <- selected_genes_subset[shuffled_test_indices, ]
  cv_test_labels <- status_Trainshuffled[shuffled_test_indices]
  
  cv_train <- selected_genes_subset[-cv_indices, ]
  cv_train_labels <- status_Trainshuffled[-cv_indices]
  
  # BEGIN UNDERSAMPLING
  # Identify indices of class 1 and 0
  indices_class1 <- which(cv_train_labels == 1)
  indices_class0 <- which(cv_train_labels == 0)
  
  # Count of the minority class
  min_class_count <- min(length(indices_class1), length(indices_class0))
  
  # Sample from the majority class to match the count of the minority class
  undersample_indices_class1 <- sample(indices_class1, min_class_count)
  
  # Combine the undersampled majority class indices with minority class indices
  balanced_indices <- c(undersample_indices_class1, indices_class0)
  
  # Extract balanced data
  cv_train <- cv_train[balanced_indices, ]
  cv_train_labels <- cv_train_labels[balanced_indices]
  # END UNDERSAMPLING
  
  cv_rf <- randomForest(x = cv_train, y = cv_train_labels, mtry = best_mtry, importance = TRUE)
  
  # Predict the probabilities, not the class labels
  cv_predictions <- predict(cv_rf, newdata = cv_test, type = "response")
  
  # Calculate ROC curve for this fold
  roc_data <- calculate_roc(cv_predictions, cv_test_labels)
  
  # Calculate the AUC from the manually calculated ROC data
  auc_val <- calculate_auc(roc_data)
  
  # Use the pROC library to get the CI for the AUC
  roc_obj <- roc(cv_test_labels, cv_predictions)
  ci_vals <- ci.auc(roc_obj, conf.level = 0.95)
  
  # Manually plot the ROC curve with Sensitivity against 1-Specificity
  plot(roc_data$fpr, roc_data$tpr, type="l", col="red", 
       xlab="False Positive Rate (1-Specificity)", 
       ylab="True Positive Rate (Sensitivity)", 
       main = paste("ROC Curve for Fold", k))
  abline(0, 1, lty=2, col="gray")  # Adds a diagonal line for reference
  
  # Display the AUC and its CI on the plot
  text(x = 0.5, y = 0.5, labels = paste("AUC for Fold", k, ": ", round(auc_val, 6), "\n95% CI: [", round(ci_vals[1], 6), ", ", round(ci_vals[2], 6), "]"), adj = c(0.5,0.5))
  
  
  '# Plot ROC curve for this fold
  plot(roc_obj, main = paste("ROC Curve for Fold", k), col = "red")
  s <- ci.se(roc_obj)
  lines(s, type = "b", pch = 19, lty = 2, col = "red")
  text(x = 0.5, y = 0.5, labels = paste("AUC for Fold", k, ": ", round(roc_obj$auc, 6), "\n95% CI: [", round(ci_vals[2], 6), ", ", round(ci_vals[4], 6), "]")) # Moved AUC value to the center
  '
  
  # Calculate accuracy using confusionMatrix
  predicted_classes <- ifelse(cv_predictions >= 0.5, 1, 0)
  cm <- confusionMatrix(as.factor(predicted_classes), as.factor(cv_test_labels))
  fold_accuracy <- cm$overall["Accuracy"]
  
  # Store the ROC, accuracy, and gene importances for this fold
  all_fold_roc[[k]] <- roc_data
  all_fold_accuracies[k] <- fold_accuracy
  all_gene_importances[[k]] <- importance(cv_rf)
  
  # Print fold results
  cat("Fold", k, "\n")
  cat("Predictions: ", cv_predictions, "\n")
  cat("Actual: ", cv_test_labels, "\n")
  cat("Accuracy for Fold", k, ": ", fold_accuracy * 100, "%\n\n")
}

# Print accuracy for each fold
for (k in 1:num_folds) {
  cat("Accuracy for Fold", k, ": ", all_fold_accuracies[k] * 100, "%\n")
}
par(mfrow=c(1,1)) # Arrange plots in a 1x1 grid

#----------------------------------------------------

rf_model_optimal <- randomForest(x = selected_genes_subset, y = status_Trainshuffled, mtry = best_mtry, importance = TRUE)

# Predict on the Training Dataset:
predictions <- predict(rf_model_optimal, newdata = train_subset, type = "response")

# Calculate ROC, AUC and its 95% CI:
roc_obj <- roc(status_Trainshuffled, predictions)
auc_val <- auc(roc_obj)
ci_vals <- ci(roc_obj, method = "delong")

# Plot the ROC curve, display AUC and its CI:

plot(roc_obj, col = "red", main = "ROC Curve for Training Dataset")
text(x = 0.5, y = 0.5, labels = paste("AUC: ", round(auc_val, 6), 
                                      "\n95% CI: (", round(ci_vals[1], 6), 
                                      "-", round(ci_vals[3], 6), ")", sep=""), col = "red")

# Calculate accuracy on the Training Dataset using confusionMatrix
predicted_classes <- ifelse(predictions >= 0.5, 1, 0)
cm <- confusionMatrix(as.factor(predicted_classes), as.factor(status_Trainshuffled))
accuracy <- cm$overall["Accuracy"]
cat("Accuracy for the Training Dataset: ", accuracy * 100, "%\n")

#----------------------------------------------------

test_predictions_prob <- predict(rf_model_optimal, newdata = test_subset, type = "response")

# Calculate ROC, AUC and its 95% CI:
roc_obj <- roc(status_Testshuffled, test_predictions_prob)
auc_val <- auc(roc_obj)
ci_vals <- ci(roc_obj, method = "delong")

# Plot the ROC curve, display AUC and its CI:
plot(roc_obj, col = "blue", main = "ROC Curve for Testing Dataset")
text(x = 0.5, y = 0.5, labels = paste("AUC: ", round(auc_val, 6), 
                                      "\n95% CI: (", round(ci_vals[1], 6), 
                                      "-", round(ci_vals[3], 6), ")", sep=""), col = "blue")

# Calculate accuracy on the Testing Dataset using confusionMatrix
predicted_classes <- ifelse(test_predictions_prob >= 0.5, 1, 0)
cm <- confusionMatrix(as.factor(predicted_classes), as.factor(status_Testshuffled))
accuracy <- cm$overall["Accuracy"]
cat("Accuracy for the Testing Dataset: ", accuracy * 100, "%\n")

#----------------------------------------------------

# Create two subsets:
# 1. Only the gene DUSP1
only_dusp1 <- selected_genes_subset[, colnames(selected_genes_subset) == "DUSP1"]

# 2. All genes EXCEPT DUSP1
without_dusp1 <- selected_genes_subset[, !(colnames(selected_genes_subset) %in% c("DUSP1"))]

# Use the existing grid search for mtry and training setup:

train_control <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
tuneGrid_dusp1 <- expand.grid(.mtry=1)
tuneGrid_without_dusp1 <- expand.grid(.mtry=1:ncol(without_dusp1))
only_dusp1_df <- data.frame(DUSP1 = only_dusp1)

set.seed(123)
#rf_gridsearch_dusp1 <- train(only_dusp1_df, status_Trainshuffled, method="rf", tuneGrid=tuneGrid_dusp1, trControl=train_control)
#best_mtry_dusp1 <- rf_gridsearch_dusp1$bestTune$mtry

#rf_gridsearch_without_dusp1 <- train(without_dusp1, status_Trainshuffled, method="rf", tuneGrid=tuneGrid_without_dusp1, trControl=train_control)
#best_mtry_without_dusp1 <- rf_gridsearch_without_dusp1$bestTune$mtry

# Next, train the optimal RF on the entire dataset for both data subsets:
only_dusp1 <- as.data.frame(only_dusp1)

rf_model_optimal_dusp1 <- randomForest(x = only_dusp1, y = status_Trainshuffled, importance = TRUE)
rf_model_optimal_without_dusp1 <- randomForest(x = without_dusp1, y = status_Trainshuffled, importance = TRUE)

'# Predict on the Training Dataset for both:
predictions_dusp1 <- predict(rf_model_optimal_dusp1, newdata = only_dusp1, type = "response")
predictions_without_dusp1 <- predict(rf_model_optimal_without_dusp1, newdata = without_dusp1, type = "response")'

# Now, predict on the Testing Dataset:
test_only_dusp1 <- test_subset[, colnames(test_subset) == "DUSP1"]
test_only_dusp1 <- as.data.frame(test_only_dusp1)

test_without_dusp1 <- test_subset[, !(colnames(test_subset) %in% c("DUSP1"))]
test_without_dusp1 <- as.data.frame(test_without_dusp1)

colnames(test_only_dusp1) <- colnames(only_dusp1)

test_predictions_prob_dusp1 <- predict(rf_model_optimal_dusp1, newdata = test_only_dusp1, type = "response")
test_predictions_prob_without_dusp1 <- predict(rf_model_optimal_without_dusp1, newdata = test_without_dusp1, type = "response")

# Calculate ROC for both scenarios
roc_obj_dusp1 <- roc(status_Testshuffled, test_predictions_prob_dusp1)
roc_obj_without_dusp1 <- roc(status_Testshuffled, test_predictions_prob_without_dusp1)

# Plot the ROCs
plot(roc_obj_dusp1, main="ROC Curves for RF Models", col="darkcyan", lwd=2)
lines(roc_obj_without_dusp1, col="lightblue", lwd=2)
legend("bottomright", inset=c(0.05,0.1), legend=c("Only DUSP1", "Without DUSP1"), col=c("darkcyan", "lightblue"), lwd=2)

# Calculate AUC values
auc_dusp1 <- auc(roc_obj_dusp1)
auc_without_dusp1 <- auc(roc_obj_without_dusp1)

# Add AUC values to the plot
usr <- par("usr")
x_pos <- usr[2] - 0.1*(usr[2]-usr[1])  # 10% from the right edge of the plot
y_pos_dusp1 <- usr[3] + 0.6*(usr[4]-usr[3])  # 90% from the bottom of the plot for GADD45B
y_pos_without_dusp1 <- y_pos_dusp1 - 0.05*(usr[4]-usr[3])  # A bit below the first AUC value for Without GADD45B

text(x=x_pos, y=y_pos_dusp1, labels=sprintf("AUC (Only DUSP1) = %.2f", auc_dusp1), col="darkcyan", adj=c(1,0))
text(x=x_pos, y=y_pos_without_dusp1, labels=sprintf("AUC (Without DUSP1) = %.2f", auc_without_dusp1), col="lightblue", adj=c(1,0))

#----------------------------------------------------
library(hugene21sttranscriptcluster.db)

#External validation of the random forest model

#GSE180394

# Create an empty list to store the datasets, expression data, and labels
datasets <- list()

# Fetch the dataset from GEO
dataset <- getGEO("GSE180394", GSEMatrix = TRUE)

# Extract the expression data table
expression_dataGSE6280 <- dataset$GSE180394_series_matrix.txt.gz@assayData$exprs

#Labels for GSE62792
# Extract the source name column
source_name <- dataset$GSE180394_series_matrix.txt.gz@phenoData@data$characteristics_ch1

# Initialize an empty vector to store the labels
labelsGSE6280 <- numeric(length(source_name))

# Assign labels based on source name
for (i in 1:length(source_name)) {
  if (source_name[i] == "sample group: 2' FSGS" ||
      source_name[i] == "sample group: chronic Glomerulonephritis (GN) with infiltration by CLL" ||
      source_name[i] == "sample group: DN" ||
      source_name[i] == "sample group: FGGS" ||
      source_name[i] == "sample group: FSGS" ||
      source_name[i] == "sample group: Hydronephrosis" ||
      source_name[i] == "sample group: IgAN" ||
      source_name[i] == "sample group: Interstitial nephritis" ||
      source_name[i] == "sample group: Hypertensive Nephrosclerosis" ||
      source_name[i] == "sample group: Light-Chain Deposit Disease (IgG lambda)" ||
      source_name[i] == "sample group: LN-WHO III" ||
      source_name[i] == "sample group: LN-WHO III+V" ||
      source_name[i] == "sample group: LN-WHO IV" ||
      source_name[i] == "sample group: LN-WHO IV+V" ||
      source_name[i] == "sample group: LN-WHO V" ||
      source_name[i] == "sample group: LN-WHO-I/II" ||
      source_name[i] == "sample group: MCD" ||
      source_name[i] == "sample group: MN" ||
      source_name[i] == "sample group: CKD with mod-severe Interstitial fibrosis" ||
      source_name[i] == "sample group: Hypertensive Nephrosclerosis" ||
      source_name[i] == "sample group: Thin-BMD" ||
      source_name[i] == "sample group: Unaffected parts of Tumor Nephrectomy") {
    labelsGSE6280[i] <- 1  # Positive for CKD
  } else if (source_name[i] == "sample group: Living donor" ) {
    labelsGSE6280[i] <- 0  # Control
  } else {
    labelsGSE6280[i] <- NA  # Unknown label
  }
}

# Fetch the annotation dataset from GEO
annotation_dataset <- getGEO("GPL19983", GSEMatrix = TRUE)
# Extract the mapping from probe ID to ENTREZ_GENE_ID
probe_to_entrez <- annotation_dataset@dataTable@table$ENTREZ_GENE_ID
names(probe_to_entrez) <- annotation_dataset@dataTable@table$ID

# Map your expression data's probe IDs to Entrez IDs
entrez_ids <- probe_to_entrez[rownames(expression_dataGSE6280)]
valid_indices <- which(!is.na(entrez_ids))
expression_dataGSE6280 <- expression_dataGSE6280[valid_indices, ]
entrez_ids <- entrez_ids[valid_indices]
entrez_ids <- as.character(entrez_ids)

# Convert Entrez IDs to gene symbols using org.Hs.eg.db
gene_symbols <- mapIds(org.Hs.eg.db,
                       keys = entrez_ids,
                       keytype = "ENTREZID",
                       column = "SYMBOL",
                       multiVals = "first")

# Update the row names of the expression matrix to the corresponding gene symbols
rownames(expression_dataGSE6280) <- gene_symbols

# 1. Check for duplicated genes and remove them
duplicated_genes <- which(duplicated(rownames(expression_dataGSE6280)))
expression_dataGSE6280 <- expression_dataGSE6280[-duplicated_genes, ]

# 2. Combine all_expression_data and expression_dataGSE62792 (using common genes)
common_genes <- intersect(rownames(all_expression_data), rownames(expression_dataGSE6280))
all_expression_subset <- all_expression_data[common_genes, ] # Extract the expression data only for the common genes
expression_dataGSE6280_subset <- expression_dataGSE6280[common_genes, ]
combined_expression <- cbind(all_expression_subset, expression_dataGSE6280_subset) # Combine the datasets column-wise


# 3. Combine all_labels and labelsGSE62792
combined_labels <- c(all_labels, labelsGSE6280)

# Before splitting, shuffle the dataset and the labels
set.seed(123)  # setting a seed ensures reproducibility
shuffle_order <- sample(1:ncol(combined_expression))
combined_expression <- combined_expression[, shuffle_order]
combined_labels <- combined_labels[shuffle_order]

# 4. Append the labels to the end of the patient names (columns) "_1" for positive" and "_0" for control
colnames(combined_expression) <- paste0(colnames(combined_expression), "_", combined_labels)

# 5. Normalize the dataset (while retaining the rows and column names)
normalized_expression <- normalizeBetweenArrays(combined_expression, method="quantile")

# 6. Log2 transform the data
log2_transformed_expression <- log2(normalized_expression + 1) # adding 1 to avoid log(0)

# 7. Split the data into two batches
batch <- c(rep(1, ncol(all_expression_data)), rep(2, ncol(expression_dataGSE6280_subset)))

# 8. Perform ComBat to correct for batch effects
combat_corrected_expression <- ComBat(dat=log2_transformed_expression, batch=batch, mod=NULL, par.prior=TRUE, prior.plots=FALSE)

# Patient names from the external dataset
external_patients <- c(
  "GSM5607814", "GSM5607815", "GSM5607816", "GSM5607817", "GSM5607818", "GSM5607819", "GSM5607820", "GSM5607821", 
  "GSM5607822", "GSM5607823", "GSM5607824", "GSM5607825", "GSM5607826", "GSM5607827", "GSM5607828", "GSM5607829", 
  "GSM5607830", "GSM5607831", "GSM5607832", "GSM5607833", "GSM5607834", "GSM5607835", "GSM5607836", "GSM5607837", 
  "GSM5607838", "GSM5607839", "GSM5607840", "GSM5607841", "GSM5607842", "GSM5607843", "GSM5607844", "GSM5607845", 
  "GSM5607846", "GSM5607847", "GSM5607848", "GSM5607849", "GSM5607850", "GSM5607851", "GSM5607852", "GSM5607853", 
  "GSM5607854", "GSM5607855", "GSM5607856", "GSM5607857", "GSM5607858", "GSM5607859", "GSM5607860", "GSM5607861", 
  "GSM5607862", "GSM5607863", "GSM5607864", "GSM5607865", "GSM5607866", "GSM5607867", "GSM5607868", "GSM5607869", 
  "GSM5607870", "GSM5607871", "GSM5607872"
)

# Identify columns that match the external patient names
matching_cols <- grep(paste(external_patients, collapse="|"), colnames(combat_corrected_expression))
#matching_cols <- grep(paste(external_patients, collapse="|"), colnames(log2_transformed_expression))


# Extract only the external dataset using the identified columns
external_data <- combat_corrected_expression[, matching_cols]
#external_data <- log2_transformed_expression[, matching_cols]
external_labels <- as.numeric(gsub(".*_", "", colnames(external_data)))

#----------------------------------------------------

# 1. Melt the dataframes for ggplot

# Before ComBat Correction
melted_data_before <- melt(log2_transformed_expression)
melted_data_before$Correction <- 'After'

# After ComBat Correction
melted_data_after <- melt(combat_corrected_expression)
melted_data_after$Correction <- 'Before'

# 2. Combine melted dataframes
all_melted_data <- rbind(melted_data_before, melted_data_after)

# 3. Assign dataset identities based on column names of expression data
# Assuming you have two datasets: all_expression_data and expression_dataGSE6280_subset
dataset_assignments <- c(rep("Combined Dataset", ncol(all_expression_subset)), rep("External Dataset (GSE180394)", ncol(expression_dataGSE6280_subset)))

names(dataset_assignments) <- colnames(combined_expression)

# 4. Add dataset information to the melted data
all_melted_data$Dataset <- dataset_assignments[all_melted_data$Var2]

# 5. Create the box plot with the dataset as both the fill and color variable
ggplot(all_melted_data, aes(x = Var2, y = value, fill = Dataset, color = Dataset)) +
  geom_boxplot(outlier.shape = NA) + # hide outliers for cleaner look
  facet_wrap(~Correction, ncol = 1, scales = "free_x") +
  coord_cartesian(ylim = c(1, 5.25)) +  # This line zooms in on the y-axis without removing data outside this range
  theme_minimal() +
  theme(axis.text.x = element_blank(), # remove x-axis labels
        axis.ticks.x = element_blank()) + # remove x-axis ticks
  labs(title = "Gene Distribution Before and After Batch Effect Correction",
       y = "Expression Value",
       x = "Sample") +
  scale_fill_manual(values = c("Combined Dataset" = "dodgerblue", "External Dataset (GSE180394)" = "olivedrab3")) + 
  scale_color_manual(values = c("Combined Dataset" = "dodgerblue", "External Dataset (GSE180394)" = "olivedrab3"))   


#----------------------------------------------------

# 1.1 Incorporate important genes
selected_genes_subset_external <- external_data[importantgenes, ]
selected_genes_subset_external <- t(selected_genes_subset_external)

# 2. Grid search for the mtry parameter using repeated 10-fold CV
train_control <- trainControl(method = "repeatedcv", number = 10, repeats = 5)

# Define the search grid
tuneGrid <- expand.grid(.mtry = 1:length(importantgenes))

# Perform the grid search for external data
set.seed(123)
rf_gridsearch_external <- train(selected_genes_subset_external, external_labels, method = "rf", tuneGrid = tuneGrid, trControl = train_control)

# Extract the best mtry
best_mtry <- rf_gridsearch_external$bestTune$mtry

#----------------------------------------------------
#Shuffle both the data and labels
set.seed(47)
shuffle_order <- sample(1:ncol(external_data))
external_data <- external_data[, shuffle_order]
external_labels <- external_labels[shuffle_order]

# Identify the common genes between external_data_post_combat and selected_genes_subset
common_genes <- intersect(rownames(external_data), colnames(selected_genes_subset))

# Subset the external dataset to only those genes
external_subset <- external_data[common_genes, ]

# Transpose the external_subset so genes are now columns
external_subset <- t(external_subset)


# Set up training control for 10-fold cross-validation
control <- trainControl(method="cv", number=5)

# Train the model with adjusted class weights
model <- train(selected_genes_subset, 
               status_Trainshuffled, 
               method="rf", 
               trControl=control,
               importance=TRUE,
               tuneGrid=data.frame(.mtry=seq(2, sqrt(ncol(selected_genes_subset)), by=2))) 


predictions <- predict(model, newdata=external_subset)

# Scale the predictions to the [0,1] range
#scaled_predictions <- (predictions - min(predictions)) / (max(predictions) - min(predictions))

#predicted_class <- ifelse(scaled_predictions > 0.51, 1, 0)
predicted_class <- ifelse(predictions > 0.65, 1, 0)


# Create a data frame to see prediction vs actual label per patient
results <- data.frame(
  Patient = rownames(external_subset),
  Predicted = predictions,
  Class = predicted_class,
  Actual = external_labels
)

print(results)
predicted_class <- factor(predicted_class)
external_labels <- factor(external_labels)

cm <- confusionMatrix(predicted_class, external_labels)
# Extracting the accuracy from the confusion matrix:
accuracy <- cm$overall['Accuracy']
print(paste("Accuracy:", round(accuracy * 100, 2), "%"))

# Extracting the Sensitivity (True Positive Rate) from the confusion matrix:
sensitivity <- cm$byClass['Sensitivity']
print(paste("Sensitivity:", round(sensitivity * 100, 2), "%"))

#roc_external <- roc(external_labels, scaled_predictions)
roc_external <- roc(external_labels, predictions)
auc_external <- auc(roc_external)

# Calculate 95% CI for the AUC
ci <- ci.auc(roc_external)

# Plotting ROC Curve
plot(roc_external, main = "ROC Curve for External Dataset", col = "gold")
# Placing the AUC and CI at the center of the plot
mid_x <- 0.5
mid_y <- 0.5
text(x = mid_x, y = mid_y, labels = paste("AUC: ", round(auc_external, 4)), col = "gold")
text(x = mid_x, y = mid_y - 0.05, labels = paste("95% CI: [", round(ci[1], 4), ", ", round(ci[2], 4), "]", sep = ""), col = "gold")

#----------------------------------------------------