setwd("D:/浙大海宁学习大一/DST2/ICA_mysql/BRCA")
library(dplyr)
library(data.table)
expr_matrix <- fread("D:/浙大海宁学习大一/DST2/ICA_mysql/BRCA/原始数据/GDC-PANCAN.htseq_fpkm.tsv",sep="\t",header = TRUE)
dim(expr_matrix)
colnames(expr_matrix)


gene_ids <- expr_matrix[[1]]
expr_numeric <- as.matrix(expr_matrix[, -1])
rownames(expr_numeric) <- gene_ids

candidate_genes <- fread("D:/浙大海宁学习大一/DST2/ICA_mysql/BRCA/receptor+secreted_using.csv", header = TRUE)
target_gene_ids <- candidate_genes[[1]]

#change
sample_info <- fread("D:/浙大海宁学习大一/DST2/ICA_mysql/BRCA/GDC_pancancer_usingCOAD.csv", header = TRUE)

sample_info$type_code <- as.numeric(substr(sample_info$sample, 14, 15))
tumor_samples <- sample_info$sample[sample_info$type_code >= 1 & sample_info$type_code <= 9]
normal_samples <- sample_info$sample[sample_info$type_code >= 10 & sample_info$type_code <= 19]

# ==================== part 2 ====================
calculate_gene_score <- function(gene_id, expr_mat, tumor_samps, normal_samps) {
  # check1: confirm gene in rownames
  if (!gene_id %in% rownames(expr_mat)) {
    return(data.frame(
      gene_id = gene_id,
      error = "Gene not in matrix rownames",
      stringsAsFactors = FALSE
    ))
  }
  
  # 2nd check: confirm sample exists
  tumor_samps <- intersect(tumor_samps, colnames(expr_mat))
  normal_samps <- intersect(normal_samps, colnames(expr_mat))
  
  if (length(tumor_samps) < 3 | length(normal_samps) < 3) {
    return(data.frame(
      gene_id = gene_id,
      error = "Not enough samples after intersection",
      stringsAsFactors = FALSE
    ))
  }
  tumor_expr <- as.numeric(expr_mat[gene_id, tumor_samps])
  normal_expr <- as.numeric(expr_mat[gene_id, normal_samps])
  
  tumor_expr <- tumor_expr[!is.na(tumor_expr)]
  normal_expr <- normal_expr[!is.na(normal_expr)]
  
  if (length(tumor_expr) < 3 | length(normal_expr) < 3) {
    return(data.frame(
      gene_id = gene_id,
      median_tumor = median(tumor_expr, na.rm=TRUE),#add
      median_normal = median(normal_expr, na.rm=TRUE),#add
      log2fc = NA,
      p_value = NA
    ))
  }
  
  tumor_expr_log <- log2(tumor_expr + 1)
  normal_expr_log <- log2(normal_expr + 1)
  
  median_tumor <- median(tumor_expr_log)
  median_normal <- median(normal_expr_log)
  log2fc <- median_tumor - median_normal#改
  
  p_val <- wilcox.test(tumor_expr, normal_expr)$p.value
  #evidence_score <- -log10(p_val + 1e-10) * abs(log2fc)
  
  return(data.frame(
    gene_id = gene_id,
    median_tumor = median_tumor,
    median_normal = median_normal,
    log2fc = log2fc,
    p_value = p_val
  ))
}

# ==================== Part 3 =================
rownames(expr_numeric) <- sub("\\..*", "", rownames(expr_numeric))
target_gene_ids <- sub("\\..*", "", target_gene_ids) 

valid_genes <- intersect(target_gene_ids, rownames(expr_numeric))

results_list <- lapply(valid_genes, function(gene) {
  calculate_gene_score(gene, expr_numeric, tumor_samples, normal_samples)
})

results_df <- bind_rows(results_list)

#check
colnames(results_df)
length(valid_genes)
head(rownames(expr_numeric),5)
head(target_gene_ids,5)
head(results_df)
length(expr_numeric)


# ==================== Part 4 : Filtering ====================

results_clean <- results_df[complete.cases(results_df), ]
results_clean$adj_p_value <- p.adjust(results_clean$p_value, method = "fdr")

significant_results <- results_clean %>%
  filter(log2fc > 0.585, adj_p_value < 0.05) %>%
  arrange(desc(log2fc))

# ==================== Output ====================
significant_results <- significant_results %>%
  arrange(desc(log2fc)) %>%
  distinct(gene_id, .keep_all = TRUE)

gene_annotation <- fread("D:/浙大海宁学习大一/DST2/ICA_mysql/BRCA/receptor+secreted_using.csv")
colnames(gene_annotation) <- c("ensembl_id", "type", "gene_symbol", "uniprot_id")

final_results <- significant_results %>%
  left_join(gene_annotation, by = c("gene_id" = "ensembl_id")) %>%
  select(gene_id, gene_symbol, uniprot_id, type, everything())

write.csv(final_results, "D:/浙大海宁学习大一/DST2/ICA_mysql/BRCA/Results/COAD_final_tumor_normal.csv", row.names = FALSE)

if ("type" %in% colnames(final_results)) {
  secreted_genes <- final_results %>% filter(type == "secreted")
  receptor_genes <- final_results %>% filter(type == "receptor")
  write.csv(secreted_genes, "secreted_genes_ranked.csv", row.names = FALSE)
  write.csv(receptor_genes, "receptor_genes_ranked.csv", row.names = FALSE)
}