setwd("D:/浙大海宁学习大一/DST2/ICA_mysql/BRCA")
library(dplyr)
library(data.table)

df1 <- fread("D:/浙大海宁学习大一/DST2/ICA_mysql/BRCA/Results/LUAD_final_tumor_normal.csv", header = TRUE)
df2 <- fread("D:/浙大海宁学习大一/DST2/ICA_mysql/BRCA/Results/BRCA_final_tumor_normal.csv", header = TRUE)
df3 <- fread("D:/浙大海宁学习大一/DST2/ICA_mysql/BRCA/Results/COAD_final_tumor_normal.csv", header = TRUE)

  colnames(df1) <- tolower(colnames(df1))
  colnames(df2) <- tolower(colnames(df2))
  colnames(df3) <- tolower(colnames(df3))
  
  genes1 <- unique(df1[, c("gene_id", "gene_symbol","uniprot_id","type")])
  genes2 <- unique(df2[, c("gene_id", "gene_symbol","uniprot_id","type")])
  genes3 <- unique(df3[, c("gene_id", "gene_symbol","uniprot_id","type")])
  
  # merge
  all_genes <- unique(rbind(genes1, genes2, genes3))
  
  # count 1
  gene_counts <- data.frame(
    gene_id = all_genes$gene_id,
    gene_symbol = all_genes$gene_symbol,
    gene_type = all_genes$type,
    uniprot_id = all_genes$uniprot_id,
    count = 0,
    stringsAsFactors = FALSE
  )
  
  # count 2
  for(i in 1:nrow(gene_counts)) {
    ensembl <- gene_counts$gene_id[i]
    symbol <- gene_counts$gene_symbol[i]
    type <-gene_counts$gene_type[i]
    uniprot <- gene_counts$uniprot_id[i]
    
    count <- 0
    if(any(df1$gene_id == ensembl | df1$gene_symbol == symbol)) count <- count + 1
    if(any(df2$gene_id == ensembl | df2$gene_symbol == symbol)) count <- count + 1
    if(any(df3$gene_id == ensembl | df3$gene_symbol == symbol)) count <- count + 1
    
    gene_counts$count[i] <- count
  }
  
  gene_counts <- gene_counts[order(-gene_counts$count, gene_counts$gene_symbol), ]
  
  gene_counts$in_files <- sapply(1:nrow(gene_counts), function(i) {
    ensembl <- gene_counts$gene_id[i]
    symbol <- gene_counts$gene_symbol[i]
    type <-gene_counts$gene_type[i]
    uniprot <- gene_counts$uniprot_id[i]
    
    files <- c()
    if(any(df1$gene_id == ensembl | df1$gene_symbol == symbol)) files <- c(files, "LUAD")
    if(any(df2$gene_id == ensembl | df2$gene_symbol == symbol)) files <- c(files, "COAD")
    if(any(df3$gene_id == ensembl | df3$gene_symbol == symbol)) files <- c(files, "BRCA")
    paste(files, collapse = ", ")
  })
  
  write.csv(gene_counts,"D:/浙大海宁学习大一/DST2/ICA_mysql/BRCA/L_C_B.csv",row.names = FALSE)
getwd()
  
}