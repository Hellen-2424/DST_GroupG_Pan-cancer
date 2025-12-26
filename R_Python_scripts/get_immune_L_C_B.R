library(data.table)

# define immune category
IMMUNE_GENE_CATEGORIES <- list(
  chemokine_ligand = c("CXCL", "CCL", "XCL", "CX3CL"),
  chemokine_receptor = c("CXCR", "CCR", "XCR", "CX3CR"),
  cytokine_ligand = c("IL", "TNF", "IFN", "TGFB", "CSF", "KITLG"),
  cytokine_receptor = c("IL", "TNFR", "IFNR", "TGFBR", "CSFR", "KIT"),
  checkpoint = c("PDCD1", "CD274", "PDCD1LG2", "CTLA4", "LAG3", "TIGIT", "HAVCR2"),
  immunoglobulin = c("IGH", "IGK", "IGL", "IGV", "IGJ"),
  fc_receptor = c("FCRL", "FCGR", "FCER", "FCMR"),
  tnfrsf = c("TNFRSF"),
  other_immune_receptor = c("TLR", "NLR", "CLEC", "CD", "HLA", "KIR", "LILR"),
  antigen_presentation = c("HLA", "B2M", "TAP", "CALR")
)

df1 <- fread("D:/浙大海宁学习大一/DST2/ICA_mysql/BRCA/L_C_B.csv", header = TRUE)

get_immune_category <- function(gene) {
  if(gene %in% IMMUNE_GENE_CATEGORIES$checkpoint) return("checkpoint")
  
  for(cat in names(IMMUNE_GENE_CATEGORIES)) {
    if(cat == "checkpoint") next
    for(prefix in IMMUNE_GENE_CATEGORIES[[cat]]) {
      if(grepl(paste0("^", prefix), gene)) return(cat)
    }
  }
  return(NA)
}

df1$immune_category <- sapply(df1$gene_symbol, get_immune_category)

write.csv(df1, "D:/浙大海宁学习大一/DST2/ICA_mysql/BRCA/immune_L_C_B.csv", row.names = FALSE)


