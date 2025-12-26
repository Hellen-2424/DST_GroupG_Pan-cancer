import pandas as pd

reactome_file = "Ensembl2Reactome_All_Levels.txt"
immune_file = "immune_L_C_B.csv"

colnames = ["ensembl_id", "reactome_id", "pathway_name", "evidence", "species"]
reactome = pd.read_csv(reactome_file, sep="\t", header=None, names=colnames)

reactome = reactome[reactome["species"] == "Homo sapiens"]

immune = pd.read_csv(immune_file)
immune["gene_id"] = immune["gene_id"].astype(str).str.strip()
reactome["ensembl_id"] = reactome["ensembl_id"].astype(str).str.strip()

immune_ids = set(immune["gene_id"].unique())
reactome_immune = reactome[reactome["ensembl_id"].isin(immune_ids)]

pathways = reactome_immune[["pathway_name"]].drop_duplicates()

def categorize_pathway(name):
    n = str(name).lower()
    if any(x in n for x in ["immune system", "chemokine", "interleukin", "interferon", "cytokine", "adaptive immune", "innate immune", "b cell", "t cell", "lymphocyte"]):
        return "Immune system"
    if any(x in n for x in ["signal transduction", "signaling by"]):
        return "Signal transduction"
    return "Other"

pathways["category"] = pathways["pathway_name"].apply(categorize_pathway)
pathways["source"] = "Reactome"
pathways = pathways.reset_index(drop=True)
pathways["pathway_id"] = pathways.index + 1

pathways = pathways[["pathway_id", "pathway_name", "category", "source"]]
pathways.to_csv("pathway.csv", index=False)
