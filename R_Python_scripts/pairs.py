import pandas as pd
import numpy as np

print("The biological errors in the ligand-receptor pairing table are being corrected...")

# --------------------------------
# 1. Read the original pairing table
# --------------------------------
try:
    pairs_df = pd.read_csv("simple_ligand_receptor_pairs.csv")
    print(f"Read the original pairing table, which consists of {len(pairs_df)} rows")
except FileNotFoundError:
    print("The original pairing table was not found. A new pairing table was created based on biological knowledge")
    
    # Create the correct pairing based on the immune gene data you provide
    # First, read the original genetic data
    gene_data = pd.read_csv("immune_L_C_B.csv")
    
    # Create the correct ligand-receptor pair (based on biological literature)
    correct_pairs = [
        # =================== Chemical factor ligand-receptor pairs ===================
        ("CXCL10", "chemokine_ligand", "CXCR3", "chemokine_receptor", "chemokine"),
        ("CXCL11", "chemokine_ligand", "CXCR3", "chemokine_receptor", "chemokine"),
        ("CXCL13", "chemokine_ligand", "CXCR5", "chemokine_receptor", "chemokine"),
        ("CCL11", "chemokine_ligand", "CCR3", "chemokine_receptor", "chemokine"),
        ("CXCL1", "chemokine_ligand", "CXCR2", "chemokine_receptor", "chemokine"),
        ("CXCL2", "chemokine_ligand", "CXCR2", "chemokine_receptor", "chemokine"),
        ("CXCL3", "chemokine_ligand", "CXCR2", "chemokine_receptor", "chemokine"),
        ("CXCL5", "chemokine_ligand", "CXCR2", "chemokine_receptor", "chemokine"),
        ("CXCL8", "chemokine_ligand", "CXCR1", "chemokine_receptor", "chemokine"),
        ("CXCL8", "chemokine_ligand", "CXCR2", "chemokine_receptor", "chemokine"),
        ("CXCL9", "chemokine_ligand", "CXCR3", "chemokine_receptor", "chemokine"),
        
        # =================== TNF Family ligand-receptor pairs ===================
        # Ligand: Starting with "TNFSF", Receptor: Starting with "TNFRSF
        ("TNFSF12", "cytokine_ligand", "TNFRSF12A", "cytokine_receptor", "TNF_family"),
        ("TNFSF13B", "cytokine_ligand", "TNFRSF17", "cytokine_receptor", "TNF_family"),
        ("TNFSF15", "cytokine_ligand", "TNFRSF25", "cytokine_receptor", "TNF_family"),
        ("TNFSF4", "cytokine_ligand", "TNFRSF4", "cytokine_receptor", "TNF_family"),
        ("TNFSF9", "cytokine_ligand", "TNFRSF9", "cytokine_receptor", "TNF_family"),
        ("TNFSF18", "cytokine_ligand", "TNFRSF18", "cytokine_receptor", "TNF_family"),
        
        # =================== 细Cytokine ligand-receptor pairs ===================
        ("IL11", "cytokine_ligand", "IL11RA", "cytokine_receptor", "cytokine"),
        ("IL23A", "cytokine_ligand", "IL23R", "cytokine_receptor", "cytokine"),
        ("IL1RN", "cytokine_ligand", "IL1R1", "cytokine_receptor", "cytokine"),
        ("IL37", "cytokine_ligand", "IL18R1", "cytokine_receptor", "cytokine"),
        
        # =================== Fc ===================
        # FCGR1A is a receptor, and the ligand is immunoglobulin (specific genes are required)
        ("IGHG1", "immunoglobulin", "FCGR1A", "fc_receptor", "Fc_receptor"),
        ("IGHG2", "immunoglobulin", "FCGR1A", "fc_receptor", "Fc_receptor"),
        ("IGHG3", "immunoglobulin", "FCGR1A", "fc_receptor", "Fc_receptor"),
        ("IGHG4", "immunoglobulin", "FCGR1A", "fc_receptor", "Fc_receptor"),
        ("IGHG1", "immunoglobulin", "FCRL5", "fc_receptor", "Fc_receptor_like"),
        
        # =================== B-cell receptor complex ===================
        ("CD19", "other_immune_receptor", "CD81", "other_immune_receptor", "B_cell_receptor_complex"),
        ("CD79A", "other_immune_receptor", "CD79B", "other_immune_receptor", "B_cell_receptor_complex"),
        
        # =================== immune ===================
        ("CD274", "immune_checkpoint", "PDCD1", "immune_checkpoint", "immune_checkpoint"),
        ("CD80", "immune_checkpoint", "CD28", "immune_checkpoint", "co-stimulation"),
        ("CD86", "immune_checkpoint", "CD28", "immune_checkpoint", "co-stimulation"),
    ]
    
    # shift to DataFrame
    pairs_df = pd.DataFrame(correct_pairs, 
                           columns=['ligand', 'ligand_category', 'receptor', 'receptor_category', 'pair_type'])
    
    # Mark which genes are in our data
    gene_in_data = set(gene_data['gene_symbol'])
    pairs_df['ligand_in_data'] = pairs_df['ligand'].isin(gene_in_data)
    pairs_df['receptor_in_data'] = pairs_df['receptor'].isin(gene_in_data)
    pairs_df['both_in_data'] = pairs_df['ligand_in_data'] & pairs_df['receptor_in_data']

# --------------------------------
# 2. Correct biological errors
# --------------------------------
print("\nCorrecting biological errors.")

def fix_biological_errors(df):
    """Correcting biological errors in the ligand-receptor relationship"""
    
    # copy data
    corrected_df = df.copy()
    
    # Correction rule dictionary: Requires swapping the pairing of ligand and receptor
    swap_pairs = {
        # If the ligand starts with "TNFRSF" and the receptor starts with "TNFSF", they should be swapped.
        # Because TNFRSF is the receptor and TNFSF is the ligand.
    }

# A simpler method: Directly identify and correct the incorrect pairings
    corrections_made = 0
    
    for idx, row in corrected_df.iterrows():
        ligand = row['ligand']
        receptor = row['receptor']
        
        # Rule 1: If the ligand starts with "TNFRSF" and the receptor starts with "TNFSF", they should be swapped.
        if ligand.startswith('TNFRSF') and receptor.startswith('TNFSF'):
            # Exchange ligands and receptors
            corrected_df.at[idx, 'ligand'] = receptor
            corrected_df.at[idx, 'receptor'] = ligand
            corrected_df.at[idx, 'ligand_category'] = row['receptor_category']
            corrected_df.at[idx, 'receptor_category'] = row['ligand_category']
            
            corrections_made += 1
            print(f"Correct the {idx+1}th line: {ligand}↔{receptor} role exchange")
        
        # Rule 2: Ensure ligand and receptor categories are reasonable
        elif ligand.startswith('TNFRSF') and row['ligand_category'] == 'cytokine_ligand':
            corrected_df.at[idx, 'ligand_category'] = 'cytokine_receptor'
            corrections_made += 1
            print(f"Correct the {idx+1}th line: {ligand} change the category to cytokine_receptor")
        
        elif receptor.startswith('TNFSF') and row['receptor_category'] == 'cytokine_receptor':
            corrected_df.at[idx, 'receptor_category'] = 'cytokine_ligand'
            corrections_made += 1
            print(f"Correct the {idx+1}th line : {receptor} change the category to cytokine_ligand")
    
    print(f"Total {corrections_made} biological errors have been corrected.")
    return corrected_df

# Apply corrections
corrected_pairs = fix_biological_errors(pairs_df)

# --------------------------------
# 3. Create a more accurate classification system
# --------------------------------
print("\nOptimize the classification system...")

# Define a more detailed classification of immune molecules
category_mapping = {
    'chemokine_ligand': 'Chemokine Ligand',
    'chemokine_receptor': 'Chemokine Receptor',
    
    'cytokine_ligand': 'Cytokine Ligand',
    'cytokine_receptor': 'Cytokine Receptor',
    
    'TNF_ligand': 'TNF Superfamily Ligand',
    'TNF_receptor': 'TNF Superfamily Receptor',
    
    'fc_receptor': 'Fc Receptor',
    'immunoglobulin': 'Immunoglobulin',
    'Fc_receptor_like': 'Fc Receptor-like',
    
    'other_immune_receptor': 'Other Immune Receptor',
    'immune_checkpoint': 'Immune Checkpoint',
    'Unknown': 'Unknown'
}

# Apply more detailed classification
corrected_pairs['ligand_category_detailed'] = corrected_pairs['ligand_category'].map(
    lambda x: category_mapping.get(x, x)
)
corrected_pairs['receptor_category_detailed'] = corrected_pairs['receptor_category'].map(
    lambda x: category_mapping.get(x, x)
)

# --------------------------------
# 4. Add biological annotations
# --------------------------------
print("Add biological annotations...")

# Add functional annotations
def add_functional_annotation(row):
    ligand = row['ligand']
    receptor = row['receptor']
    pair_type = row['pair_type']
    
    functional_annotations = {
        'chemokine': 'Immune cell chemotaxis and migration',
        'TNF_family': 'Cell survival, apoptosis and inflammation regulation',
        'cytokine': 'Immune cell activation and differentiation signals',  
        'Fc_receptor': 'Antibody-mediated immune response',
        'Fc_receptor_like': 'Regulatory B cell function',
        'B_cell_receptor_complex': 'B cell antigen recognition and activation',
        'immune_checkpoint': 'Regulation and inhibition of immune response'
    }
    
    return functional_annotations.get(pair_type, 'Other immune functions')

corrected_pairs['biological_function'] = corrected_pairs.apply(add_functional_annotation, axis=1)

# Add Channel information
def add_pathway_info(row):
    ligand = row['ligand']
    
    if ligand.startswith('CXCL'):
        return 'Chemokine Signaling Pathway'
    elif ligand.startswith('TNF'):
        return 'TNF Signaling Pathway'
    elif ligand.startswith('IL'):
        return 'Cytokine-Cytokine Receptor Interaction'
    elif 'IGH' in ligand or 'IGK' in ligand or 'IGL' in ligand:
        return 'B Cell Receptor Signaling'
    elif 'FCGR' in ligand or 'FCRL' in ligand:
        return 'Fc Receptor Signaling'
    else:
        return 'Other Immune Pathway'

corrected_pairs['pathway'] = corrected_pairs.apply(add_pathway_info, axis=1)

# --------------------------------
# 5. Filter the pairings related to our data
# --------------------------------
print("Filter the pairings related to our data...")

# Only retain at least one pairing in our genetic data
if 'ligand_in_data' in corrected_pairs.columns:
    relevant_pairs = corrected_pairs[
        corrected_pairs['ligand_in_data'] | corrected_pairs['receptor_in_data']
    ]
else:
    # If these columns are missing, create them
    gene_data = pd.read_csv("immune_L_C_B.csv")
    gene_in_data = set(gene_data['gene_symbol'])
    
    corrected_pairs['ligand_in_data'] = corrected_pairs['ligand'].isin(gene_in_data)
    corrected_pairs['receptor_in_data'] = corrected_pairs['receptor'].isin(gene_in_data)
    corrected_pairs['both_in_data'] = corrected_pairs['ligand_in_data'] & corrected_pairs['receptor_in_data']
    
    relevant_pairs = corrected_pairs[
        corrected_pairs['ligand_in_data'] | corrected_pairs['receptor_in_data']
    ]

print(f"Filter out {len(relevant_pairs)} related pairs")

# --------------------------------
# 6. Save the corrected result
# --------------------------------
print("\nSave the corrected resul...")

# Create an output directory
import os
output_dir = "corrected_results"
os.makedirs(output_dir, exist_ok=True)

# Save the complete correction table
full_output = os.path.join(output_dir, "corrected_ligand_receptor_pairs_full.csv")
corrected_pairs.to_csv(full_output, index=False, encoding='utf-8-sig')
print(f"The complete correction table has been saved: {full_output}")

# Save the relevant pairing table
relevant_output = os.path.join(output_dir, "corrected_ligand_receptor_pairs_relevant.csv")
relevant_pairs.to_csv(relevant_output, index=False, encoding='utf-8-sig')
print(f"The relevant pairing table has been saved: {relevant_output}")

# Statistics by category
category_stats = relevant_pairs.groupby(['ligand_category_detailed', 'receptor_category_detailed']).size().reset_index(name='pair_count')
category_stats = category_stats.sort_values('pair_count', ascending=False)

stats_output = os.path.join(output_dir, "pair_category_stats.csv")
category_stats.to_csv(stats_output, index=False, encoding='utf-8-sig')
print(f"Category statistics have been saved : {stats_output}")

# --------------------------------
# 7. Generate a correction report
# --------------------------------
print("\nGenerate a correction report...")

report_content = f"""
{'='*80}
Ligand-Receptor Pairing Table Biological Correction Report
{'='*80}

Correction time: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}

[1] Correction Overview
    Original pairing count: {len(pairs_df)}
    Corrected pairing count: {len(corrected_pairs)}
    Relevant pair count: {len(relevant_pairs)}
Pairing present in both datasets: {len(relevant_pairs[relevant_pairs['both_in_data']])}

[2] Main corrections:
    1. Corrected the issue of reversed roles of TNF family ligands and receptors
    2. Optimized the molecular classification system
    3. Added biological function annotations
    4. Added pathway information

[3] Key Pairing Type Statistics
"""
# Statistical pairing types
pair_type_stats = relevant_pairs['pair_type'].value_counts()
for pair_type, count in pair_type_stats.items():
    report_content += f"    {pair_type}: {count} 对\n"

# Show the important pairings where both parties are in the data
both_in_data = relevant_pairs[relevant_pairs['both_in_data']]
if len(both_in_data) > 0:
    report_content += f"""
[4] Important ppairs in both parties' data (a total of {len(both_in_data)} pairs)
"""
    for idx, row in both_in_data.head(10).iterrows():
        report_content += f"    {row['ligand']} → {row['receptor']} ({row['biological_function']})\n"


report_content += f"""


# Save the report
report_file = os.path.join(output_dir, "correction_report.txt")
with open(report_file, 'w', encoding='utf-8') as f:
    f.write(report_content)

print(f"The correction report has been saved: {report_file}")

# --------------------------------
# 8. show important results
# --------------------------------

print("\n" + "="*60)
print("Summary of Important Revision results")
print("="*60)

print(f"\n statistical result:")
print(f"  Overall correction pairing: {len(corrected_pairs)}")
print(f"  Related to your data: {len(relevant_pairs)}")
print(f"  Both sides are in the data: {len(both_in_data)}")

print(f"\n Key immune pathways:")
pathway_stats = relevant_pairs['pathway'].value_counts()
for pathway, count in pathway_stats.head(5).items():
    print(f"  • {pathway}: {count} 对")

print(f"\n Recommended pairings for focused research:")
if len(both_in_data) > 0:
    for idx, row in both_in_data.head(3).iterrows():
        print(f"  • {row['ligand']} → {row['receptor']} ({row['biological_function']})")

print(f"\n Output file location: {output_dir}/")
print("   - corrected_ligand_receptor_pairs_*.csv")
print("   - pair_category_stats.csv")
print("   - correction_report.txt")

print("\n" + "="*60)
print("Biological error correction completed!")
print("="*60)



