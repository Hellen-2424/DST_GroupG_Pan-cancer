import pandas as pd
import numpy as np

def pan_cancer_data_filter(input_file, output_file):
    """
    Pan-cancer research data screening:
    1. Retain one tumor sample and one normal sample per case
    2. Filter out records with missing core fields
    3. Exclude cancer types with a sample size <50
    4. Mandatory retention of records with follow-up/recurrence data
    :param input_file: The integrated TCGA Excel file
    :param output_file: The path to save the filtered results
    """
    # 1. Read the data and standardize the missing value indicators
    df = pd.read_excel(input_file, engine="openpyxl")
    df = df.replace('--', pd.NA)  # Unified TCGA Missing Value Identification
    print(f"Total number of records of raw data：{len(df)}")
    print(f"Original number of cancer types involved：{df['cases.disease_type'].nunique()}")

    # 2. Step 1: Label the sample type（tumor/normal）
    # Adapt to the naming rules of TCGA sample types (adjustable according to actual field values)
    def label_sample_type(sample_type):
        if pd.isna(sample_type):
            return np.nan
        # TCGA sample type keyword matching（Tumor/Normal）
        if "Tumor" in str(sample_type) or "tumor" in str(sample_type):
            return "tumor"
        elif "Normal" in str(sample_type) or "normal" in str(sample_type):
            return "normal"
        else:
            return np.nan  # Non target sample type
    
    df["sample_category"] = df["samples.sample_type"].apply(label_sample_type)
    # Remove non tumor/normal samples
    df = df.dropna(subset=["sample_category"])
    print(f"Record number after filtering non tumor/normal samples：{len(df)}")

    # 3. Step 2: Keep records of non missing core fields
    core_clinical_fields = [
        "cases.disease_type",       # Cancer type (pan cancer core)
        "demographic.gender",       # Gender
        "demographic.age_at_index", # Diagnosis age
        "cases.submitter_id",       # Patient ID
        "samples.sample_id",        # Sample ID
        "follow_ups.days_to_follow_up",          # Follow-up days (mandatory)
        "follow_ups.progression_or_recurrence"   # Recurrence status (mandatory)
    ]
    df = df.dropna(subset=core_clinical_fields, how="any")
    print(f"Number of non missing records for core fields (including follow-up recurrence)：{len(df)}")

    # 4. Step 3: Reserve 1 tumor and 1 normal sample for each case (patient ID)
    def get_one_tumor_one_normal(group):
        """Grouping function: Each patient only retains 1 tumor and 1 normal sample"""
        tumor_samples = group[group["sample_category"] == "tumor"].head(1)
        normal_samples = group[group["sample_category"] == "normal"].head(1)
        return pd.concat([tumor_samples, normal_samples], ignore_index=True)
    
    # Extract 1 pair of samples grouped by patient ID
    df_case_samples = df.groupby("cases.submitter_id").apply(get_one_tumor_one_normal).reset_index(drop=True)
    # Filter: Only retain patients with both tumor and normal (at least 2 records per group)
    case_sample_count = df_case_samples["cases.submitter_id"].value_counts()
    valid_cases = case_sample_count[case_sample_count >= 2].index
    df_case_samples = df_case_samples[df_case_samples["cases.submitter_id"].isin(valid_cases)]
    print(f"Record the number of cases after retaining 1 tumor and 1 normal for each case：{len(df_case_samples)}")
    print(f"Effective number of patients：{len(valid_cases)}")

    # 5. Step 4: Screen out cancer types with a sample size<50 (to ensure statistical validity)
    cancer_sample_count = df_case_samples["cases.disease_type"].value_counts()
    valid_cancers = cancer_sample_count[cancer_sample_count >= 50].index.tolist()
    df_final = df_case_samples[df_case_samples["cases.disease_type"].isin(valid_cancers)]
    print(f"Number of cancer types with a sample size of ≥ 50 retained：{len(valid_cancers)}")
    print(f"Number of records after final screening：{len(df_final)}")

    # 6. Reset the index and save the results
    df_final = df_final.reset_index(drop=True)
    df_final.to_excel(output_file, index=False, engine="openpyxl")
    
    # Overview of filtered data output
    print("\n===== Overview of data after pan cancer screening =====")
    print(f"Final number of patients：{df_final['cases.submitter_id'].nunique()}")
    print(f"Final number of cancer types：{df_final['cases.disease_type'].nunique()}")
    print("\nSample size distribution of various cancer types：")
    print(df_final["cases.disease_type"].value_counts())
    print(f"\nThe filtering results have been saved to：{output_file}")
    return df_final

# -------------------------- run setup --------------------------
if __name__ == "__main__":
    # Replace with actual file path
    INPUT_FILE = "tcga_integrated_data.xlsx"   # Input integrated data
    OUTPUT_FILE = "pan_cancer_filtered_data.xlsx"  # Output filtered data
    # Execute filtering function
    filtered_df = pan_cancer_data_filter(INPUT_FILE, OUTPUT_FILE)