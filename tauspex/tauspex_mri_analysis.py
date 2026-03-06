import os
import itertools
import numpy as np
import pandas as pd
from os.path import join
import sys

# Append anapyze path
sys.path.append(r"C:\Users\ovourkas\Documents\anapyze-main\src")
from anapyze.analysis import two_samples

# -------------------- PATHS --------------------

db_xlsx = r"C:\Users\ovourkas\Documents\Projects\FDG_vs_MRI_TauPETStaging\Final_Database.xlsx"
dir_patients = r"C:\Users\ovourkas\Documents\Projects\FDG_vs_MRI_TauPETStaging\nifti_MRI_TauPET_Staging"
stat_analysis_dir = r"C:\Users\ovourkas\Documents\Statistical_Analyses_ANOVA_agecovariate"
atlas_nii = r"C:\Users\ovourkas\Documents\Atlases_Templates\Harvard_Oxford_atlas\cortex_Harvard_Oxford_cat12_1.5mm_LZ.nii"

# -------------------- COVARIATES --------------------
# Sex DISCARDED
covar_cols = ["Age_MRI", "Interval_MRI_FDG"]
covar_names = ["Age", "MRI_FDG_Interval"]

# -------------------- HELPER FUNCTION --------------------

def collect_images(df_group, covar_cols):
    """Collect MRI file paths and multiple covariates for a group."""
    images = []
    covars = []

    for _, row in df_group.iterrows():
        folder_name = row["FOLDER_NAME"]
        folder = join(dir_patients, folder_name, "mri")

        if not os.path.exists(folder):
            print(f"Warning: MRI folder not found: {folder}")
            continue

        nii_files = [
            f for f in os.listdir(folder)
            if f.endswith((".nii", ".nii.gz")) and folder_name in f
        ]

        if not nii_files:
            print(f"Warning: No MRI file found matching {folder_name} in {folder}")
            continue

        images.append(join(folder, nii_files[0]))
        covars.append([row[col] for col in covar_cols])

    return images, np.array(covars, dtype=float)

# -------------------- MAIN --------------------

def main():

    # Load database
    df = pd.read_excel(db_xlsx)

    # Assign folder names
    df["FOLDER_NAME"] = df["Subject_ID"]

    # Save for traceability
    df.to_csv(
        r"C:\Users\ovourkas\Documents\anapyze-main\tau_pet_classified_fixed_subjids.csv",
        index=False
    )

    print("Database loaded and folder names assigned.")

    # Tau stages
    stages = [0, 1, 2, 3]

    # Iterate over all pairwise comparisons
    for pair in itertools.combinations(stages, 2):

        try:
            print(f"\n=== Running analysis for Stage {pair[0]} vs Stage {pair[1]} ===")

            output_dir = join(
                stat_analysis_dir,
                f"ROIBased_MRI_Stage{pair[0]}_Stage{pair[1]}"
            )
            os.makedirs(output_dir, exist_ok=True)

            # -------------------- GROUP SELECTION --------------------

            df_group_1 = df[df["TAU_SPEX_New_Stage"] == pair[0]]
            if pair[0] == 0:
                df_group_1 = df_group_1[df_group_1["Cognitive_Status"] == 0]

            df_group_2 = df[df["TAU_SPEX_New_Stage"] == pair[1]]

            # Drop missing data
            df_group_1 = df_group_1.dropna(subset=["FOLDER_NAME"] + covar_cols)
            df_group_2 = df_group_2.dropna(subset=["FOLDER_NAME"] + covar_cols)

            # -------------------- IMAGE + COVARIATE COLLECTION --------------------

            images_group_1, covars_group_1 = collect_images(df_group_1, covar_cols)
            images_group_2, covars_group_2 = collect_images(df_group_2, covar_cols)

            print(f"Group {pair[0]}: {len(images_group_1)} subjects")
            print(f"Group {pair[1]}: {len(images_group_2)} subjects")

            # Sanity checks
            if len(images_group_1) < 2 or len(images_group_2) < 2:
                print("Skipping comparison: not enough subjects.")
                continue

            print("Covariate shapes:",
                  covars_group_1.shape,
                  covars_group_2.shape)

            # -------------------- STATISTICAL ANALYSIS --------------------

            two_samples.run_2sample_anova_with_multiple_covariates_atlas(
                images_group_1,
                images_group_2,
                covars_group_1,
                covars_group_2,
                atlas_nii,
                output_dir,
                operation="mean",
                covar_names=covar_names
            )

        except Exception as e:
            print(f"ERROR for stage pair {pair}: {e}")

# -------------------- ENTRY POINT --------------------

if __name__ == "__main__":
    main()