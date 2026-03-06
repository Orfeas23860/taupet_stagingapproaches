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
dir_patients_fdg = r"C:\Users\ovourkas\Documents\Projects\FDG_vs_MRI_TauPETStaging\nifti_FDG_TauPET_Staging"
stat_analysis_dir = r"C:\Users\ovourkas\Documents\Statistical_Analyses_ANOVA_agecovariate_FDG"
atlas_nii = r"C:\Users\ovourkas\Documents\Atlases_Templates\Harvard_Oxford_atlas\cortex_Harvard_Oxford_cat12_1.5mm_LZ.nii"

# -------------------- COVARIATES --------------------
# Sex discarded, same second covariate as MRI

covar_cols = ["Age_FDG", "Interval_MRI_FDG"]
covar_names = ["Age", "MRI_FDG_Interval"]

# -------------------- HELPER FUNCTION --------------------

def collect_fdg_images(df_group, covar_cols):
    """
    Collect FDG file paths and multiple covariates for a group.
    Looks for 'swfdg_normhist.nii' inside each subject's folder.
    """
    images = []
    covars = []

    for _, row in df_group.iterrows():
        folder_name = row["FOLDER_NAME"]
        folder = join(dir_patients_fdg, folder_name)
        fdg_file = join(folder, "swfdg_normhist.nii")

        if not os.path.exists(fdg_file):
            print(f"Warning: FDG file not found for {folder_name}: {fdg_file}")
            continue

        images.append(fdg_file)
        covars.append([row[col] for col in covar_cols])

    return images, np.array(covars, dtype=float)

# -------------------- MAIN --------------------

def main_fdg():

    # Load spreadsheet
    df = pd.read_excel(db_xlsx)

    # Assign folder names
    df["FOLDER_NAME"] = df["Subject_ID"]

    # Save for traceability
    df.to_csv(
        r"C:\Users\ovourkas\Documents\anapyze-main\tau_pet_classified_fixed_subjids_fdg.csv",
        index=False
    )

    print("Database loaded and folder names assigned for FDG.")

    # Tau stages
    stages = [0, 1, 2, 3]

    # Iterate over stage pairs
    for pair in itertools.combinations(stages, 2):

        try:
            print(f"\n=== Running FDG analysis for Stage {pair[0]} vs Stage {pair[1]} ===")

            output_dir = join(
                stat_analysis_dir,
                f"ROIBased_FDG_Stage{pair[0]}_Stage{pair[1]}"
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

            images_group_1, covars_group_1 = collect_fdg_images(df_group_1, covar_cols)
            images_group_2, covars_group_2 = collect_fdg_images(df_group_2, covar_cols)

            print(f"Group {pair[0]}: {len(images_group_1)} FDG subjects")
            print(f"Group {pair[1]}: {len(images_group_2)} FDG subjects")

            # Sanity checks
            if len(images_group_1) < 2 or len(images_group_2) < 2:
                print("Skipping comparison: not enough subjects.")
                continue

            print(
                "Covariate shapes:",
                covars_group_1.shape,
                covars_group_2.shape
            )

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
            print(f"ERROR for FDG stage pair {pair}: {e}")

# -------------------- ENTRY POINT --------------------

if __name__ == "__main__":
    main_fdg()
