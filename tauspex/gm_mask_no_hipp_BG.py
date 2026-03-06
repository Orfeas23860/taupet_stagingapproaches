import os
import nibabel as nib
import numpy as np
from nilearn.datasets import fetch_atlas_harvard_oxford
from nilearn.image import resample_to_img

# ----------------------------
# 1️⃣ Paths
# ----------------------------
gm_path = r"C:\Users\ovourkas\Documents\Atlases_Templates\gm_cat12\gm_cat12.nii"
output_path = r"C:\Users\ovourkas\Documents\Atlases_Templates\gm_cat12\gm_no_hippo_BG.nii"

# ----------------------------
# 2️⃣ Load CAT12 GM mask
# ----------------------------
gm_img = nib.load(gm_path)
gm_data = gm_img.get_fdata()
gm_data = np.squeeze(gm_data)  # remove singleton dimensions if any

# ----------------------------
# 3️⃣ Download subcortical Harvard-Oxford atlas
# ----------------------------
sub_atlas = fetch_atlas_harvard_oxford('sub-maxprob-thr50-1mm')

# In recent Nilearn versions, sub_atlas.maps is already a Nifti1Image
atlas_img = sub_atlas.maps
atlas_data = atlas_img.get_fdata()
labels = sub_atlas.labels

print("Atlas loaded with shape:", atlas_data.shape)
print("Atlas labels:", labels)

# ----------------------------
# 4️⃣ Identify hippocampus + basal ganglia IDs
# ----------------------------
exclude_labels = [
    'Left Hippocampus', 'Right Hippocampus',
    'Left Caudate', 'Right Caudate',
    'Left Putamen', 'Right Putamen',
    'Left Pallidum', 'Right Pallidum'
]

# IDs are 1-based according to the label order
exclude_ids = [labels.index(lbl) + 1 for lbl in exclude_labels]

# Build exclusion mask
exclude_mask = np.isin(atlas_data, exclude_ids).astype(np.uint8)

# ----------------------------
# 5️⃣ Resample exclusion mask to GM mask space
# ----------------------------
exclude_mask_resampled = resample_to_img(
    nib.Nifti1Image(exclude_mask, atlas_img.affine),
    gm_img,
    interpolation='nearest'
)
exclude_mask_data = np.squeeze(exclude_mask_resampled.get_fdata()).astype(bool)

# ----------------------------
# 6️⃣ Apply exclusion mask to GM
# ----------------------------
gm_clean_data = gm_data * (~exclude_mask_data)

# ----------------------------
# 7️⃣ Save cleaned GM mask
# ----------------------------
gm_clean_img = nib.Nifti1Image(gm_clean_data, gm_img.affine, gm_img.header)
nib.save(gm_clean_img, output_path)

print(f"Cleaned GM mask saved: {output_path}")
print("You can now open it in SPM, MRIcron, or FSLeyes to inspect it.")
