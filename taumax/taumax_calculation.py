
# Tau-MaX COMPUTATION 

import os
import numpy as np
import nibabel as nib
import pandas as pd
import matplotlib.pyplot as plt
from nilearn.image import resample_to_img
from sklearn.mixture import GaussianMixture

base_dir = r"C:\Users\ovourkas\Documents\Projects\FDG_vs_MRI_TauPETStaging\nifti_TAU_TauPET_Staging"
atlas_path = r"C:\Users\ovourkas\Documents\Atlases_Templates\Harvard_Oxford_atlas\cortex_Harvard_Oxford_cat12_1.5mm_LZ.nii"
gm_mask_path = r"C:\Users\ovourkas\Documents\Atlases_Templates\gm_cat12\gm_cat12.nii"

tau_files = sorted([
    os.path.join(r, f)
    for r, _, fs in os.walk(base_dir)
    for f in fs if f.endswith((".nii", ".nii.gz"))
])

n_subj = len(tau_files)
if n_subj == 0:
    raise RuntimeError("No Tau-PET files found")

print(f"Found {n_subj} Tau-PET scans")

ref_img = nib.load(tau_files[0])

atlas = resample_to_img(
    atlas_path, ref_img,
    interpolation="nearest",
    force_resample=True
).get_fdata().astype(int)

gm_mask = resample_to_img(
    gm_mask_path, ref_img,
    interpolation="nearest",
    force_resample=True
).get_fdata()

if gm_mask.ndim == 4:
    gm_mask = np.squeeze(gm_mask)

gm_mask = gm_mask > 0.2

roi_labels = np.unique(atlas)
roi_labels = roi_labels[roi_labels > 0]
n_rois = len(roi_labels)

print(f"Total cortical ROIs: {n_rois}")

voxel_volume = np.prod(ref_img.header.get_zooms())

roi_volumes = np.array([
    np.sum((atlas == roi) & gm_mask) * voxel_volume
    for roi in roi_labels
])

relative_roi_volumes = roi_volumes / roi_volumes.sum()

roi_suvr = np.zeros((n_subj, n_rois))

for s, fpath in enumerate(tau_files):

    img = nib.load(fpath)
    data = resample_to_img(
        img, ref_img,
        interpolation="linear",
        force_resample=True
    ).get_fdata()

    if data.ndim == 4:
        data = np.squeeze(data)

    for r, roi in enumerate(roi_labels):
        mask = (atlas == roi) & gm_mask
        roi_suvr[s, r] = np.nanmean(data[mask]) if mask.any() else np.nan

valid_roi_mask = np.ones(n_rois, dtype=bool)
TPI = np.zeros_like(roi_suvr)
tau_positive = np.zeros_like(roi_suvr, dtype=bool)

for r in range(n_rois):

    vals = roi_suvr[:, r]
    vals = vals[~np.isnan(vals)].reshape(-1, 1)

    if len(vals) < 10:
        valid_roi_mask[r] = False
        continue

    gmm_2 = GaussianMixture(n_components=2, n_init=20, random_state=0)
    gmm_2.fit(vals)
    bic_2 = gmm_2.bic(vals)

    gmm_1 = GaussianMixture(n_components=1, n_init=10, random_state=0)
    gmm_1.fit(vals)
    bic_1 = gmm_1.bic(vals)

    if bic_1 < bic_2:
        valid_roi_mask[r] = False
        continue

    means = gmm_2.means_.flatten()
    covs = gmm_2.covariances_.flatten()

    normal_idx = np.argmin(means)
    path_idx = np.argmax(means)

    mu_norm = means[normal_idx]
    sigma_norm = np.sqrt(covs[normal_idx])

    tau_threshold = mu_norm + 2.0 * sigma_norm

    post = gmm_2.predict_proba(roi_suvr[:, r].reshape(-1,1))
    TPI[:, r] = 100 * (post[:, path_idx] - post[:, normal_idx])

    tau_positive[:, r] = roi_suvr[:, r] > tau_threshold

roi_suvr = roi_suvr[:, valid_roi_mask]
TPI = TPI[:, valid_roi_mask]
tau_positive = tau_positive[:, valid_roi_mask]
relative_roi_volumes = relative_roi_volumes[valid_roi_mask]

print(f"ROIs after exclusion: {roi_suvr.shape[1]}")

mean_tau_suvr = np.zeros(n_subj)
tau_extent = np.zeros(n_subj)
tau_max = np.zeros(n_subj)

for s in range(n_subj):

    mean_tau_suvr[s] = np.sum(
        roi_suvr[s] * relative_roi_volumes
    )

    tau_extent[s] = (
        np.sum(relative_roi_volumes[tau_positive[s]]) * 100
    )

    tau_max[s] = np.sum(
        TPI[s, tau_positive[s]] *
        relative_roi_volumes[tau_positive[s]]
    )


results_df = pd.DataFrame({
    "Subject": [os.path.basename(f) for f in tau_files],
    "Mean_Tau_SUVR": mean_tau_suvr,
    "Tau_Extent_percent": tau_extent,
    "Tau_MaX": tau_max
})

excel_path = os.path.join(base_dir, "TauMeasures_BrownEtAl.xlsx")
results_df.to_excel(excel_path, index=False)

print(f"Results saved to {excel_path}")

plt.figure(figsize=(6,5))
plt.scatter(tau_extent, tau_max, alpha=0.7)
plt.xlabel("Tau Extent (%)")
plt.ylabel("Tau-MaX")
plt.title("Validation: Tau-MaX vs Tau Extent")
plt.grid(alpha=0.3)
plt.tight_layout()
plt.show()

plt.figure(figsize=(6,5))
plt.scatter(mean_tau_suvr, tau_max, alpha=0.7)
plt.xlabel("Mean Tau SUVR")
plt.ylabel("Tau-MaX")
plt.title("Validation: Tau-MaX vs Mean Tau SUVR")
plt.grid(alpha=0.3)
plt.tight_layout()
plt.show()
