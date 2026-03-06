% Main directory with subject folders
main_dir = 'C:\Users\ovourkas\Documents\Projects\FDG_vs_MRI_TauPETStaging\Staging_Approaches\InVivo_Braak_Staging\CenTauR_Scaling_Tau_PET\Tau_PET_data';

% Path to Excel file with reference values
ref_excel_path = 'C:\Users\ovourkas\Documents\Projects\FDG_vs_MRI_TauPETStaging\Staging_Approaches\InVivo_Braak_Staging\CenTauR_Scaling_Tau_PET\Cerebellum_TauSignal.xlsx';
ref_table = readtable(ref_excel_path);

% GM mask template path (fixed for all subjects)
gm_path = 'C:\Users\ovourkas\Documents\Atlases_Templates\gm_cat12\gm_no_hippo_BG.nii';

% List subject folders (skip '.', '..', hidden)
subject_folders = dir(main_dir);
subject_folders = subject_folders([subject_folders.isdir] & ~startsWith({subject_folders.name}, '.'));

% Prepare storage
subject_names = {};
tau_spex_values = [];

% Initialize SPM
spm('defaults', 'FMRI');
spm_jobman('initcfg');

for i = 1:length(subject_folders)
    subject_name = subject_folders(i).name;
    subject_dir = fullfile(main_dir, subject_name);

    try
        % Get reference mean from Excel for this subject
        ref_row = find(strcmp(ref_table{:,1}, subject_name));
        if isempty(ref_row)
            error('Reference region value not found for subject %s in Excel.', subject_name);
        end
        ref_mean = ref_table{ref_row, 2};

        % Find PET file (wr*.nii)
        pet_file = dir(fullfile(subject_dir, 'wr*.nii'));
        if isempty(pet_file)
            error('No PET file found.');
        end
        pet_path = fullfile(subject_dir, pet_file(1).name);

        % ---- Coregistration Reslice GM mask template onto PET space ----
        
        % Parse GM mask template file name
        [~, gm_name, ext] = fileparts(gm_path);

        % Prepare SPM coregister reslice batch (only reslicing)
        matlabbatch = [];
        matlabbatch{1}.spm.spatial.coreg.write.ref = {pet_path};
        matlabbatch{1}.spm.spatial.coreg.write.source = {gm_path};
        matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 0; % nearest neighbor
        matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
        matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
        matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'r';

        % Run reslice job
        spm_jobman('run', matlabbatch);

        % Copy resliced GM mask from template folder to subject folder
        resliced_gm_source = fullfile(fileparts(gm_path), ['r' gm_name ext]);
        resliced_gm_dest = fullfile(subject_dir, ['r' gm_name ext]);
        copyfile(resliced_gm_source, resliced_gm_dest);

        % Load PET and resliced GM mask volumes from subject folder
        pet_vol = spm_read_vols(spm_vol(pet_path));
        gm_vol = spm_read_vols(spm_vol(resliced_gm_dest));

        % Threshold PET in GM voxels
        gm_mask = gm_vol > 0.5;  % threshold probability map at 0.5 to binarize
        threshold = 1.65 * ref_mean;

        pet_gm = pet_vol;
        pet_gm(~gm_mask) = 0;

        above_thresh = pet_gm > threshold;

        % Cluster size thresholding (≥ 50 voxels)
        CC = bwconncomp(above_thresh, 26);
        cluster_sizes = cellfun(@numel, CC.PixelIdxList);
        large_clusters_idx = find(cluster_sizes >= 50);

        large_clusters_mask = false(size(above_thresh));
        for k = 1:length(large_clusters_idx)
            large_clusters_mask(CC.PixelIdxList{large_clusters_idx(k)}) = true;
        end

        % Calculate TAU-SPEX %
        num_gm_voxels = sum(gm_mask(:));
        num_above_thresh = sum(large_clusters_mask(:));
        TAU_SPEX_percent = 100 * num_above_thresh / num_gm_voxels;

        % Store results
        subject_names{end+1,1} = subject_name;
        tau_spex_values(end+1,1) = TAU_SPEX_percent;

        fprintf('Processed %s: TAU-SPEX = %.2f%%\n', subject_name, TAU_SPEX_percent);

    catch ME
        fprintf('Error processing %s: %s\n', subject_name, ME.message);
    end
end

% Save results to Excel
results_table = table(subject_names, tau_spex_values, 'VariableNames', {'Subject', 'TAU_SPEX_Percent'});
output_excel = fullfile(main_dir, 'TAU_SPEX__newresults.xlsx');
writetable(results_table, output_excel);

fprintf('Done. Results saved to %s\n', output_excel);
