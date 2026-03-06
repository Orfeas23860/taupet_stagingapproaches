addpath('C:\Users\ovourkas\Documents\CenTauR_Pipelines_ROIs\ROIs\CenTauR_to_MCALT')

% Define input directories
pet_dir = 'C:\Users\ovourkas\Documents\Projects\FDG_vs_MRI_VisualTauRead\nifti_TAU_followup';   
t1_dir = 'C:\Users\ovourkas\Documents\Projects\FDG_vs_MRI_VisualTauRead\nifti_MRI_followup';     

% Get list of subfolders (subjects) in PET and T1 folders
pet_folders_struct = dir(pet_dir);
pet_folders = {pet_folders_struct([pet_folders_struct.isdir] & ~startsWith({pet_folders_struct.name}, '.')).name};

t1_folders_struct = dir(t1_dir);
t1_folders = {t1_folders_struct([t1_folders_struct.isdir] & ~startsWith({t1_folders_struct.name}, '.')).name};

nSubjects = length(pet_folders);

% Output base directory
outDir = 'C:\Users\ovourkas\Documents\Projects\FDG_vs_MRI_VisualTauRead\Output_Processing';

if ~exist(outDir,'dir')
    mkdir(outDir);
end

for i = 412:nSubjects
    % Get full paths to subject's PET and T1 subfolders
    pet_path = fullfile(pet_dir, pet_folders{i});
    t1_path = fullfile(t1_dir, t1_folders{i});
    
    % Find first .nii or .nii.gz file in each folder
    pet_file = dir(fullfile(pet_path, '*.nii*'));
    t1_file = dir(fullfile(t1_path, '*.nii*'));
    
    pet_nii = fullfile(pet_path, pet_file(1).name);
    t1_nii = fullfile(t1_path, t1_file(1).name);
    
    % Create subject-specific output folder to avoid overwrites
    subject_outDir = fullfile(outDir, pet_folders{i});
    if ~exist(subject_outDir, 'dir')
        mkdir(subject_outDir);
    end
    
    % Copy files to subject output folder
    [~, t1name, ext1] = fileparts(t1_nii);
    [~, petname, ext2] = fileparts(pet_nii);
    
    copyfile(t1_nii, fullfile(subject_outDir, [t1name ext1]));
    copyfile(pet_nii, fullfile(subject_outDir, [petname ext2]));
    
    % Files to fix origin (example placeholders)
    files_to_fix = { ...
        fullfile(subject_outDir, [t1name ext1]), ...
        fullfile(subject_outDir, [petname ext2]) ...
    };

for i = 1:numel(files_to_fix)
    V = spm_vol(files_to_fix{i});
    img = spm_read_vols(V);
    center_voxel = (V.dim + 1) / 2;         
    origin_world = -V.mat(1:3,1:3) * center_voxel';
    V.mat(1:3,4) = origin_world;
    [folder, name, ext] = fileparts(files_to_fix{i});
    V.fname = fullfile(folder, [name '_s' ext]);
    spm_write_vol(V, img);
    fprintf('Saved origin-centered copy: %s\n', V.fname);
end

ROIs = struct('name', {}, 'mask', {});
ROIs(1).name = 'Universal';
ROIs(1).mask = which('MCALT_CenTauR.nii');
ROIs(2).name = 'MesialTemporal';
ROIs(2).mask = which('MCALT_Mesial_CenTauR.nii');
ROIs(3).name = 'MetaTemporal';
ROIs(3).mask = which('MCALT_Meta_CenTauR.nii');
ROIs(4).name = 'TemporoParietal';
ROIs(4).mask = which('MCALT_TP_CenTauR.nii');
ROIs(5).name = 'Frontal';
ROIs(5).mask = which('MCALT_Frontal_CenTauR.nii');

centaur_constants = which('centaur_constants.csv');
if(isempty(centaur_constants))
    error('Centaur table of constants was not found in your matlab path.')
end
c = table2struct(readtable(centaur_constants));

pipeline_constants = which('Pipeline_constants_SPM12.csv');
if(isempty(pipeline_constants))
    error('Pipeline table of constants was not found in your matlab path.')
end
cPipe = table2struct(readtable(pipeline_constants));

ref_mask = which('MCALT_voi_CerebGry_tau_2mm.nii');
if(isempty(ref_mask))
    error('Could not find reference mask file');
end

MCALT_tpm = which('MCALT_tpm.nii');
if(isempty(MCALT_tpm))
    error('Could not find MCALT_tpm.nii');
end

MCALT_T1 = which('MCALT_T1.nii');
if(isempty(MCALT_T1))
    error('Could not find MCALT_T1.nii');
end

tracer = 'FTP';

supported_tracers = {'RO948','FTP','MK6240','GTP1','PM-PBB3','PI2620'};
if ~ismember(tracer, supported_tracers)
    error(['''' tracer '''' ' is not a supported tracer. Supported tracers: ' strjoin(supported_tracers, ', ')]);
end

owd = pwd;
cd(outDir);
cleanupVar = onCleanup(@()cd(owd));

disp('Running SPM12 Centaur pipeline');

clear matlabbatch

%% coreg PET to T1 (and reslice it)

matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {fullfile(subject_outDir, [t1name '_s' ext1])}; % Use origin-centered T1
matlabbatch{1}.spm.spatial.coreg.estwrite.source = {fullfile(subject_outDir, [petname '_s' ext2])}; % Use origin-centered PET
matlabbatch{1}.spm.spatial.coreg.estwrite.other = {''};
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = ...
   [0.0200 0.0200 0.0200 0.0010 0.0010 0.0010 0.0100 0.0100 0.0100 0.0010 0.0010 0.0010];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';

spm_jobman('run', matlabbatch)

clear matlabbatch

%% segment T1

matlabbatch = [];
matlabbatch{1}.spm.spatial.preproc.channel.vols = {fullfile(subject_outDir, [t1name '_s' ext1])};
matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001;
matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
matlabbatch{1}.spm.spatial.preproc.channel.write = [0 1];

% Template path
tpm = fullfile(spm('Dir'),'tpm','TPM.nii');
for i=1:6
    matlabbatch{1}.spm.spatial.preproc.tissue(i).tpm = {[tpm ',' num2str(i)]};
end

% Tissue settings
matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1; matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 0]; matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 1; matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 0]; matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2; matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 0]; matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3; matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [0 0]; matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4; matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [0 0]; matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2; matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [0 0]; matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];

% Warping and cleanup
matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
matlabbatch{1}.spm.spatial.preproc.warp.write = [1 1];

% Initialize with affine to MCALT template
disp('Initializing subject->template registration with spm_coreg');
t1_file = fullfile(subject_outDir, [t1name '_s' ext1]);
regInit = spm_coreg(t1_file, MCALT_T1);
segAffine = spm_matrix(regInit);
matlabbatch{1}.spm.spatial.preproc.warp.Affine = segAffine;

% Run
spm_jobman('initcfg');
spm_jobman('run', matlabbatch);

clear matlabbatch

%% Transform ROIs to T1 space

sn_mat = dir(fullfile(subject_outDir, 'iy_*.nii'));
sn_mat = fullfile(subject_outDir, sn_mat(1).name);

all_mask_paths = [{ROIs.mask}, {ref_mask}]'; 

    matlabbatch{1}.spm.util.defs.comp{1}.def = {sn_mat};
    matlabbatch{1}.spm.util.defs.out{1}.pull.fnames = all_mask_paths;
    matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.saveusr = {char(subject_outDir)};
    matlabbatch{1}.spm.util.defs.out{1}.pull.interp = -1;
    matlabbatch{1}.spm.util.defs.out{1}.pull.mask = 0;
    matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
    matlabbatch{1}.spm.util.defs.out{1}.pull.prefix = 'w';

    spm_jobman('run', matlabbatch);
    clear matlabbatch

%%Transform ROIs to T1 space

ref_name = 'wMCALT_voi_CerebGry_tau_2mm.nii';

roi_files = dir(fullfile(subject_outDir, 'w*.nii'));

% Exclude reference region from ROIs
roi_files = roi_files(~strcmp({roi_files.name}, ref_name));

% Load ROI masks
for i = 1:numel(roi_files)
    roi_path = fullfile(subject_outDir, roi_files(i).name);
    ROIs(i).mask = roi_path;
    ROIs(i).name = roi_files(i).name;
    ROIs(i).hdr = spm_vol(roi_path);
    ROIs(i).vol = spm_read_vols(ROIs(i).hdr);
end

% Load reference mask separately
wRef = fullfile(subject_outDir, ref_name);
roi_ref_hdr = spm_vol(wRef);
roi_ref_vol = spm_read_vols(roi_ref_hdr);

%% Quantify

rpet_file = dir(fullfile(subject_outDir, 'r*.nii'));
rpet_coreg = fullfile(subject_outDir, rpet_file(1).name);

rPET_hdr = spm_vol(rpet_coreg);
rPET_vol = spm_read_vols(rPET_hdr);
roi_ref_mean = mean(rPET_vol(roi_ref_vol==1));

% Unlike the Centiloid standard SPM8 pipeline, we are operating in T1 image native space. We use the segmentations to
% remove voxels considered mostly CSF. This is sometimes called "tissue sharpening".

c1_file = dir(fullfile(subject_outDir, 'c1*.nii'));
c1_hdr = spm_vol(fullfile(subject_outDir, c1_file(1).name));
c1_vol = spm_read_vols(c1_hdr);

c2_file = dir(fullfile(subject_outDir, 'c2*.nii'));
c2_hdr = spm_vol(fullfile(subject_outDir, c2_file(1).name));
c2_vol = spm_read_vols(c2_hdr);

for i = 1:numel(ROIs)
    switch ROIs(i).name
        case 'wMCALT_CenTauR.nii'
            ROIs(i).name = 'Universal';
        case 'wMCALT_Mesial_CenTauR.nii'
            ROIs(i).name = 'MesialTemporal';
        case 'wMCALT_Meta_CenTauR.nii'
            ROIs(i).name = 'MetaTemporal';
        case 'wMCALT_TP_CenTauR.nii'
            ROIs(i).name = 'TemporoParietal';
        case 'wMCALT_Frontal_CenTauR.nii'
            ROIs(i).name = 'Frontal';
    end
end

CentaurZ_c = c(strcmp({c.Units},'CentaurZ') & strcmp({c.Tracer},tracer));
Pipe_c = cPipe;

for i=1:numel(ROIs)
    ROIs(i).mean = mean(rPET_vol(ROIs(i).vol==1 & (c1_vol + c2_vol) >= 0.5));
    ROIs(i).SUVR = ROIs(i).mean / roi_ref_mean;
    fprintf('  %s * %f + %f for SPM12 to SPM8, * %f + %f for SPM8 to CentaurZ (%s)\n',ROIs(i).name,Pipe_c(strcmp({Pipe_c.Region},ROIs(i).name),:).slope,Pipe_c(strcmp({Pipe_c.Region},ROIs(i).name),:).inter,CentaurZ_c.([ROIs(i).name '_slope']),CentaurZ_c.([ROIs(i).name '_inter']),tracer);
    ROIs(i).SUVR_SPM8 = ROIs(i).SUVR * Pipe_c(strcmp({Pipe_c.Region},ROIs(i).name),:).slope + Pipe_c(strcmp({Pipe_c.Region},ROIs(i).name),:).inter;
    ROIs(i).CentaurZ = ROIs(i).SUVR_SPM8 * CentaurZ_c.([ROIs(i).name '_slope']) + CentaurZ_c.([ROIs(i).name '_inter']);

end

ROIs = rmfield(ROIs,'hdr');
ROIs = rmfield(ROIs,'vol');
ROIs = rmfield(ROIs,'mask');
ROIs = rmfield(ROIs,'mean');
ROIs = struct2table(ROIs);
output_file = fullfile(subject_outDir, 'ROIs_metrics.xlsx');
writetable(ROIs, output_file);

clear matlabbatch

end

