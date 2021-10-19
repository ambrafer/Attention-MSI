%% 1st level model specification, estimation and contrasts manager
% main experiment
% model specification based on stimuli presentation

%% Set parameters

clear
close all
clc

% steps to run (0 = no; 1 = yes)
modelspec = 1; % model specification
estimate = 1;  % SPM estimation
contrasts = 1; % contrasts manager
anova = 1; % run contasts for 2nd level anova

% exp info
PID = 'sub-MA01';
numDynamics_ve = 276;
numCond = 38;

% type of images for analysis
norm = 1; % 0 = no; 1 = yes
smooth = 8; % 0 3 6 8

%% Create image lists

% GLM info
glm = fullfile('E:\Data\MAMSI_MRI',PID,'derivatives\ve\GLM_info');
cd(glm);
glmlist = dir([PID '_task_ve_stim_GLM*.mat']);
glmnames = {glmlist.name}';

% scans
root = fullfile('E:\Data\MAMSI_MRI',PID,'derivatives\prep');
cd(root);

if smooth ~= 0 && norm
    runlist_ve = dir(['s' num2str(smooth) 'wau' PID '_task_ve*.nii']);
elseif smooth ~= 0 && ~norm
    runlist_ve = dir(['s' num2str(smooth) 'au' PID '_task_ve*.nii']);
elseif smooth == 0 && ~norm
    runlist_ve = dir(['au' PID '_task_ve*.nii']);
end

numruns_ve = length({runlist_ve.name});
filenames_ve = {runlist_ve.name}';
disp([num2str(numruns_ve) ' files found for 1st level fMRI model specification']);

% realignment parameters
rplist = dir(['rp_' PID '_task_ve_bold*.txt']);
rpnames = {rplist.name}';

if norm
    matlabbatch_modelspec{1}.spm.stats.fmri_spec.dir = {fullfile('E:\Data\MAMSI_MRI',PID,['derivatives\ve\anl_ve_smooth' num2str(smooth)])};
else
    matlabbatch_modelspec{1}.spm.stats.fmri_spec.dir = {fullfile('E:\Data\MAMSI_MRI',PID,['derivatives\ve\anl_ve_native_smooth' num2str(smooth)])};
end

anldir = matlabbatch_modelspec{1}.spm.stats.fmri_spec.dir{:};

matlabbatch_modelspec{1}.spm.stats.fmri_spec.timing.units = 'secs';
matlabbatch_modelspec{1}.spm.stats.fmri_spec.timing.RT = 2.8;
matlabbatch_modelspec{1}.spm.stats.fmri_spec.timing.fmri_t = 38;
matlabbatch_modelspec{1}.spm.stats.fmri_spec.timing.fmri_t0 = 19;

for iRuns = 1:numruns_ve
    for iDynams = 1:numDynamics_ve
        
        imageList_ve{iRuns,iDynams} = [fullfile(root,filenames_ve{iRuns}) ',' num2str(iDynams)];
        
    end
end

for iRuns = 1:numruns_ve
    matlabbatch_modelspec{1}.spm.stats.fmri_spec.sess(iRuns).scans = ...
        imageList_ve(iRuns,:)';
    
    matlabbatch_modelspec{1}.spm.stats.fmri_spec.sess(iRuns).cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
    matlabbatch_modelspec{1}.spm.stats.fmri_spec.sess(iRuns).multi = {fullfile('E:\Data\MAMSI_MRI', PID, 'derivatives\ve\GLM_info', glmnames{iRuns})};
    matlabbatch_modelspec{1}.spm.stats.fmri_spec.sess(iRuns).regress = struct('name', {}, 'val', {});
    matlabbatch_modelspec{1}.spm.stats.fmri_spec.sess(iRuns).multi_reg = {fullfile(root, rpnames{iRuns})};
    matlabbatch_modelspec{1}.spm.stats.fmri_spec.sess(iRuns).hpf = 128;
end

matlabbatch_modelspec{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
matlabbatch_modelspec{1}.spm.stats.fmri_spec.bases.hrf.derivs = [1 0];
matlabbatch_modelspec{1}.spm.stats.fmri_spec.volt = 1;
matlabbatch_modelspec{1}.spm.stats.fmri_spec.global = 'None';
if smooth ~= 0
    matlabbatch_modelspec{1}.spm.stats.fmri_spec.mthresh = 0.8;
    matlabbatch_modelspec{1}.spm.stats.fmri_spec.mask = {''};
else
    % unsmoothed images with smoothed mask (for MVPA)
    matlabbatch_modelspec{1}.spm.stats.fmri_spec.mthresh = 0;
    matlabbatch_modelspec{1}.spm.stats.fmri_spec.mask = {fullfile('E:\Data\MAMSI_MRI', PID, 'derivatives\ve\anl_ve_native_smooth8\mask.nii,1')};
end
matlabbatch_modelspec{1}.spm.stats.fmri_spec.cvi = 'AR(1)';

%% Run 1st level model specification
if modelspec
    spm_jobman('run', matlabbatch_modelspec);
end

%% Estimate SPM

cd(anldir)
load SPM;

if estimate
    spm_spm(SPM);
end

%% Run contrasts (if model is specified for univariate analysis)

if smooth ~= 0 && norm    
    % load contrasts vectors
    cd('E:\Data\MAMSI_MRI\contrasts');    
    % replicate and scale across sessions 'replsc', or just replicate 'repl'?
    sessScale = 'repl';    
    load ve_stim_contrasts_2lev_anova_sessions_14;
    
    for icon = 1:numCond-2 % no Attentional Cues? put -2
        
        matlabbatch_contrast{1}.spm.stats.con.spmmat = {fullfile(...
            anldir,'SPM.mat')};
        matlabbatch_contrast{1}.spm.stats.con.consess{icon}.tcon.name = names{icon};
        matlabbatch_contrast{1}.spm.stats.con.consess{icon}.tcon.weights = ...
            con_tot{icon}(1:numCond*2); % account for temporal derivative
        matlabbatch_contrast{1}.spm.stats.con.consess{icon}.tcon.sessrep = sessScale;
        
    end
    
    matlabbatch_contrast{1}.spm.stats.con.delete = 1;
    
    if contrasts
        spm_jobman('run', matlabbatch_contrast);
    end    
    cd(anldir);    
end
