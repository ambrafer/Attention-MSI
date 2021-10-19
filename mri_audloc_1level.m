%% 1st level model specification, estimation and contrasts manager
% auditory localiser
% model specification based on stimuli presentation

%% Set parameters

clear
close all
clc

addpath(genpath('E:\AMBRA\UoB\Programs\spm12'));

% steps to run (0 = no; 1 = yes)
modelspec = 1; % model specification
estimate = 1;  % SPM estimation
contrasts = 1; % contrasts manager
anova = 1; % run contasts for 2nd level anova

% exp info
PID_list = {'sub-MA01'};
numDynamics = 128;
numCond = 3;

% type of images for analysis
norm = 1; % 0 = no; 1 = yes
smooth = 8; % 0 3 6 8

%% Create image lists

for isubj=1:length(PID_list)
    
    clearvars -except isubj PID_list modelspec estimate contrasts anova numDynamics numCond norm smooth
    
    PID=PID_list{isubj};
    
    % GLM info
    glm = fullfile('E:\AMBRA\UoB\Data\MAMSI_MRI',PID,'derivatives\audloc\GLM_info');
    cd(glm);
    glmlist = dir([PID '_task_aloc_GLM*.mat']);
    glmnames = {glmlist.name}';
    
    % scans
    root = fullfile('E:\AMBRA\UoB\Data\MAMSI_MRI',PID,'derivatives\prep');
    cd(root);
    
    if smooth ~= 0 && norm
        runlist = dir(['s' num2str(smooth) 'wau' PID '_task_aloc*.nii']);
    elseif smooth ~= 0 && ~norm
        runlist = dir(['s' num2str(smooth) 'au' PID '_task_aloc*.nii']);
    elseif smooth == 0 && ~norm
        runlist = dir(['au' PID '_task_aloc*.nii']);
    end
    
    numruns = length({runlist.name});
    filenames = {runlist.name}';
    disp([num2str(numruns) ' files found for 1st level fMRI model specification']);
    
    % realignment parameters
    rplist = dir(['rp_' PID '_task_aloc_bold*.txt']);
    rpnames = {rplist.name}';
    
    if norm
        matlabbatch_modelspec{1}.spm.stats.fmri_spec.dir = {fullfile('E:\AMBRA\UoB\Data\MAMSI_MRI',PID,['derivatives\audloc\anl_audloc_smooth' num2str(smooth)])};
    else
        matlabbatch_modelspec{1}.spm.stats.fmri_spec.dir = {fullfile('E:\AMBRA\UoB\Data\MAMSI_MRI',PID,['derivatives\audloc\anl_audloc_native_smooth' num2str(smooth)])};
    end
    
    anldir = matlabbatch_modelspec{1}.spm.stats.fmri_spec.dir{:};
    
    matlabbatch_modelspec{1}.spm.stats.fmri_spec.timing.units = 'secs';
    matlabbatch_modelspec{1}.spm.stats.fmri_spec.timing.RT = 2.8;
    matlabbatch_modelspec{1}.spm.stats.fmri_spec.timing.fmri_t = 38;
    matlabbatch_modelspec{1}.spm.stats.fmri_spec.timing.fmri_t0 = 19;
    
    for iRuns = 1:numruns
        for iDynams = 1:numDynamics
            
            imageList{iRuns,iDynams} = [fullfile(root,filenames{iRuns}) ',' num2str(iDynams)];
            
        end
    end
    
    for iRuns = 1:numruns
        matlabbatch_modelspec{1}.spm.stats.fmri_spec.sess(iRuns).scans = ...
            imageList(iRuns,:)';
        
        matlabbatch_modelspec{1}.spm.stats.fmri_spec.sess(iRuns).cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
        matlabbatch_modelspec{1}.spm.stats.fmri_spec.sess(iRuns).multi = {fullfile('E:\AMBRA\UoB\Data\MAMSI_MRI', PID, 'derivatives\audloc\GLM_info', glmnames{iRuns})};
        matlabbatch_modelspec{1}.spm.stats.fmri_spec.sess(iRuns).regress = struct('name', {}, 'val', {});
        matlabbatch_modelspec{1}.spm.stats.fmri_spec.sess(iRuns).multi_reg = {fullfile(root, rpnames{iRuns})};
        matlabbatch_modelspec{1}.spm.stats.fmri_spec.sess(iRuns).hpf = 128;
    end
    
    matlabbatch_modelspec{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
    matlabbatch_modelspec{1}.spm.stats.fmri_spec.bases.hrf.derivs = [1 0];
    matlabbatch_modelspec{1}.spm.stats.fmri_spec.volt = 1;
    matlabbatch_modelspec{1}.spm.stats.fmri_spec.global = 'None';
    matlabbatch_modelspec{1}.spm.stats.fmri_spec.mthresh = 0.8;
    matlabbatch_modelspec{1}.spm.stats.fmri_spec.mask = {''};
    matlabbatch_modelspec{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
    
    %% Run 1st level model specification
    if modelspec
        spm_jobman('run', matlabbatch_modelspec);
    end
    
    %% Estimate SPM
    if estimate
        cd(anldir);
        load SPM;
        spm_spm(SPM);
    end
    
    %% Run contrasts
    
    % load contrasts vectors
    cd('E:\AMBRA\UoB\Data\MAMSI_MRI\contrasts');    
    % replicate and scale across sessions 'replsc', or just replicate 'repl'?
    sessScale = 'repl';    
    load audloc_contrasts_2lev_anova_sessions_4; % does not really matter which one we select, we will anyway just keep the first run (and use 'repl')
    
    for icon = 1:numCond        
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
end
