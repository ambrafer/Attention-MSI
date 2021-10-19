%% 2nd level model specification, estimation and contrasts manager
% main experiment: 3 A location x 3 V locations x 2 Attention x 2 Response

clear
close all
clc

addpath(genpath('E:\AMBRA\UoB\Programs\spm12\'));

% steps to run (0 = no; 1 = yes)
modelspec = 1; % model specification
estimate = 1;  % SPM estimation
contrasts = 1; % contrasts manager

% number of contrasts
numcon = {'01'; '02'; '03'; '04'; '05'; '06'; '07'; '08'; '09'; '10'; '11'; '12'; ...
    '13'; '14'; '15'; '16'; '17'; '18'; '19'; '20'; '21'; '22'; '23'; '24'; '25'; ...
    '26'; '27'; '28'; '29'; '30'; '31'; '32'; '33'; '34'; '35'; '36'};
numcontot = size(numcon,1);

for icon = 1:numcontot
    matlabbatch_modelspec{1}.spm.stats.factorial_design.dir = {'E:\AMBRA\UoB\Data\MAMSI_MRI\group\main'};
    matlabbatch_modelspec{1}.spm.stats.factorial_design.des.anova.icell(icon).scans = {
        ['E:\AMBRA\UoB\Data\MAMSI_MRI\sub-MA01\derivatives\ve\anl_ve_stim_smooth8\con_2lev_anova\con_00' numcon{icon} '.nii,1']
        ['E:\AMBRA\UoB\Data\MAMSI_MRI\sub-MA02\derivatives\ve\anl_ve_stim_smooth8\con_2lev_anova\con_00' numcon{icon} '.nii,1']
        ['E:\AMBRA\UoB\Data\MAMSI_MRI\sub-MA03\derivatives\ve\anl_ve_stim_smooth8\con_2lev_anova\con_00' numcon{icon} '.nii,1']
        ['E:\AMBRA\UoB\Data\MAMSI_MRI\sub-MA04\derivatives\ve\anl_ve_stim_smooth8\con_2lev_anova\con_00' numcon{icon} '.nii,1']
        ['E:\AMBRA\UoB\Data\MAMSI_MRI\sub-MA05\derivatives\ve\anl_ve_stim_smooth8\con_2lev_anova\con_00' numcon{icon} '.nii,1']
        ['E:\AMBRA\UoB\Data\MAMSI_MRI\sub-MA06\derivatives\ve\anl_ve_stim_smooth8\con_2lev_anova\con_00' numcon{icon} '.nii,1']
        ['E:\AMBRA\UoB\Data\MAMSI_MRI\sub-MA07\derivatives\ve\anl_ve_stim_smooth8\con_2lev_anova\con_00' numcon{icon} '.nii,1']
        ['E:\AMBRA\UoB\Data\MAMSI_MRI\sub-MA08\derivatives\ve\anl_ve_stim_smooth8\con_2lev_anova\con_00' numcon{icon} '.nii,1']
        ['E:\AMBRA\UoB\Data\MAMSI_MRI\sub-MA08\derivatives\ve\anl_ve_stim_smooth8\con_2lev_anova\con_00' numcon{icon} '.nii,1']
        ['E:\AMBRA\UoB\Data\MAMSI_MRI\sub-MA010\derivatives\ve\anl_ve_stim_smooth8\con_2lev_anova\con_00' numcon{icon} '.nii,1']
        ['E:\AMBRA\UoB\Data\MAMSI_MRI\sub-MA011\derivatives\ve\anl_ve_stim_smooth8\con_2lev_anova\con_00' numcon{icon} '.nii,1']
        ['E:\AMBRA\UoB\Data\MAMSI_MRI\sub-MA012\derivatives\ve\anl_ve_stim_smooth8\con_2lev_anova\con_00' numcon{icon} '.nii,1']
        };
    matlabbatch_modelspec{1}.spm.stats.factorial_design.des.anova.dept = 1;
    matlabbatch_modelspec{1}.spm.stats.factorial_design.des.anova.variance = 1;
    matlabbatch_modelspec{1}.spm.stats.factorial_design.des.anova.gmsca = 0;
    matlabbatch_modelspec{1}.spm.stats.factorial_design.des.anova.ancova = 0;
    matlabbatch_modelspec{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch_modelspec{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch_modelspec{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
    matlabbatch_modelspec{1}.spm.stats.factorial_design.masking.im = 1;
    matlabbatch_modelspec{1}.spm.stats.factorial_design.masking.em = {''};
    matlabbatch_modelspec{1}.spm.stats.factorial_design.globalc.g_omit = 1;
    matlabbatch_modelspec{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
    matlabbatch_modelspec{1}.spm.stats.factorial_design.globalm.glonorm = 1;
end
anldir = matlabbatch_modelspec{1}.spm.stats.factorial_design.dir{:};

%% RUN

if modelspec==1
    % specify model
    spm_jobman('run', matlabbatch_modelspec(1));
end

if estimate==1
    % estimate SPM
    cd(anldir);
    load SPM;
    spm_spm(SPM);
end

% define contrasts

load('E:\AMBRA\UoB\Data\MAMSI_MRI\contrasts\ve_stim_contrasts_2lev_anova_sessions_14');
names = names(1:36); % discard 2 attentional cues

% contrasts names
con(:,1) = {'baseline'; 'all'; 'AttA_baseline'; 'baseline_AttA'; 'AttV_baseline'; 'baseline_AttV'; 'AttA_AttV'; 'AttV_AttA'; ...
    
    'RespA_baseline'; 'baseline_RespA'; 'RespV_baseline'; 'baseline_RespV'; 'RespA_RespV'; 'RespV_RespA'; ...
    
    'Invalid_Valid'; 'Valid_Invalid'; ...
    
    'AttARespA_AttARespV'; 'AttARespA_AttVRespA'; 'AttARespA_AttVRespV'; ...
    'AttARespV_AttARespA'; 'AttARespV_AttVRespA'; 'AttARespV_AttVRespV'; ...
    'AttVRespA_AttARespA'; 'AttVRespA_AttARespV'; 'AttVRespA_AttVRespV'; ...
    'AttVRespV_AttARespA'; 'AttVRespV_AttARespV'; 'AttVRespV_AttVRespA'; ...

    'VisualL_VisualR'; 'VisualR_VisualL'; 'AuditoryL_AuditoryR'; 'AuditoryR_AuditoryL'; ...
    
    'Congruent_baseline'; 'Incongruent_baseline'; 'Incongruent_Congruent'; 'Congruent_Incongruent'; ...
    
    'Incongruent_Congruent_repA_repV'; 'Incongruent_Congruent_repV_repA'; ...
    
    'AL_AR_AttA'; 'AL_AR_AttV'; 'AL_AR_RespA'; 'AL_AR_RespV'; ...
    'AR_AL_AttA'; 'AR_AL_AttV'; 'AR_AL_RespA'; 'AR_AL_RespV'; ...
    
    'AL_AR_AttA_AttV'; 'AL_AR_AttV_AttA'; 'AL_AR_RespA_RespV'; 'AL_AR_RespV_RespA'; ...
    'AR_AL_AttA_AttV'; 'AR_AL_AttV_AttA'; 'AR_AL_RespA_RespV'; 'AR_AL_RespV_RespA'; ...
    
    'AL_AR_AttARespA'; 'AL_AR_AttARespV'; 'AL_AR_AttVRespA'; 'AL_AR_AttVRespV'; ...
    'AR_AL_AttARespA'; 'AR_AL_AttARespV'; 'AR_AL_AttVRespA'; 'AR_AL_AttVRespV'; ...
    
    'AL_AR_AttARespA_AttVRespA'; 'AL_AR_AttVRespA_AttARespA'; 'AL_AR_AttVRespV_AttARespV'; 'AL_AR_AttARespV_AttVRespV'; ...
    'AR_AL_AttARespA_AttVRespA'; 'AR_AL_AttVRespA_AttARespA'; 'AR_AL_AttVRespV_AttARespV'; 'AR_AL_AttARespV_AttVRespV'};

% contrasts values
baseline = eye(36);

all = ones(1,36);

AttA = regexp(names,'AttA_', 'once');
AttA = double(cellfun(@(x) ~isempty(x),AttA));

AttV = regexp(names,'AttV_', 'once');
AttV = double(cellfun(@(x) ~isempty(x),AttV));

RespA = regexp(names,'RespA', 'once');
RespA = double(cellfun(@(x) ~isempty(x),RespA));

RespV = regexp(names,'RespV', 'once');
RespV = double(cellfun(@(x) ~isempty(x),RespV));

AttARespA = regexp(names,'AttA_RespA', 'once');
AttARespA = double(cellfun(@(x) ~isempty(x),AttARespA));

AttARespV = regexp(names,'AttA_RespV', 'once');
AttARespV = double(cellfun(@(x) ~isempty(x),AttARespV));

AttVRespA = regexp(names,'AttV_RespA', 'once');
AttVRespA = double(cellfun(@(x) ~isempty(x),AttVRespA));

AttVRespV = regexp(names,'AttV_RespV', 'once');
AttVRespV = double(cellfun(@(x) ~isempty(x),AttVRespV));

Valid = double(or(AttARespA,AttVRespV));

Invalid = double(or(AttARespV,AttVRespA));

VL = regexp(names,'VL', 'once');
VL = double(cellfun(@(x) ~isempty(x),VL));

VC = regexp(names,'VC', 'once');
VC = double(cellfun(@(x) ~isempty(x),VC));

VR = regexp(names,'VR', 'once');
VR = double(cellfun(@(x) ~isempty(x),VR));

AL = regexp(names,'AL', 'once');
AL = double(cellfun(@(x) ~isempty(x),AL));

AC = regexp(names,'AC', 'once');
AC = double(cellfun(@(x) ~isempty(x),AC));

AR = regexp(names,'AR', 'once');
AR = double(cellfun(@(x) ~isempty(x),AR));

Congruent = repmat([1 0 0 0 1 0 0 0 1],1,4);
Congruent_repA = repmat([1 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0 0 0],1,2);
Congruent_repV = repmat([0 0 0 0 0 0 0 0 0 1 0 0 0 1 0 0 0 1],1,2);

Incongruent = double(~Congruent);
Incongruent_repA = repmat([0 1 1 1 0 1 1 1 0 0 0 0 0 0 0 0 0 0],1,2);
Incongruent_repV = repmat([0 0 0 0 0 0 0 0 0 0 1 1 1 0 1 1 1 0],1,2);

AL_AttA = and(AL, AttA);

AL_AttV = and(AL, AttV);

AR_AttA = and(AR, AttA);

AR_AttV = and(AR, AttV);

AL_RespA = and(AL, RespA);

AL_RespV = and(AL, RespV);

AR_RespA = and(AR, RespA);

AR_RespV = and(AR, RespV);

% main effect of Attention and Response
AL_AR_AttA = AL_AttA-AR_AttA;
AL_AR_AttV = AL_AttV-AR_AttV;
AL_AR_RespA = AL_RespA-AR_RespA;
AL_AR_RespV = AL_RespV-AR_RespV;

AR_AL_AttA = AR_AttA-AL_AttA;
AR_AL_AttV = AR_AttV-AL_AttV;
AR_AL_RespA = AR_RespA-AL_RespA;
AR_AL_RespV = AR_RespV-AL_RespV;

AL_AR_AttA_AttV = AL_AR_AttA-AL_AR_AttV;
AL_AR_AttV_AttA = AL_AR_AttV-AL_AR_AttA;
AL_AR_RespA_RespV = AL_AR_RespA-AL_AR_RespV;
AL_AR_RespV_RespA = AL_AR_RespV-AL_AR_RespA;

AR_AL_AttA_AttV = AR_AL_AttA-AR_AL_AttV;
AR_AL_AttV_AttA = AR_AL_AttV-AR_AL_AttA;
AR_AL_RespA_RespV = AR_AL_RespA-AR_AL_RespV;
AR_AL_RespV_RespA = AR_AL_RespV-AR_AL_RespA;

% interaction Attention and Response (validity)
AL_AttARespA = and(AL_AttA,AL_RespA);
AL_AttARespV = and(AL_AttA,AL_RespV);
AL_AttVRespA = and(AL_AttV,AL_RespA);
AL_AttVRespV = and(AL_AttV,AL_RespV);

AR_AttARespA = and(AR_AttA,AR_RespA);
AR_AttARespV = and(AR_AttA,AR_RespV);
AR_AttVRespA = and(AR_AttV,AR_RespA);
AR_AttVRespV = and(AR_AttV,AR_RespV);

AL_AR_AttARespA = AL_AttARespA-AR_AttARespA;
AL_AR_AttARespV = AL_AttARespV-AR_AttARespV;
AL_AR_AttVRespA = AL_AttVRespA-AR_AttVRespA;
AL_AR_AttVRespV = AL_AttVRespV-AR_AttVRespV;

AL_AR_AttARespA_AttVRespA = AL_AR_AttARespA-AL_AR_AttVRespA;
AL_AR_AttVRespA_AttARespA = AL_AR_AttVRespA-AL_AR_AttARespA;
AL_AR_AttVRespV_AttARespV = AL_AR_AttVRespV-AL_AR_AttARespV;
AL_AR_AttARespV_AttVRespV = AL_AR_AttARespV-AL_AR_AttVRespV;

AR_AL_AttARespA = AR_AttARespA-AL_AttARespA;
AR_AL_AttARespV = AR_AttARespV-AL_AttARespV;
AR_AL_AttVRespA = AR_AttVRespA-AL_AttVRespA;
AR_AL_AttVRespV = AR_AttVRespV-AL_AttVRespV;

AR_AL_AttARespA_AttVRespA = AR_AL_AttARespA-AR_AL_AttVRespA;
AR_AL_AttVRespA_AttARespA = AR_AL_AttVRespA-AR_AL_AttARespA;
AR_AL_AttVRespV_AttARespV = AR_AL_AttVRespV-AR_AL_AttARespV;
AR_AL_AttARespV_AttVRespV = AR_AL_AttARespV-AR_AL_AttVRespV;

con(:,2) = {baseline; all; AttA; -AttA; AttV; -AttV; AttA - AttV; AttV - AttA; ...
    
    RespA; -RespA; RespV; -RespV; RespA - RespV; RespV - RespA; ...
    
    Invalid - Valid; Valid - Invalid; ...
    
    AttARespA - AttARespV; AttARespA - AttVRespA; AttARespA - AttVRespV; ...
    AttARespV - AttARespA; AttARespV - AttVRespA; AttARespV - AttVRespV; ...
    AttVRespA - AttARespA; AttVRespA - AttARespV; AttVRespA - AttVRespV; ...
    AttVRespV - AttARespA; AttVRespV - AttARespV; AttVRespV - AttVRespA; ...
    
    VL - VR; VR - VL; AL - AR; AR - AL; ...
    
    Congruent; Incongruent; 0.5*Incongruent - Congruent; Congruent - 0.5*Incongruent; ...
    
    (0.5*Incongruent_repA - Congruent_repA) - (0.5*Incongruent_repV - Congruent_repV); ...
    (0.5*Incongruent_repV - Congruent_repV) - (0.5*Incongruent_repA - Congruent_repA); ...
    
    AL_AR_AttA; AL_AR_AttV; AL_AR_RespA; AL_AR_RespV; ...
    AR_AL_AttA; AR_AL_AttV; AR_AL_RespA; AR_AL_RespV; ...
    
    AL_AR_AttA_AttV; AL_AR_AttV_AttA; AL_AR_RespA_RespV; AL_AR_RespV_RespA; ...
    AR_AL_AttA_AttV; AR_AL_AttV_AttA; AR_AL_RespA_RespV; AR_AL_RespV_RespA; ...
    
    AL_AR_AttARespA; AL_AR_AttARespV; AL_AR_AttVRespA; AL_AR_AttVRespV; ...
    AR_AL_AttARespA; AR_AL_AttARespV; AR_AL_AttVRespA; AR_AL_AttVRespV; ...
    
    AL_AR_AttARespA_AttVRespA; AL_AR_AttVRespA_AttARespA; AL_AR_AttVRespV_AttARespV; AL_AR_AttARespV_AttVRespV; ...
    AR_AL_AttARespA_AttVRespA; AR_AL_AttVRespA_AttARespA; AR_AL_AttVRespV_AttARespV; AR_AL_AttARespV_AttVRespV};  

matlabbatch_contrast{1}.spm.stats.con.spmmat = {fullfile(...
    anldir,'SPM.mat')};

for icon = 1:size(con,1)
    
    if icon == 1
        matlabbatch_contrast{1}.spm.stats.con.consess{icon}.fcon.name = con{icon,1};
        matlabbatch_contrast{1}.spm.stats.con.consess{icon}.fcon.weights = con{icon,2};
        matlabbatch_contrast{1}.spm.stats.con.consess{icon}.fcon.sessrep = 'none';        
    else        
        matlabbatch_contrast{1}.spm.stats.con.consess{icon}.tcon.name = con{icon,1};
        matlabbatch_contrast{1}.spm.stats.con.consess{icon}.tcon.weights = con{icon,2};
        matlabbatch_contrast{1}.spm.stats.con.consess{icon}.tcon.sessrep = 'none';
        
    end
end
matlabbatch_contrast{1}.spm.stats.con.delete = 1;

if contrasts
    spm_jobman('run', matlabbatch_contrast);
end

cd('E:\AMBRA\UoB\Data\MAMSI_MRI\group\main');