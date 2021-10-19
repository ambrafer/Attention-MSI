%% 2nd level model specification, estimation and contrasts manager
% audloc experiment: 3 A location

clear
close all
clc

% steps to run (0 = no; 1 = yes)
modelspec = 1; % model specification
estimate = 1;  % SPM estimation
contrasts = 1; % contrasts manager

% number of contrasts
numcon = {'01'; '02'; '03'};
numcontot = size(numcon,1);

for icon = 1:numcontot
    
    matlabbatch_modelspec{1}.spm.stats.factorial_design.dir = {'E:\Data\MAMSI_MRI\group\audloc'};
    matlabbatch_modelspec{1}.spm.stats.factorial_design.des.anova.icell(icon).scans = {
        ['E:\Data\MAMSI_MRI\sub-MA01\derivatives\audloc\anl_audloc_smooth8\con_00' numcon{icon} '.nii,1']
        ['E:\Data\MAMSI_MRI\sub-MA02\derivatives\audloc\anl_audloc_smooth8\con_00' numcon{icon} '.nii,1']
        ['E:\Data\MAMSI_MRI\sub-MA03\derivatives\audloc\anl_audloc_smooth8\con_00' numcon{icon} '.nii,1']
        ['E:\Data\MAMSI_MRI\sub-MA04\derivatives\audloc\anl_audloc_smooth8\con_00' numcon{icon} '.nii,1']
        ['E:\Data\MAMSI_MRI\sub-MA05\derivatives\audloc\anl_audloc_smooth8\con_00' numcon{icon} '.nii,1']
        ['E:\Data\MAMSI_MRI\sub-MA06\derivatives\audloc\anl_audloc_smooth8\con_00' numcon{icon} '.nii,1']
        ['E:\Data\MAMSI_MRI\sub-MA07\derivatives\audloc\anl_audloc_smooth8\con_00' numcon{icon} '.nii,1']
        ['E:\Data\MAMSI_MRI\sub-MA08\derivatives\audloc\anl_audloc_smooth8\con_00' numcon{icon} '.nii,1']
        ['E:\Data\MAMSI_MRI\sub-MA09\derivatives\audloc\anl_audloc_smooth8\con_00' numcon{icon} '.nii,1']
        ['E:\Data\MAMSI_MRI\sub-MA010\derivatives\audloc\anl_audloc_smooth8\con_00' numcon{icon} '.nii,1']
        ['E:\Data\MAMSI_MRI\sub-MA011\derivatives\audloc\anl_audloc_smooth8\con_00' numcon{icon} '.nii,1']
        ['E:\Data\MAMSI_MRI\sub-MA012\derivatives\audloc\anl_audloc_smooth8\con_00' numcon{icon} '.nii,1']
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

if modelspec
    % specify model
    spm_jobman('run', matlabbatch_modelspec(1));
end

if estimate    
    % estimate SPM
    cd(anldir);
    load SPM;
    spm_spm(SPM);
end

% define contrasts
load('E:\Data\MAMSI_MRI\contrasts\audloc_contrasts_2lev_anova_sessions_4');

% contrasts names
con(:,1) = {'baseline'; 'all'; 'L>R'; 'R>L'; 'L>C'; 'C>L'; 'R>C'; 'C>R'};
    
% contrasts values
baseline = eye(3);
all = ones(1,3);
L = regexp(names,'L', 'once');
L = double(cellfun(@(x) ~isempty(x),L));
C = regexp(names,'C', 'once');
C = double(cellfun(@(x) ~isempty(x),C));
R = regexp(names,'R', 'once');
R = double(cellfun(@(x) ~isempty(x),R));

con(:,2) = {baseline; all; L - R; R - L; L - C; C - L; R - C; C - R};

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

cd('E:\Data\MAMSI_MRI\group\audloc');
