%% SPM preprocessing job script
%(Always run this whole section first, then choose steps at the bottom)

clear
close all
clc

PID = 'sub-MA01';

root = fullfile('E:\Data\MAMSI_MRI',PID,'derivatives\prep');

numDynamics_audloc = 128;
numRuns_audloc = 2;

numDynamics_ve = 276;
numRuns_ve = 14;

cd(root);

runlist_audloc = dir([PID '_task_aloc*.nii']);
numruns_audloc = length({runlist_audloc.name});
filenames_audloc = {runlist_audloc.name}';
disp([num2str(numruns_audloc) ' files found for audloc preprocessing']);

runlist_ve = dir([PID '_task_ve*.nii']);
numruns_ve = length({runlist_ve.name});
filenames_ve = {runlist_ve.name}';
disp([num2str(numruns_ve) ' files found for ve preprocessing']);

for iRuns = 1:numruns_audloc
    for iDynams = 1:numDynamics_audloc
        
        imageList_audloc{iRuns,iDynams} = [fullfile(root,filenames_audloc{iRuns}) ',' num2str(iDynams)];
        
    end
end

for iRuns = 1:numruns_ve
    for iDynams = 1:numDynamics_ve
        
        imageList_ve{iRuns,iDynams} = [fullfile(root,filenames_ve{iRuns}) ',' num2str(iDynams)];
        
    end
end

anatfilename = fullfile(root,[PID '_T1w_orig.nii,1']);

% REALIGN AND UNWARP

for iRuns = 1:numruns_audloc
    matlabbatch_realignunwarp{1}.spm.spatial.realignunwarp.data(iRuns).scans = ...
        imageList_audloc(iRuns,:)';
    
    matlabbatch_realignunwarp{1}.spm.spatial.realignunwarp.data(iRuns).pmscan = '';
end

for iRuns = numruns_audloc+1:numruns_ve+numruns_audloc
    matlabbatch_realignunwarp{1}.spm.spatial.realignunwarp.data(iRuns).scans = ...
        imageList_ve(iRuns-numruns_audloc,:)';
    
    matlabbatch_realignunwarp{1}.spm.spatial.realignunwarp.data(iRuns).pmscan = '';
end

matlabbatch_realignunwarp{1}.spm.spatial.realignunwarp.eoptions.quality = 1;
matlabbatch_realignunwarp{1}.spm.spatial.realignunwarp.eoptions.sep = 4;
matlabbatch_realignunwarp{1}.spm.spatial.realignunwarp.eoptions.fwhm = 5;
matlabbatch_realignunwarp{1}.spm.spatial.realignunwarp.eoptions.rtm = 0;
matlabbatch_realignunwarp{1}.spm.spatial.realignunwarp.eoptions.einterp = 2;
matlabbatch_realignunwarp{1}.spm.spatial.realignunwarp.eoptions.ewrap = [0 0 0];
matlabbatch_realignunwarp{1}.spm.spatial.realignunwarp.eoptions.weight = '';
matlabbatch_realignunwarp{1}.spm.spatial.realignunwarp.uweoptions.basfcn = [12 12];
matlabbatch_realignunwarp{1}.spm.spatial.realignunwarp.uweoptions.regorder = 1;
matlabbatch_realignunwarp{1}.spm.spatial.realignunwarp.uweoptions.lambda = 100000;
matlabbatch_realignunwarp{1}.spm.spatial.realignunwarp.uweoptions.jm = 0;
matlabbatch_realignunwarp{1}.spm.spatial.realignunwarp.uweoptions.fot = [4 5];
matlabbatch_realignunwarp{1}.spm.spatial.realignunwarp.uweoptions.sot = [];
matlabbatch_realignunwarp{1}.spm.spatial.realignunwarp.uweoptions.uwfwhm = 4;
matlabbatch_realignunwarp{1}.spm.spatial.realignunwarp.uweoptions.rem = 1;
matlabbatch_realignunwarp{1}.spm.spatial.realignunwarp.uweoptions.noi = 5;
matlabbatch_realignunwarp{1}.spm.spatial.realignunwarp.uweoptions.expround = 'Average';
matlabbatch_realignunwarp{1}.spm.spatial.realignunwarp.uwroptions.uwwhich = [2 1];
matlabbatch_realignunwarp{1}.spm.spatial.realignunwarp.uwroptions.rinterp = 4;
matlabbatch_realignunwarp{1}.spm.spatial.realignunwarp.uwroptions.wrap = [0 0 0];
matlabbatch_realignunwarp{1}.spm.spatial.realignunwarp.uwroptions.mask = 1;
matlabbatch_realignunwarp{1}.spm.spatial.realignunwarp.uwroptions.prefix = 'u';

for iRuns = 1:numruns_audloc
    for iDynams = 1:numDynamics_audloc
        [path, name, ext] = fileparts(imageList_audloc{iRuns,iDynams});
        imageList_audloc{iRuns,iDynams} = fullfile(path,['u' name ext]);
        clear path name ext
    end
end

for iRuns = 1:numruns_ve
    for iDynams = 1:numDynamics_ve
        [path, name, ext] = fileparts(imageList_ve{iRuns,iDynams});
        imageList_ve{iRuns,iDynams} = fullfile(path,['u' name ext]);
        clear path name ext
    end
end

[path, name, ext] = fileparts(imageList_audloc{1,1});
meanEPIfilename = fullfile(path,['mean' name ext]);
pathreal = path;
namereal = name;
extreal = ext;

% SLICE TIMING CORRECTION

for iRuns = 1:numruns_audloc
    matlabbatch_slicetiming{1}.spm.temporal.st.scans{iRuns} = imageList_audloc(iRuns,:)';
end

for iRuns = numruns_audloc+1:numruns_ve+numruns_audloc
    matlabbatch_slicetiming{1}.spm.temporal.st.scans{iRuns} = imageList_ve(iRuns-numruns_audloc,:)';
end

matlabbatch_slicetiming{1}.spm.temporal.st.nslices = 38;
matlabbatch_slicetiming{1}.spm.temporal.st.tr = 2.8;
matlabbatch_slicetiming{1}.spm.temporal.st.ta = 2.72631578947368;
matlabbatch_slicetiming{1}.spm.temporal.st.so = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38];
matlabbatch_slicetiming{1}.spm.temporal.st.refslice = 19;
matlabbatch_slicetiming{1}.spm.temporal.st.prefix = 'a';

for iRuns = 1:numruns_audloc
    for iDynams = 1:numDynamics_audloc
        [path, name, ext] = fileparts(imageList_audloc{iRuns,iDynams});
        imageList_audloc{iRuns,iDynams} = fullfile(path,['a' name ext]);
        clear path name ext
    end
end

for iRuns = 1:numruns_ve
    for iDynams = 1:numDynamics_ve
        [path, name, ext] = fileparts(imageList_ve{iRuns,iDynams});
        imageList_ve{iRuns,iDynams} = fullfile(path,['a' name ext]);
        clear path name ext
    end
end

% COREGISTER

matlabbatch_coreg{1}.spm.spatial.coreg.estimate.ref = {meanEPIfilename};
matlabbatch_coreg{1}.spm.spatial.coreg.estimate.source = {anatfilename};
matlabbatch_coreg{1}.spm.spatial.coreg.estimate.other = {''};
matlabbatch_coreg{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
matlabbatch_coreg{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
matlabbatch_coreg{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch_coreg{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];


% SEGMENT

matlabbatch_segment{1}.spm.spatial.preproc.channel.vols = {anatfilename};
matlabbatch_segment{1}.spm.spatial.preproc.channel.biasreg = 0.001;
matlabbatch_segment{1}.spm.spatial.preproc.channel.biasfwhm = 60;
matlabbatch_segment{1}.spm.spatial.preproc.channel.write = [0 1];
matlabbatch_segment{1}.spm.spatial.preproc.tissue(1).tpm = {'C:\toolbox\spm12\tpm\TPM.nii,1'};
matlabbatch_segment{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
matlabbatch_segment{1}.spm.spatial.preproc.tissue(1).native = [1 0];
matlabbatch_segment{1}.spm.spatial.preproc.tissue(1).warped = [0 0];
matlabbatch_segment{1}.spm.spatial.preproc.tissue(2).tpm = {'C:\toolbox\spm12\tpm\TPM.nii,2'};
matlabbatch_segment{1}.spm.spatial.preproc.tissue(2).ngaus = 1;
matlabbatch_segment{1}.spm.spatial.preproc.tissue(2).native = [1 0];
matlabbatch_segment{1}.spm.spatial.preproc.tissue(2).warped = [0 0];
matlabbatch_segment{1}.spm.spatial.preproc.tissue(3).tpm = {'C:\toolbox\spm12\tpm\TPM.nii,3'};
matlabbatch_segment{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
matlabbatch_segment{1}.spm.spatial.preproc.tissue(3).native = [1 0];
matlabbatch_segment{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
matlabbatch_segment{1}.spm.spatial.preproc.tissue(4).tpm = {'C:\toolbox\spm12\tpm\TPM.nii,4'};
matlabbatch_segment{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
matlabbatch_segment{1}.spm.spatial.preproc.tissue(4).native = [1 0];
matlabbatch_segment{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
matlabbatch_segment{1}.spm.spatial.preproc.tissue(5).tpm = {'C:\toolbox\spm12\tpm\TPM.nii,5'};
matlabbatch_segment{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
matlabbatch_segment{1}.spm.spatial.preproc.tissue(5).native = [1 0];
matlabbatch_segment{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
matlabbatch_segment{1}.spm.spatial.preproc.tissue(6).tpm = {'C:\toolbox\spm12\tpm\TPM.nii,6'};
matlabbatch_segment{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
matlabbatch_segment{1}.spm.spatial.preproc.tissue(6).native = [0 0];
matlabbatch_segment{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
matlabbatch_segment{1}.spm.spatial.preproc.warp.mrf = 1;
matlabbatch_segment{1}.spm.spatial.preproc.warp.cleanup = 1;
matlabbatch_segment{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
matlabbatch_segment{1}.spm.spatial.preproc.warp.affreg = 'mni';
matlabbatch_segment{1}.spm.spatial.preproc.warp.fwhm = 0;
matlabbatch_segment{1}.spm.spatial.preproc.warp.samp = 3;
matlabbatch_segment{1}.spm.spatial.preproc.warp.write = [1 1];

[path, name, ext] = fileparts(anatfilename);
biascorrectedanatfilename = fullfile(path,['m' name ext]);
ext = ext(1:4);
forwarddeformationfilename = fullfile(path,['y_' name ext]);

% NORMALISE FUNCTIONALS (Convert EPI from native to MNI space)
filestowrite = cell(numDynamics_audloc*numRuns_audloc+numDynamics_ve*numRuns_ve,1);
filestowrite(1:numDynamics_audloc*numRuns_audloc) = imageList_audloc(:);
filestowrite(numDynamics_audloc*numRuns_audloc+1:numDynamics_audloc*numRuns_audloc+numDynamics_ve*numRuns_ve) = imageList_ve(:);
filestowrite(end+1) = {meanEPIfilename};

matlabbatch_normalise{1}.spm.spatial.normalise.write.subj.def = {forwarddeformationfilename};
matlabbatch_normalise{1}.spm.spatial.normalise.write.subj.resample = filestowrite;
matlabbatch_normalise{1}.spm.spatial.normalise.write.woptions.bb = nan(2,3);
matlabbatch_normalise{1}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
matlabbatch_normalise{1}.spm.spatial.normalise.write.woptions.interp = 4;
matlabbatch_normalise{1}.spm.spatial.normalise.write.woptions.prefix = 'w';

for iRuns = 1:numruns_audloc
    for iDynams = 1:numDynamics_audloc
        [path, name, ext] = fileparts(imageList_audloc{iRuns,iDynams});
        imageList_audloc{iRuns,iDynams} = fullfile(path,['w' name ext]);
        clear path name ext
    end
end

for iRuns = 1:numruns_ve
    for iDynams = 1:numDynamics_ve
        [path, name, ext] = fileparts(imageList_ve{iRuns,iDynams});
        imageList_ve{iRuns,iDynams} = fullfile(path,['w' name ext]);
        clear path name ext
    end
end
meannormEPIfilename = fullfile(pathreal,['wmean' namereal extreal]);

% NORMALISE ANATOMICAL (Convert anatomical from native to MNI space)

matlabbatch_norm_anat{1}.spm.spatial.normalise.write.subj.def = {forwarddeformationfilename};
matlabbatch_norm_anat{1}.spm.spatial.normalise.write.subj.resample = {biascorrectedanatfilename};
matlabbatch_norm_anat{1}.spm.spatial.normalise.write.woptions.bb = nan(2,3);
matlabbatch_norm_anat{1}.spm.spatial.normalise.write.woptions.vox = [1 1 1];
matlabbatch_norm_anat{1}.spm.spatial.normalise.write.woptions.interp = 4;
matlabbatch_norm_anat{1}.spm.spatial.normalise.write.woptions.prefix = 'w';

[path, name, ext] = fileparts(biascorrectedanatfilename);
normalisedanatfilename = fullfile(path,['w' name ext]);

filestowritenorm = cell(numDynamics_audloc*numRuns_audloc+numDynamics_ve*numRuns_ve,1);
filestowritenorm(1:numDynamics_audloc*numRuns_audloc) = imageList_audloc(:);
filestowritenorm(numDynamics_audloc*numRuns_audloc+1:numDynamics_audloc*numRuns_audloc+numDynamics_ve*numRuns_ve) = imageList_ve(:);
filestowritenorm(end+1) = {meannormEPIfilename};

% SMOOTH NATIVE IMAGES TO 3MM FOR MVPA

matlabbatch_smooth_native3{1}.spm.spatial.smooth.data = filestowrite(:);
matlabbatch_smooth_native3{1}.spm.spatial.smooth.fwhm = [3 3 3];
matlabbatch_smooth_native3{1}.spm.spatial.smooth.dtype = 0;
matlabbatch_smooth_native3{1}.spm.spatial.smooth.im = 0;
matlabbatch_smooth_native3{1}.spm.spatial.smooth.prefix = 's3';

% SMOOTH NORMALISED IMAGES TO 8MM FOR GLM

matlabbatch_smooth_normalised8{1}.spm.spatial.smooth.data = filestowritenorm(:);
matlabbatch_smooth_normalised8{1}.spm.spatial.smooth.fwhm = [8 8 8];
matlabbatch_smooth_normalised8{1}.spm.spatial.smooth.dtype = 0;
matlabbatch_smooth_normalised8{1}.spm.spatial.smooth.im = 0;
matlabbatch_smooth_normalised8{1}.spm.spatial.smooth.prefix = 's8';


%% Run the steps (do these selectively as required)

spm_jobman('run', matlabbatch_realignunwarp);
spm_jobman('run', matlabbatch_slicetiming);
spm_jobman('run', matlabbatch_coreg);
spm_jobman('run', matlabbatch_segment);
spm_jobman('run', matlabbatch_normalise);
spm_jobman('run', matlabbatch_norm_anat);
spm_jobman('run', matlabbatch_smooth_native3);
spm_jobman('run', matlabbatch_smooth_normalised8);
spm_surf(['c1' PID '_T1w_orig.nii';'c2' PID '_T1w_orig.nii'],3);

%% Plot realignment parameters

rplist = dir('rp*.txt');
rpnum = size(rplist,1);

for irp = 1: rpnum
    figure;
    
    rp = dlmread(rplist(irp).name);
    
    % translations
    subplot(2,1,1)
    plot(rp(:,1:3))
    
    if irp < 3
        title(['rp ' PID ' task ' rplist(irp).name(end-15:end-12) ' bold ' rplist(irp).name(end-5:end-4)],...
            'FontName','Helvetica','FontSize', 15)
    else
        title(['rp ' PID ' task ' rplist(irp).name(end-13:end-12) ' bold ' rplist(irp).name(end-5:end-4)],...
            'FontName','Helvetica','FontSize', 15)
    end
    
    % add grid
    grid on
    set(gca,'GridLineStyle',':');
    
    % rotations
    subplot(2,1,2)
    plot(rp(:,4:6))
    
    % add grid
    grid on
    set(gca,'GridLineStyle',':');
    
    saveas(gcf,fullfile(root, [rplist(irp).name(1:end-4) '.emf']));
    
end
