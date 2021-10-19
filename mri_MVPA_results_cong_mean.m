%% MVPA analysis: congruent responses
% compute mean decoded congruent spatial locations across participants
% take data from output of mri_MVPA_results_cong.m

clear;
close all;
clc;

% Exp info
runNr = 7;
scaling = 'min0max1_concat_scale';
save_notes = 'results';

% Subject ID
subjID_list = {'sub-MA01'};
% rois list
roi_list = {'V1-3';'IPS0-2';'IPS3-4_SPL1';'TE1.0-1.1';'PT'};

mean_roi = nan(4,3,length(roi_list),length(subjID_list)); % 4 AttxResp x 3 AVcomb
mean_tot_roi = nan(1,3,length(roi_list),length(subjID_list)); % 3 AVcomb

for iroi = 1:length(roi_list)
    
    ROI = roi_list{iroi};
    
    for isubj = 1:length(subjID_list)        
        subjID = subjID_list{isubj};        
        % Select the data folders of the subject
        dataPath = fullfile('E:\AMBRA\UoB\Data\MAMSI_MRI', subjID, 'derivatives\ve\MVPA', save_notes, scaling);
        cd(dataPath);        
        % Get list of .mat files        
        Files_con = dir('res_cong_cong_pred_labels.mat');
        load(Files_con.name);
        res_con = res;        
        AV_con = AV;
        clear res       
        mean_roi(:,:,iroi,isubj)=res_con(iroi).mean_bycond;
        mean_tot_roi(:,:,iroi,isubj)=res_con(iroi).mean;
    end    
    clear res_con
end

mean_group_mean=mean(mean_roi(:,:,:,:),4);
mean_tot=[];
for isubj=1:length(subjID_list)
    mean_tot=[mean_tot;mean_tot_roi(:,:,:,isubj)];
end
mean_tot=mean(mean_tot);
save(fullfile('E:\AMBRA\UoB\Data\MAMSI_MRI\group\MVPA',scaling,'cong_cong_mean_tot_bycond'),'mean_tot','mean_group_mean');
