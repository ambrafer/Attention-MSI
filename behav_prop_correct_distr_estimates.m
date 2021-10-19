%% Compute proportions of correct responses and plot distributions of spatial estimates in the fMRI experiment

close all;
clear;
clc;

%% general settings

% 9 AV locations in factorial design
ni = 3; % A loc
nj = 3; % V loc
nCond = 4; % 2 Att x 2 Resp
savePath = 'E:\AMBRA\UoB\Data\MAMSI_MRI\group\behav';
subjList={'sub-MA01'};
nSubj=length(subjList);

hist_mat=nan(9,3,nCond,nSubj);
prop_mat=nan(nj,nCond,nSubj);
stats_mat=nan(nSubj,12);

for isubj=1:nSubj
    
    % Select the data folder of the subject (subj_mat)
    dataPath = ['E:\AMBRA\UoB\Data\MAMSI_MRI\' subjList{isubj} '\behav\scanner'];
    
    % Load current subject dataset
    Files = dir(fullfile(dataPath,'*Exp_All_Sessions*.mat'));
    load(fullfile(dataPath,Files.name), 'tdata');
    
    ALVL_AA=tdata.Response(strcmp(tdata.ResponseModality,'Aud') & strcmp(tdata.AttentionModality,'Aud') & ...
        tdata.TargetLoc==-9 & tdata.NonTargetLoc==-9 & tdata.IncorResp==0 & tdata.MissedResp==0 & tdata.AntResp==0,:);
    ALVC_AA=tdata.Response(strcmp(tdata.ResponseModality,'Aud') & strcmp(tdata.AttentionModality,'Aud') & ...
        tdata.TargetLoc==-9 & tdata.NonTargetLoc==0 & tdata.IncorResp==0 & tdata.MissedResp==0 & tdata.AntResp==0,:);
    ALVR_AA=tdata.Response(strcmp(tdata.ResponseModality,'Aud') & strcmp(tdata.AttentionModality,'Aud') & ...
        tdata.TargetLoc==-9 & tdata.NonTargetLoc==9 & tdata.IncorResp==0 & tdata.MissedResp==0 & tdata.AntResp==0,:);
    ACVL_AA=tdata.Response(strcmp(tdata.ResponseModality,'Aud') & strcmp(tdata.AttentionModality,'Aud') & ...
        tdata.TargetLoc==0 & tdata.NonTargetLoc==-9 & tdata.IncorResp==0 & tdata.MissedResp==0 & tdata.AntResp==0,:);
    ACVC_AA=tdata.Response(strcmp(tdata.ResponseModality,'Aud') & strcmp(tdata.AttentionModality,'Aud') & ...
        tdata.TargetLoc==0 & tdata.NonTargetLoc==0 & tdata.IncorResp==0 & tdata.MissedResp==0 & tdata.AntResp==0,:);
    ACVR_AA=tdata.Response(strcmp(tdata.ResponseModality,'Aud') & strcmp(tdata.AttentionModality,'Aud') & ...
        tdata.TargetLoc==0 & tdata.NonTargetLoc==9 & tdata.IncorResp==0 & tdata.MissedResp==0 & tdata.AntResp==0,:);
    ARVL_AA=tdata.Response(strcmp(tdata.ResponseModality,'Aud') & strcmp(tdata.AttentionModality,'Aud') & ...
        tdata.TargetLoc==9 & tdata.NonTargetLoc==-9 & tdata.IncorResp==0 & tdata.MissedResp==0 & tdata.AntResp==0,:);
    ARVC_AA=tdata.Response(strcmp(tdata.ResponseModality,'Aud') & strcmp(tdata.AttentionModality,'Aud') & ...
        tdata.TargetLoc==9 & tdata.NonTargetLoc==0 & tdata.IncorResp==0 & tdata.MissedResp==0 & tdata.AntResp==0,:);
    ARVR_AA=tdata.Response(strcmp(tdata.ResponseModality,'Aud') & strcmp(tdata.AttentionModality,'Aud') & ...
        tdata.TargetLoc==9 & tdata.NonTargetLoc==9 & tdata.IncorResp==0 & tdata.MissedResp==0 & tdata.AntResp==0,:);
    
    ALVL_AV=tdata.Response(strcmp(tdata.ResponseModality,'Vis') & strcmp(tdata.AttentionModality,'Aud') & ...
        tdata.TargetLoc==-9 & tdata.NonTargetLoc==-9 & tdata.IncorResp==0 & tdata.MissedResp==0 & tdata.AntResp==0,:);
    ALVC_AV=tdata.Response(strcmp(tdata.ResponseModality,'Vis') & strcmp(tdata.AttentionModality,'Aud') & ...
        tdata.TargetLoc==0 & tdata.NonTargetLoc==-9 & tdata.IncorResp==0 & tdata.MissedResp==0 & tdata.AntResp==0,:);
    ALVR_AV=tdata.Response(strcmp(tdata.ResponseModality,'Vis') & strcmp(tdata.AttentionModality,'Aud') & ...
        tdata.TargetLoc==9 & tdata.NonTargetLoc==-9 & tdata.IncorResp==0 & tdata.MissedResp==0 & tdata.AntResp==0,:);
    ACVL_AV=tdata.Response(strcmp(tdata.ResponseModality,'Vis') & strcmp(tdata.AttentionModality,'Aud') & ...
        tdata.TargetLoc==-9 & tdata.NonTargetLoc==0 & tdata.IncorResp==0 & tdata.MissedResp==0 & tdata.AntResp==0,:);
    ACVC_AV=tdata.Response(strcmp(tdata.ResponseModality,'Vis') & strcmp(tdata.AttentionModality,'Aud') & ...
        tdata.TargetLoc==0 & tdata.NonTargetLoc==0 & tdata.IncorResp==0 & tdata.MissedResp==0 & tdata.AntResp==0,:);
    ACVR_AV=tdata.Response(strcmp(tdata.ResponseModality,'Vis') & strcmp(tdata.AttentionModality,'Aud') & ...
        tdata.TargetLoc==9 & tdata.NonTargetLoc==0 & tdata.IncorResp==0 & tdata.MissedResp==0 & tdata.AntResp==0,:);
    ARVL_AV=tdata.Response(strcmp(tdata.ResponseModality,'Vis') & strcmp(tdata.AttentionModality,'Aud') & ...
        tdata.TargetLoc==-9 & tdata.NonTargetLoc==9 & tdata.IncorResp==0 & tdata.MissedResp==0 & tdata.AntResp==0,:);
    ARVC_AV=tdata.Response(strcmp(tdata.ResponseModality,'Vis') & strcmp(tdata.AttentionModality,'Aud') & ...
        tdata.TargetLoc==0 & tdata.NonTargetLoc==9 & tdata.IncorResp==0 & tdata.MissedResp==0 & tdata.AntResp==0,:);
    ARVR_AV=tdata.Response(strcmp(tdata.ResponseModality,'Vis') & strcmp(tdata.AttentionModality,'Aud') & ...
        tdata.TargetLoc==9 & tdata.NonTargetLoc==9 & tdata.IncorResp==0 & tdata.MissedResp==0 & tdata.AntResp==0,:);
    
    ALVL_VA=tdata.Response(strcmp(tdata.ResponseModality,'Aud') & strcmp(tdata.AttentionModality,'Vis') & ...
        tdata.TargetLoc==-9 & tdata.NonTargetLoc==-9 & tdata.IncorResp==0 & tdata.MissedResp==0 & tdata.AntResp==0,:);
    ALVC_VA=tdata.Response(strcmp(tdata.ResponseModality,'Aud') & strcmp(tdata.AttentionModality,'Vis') & ...
        tdata.TargetLoc==-9 & tdata.NonTargetLoc==0 & tdata.IncorResp==0 & tdata.MissedResp==0 & tdata.AntResp==0,:);
    ALVR_VA=tdata.Response(strcmp(tdata.ResponseModality,'Aud') & strcmp(tdata.AttentionModality,'Vis') & ...
        tdata.TargetLoc==-9 & tdata.NonTargetLoc==9 & tdata.IncorResp==0 & tdata.MissedResp==0 & tdata.AntResp==0,:);
    ACVL_VA=tdata.Response(strcmp(tdata.ResponseModality,'Aud') & strcmp(tdata.AttentionModality,'Vis') & ...
        tdata.TargetLoc==0 & tdata.NonTargetLoc==-9 & tdata.IncorResp==0 & tdata.MissedResp==0 & tdata.AntResp==0,:);
    ACVC_VA=tdata.Response(strcmp(tdata.ResponseModality,'Aud') & strcmp(tdata.AttentionModality,'Vis') & ...
        tdata.TargetLoc==0 & tdata.NonTargetLoc==0 & tdata.IncorResp==0 & tdata.MissedResp==0 & tdata.AntResp==0,:);
    ACVR_VA=tdata.Response(strcmp(tdata.ResponseModality,'Aud') & strcmp(tdata.AttentionModality,'Vis') & ...
        tdata.TargetLoc==0 & tdata.NonTargetLoc==9 & tdata.IncorResp==0 & tdata.MissedResp==0 & tdata.AntResp==0,:);
    ARVL_VA=tdata.Response(strcmp(tdata.ResponseModality,'Aud') & strcmp(tdata.AttentionModality,'Vis') & ...
        tdata.TargetLoc==9 & tdata.NonTargetLoc==-9 & tdata.IncorResp==0 & tdata.MissedResp==0 & tdata.AntResp==0,:);
    ARVC_VA=tdata.Response(strcmp(tdata.ResponseModality,'Aud') & strcmp(tdata.AttentionModality,'Vis') & ...
        tdata.TargetLoc==9 & tdata.NonTargetLoc==0 & tdata.IncorResp==0 & tdata.MissedResp==0 & tdata.AntResp==0,:);
    ARVR_VA=tdata.Response(strcmp(tdata.ResponseModality,'Aud') & strcmp(tdata.AttentionModality,'Vis') & ...
        tdata.TargetLoc==9 & tdata.NonTargetLoc==9 & tdata.IncorResp==0 & tdata.MissedResp==0 & tdata.AntResp==0,:);
    
    ALVL_VV=tdata.Response(strcmp(tdata.ResponseModality,'Vis') & strcmp(tdata.AttentionModality,'Vis') & ...
        tdata.TargetLoc==-9 & tdata.NonTargetLoc==-9 & tdata.IncorResp==0 & tdata.MissedResp==0 & tdata.AntResp==0,:);
    ALVC_VV=tdata.Response(strcmp(tdata.ResponseModality,'Vis') & strcmp(tdata.AttentionModality,'Vis') & ...
        tdata.TargetLoc==0 & tdata.NonTargetLoc==-9 & tdata.IncorResp==0 & tdata.MissedResp==0 & tdata.AntResp==0,:);
    ALVR_VV=tdata.Response(strcmp(tdata.ResponseModality,'Vis') & strcmp(tdata.AttentionModality,'Vis') & ...
        tdata.TargetLoc==9 & tdata.NonTargetLoc==-9 & tdata.IncorResp==0 & tdata.MissedResp==0 & tdata.AntResp==0,:);
    ACVL_VV=tdata.Response(strcmp(tdata.ResponseModality,'Vis') & strcmp(tdata.AttentionModality,'Vis') & ...
        tdata.TargetLoc==-9 & tdata.NonTargetLoc==0 & tdata.IncorResp==0 & tdata.MissedResp==0 & tdata.AntResp==0,:);
    ACVC_VV=tdata.Response(strcmp(tdata.ResponseModality,'Vis') & strcmp(tdata.AttentionModality,'Vis') & ...
        tdata.TargetLoc==0 & tdata.NonTargetLoc==0 & tdata.IncorResp==0 & tdata.MissedResp==0 & tdata.AntResp==0,:);
    ACVR_VV=tdata.Response(strcmp(tdata.ResponseModality,'Vis') & strcmp(tdata.AttentionModality,'Vis') & ...
        tdata.TargetLoc==9 & tdata.NonTargetLoc==0 & tdata.IncorResp==0 & tdata.MissedResp==0 & tdata.AntResp==0,:);
    ARVL_VV=tdata.Response(strcmp(tdata.ResponseModality,'Vis') & strcmp(tdata.AttentionModality,'Vis') & ...
        tdata.TargetLoc==-9 & tdata.NonTargetLoc==9 & tdata.IncorResp==0 & tdata.MissedResp==0 & tdata.AntResp==0,:);
    ARVC_VV=tdata.Response(strcmp(tdata.ResponseModality,'Vis') & strcmp(tdata.AttentionModality,'Vis') & ...
        tdata.TargetLoc==0 & tdata.NonTargetLoc==9 & tdata.IncorResp==0 & tdata.MissedResp==0 & tdata.AntResp==0,:);
    ARVR_VV=tdata.Response(strcmp(tdata.ResponseModality,'Vis') & strcmp(tdata.AttentionModality,'Vis') & ...
        tdata.TargetLoc==9 & tdata.NonTargetLoc==9 & tdata.IncorResp==0 & tdata.MissedResp==0 & tdata.AntResp==0,:);
    
    %% proportions of correct responses
    
    ALVL_AA_prop=sum(ALVL_AA==-9)/length(ALVL_AA);
    ALVC_AA_prop=sum(ALVC_AA==-9)/length(ALVC_AA);
    ALVR_AA_prop=sum(ALVR_AA==-9)/length(ALVR_AA);
    ACVL_AA_prop=sum(ACVL_AA==0)/length(ACVL_AA);
    ACVC_AA_prop=sum(ACVC_AA==0)/length(ACVC_AA);
    ACVR_AA_prop=sum(ACVR_AA==0)/length(ACVR_AA);
    ARVL_AA_prop=sum(ARVL_AA==9)/length(ARVL_AA);
    ARVC_AA_prop=sum(ARVC_AA==9)/length(ARVC_AA);
    ARVR_AA_prop=sum(ARVR_AA==9)/length(ARVR_AA);
    
    ALVL_AV_prop=sum(ALVL_AV==-9)/length(ALVL_AV);
    ALVC_AV_prop=sum(ALVC_AV==0)/length(ALVC_AV);
    ALVR_AV_prop=sum(ALVR_AV==9)/length(ALVR_AV);
    ACVL_AV_prop=sum(ACVL_AV==-9)/length(ACVL_AV);
    ACVC_AV_prop=sum(ACVC_AV==0)/length(ACVC_AV);
    ACVR_AV_prop=sum(ACVR_AV==9)/length(ACVR_AV);
    ARVL_AV_prop=sum(ARVL_AV==-9)/length(ARVL_AV);
    ARVC_AV_prop=sum(ARVC_AV==0)/length(ARVC_AV);
    ARVR_AV_prop=sum(ARVR_AV==9)/length(ARVR_AV);
    
    ALVL_VA_prop=sum(ALVL_VA==-9)/length(ALVL_VA);
    ALVC_VA_prop=sum(ALVC_VA==-9)/length(ALVC_VA);
    ALVR_VA_prop=sum(ALVR_VA==-9)/length(ALVR_VA);
    ACVL_VA_prop=sum(ACVL_VA==0)/length(ACVL_VA);
    ACVC_VA_prop=sum(ACVC_VA==0)/length(ACVC_VA);
    ACVR_VA_prop=sum(ACVR_VA==0)/length(ACVR_VA);
    ARVL_VA_prop=sum(ARVL_VA==9)/length(ARVL_VA);
    ARVC_VA_prop=sum(ARVC_VA==9)/length(ARVC_VA);
    ARVR_VA_prop=sum(ARVR_VA==9)/length(ARVR_VA);
    
    ALVL_VV_prop=sum(ALVL_VV==-9)/length(ALVL_VV);
    ALVC_VV_prop=sum(ALVC_VV==0)/length(ALVC_VV);
    ALVR_VV_prop=sum(ALVR_VV==9)/length(ALVR_VV);
    ACVL_VV_prop=sum(ACVL_VV==-9)/length(ACVL_VV);
    ACVC_VV_prop=sum(ACVC_VV==0)/length(ACVC_VV);
    ACVR_VV_prop=sum(ACVR_VV==9)/length(ACVR_VV);
    ARVL_VV_prop=sum(ARVL_VV==-9)/length(ARVL_VV);
    ARVC_VV_prop=sum(ARVC_VV==0)/length(ARVC_VV);
    ARVR_VV_prop=sum(ARVR_VV==9)/length(ARVR_VV);
    
    %% proportions of correct responses as a function of AV disparity
    
    prop_mat(:,1,isubj)=...
        [mean([ALVL_AA_prop ACVC_AA_prop ARVR_AA_prop]),...
        mean([ALVC_AA_prop ACVL_AA_prop ARVC_AA_prop ACVR_AA_prop]),...
        mean([ARVL_AA_prop ALVR_AA_prop])];
    
    prop_mat(:,2,isubj)=...
        [mean([ALVL_AV_prop ACVC_AV_prop ARVR_AV_prop]),...
        mean([ALVC_AV_prop ACVL_AV_prop ARVC_AV_prop ACVR_AV_prop]),...
        mean([ARVL_AV_prop ALVR_AV_prop])];
    
    prop_mat(:,3,isubj)=...
        [mean([ALVL_VA_prop ACVC_VA_prop ARVR_VA_prop]),...
        mean([ALVC_VA_prop ACVL_VA_prop ARVC_VA_prop ACVR_VA_prop]),...
        mean([ARVL_VA_prop ALVR_VA_prop])];
    
    prop_mat(:,4,isubj)=...
        [mean([ALVL_VV_prop ACVC_VV_prop ARVR_VV_prop]),...
        mean([ALVC_VV_prop ACVL_VV_prop ARVC_VV_prop ACVR_VV_prop]),...
        mean([ARVL_VV_prop ALVR_VV_prop])];
    
    stats_mat(isubj,:)=[mean([ALVL_AA_prop ACVC_AA_prop ARVR_AA_prop]),...
        mean([ALVC_AA_prop ACVL_AA_prop ARVC_AA_prop ACVR_AA_prop]),...
        mean([ARVL_AA_prop ALVR_AA_prop]),...
        mean([ALVL_VA_prop ACVC_VA_prop ARVR_VA_prop]),...
        mean([ALVC_VA_prop ACVL_VA_prop ARVC_VA_prop ACVR_VA_prop]),...
        mean([ARVL_VA_prop ALVR_VA_prop])...
        mean([ALVL_AV_prop ACVC_AV_prop ARVR_AV_prop]),...
        mean([ALVC_AV_prop ACVL_AV_prop ARVC_AV_prop ACVR_AV_prop]),...
        mean([ARVL_AV_prop ALVR_AV_prop])...
        mean([ALVL_VV_prop ACVC_VV_prop ARVR_VV_prop]),...
        mean([ALVC_VV_prop ACVL_VV_prop ARVC_VV_prop ACVR_VV_prop]),...
        mean([ARVL_VV_prop ALVR_VV_prop])];
    
    %% distributions for histograms
    loc = [-9 0 9];
    
    for iloc=1:3
        ALVL_AA_hist(iloc)=sum(ALVL_AA==loc(iloc))/length(ALVL_AA);
        ALVC_AA_hist(iloc)=sum(ALVC_AA==loc(iloc))/length(ALVC_AA);
        ALVR_AA_hist(iloc)=sum(ALVR_AA==loc(iloc))/length(ALVR_AA);
        ACVL_AA_hist(iloc)=sum(ACVL_AA==loc(iloc))/length(ACVL_AA);
        ACVC_AA_hist(iloc)=sum(ACVC_AA==loc(iloc))/length(ACVC_AA);
        ACVR_AA_hist(iloc)=sum(ACVR_AA==loc(iloc))/length(ACVR_AA);
        ARVL_AA_hist(iloc)=sum(ARVL_AA==loc(iloc))/length(ARVL_AA);
        ARVC_AA_hist(iloc)=sum(ARVC_AA==loc(iloc))/length(ARVC_AA);
        ARVR_AA_hist(iloc)=sum(ARVR_AA==loc(iloc))/length(ARVR_AA);
        
        ALVL_AV_hist(iloc)=sum(ALVL_AV==loc(iloc))/length(ALVL_AV);
        ALVC_AV_hist(iloc)=sum(ALVC_AV==loc(iloc))/length(ALVC_AV);
        ALVR_AV_hist(iloc)=sum(ALVR_AV==loc(iloc))/length(ALVR_AV);
        ACVL_AV_hist(iloc)=sum(ACVL_AV==loc(iloc))/length(ACVL_AV);
        ACVC_AV_hist(iloc)=sum(ACVC_AV==loc(iloc))/length(ACVC_AV);
        ACVR_AV_hist(iloc)=sum(ACVR_AV==loc(iloc))/length(ACVR_AV);
        ARVL_AV_hist(iloc)=sum(ARVL_AV==loc(iloc))/length(ARVL_AV);
        ARVC_AV_hist(iloc)=sum(ARVC_AV==loc(iloc))/length(ARVC_AV);
        ARVR_AV_hist(iloc)=sum(ARVR_AV==loc(iloc))/length(ARVR_AV);
        
        ALVL_VA_hist(iloc)=sum(ALVL_VA==loc(iloc))/length(ALVL_VA);
        ALVC_VA_hist(iloc)=sum(ALVC_VA==loc(iloc))/length(ALVC_VA);
        ALVR_VA_hist(iloc)=sum(ALVR_VA==loc(iloc))/length(ALVR_VA);
        ACVL_VA_hist(iloc)=sum(ACVL_VA==loc(iloc))/length(ACVL_VA);
        ACVC_VA_hist(iloc)=sum(ACVC_VA==loc(iloc))/length(ACVC_VA);
        ACVR_VA_hist(iloc)=sum(ACVR_VA==loc(iloc))/length(ACVR_VA);
        ARVL_VA_hist(iloc)=sum(ARVL_VA==loc(iloc))/length(ARVL_VA);
        ARVC_VA_hist(iloc)=sum(ARVC_VA==loc(iloc))/length(ARVC_VA);
        ARVR_VA_hist(iloc)=sum(ARVR_VA==loc(iloc))/length(ARVR_VA);
        
        ALVL_VV_hist(iloc)=sum(ALVL_VV==loc(iloc))/length(ALVL_VV);
        ALVC_VV_hist(iloc)=sum(ALVC_VV==loc(iloc))/length(ALVC_VV);
        ALVR_VV_hist(iloc)=sum(ALVR_VV==loc(iloc))/length(ALVR_VV);
        ACVL_VV_hist(iloc)=sum(ACVL_VV==loc(iloc))/length(ACVL_VV);
        ACVC_VV_hist(iloc)=sum(ACVC_VV==loc(iloc))/length(ACVC_VV);
        ACVR_VV_hist(iloc)=sum(ACVR_VV==loc(iloc))/length(ACVR_VV);
        ARVL_VV_hist(iloc)=sum(ARVL_VV==loc(iloc))/length(ARVL_VV);
        ARVC_VV_hist(iloc)=sum(ARVC_VV==loc(iloc))/length(ARVC_VV);
        ARVR_VV_hist(iloc)=sum(ARVR_VV==loc(iloc))/length(ARVR_VV);
    end
    
    %% distributions for histogams as a function of stimuli location
    
    hist_mat(:,:,1,isubj)=...
        [ALVL_AA_hist;ACVL_AA_hist;ARVL_AA_hist;...
        ALVC_AA_hist;ACVC_AA_hist;ARVC_AA_hist;...
        ALVR_AA_hist;ACVR_AA_hist;ARVR_AA_hist];
    
    hist_mat(:,:,2,isubj)=...
        [ALVL_VA_hist;ACVL_VA_hist;ARVL_VA_hist;...
        ALVC_VA_hist;ACVC_VA_hist;ARVC_VA_hist;...
        ALVR_VA_hist;ACVR_VA_hist;ARVR_VA_hist];
    
    hist_mat(:,:,3,isubj)=...
        [ALVL_AV_hist;ACVL_AV_hist;ARVL_AV_hist;...
        ALVC_AV_hist;ACVC_AV_hist;ARVC_AV_hist;...
        ALVR_AV_hist;ACVR_AV_hist;ARVR_AV_hist];
    
    hist_mat(:,:,4,isubj)=...
        [ALVL_VV_hist;ACVL_VV_hist;ARVL_VV_hist;...
        ALVC_VV_hist;ACVC_VV_hist;ARVC_VV_hist;...
        ALVR_VV_hist;ACVR_VV_hist;ARVR_VV_hist];
    
    %% BCI predictions
    % Get a list of .mat files to analyse
    FilesVE = dir(fullfile(dataPath_bci,'group_bciSimulations_attmod*.mat'));
    load(fullfile(dataPath_bci,FilesVE.name),'*best_results*');
    group_best_results=group_best_results_mod;
    
    data_behA=[];
    data_behV=[];
    data_predA=[];
    data_predV=[];
    
    for icond=1:18
        data_behA=[data_behA;group_best_results{isubj}.bci_details(icond).freq_dataA];
        data_behV=[data_behV;group_best_results{isubj}.bci_details(icond).freq_dataV];
        data_predA=[data_predA;group_best_results{isubj}.bci_details(icond).freq_predA];
        data_predV=[data_predV;group_best_results{isubj}.bci_details(icond).freq_predV];
    end
    hist_mat_bci(:,:,1,isubj)=data_predA(1:9,:);
    hist_mat_bci(:,:,2,isubj)=data_predA(10:18,:);
    hist_mat_bci(:,:,3,isubj)=data_predV(1:9,:);
    hist_mat_bci(:,:,4,isubj)=data_predV(10:18,:);
    
    hist_mat_beh(:,:,1,isubj)=data_behA(1:9,:);
    hist_mat_beh(:,:,2,isubj)=data_behA(10:18,:);
    hist_mat_beh(:,:,3,isubj)=data_behV(1:9,:);
    hist_mat_beh(:,:,4,isubj)=data_behV(10:18,:);

end

%% group means and std/sem

hist_mat_group_mean = mean(hist_mat,4);
hist_mat_group_std = std(hist_mat,0,4);
hist_mat_group_sem = hist_mat_group_std/sqrt(nSubj);

hist_mat_bci_group_mean = mean(hist_mat_bci,4);
hist_mat_bci_group_std = std(hist_mat_bci,0,4);
hist_mat_bci_group_sem = hist_mat_bci_group_std/sqrt(nSubj);

cols.yelAud1  = [165 0 33]/255;
cols.yelAud2  = [255 153 153]/255;
cols.blueVis1 = [153 204 255]/255;
cols.blueVis2 = [0 51 204]/255;

%% plot distributions (histograms)

positionXY = [0, 0, 500, 500];
loc = [1 2 3];

% auditory report
figure('color', [1 1 1], 'Position', positionXY);
% subplot loop (rows)
for i = 1:ni    
    % subplot loop (columns)
    for j = 1:nj
        % determine plot id and position
        curplotid = (i-1)*nj+j;
        % determine plot id and position
        hsp = subplot(ni,nj,curplotid);
        plot(loc, hist_mat_group_mean(curplotid,:,1),'Color',cols.yelAud1,'LineWidth',1.2); hold on
        plot(loc, hist_mat_bci_group_mean(curplotid,:,1),'--','Color',cols.yelAud1,'LineWidth',1.2); hold on
        plot(loc, hist_mat_group_mean(curplotid,:,2),'Color',cols.yelAud2,'LineWidth',1.2); hold on
        plot(loc, hist_mat_bci_group_mean(curplotid,:,2),'--','Color',cols.yelAud2,'LineWidth',1.2); hold on
        
        % plot settings
        xl = [0.8 3.2]; xlim(xl);
        yl = [-0.05 1.05]; ylim(yl);
        set(gca,'FontName', 'Arial');
        set(gca,'FontSize', 12);        
        ticksX = [];
        ticksY = [];
        set(gca, 'YTick', ticksY);
        set(gca, 'XTick', ticksX);
        h = gca; h.YAxis.Visible = 'off';
        set(gca,'LineWidth',1.2);
        box off
    end 
end

% visual report
figure('color', [1 1 1], 'Position', positionXY);
% subplot loop (rows)
for i = 1:ni    
    % subplot loop (columns)
    for j = 1:nj
        % determine plot id and position
        curplotid = (i-1)*nj+j;
        % determine plot id and position
        hsp = subplot(ni,nj,curplotid);
        plot(loc, hist_mat_group_mean(curplotid,:,4),'Color',cols.blueVis2,'LineWidth',1.2); hold on
        plot(loc, hist_mat_bci_group_mean(curplotid,:,4),'--','Color',cols.blueVis2,'LineWidth',1.2); hold on
        plot(loc, hist_mat_group_mean(curplotid,:,3),'Color',cols.blueVis1,'LineWidth',1.2); hold on
        plot(loc, hist_mat_bci_group_mean(curplotid,:,3),'--','Color',cols.blueVis1,'LineWidth',1.2); hold on
        
        % plot settings
        xl = [0.8 3.2]; xlim(xl);
        yl = [-0.05 1.05]; ylim(yl);
        set(gca,'FontName', 'Arial');
        set(gca,'FontSize', 12);        
        ticksX = [];
        ticksY = [];
        set(gca, 'YTick', ticksY);
        set(gca, 'XTick', ticksX);
        h = gca; h.YAxis.Visible = 'off';
        set(gca,'LineWidth',1.2);
        box off
    end  
end
