%% Compute R^2 (coefficient of determination)

clear;
clc;

% Path to BCI scripts
addpath(genpath('E:\AMBRA\UoB\Exp\misc'));

load('group_bciSimulations_attmod_best10','group_best_results_mod');
bci_att = group_best_results_mod;
clear group_best_results_mod
load('group_fusSimulations_attmod_best10','group_best_results_mod');
fus_att = group_best_results_mod;
clear group_best_results_mod
load('group_fusSimulations_repmod_best10','group_best_results_repmod');
fus_rep = group_best_results_repmod;
clear group_best_results_repmod
load('group_fusSimulations_attmod_switch_best10','group_best_results_mod_switch');
fus_attrep = group_best_results_mod_switch;
clear group_best_results_mod_switch

logLike_null = nan(length(subjID),1);
logLike_bci_att = nan(length(subjID),1);
logLike_bci_noatt = nan(length(subjID),1);
logLike_fus_att = nan(length(subjID),1);
logLike_fus_noatt = nan(length(subjID),1);
logLike_fus_rep = nan(length(subjID),1);
logLike_fus_attrep = nan(length(subjID),1);

R2 = nan(length(subjID),4);
maxR2 = nan(length(subjID),1);

for iSubj = 1:length(subjID)
    
    origData = dir('*Exp_All_Sessions*.mat');
    load(origData.name, 'tdata');
    origData = tdata;
    
    % Reorganize dataset 'tdata' so that I get the information I need for model estimation (locA, locV, respA, respV)
    % Get rid of incorrect responses (wrong keypad), missed responses (no answer provided) and anticipated responses (RT < 100ms)
    origData = origData(origData.IncorResp==0 & origData.MissedResp==0 & origData.AntResp==0,:);
    
    % Put locA and locV
    if ~sum(strcmp(origData.Properties.VarNames(:),'locA')) % if locA is not present (and locV as well)
        origData.locA = NaN(size(origData,1),1);
        origData.locV = NaN(size(origData,1),1);
        for j = 1:size(origData,1)
            % locA: when people have to locate A, locA = location of the target and locV = location of the distractor
            if strcmp(origData.ResponseModality(j), 'Aud')
                origData.locA(j) = origData.TargetLoc(j);
                origData.locV(j) = origData.NonTargetLoc(j);
                % locV: when people have to locate V, locV = location of the target and locA = location of the distractor
            elseif strcmp(origData.ResponseModality(j), 'Vis')
                origData.locV(j) = origData.TargetLoc(j);
                origData.locA(j) = origData.NonTargetLoc(j);
            end
        end
    end
    
    % Put respA and respV
    if ~sum(strcmp(origData.Properties.VarNames(:),'respA')) % if locA is not present (and respV as well)
        origData.respA = NaN(size(origData,1),1);
        origData.respV = NaN(size(origData,1),1);
        for j = 1:size(origData,1)
            % respA: when people have to locate A, respA = response
            if strcmp(origData.ResponseModality(j), 'Aud')
                origData.respA(j) = origData.Response(j);
                origData.respV(j) = NaN;
                % respV: when people have to locate V, respV = response
            elseif strcmp(origData.ResponseModality(j), 'Vis')
                origData.respV(j) = origData.Response(j);
                origData.respA(j) = NaN;
            end
        end
    end
    
    % Keep only useful data (locV, locA, respV, respA)
    bciData = origData(:,{'locV','locA','respV','respA'});
    actDataVA = dataset2table(bciData);
    responseLoc = unique(bciData.locV);
    
    logLike_null(iSubj) = fitModelNull(actDataVA,responseLoc);
    
    %% logLike other models
    
    % BCI model no attention
    if strcmp(expID, 'MRI')
        file = fullfile(fullfile('E:\AMBRA\UoB\Data\MAMSI_MRI', subjID{iSubj}, 'behav\scanner', [subjID{iSubj} '_bciSimulations_best10']));
    elseif strcmp(expID, 'behav') || strcmp(expID, 'behav_MRI')
        file = fullfile(fullfile('E:\AMBRA\UoB\Data\MAMSI_MRI', 'MAMSI_MRI_behav', subjID{iSubj}, [subjID{iSubj} '_bciSimulations_best10']));
    end
    
    load(file,'fm_bestlogLike');
    logLike_bci_noatt(iSubj) = fm_bestlogLike;
    clear fm_bestlogLike
    
    % fusion model no attention
    if strcmp(expID, 'MRI')
        file = fullfile(fullfile('E:\AMBRA\UoB\Data\MAMSI_MRI', subjID{iSubj}, 'behav\scanner', [subjID{iSubj} '_fusSimulations_best10']));
    elseif strcmp(expID, 'behav') || strcmp(expID, 'behav_MRI')
        file = fullfile(fullfile('E:\AMBRA\UoB\Data\MAMSI_MRI', 'MAMSI_MRI_behav', subjID{iSubj}, [subjID{iSubj} '_fusSimulations_best10']));
    end
    
    load(file,'fm_bestlogLike');
    logLike_fus_noatt(iSubj) = fm_bestlogLike;
    clear fm_bestlogLike
    
    logLike_bci_att(iSubj) = bci_att{iSubj}.fmin_bestlogLike;
    logLike_fus_att(iSubj) = fus_att{iSubj}.fmin_bestlogLike;
    logLike_fus_rep(iSubj) = fus_rep{iSubj}.fmin_bestlogLike;
    logLike_fus_attrep(iSubj) = fus_attrep{iSubj}.fmin_bestlogLike;
    
    %% Compute R2 for discrete responses based on Nagelkerke (1991)
    n = size(bciData,1);
    maxR2(iSubj,1) = 1-exp((2/n)*(-logLike_null(iSubj)));
    
    % we invert logLike_null and logLike_modelX because we have negative logLike
    R2(iSubj,1) = 1-exp((-2/n)*(logLike_null(iSubj)-logLike_fus_noatt(iSubj)));
    R2(iSubj,2) = 1-exp((-2/n)*(logLike_null(iSubj)-logLike_fus_att(iSubj)));
    R2(iSubj,3) = 1-exp((-2/n)*(logLike_null(iSubj)-logLike_fus_rep(iSubj)));
    R2(iSubj,4) = 1-exp((-2/n)*(logLike_null(iSubj)-logLike_fus_attrep(iSubj)));
    R2(iSubj,5) = 1-exp((-2/n)*(logLike_null(iSubj)-logLike_bci_noatt(iSubj)));
    R2(iSubj,6) = 1-exp((-2/n)*(logLike_null(iSubj)-logLike_bci_att(iSubj)));    
    
    R2(iSubj,1) = R2(iSubj,1)/maxR2(iSubj,1);
    R2(iSubj,2) = R2(iSubj,2)/maxR2(iSubj,1);
    R2(iSubj,3) = R2(iSubj,3)/maxR2(iSubj,1);
    R2(iSubj,4) = R2(iSubj,4)/maxR2(iSubj,1);
    R2(iSubj,5) = R2(iSubj,5)/maxR2(iSubj,1);
    R2(iSubj,6) = R2(iSubj,6)/maxR2(iSubj,1);
    
end % End of loop over subjects

meanR2 = mean(R2);
semR2 = std(R2)/sqrt(length(subjID));

