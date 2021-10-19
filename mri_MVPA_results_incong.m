%% Plot decoded labels of train cong/test cong as a function of attention and response
% we plot antisymmetric positions to easily check for attentional effect

clear;
close all;
clc;

% general settings
subjID_list = {'sub-MA01'};
anl_scheme = 'cong_incong';
scaling = 'min0max1_concat_scale';
anl_notes = '_7_sessions_nosmooth_tempder_concat';
% experiment info
% A = -9 -9  0  0  9  9
% V =  0  9 -9  9 -9  0
AV = [-9 -9  0  0  9  9;...
    0  9 -9  9 -9  0];
AVcomb = 6; % number of AV spatial combinations
condNr = 4; % number of AttxRep conditions
% number of MRI runs
runNr = 7;
% rois
roi = {'V1-3';'IPS0-2';'IPS3-4_SPL1';'TE1.0-1.1';'PT'};
% number of cross-validations
crossNr = 1;

for iSubj = 1:length(subjID_list)
    
    subjID = subjID_list{iSubj};
    
    res = struct;
    
    for iroi = 1:length(roi)
        
        % save current roi name
        res(iroi).roi_name = roi{iroi};
        
        % Select the data folder
        dataPath = fullfile('E:\AMBRA\UoB\Data\MAMSI_MRI', subjID, 'derivatives\ve\MVPA\', scaling, ['roi_' roi{iroi} '_train_test_' anl_scheme anl_notes]);
        cd(dataPath);
        load('res_predicted_labels.mat'); % individual values
        
        % res(iroi).pred_labels: rows = iterations (runs); columns = conditions (AV incongruent positions x Att x Resp)
        for iIters = 1:runNr
            res(iroi).pred_labels(iIters,:) = results.predicted_labels.output.predicted_labels{iIters}';
        end
        
        % summary results across runs divided by four AttxResp experimental conditions and three AV combinations
        %res(iroi).pred_labels = mean(res(iroi).pred_labels);
        res(iroi).byruncond = [];
        for i = 1:24:size(res(iroi).pred_labels,2)
            res(iroi).byruncond = [res(iroi).byruncond; res(iroi).pred_labels(:,i:i+23)];
        end
        
        res(iroi).mean_bycond = mean(res(iroi).byruncond);
        res(iroi).std_bycond = std(res(iroi).byruncond);
        res(iroi).sem_bycond = res(iroi).std_bycond/sqrt(crossNr*runNr);
        res(iroi).mean_bycond = [res(iroi).mean_bycond(:,1:6);res(iroi).mean_bycond(:,7:12);res(iroi).mean_bycond(:,13:18);res(iroi).mean_bycond(:,19:24)];
        res(iroi).std_bycond = [res(iroi).std_bycond(:,1:6);res(iroi).std_bycond(:,7:12);res(iroi).std_bycond(:,13:18);res(iroi).std_bycond(:,19:24)];
        res(iroi).sem_bycond = [res(iroi).sem_bycond(:,1:6);res(iroi).sem_bycond(:,7:12);res(iroi).sem_bycond(:,13:18);res(iroi).sem_bycond(:,19:24)];
        
        % summary results across runs and conditions divided by six AV combinations
        res(iroi).mean = mean([res(iroi).byruncond(:,1:6);res(iroi).byruncond(:,7:12);res(iroi).byruncond(:,13:18);res(iroi).byruncond(:,19:24)]);
        res(iroi).std = std([res(iroi).byruncond(:,1:6);res(iroi).byruncond(:,7:12);res(iroi).byruncond(:,13:18);res(iroi).byruncond(:,19:24)]);
        res(iroi).sem = res(iroi).std/sqrt(crossNr*runNr*condNr);
        
        disp('...');
        disp([subjID ' ' roi{iroi} ' done']);
        
    end % roi
    
    %% Save results
    
    save_notes = 'results';
    
    % Select the data folder
    savePath = fullfile('E:\AMBRA\UoB\Data\MAMSI_MRI', subjID, 'derivatives\ve\MVPA', save_notes, scaling);
    if ~isfolder(savePath)
        mkdir(savePath)
    end
    
    cd(savePath);
    save('res_cong_incong_pred_labels', 'res', 'AV');
    
    disp(' ');
    disp([subjID ' completed']);    
    disp(' ');
    
end % subj