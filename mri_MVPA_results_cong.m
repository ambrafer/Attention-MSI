%% MVPA analysis: save decoded labels of train cong/test cong as a function of attention and report

clear;
close all;
clc;

% general settings
subjID_list = {'sub-MA01'};
anl_scheme = 'cong_cong';
scaling = 'min0max1_concat_scale';
anl_notes = '_7_sessions_nosmooth_tempder_concat';
% experiment info
AV = [-9 0 9];
AVcomb = 3; % number of AV spatial combinations
condNr = 4; % number of AttxRep conditions
% number of MRI runs
runNr = 7;
% rois
roi = {'V1-3';'IPS0-2';'IPS3-4_SPL1';'TE1.0-1.1';'PT'};

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
        
        % res(iroi).pred_labels: rows = iterations (runs); columns = conditions (AV congruent positions x Att x Rep)
        for iIters = 1:runNr
            res(iroi).pred_labels(iIters,:) = results.predicted_labels.output.predicted_labels{iIters}';
        end
        % prepare IV (column 1 = constant; column 2 = true labels)
        X(:,1) = ones(runNr*size(res(iroi).pred_labels,2),1);
        X(:,2) = repmat([ones(runNr,1)*-9;zeros(runNr,1);ones(runNr,1)*9],4,1); % 4 repetitions: Att x Rep
        % prepare DV
        y = res(iroi).pred_labels(:);
        % perform stats
        % b(1)=intercept, b(2)=slope
        % bint = confidence intervals for each coefficient estimate
        % regstats = R^2 statistic, F-statistic and its p-value, estimate of error variance
        [res(iroi).b_tot,res(iroi).bint_tot,~,~,res(iroi).regstats_tot] = regress(y,X);
        res(iroi).RMSE_tot = sqrt(mean((y-X(:,2)).^2));
        [res(iroi).corr_rho_tot,res(iroi).corr_pval_tot] = corr(X(:,2),y,'Tail','right');
        clear X y
        
        % summary results across runs divided by four AttxResp experimental conditions and three AV combinations
        res(iroi).mean_bycond = mean(res(iroi).pred_labels);
        res(iroi).mean_bycond = [res(iroi).mean_bycond(:,1:3);res(iroi).mean_bycond(:,4:6);res(iroi).mean_bycond(:,7:9);res(iroi).mean_bycond(:,10:12)];
        res(iroi).std_bycond = std(res(iroi).pred_labels);
        res(iroi).std_bycond = [res(iroi).std_bycond(:,1:3);res(iroi).std_bycond(:,4:6);res(iroi).std_bycond(:,7:9);res(iroi).std_bycond(:,10:12)];
        res(iroi).sem_bycond = res(iroi).std_bycond/sqrt(runNr);
        
        % summary results across runs and conditions divided by six AV combinations
        res(iroi).mean = mean([res(iroi).pred_labels(:,1:3);res(iroi).pred_labels(:,4:6);res(iroi).pred_labels(:,7:9);res(iroi).pred_labels(:,10:12)]);
        res(iroi).std = std([res(iroi).pred_labels(:,1:3);res(iroi).pred_labels(:,4:6);res(iroi).pred_labels(:,7:9);res(iroi).pred_labels(:,10:12)]);
        res(iroi).sem = res(iroi).std/sqrt(runNr*condNr);
        
        % individual results of each AV combination divided by four AttxResp experimental conditions
        res(iroi).Y(:,1)=reshape(res(iroi).pred_labels(:,1:AVcomb),AVcomb*runNr,1);
        res(iroi).Y(:,2)=reshape(res(iroi).pred_labels(:,AVcomb+1:AVcomb*2),AVcomb*runNr,1);
        res(iroi).Y(:,3)=reshape(res(iroi).pred_labels(:,AVcomb*2+1:AVcomb*3),AVcomb*runNr,1);
        res(iroi).Y(:,4)=reshape(res(iroi).pred_labels(:,AVcomb*3+1:AVcomb*4),AVcomb*runNr,1);
        
        for i=1:4 % Att x Resp
            % prepare IV (column 1 = constant; column 2 = true labels)
            X(:,1) = ones(runNr*AVcomb,1);
            X(:,2) = [ones(runNr,1)*-9;zeros(runNr,1);ones(runNr,1)*9];
            % prepare DV
            y = res(iroi).Y(:,i);
            % perform stats
            % b(1)=intercept, b(2)=slope
            % bint = confidence intervals for each coefficient estimate
            % regstats = R^2 statistic, F-statistic and its p-value, estimate of error variance
            [res(iroi).b(:,i),res(iroi).bint{i},~,~,res(iroi).regstats(i,:)] = regress(y,X);
            res(iroi).RMSE(i,:) = sqrt(mean((y-X(:,2)).^2));
            [res(iroi).corr_rho(i,:),res(iroi).corr_pval(i,:)] = corr(X(:,2),y,'Tail','right');
            clear X y
        end
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
    save('res_cong_cong_pred_labels', 'res', 'AV');
    
    disp(' ');
    disp([subjID ' completed']);
    disp(' ');
    
end % subj