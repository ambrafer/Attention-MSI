%% MVPA analysis: computation of nWAV

clear;
close all;
clc;

% general settings
subjID_list = {'sub-MA01'};
scaling = 'min0max1_concat_scale';
save_notes = 'results';
% experiment info
% A = -9 -9  0  0  9  9
% V =  0  9 -9  9 -9  0
AV = [-9 -9  0  0  9  9;...
    0  9 -9  9 -9  0];
Resp_list = {'Aud';'Vis'};
% number of MRI runs
runNr = 7;
% sqrt of number of data points for SEM calculation
n = sqrt(runNr);
% rois
roi = {'V1-3';'IPS0-2';'IPS3-4_SPL1';'TE1.0-1.1';'PT'};

% Load mean decoded congruent spatial locations across participants
load(fullfile('E:\AMBRA\UoB\Data\MAMSI_MRI\group\MVPA',scaling,'cong_cong_mean_tot_bycond'));

%% Compute nWAV per subject

nWAV_group = struct;

for iSubj = 1:length(subjID_list)
    
    subjID = subjID_list{iSubj};
    nWAV_resp = struct;
    
    for iResp = 1:length(Resp_list)
        
        Resp = Resp_list{iResp};
        % define indices based on report modality (Aud/Vis)
        if strcmp(Resp, 'Aud')
            a = 1;
            b = 3;
            anl_aud = 1;
            anl_vis = 1;
            c_aud = (1:6);
            c_vis = (13:18);
        elseif strcmp(Resp, 'Vis')
            a = 2;
            b = 4;
            anl_aud = 1;
            anl_vis = 1;
            c_aud = (7:12);
            c_vis = (19:24);
        end
        
        % Select the data folders of the subject
        dataPath = fullfile('E:\AMBRA\UoB\Data\MAMSI_MRI', subjID, 'derivatives\ve\MVPA', save_notes, scaling);
        cd(dataPath);
        
        % Get list of .mat files
        Files_inc = dir('res_cong_incong_pred_labels.mat');
        Files_con = dir('res_cong_cong_pred_labels.mat');
        load(Files_inc.name);
        res_inc = res;
        AV_inc = AV;
        load(Files_con.name);
        res_con = res;
        AV_con = AV;
        
        nWAV = struct;
        
        for iroi =1:length(roi)
            
            % save current roi name
            nWAV(iroi).roi_name = roi{iroi};
            
            if corr_type==3 || use_zscore==1
                res_inc(iroi).pred_labels = zscore(res_inc(iroi).pred_labels);
                res_con(iroi).pred_labels = zscore(res_con(iroi).pred_labels);
            end
            
            % mean congruent responses
            nWAV(iroi).con_aud_mean = res_con(iroi).mean_bycond(a,:);
            nWAV(iroi).con_vis_mean = res_con(iroi).mean_bycond(b,:);
            nWAV(iroi).con_pool_mean = mean(res_con(iroi).mean_bycond);
            
            nWAV(iroi).AL_VC_con_aud = (res_inc(iroi).pred_labels(:,c_aud(1),anl_aud)-mean_tot(:,1,iroi))...
                /(mean_tot(:,2,iroi)-mean_tot(:,1,iroi));
            nWAV(iroi).AL_VC_con_aud_mean = mean(nWAV(iroi).AL_VC_con_aud);
            nWAV(iroi).AL_VC_con_aud_std = std(nWAV(iroi).AL_VC_con_aud);
            nWAV(iroi).AL_VC_con_aud_sem = nWAV(iroi).AL_VC_con_aud_std/n;
            
            nWAV(iroi).AL_VC_con_vis = (res_inc(iroi).pred_labels(:,c_vis(1),anl_vis)-mean_tot(:,1,iroi))...
                /(mean_tot(:,2,iroi)-mean_tot(:,1,iroi));
            nWAV(iroi).AL_VC_con_vis_mean = mean(nWAV(iroi).AL_VC_con_vis);
            nWAV(iroi).AL_VC_con_vis_std = std(nWAV(iroi).AL_VC_con_vis);
            nWAV(iroi).AL_VC_con_vis_sem = nWAV(iroi).AL_VC_con_vis_std/n;
            
            nWAV(iroi).AL_VR_con_aud = (res_inc(iroi).pred_labels(:,c_aud(2),anl_aud)-mean_tot(:,1,iroi))...
                /(mean_tot(:,3,iroi)-mean_tot(:,1,iroi));
            nWAV(iroi).AL_VR_con_aud_mean = mean(nWAV(iroi).AL_VR_con_aud);
            nWAV(iroi).AL_VR_con_aud_std = std(nWAV(iroi).AL_VR_con_aud);
            nWAV(iroi).AL_VR_con_aud_sem = nWAV(iroi).AL_VR_con_aud_std/n;
            
            nWAV(iroi).AL_VR_con_vis = (res_inc(iroi).pred_labels(:,c_vis(2),anl_vis)-mean_tot(:,1,iroi))...
                /(mean_tot(:,3,iroi)-mean_tot(:,1,iroi));
            nWAV(iroi).AL_VR_con_vis_mean = mean(nWAV(iroi).AL_VR_con_vis);
            nWAV(iroi).AL_VR_con_vis_std = std(nWAV(iroi).AL_VR_con_vis);
            nWAV(iroi).AL_VR_con_vis_sem = nWAV(iroi).AL_VR_con_vis_std/n;
            
            nWAV(iroi).AC_VL_con_aud = (res_inc(iroi).pred_labels(:,c_aud(3),anl_aud)-mean_tot(:,2,iroi))...
                /(mean_tot(:,1,iroi)-mean_tot(:,2,iroi));
            nWAV(iroi).AC_VL_con_aud_mean = mean(nWAV(iroi).AC_VL_con_aud);
            nWAV(iroi).AC_VL_con_aud_std = std(nWAV(iroi).AC_VL_con_aud);
            nWAV(iroi).AC_VL_con_aud_sem = nWAV(iroi).AC_VL_con_aud_std/n;
            
            nWAV(iroi).AC_VL_con_vis = (res_inc(iroi).pred_labels(:,c_vis(3),anl_vis)-mean_tot(:,2,iroi))...
                /(mean_tot(:,1,iroi)-mean_tot(:,2,iroi));
            nWAV(iroi).AC_VL_con_vis_mean = mean(nWAV(iroi).AC_VL_con_vis);
            nWAV(iroi).AC_VL_con_vis_std = std(nWAV(iroi).AC_VL_con_vis);
            nWAV(iroi).AC_VL_con_vis_sem = nWAV(iroi).AC_VL_con_vis_std/n;
            
            nWAV(iroi).AC_VR_con_aud = (res_inc(iroi).pred_labels(:,c_aud(4),anl_aud)-mean_tot(:,2,iroi))...
                /(mean_tot(:,3,iroi)-mean_tot(:,2,iroi));
            nWAV(iroi).AC_VR_con_aud_mean = mean(nWAV(iroi).AC_VR_con_aud);
            nWAV(iroi).AC_VR_con_aud_std = std(nWAV(iroi).AC_VR_con_aud);
            nWAV(iroi).AC_VR_con_aud_sem = nWAV(iroi).AC_VR_con_aud_std/n;
            
            nWAV(iroi).AC_VR_con_vis = (res_inc(iroi).pred_labels(:,c_vis(4),anl_vis)-mean_tot(:,2,iroi))...
                /(mean_tot(:,3,iroi)-mean_tot(:,2,iroi));
            nWAV(iroi).AC_VR_con_vis_mean = mean(nWAV(iroi).AC_VR_con_vis);
            nWAV(iroi).AC_VR_con_vis_std = std(nWAV(iroi).AC_VR_con_vis);
            nWAV(iroi).AC_VR_con_vis_sem = nWAV(iroi).AC_VR_con_vis_std/n;
            
            nWAV(iroi).AR_VL_con_aud = (res_inc(iroi).pred_labels(:,c_aud(5),anl_aud)-mean_tot(:,3,iroi))...
                /(mean_tot(:,1,iroi)-mean_tot(:,3,iroi));
            nWAV(iroi).AR_VL_con_aud_mean = mean(nWAV(iroi).AR_VL_con_aud);
            nWAV(iroi).AR_VL_con_aud_std = std(nWAV(iroi).AR_VL_con_aud);
            nWAV(iroi).AR_VL_con_aud_sem = nWAV(iroi).AR_VL_con_aud_std/n;
            
            nWAV(iroi).AR_VL_con_vis = (res_inc(iroi).pred_labels(:,c_vis(5),anl_vis)-mean_tot(:,3,iroi))...
                /(mean_tot(:,1,iroi)-mean_tot(:,3,iroi));
            nWAV(iroi).AR_VL_con_vis_mean = mean(nWAV(iroi).AR_VL_con_vis);
            nWAV(iroi).AR_VL_con_vis_std = std(nWAV(iroi).AR_VL_con_vis);
            nWAV(iroi).AR_VL_con_vis_sem = nWAV(iroi).AR_VL_con_vis_std/n;
            
            nWAV(iroi).AR_VC_con_aud = (res_inc(iroi).pred_labels(:,c_aud(6),anl_aud)-mean_tot(:,3,iroi))...
                /(mean_tot(:,2,iroi)-mean_tot(:,3,iroi));
            nWAV(iroi).AR_VC_con_aud_mean = mean(nWAV(iroi).AR_VC_con_aud);
            nWAV(iroi).AR_VC_con_aud_std = std(nWAV(iroi).AR_VC_con_aud);
            nWAV(iroi).AR_VC_con_aud_sem = nWAV(iroi).AR_VC_con_aud_std/n;
            
            nWAV(iroi).AR_VC_con_vis = (res_inc(iroi).pred_labels(:,c_vis(6),anl_vis)-mean_tot(:,3,iroi))...
                /(mean_tot(:,2,iroi)-mean_tot(:,3,iroi));
            nWAV(iroi).AR_VC_con_vis_mean = mean(nWAV(iroi).AR_VC_con_vis);
            nWAV(iroi).AR_VC_con_vis_std = std(nWAV(iroi).AR_VC_con_vis);
            nWAV(iroi).AR_VC_con_vis_sem = nWAV(iroi).AR_VC_con_vis_std/n;
            
            % pool disparity
            %attended
            nWAV(iroi).pooldisp_con_aud_mean = ...
                mean([nWAV(iroi).AC_VL_con_aud;nWAV(iroi).AC_VR_con_aud;nWAV(iroi).AL_VC_con_aud;nWAV(iroi).AR_VC_con_aud;nWAV(iroi).AL_VR_con_aud;nWAV(iroi).AR_VL_con_aud]);
            nWAV(iroi).pooldisp_con_aud_std = ...
                std([nWAV(iroi).AC_VL_con_aud;nWAV(iroi).AC_VR_con_aud;nWAV(iroi).AL_VC_con_aud;nWAV(iroi).AR_VC_con_aud;nWAV(iroi).AL_VR_con_aud;nWAV(iroi).AR_VL_con_aud]);
            nWAV(iroi).pooldisp_con_aud_sem = ...
                nWAV(iroi).pooldisp_con_aud_std/sqrt(length([nWAV(iroi).AC_VL_con_aud;nWAV(iroi).AC_VR_con_aud;nWAV(iroi).AL_VC_con_aud;nWAV(iroi).AR_VC_con_aud;nWAV(iroi).AL_VR_con_aud;nWAV(iroi).AR_VL_con_aud]));
            %unattended
            nWAV(iroi).pooldisp_con_vis_mean = ...
                mean([nWAV(iroi).AC_VL_con_vis;nWAV(iroi).AC_VR_con_vis;nWAV(iroi).AL_VC_con_vis;nWAV(iroi).AR_VC_con_vis;nWAV(iroi).AL_VR_con_vis;nWAV(iroi).AR_VL_con_vis]);
            nWAV(iroi).pooldisp_con_vis_std = ...
                std([nWAV(iroi).AC_VL_con_vis;nWAV(iroi).AC_VR_con_vis;nWAV(iroi).AL_VC_con_vis;nWAV(iroi).AR_VC_con_vis;nWAV(iroi).AL_VR_con_vis;nWAV(iroi).AR_VL_con_vis]);
            nWAV(iroi).pooldisp_con_vis_sem = ...
                nWAV(iroi).pooldisp_con_vis_std/sqrt(length([nWAV(iroi).AC_VL_con_vis;nWAV(iroi).AC_VR_con_vis;nWAV(iroi).AL_VC_con_vis;nWAV(iroi).AR_VC_con_vis;nWAV(iroi).AL_VR_con_vis;nWAV(iroi).AR_VL_con_vis]));    
        end % roi        
        nWAV_resp(iResp).resp = nWAV;        
    end % resp    
    nWAV_group(iSubj).subj = nWAV_resp;    
    disp('...');
    disp([subjID ' done']);    
end % subj

% Save
cd(['E:\AMBRA\UoB\Data\MAMSI_MRI\group\MVPA\' scaling]);
save('group_nWAV', 'nWAV_group');

%% Create matrices for permutation testing

for iroi = 1:length(roi)
    results = [];
    for iResp = 1:length(Resp_list)
        results_temp = [subj{iResp,iroi,3},subj{iResp,iroi,5},subj{iResp,iroi,4},subj{iResp,iroi,6}];
        results = cat(2,results, results_temp);
    end
    
    for iResp = 1:length(Resp_list)
        results_temp = [subj{iResp,iroi,7},subj{iResp,iroi,8}];
        results = cat(2,results, results_temp);
    end
    save([roi{iroi} '_nWAV'], 'results');    
end

%% Compute group nWAV

n = sqrt(length(subjID_list)); % sem calculation

subj = cell(length(Resp_list),length(roi),8);
roi_mean = cell(length(Resp_list),length(roi),8);
roi_std = cell(length(Resp_list),length(roi),8);
roi_sem = cell(length(Resp_list),length(roi),8);

for iResp = 1:length(Resp_list)
    
    for iroi = 1:length(roi)
        
        for iSubj = 1:length(subjID_list)
            
            subj{iResp,iroi,1} = cat(1, subj{iResp,iroi,1}, cmb_group(iSubj).subj(iResp).resp(iroi).con_aud_mean);
            subj{iResp,iroi,2} = cat(1, subj{iResp,iroi,2}, cmb_group(iSubj).subj(iResp).resp(iroi).con_vis_mean);
            
            subj{iResp,iroi,3} = cat(1, subj{iResp,iroi,3}, cmb_group(iSubj).subj(iResp).resp(iroi).low_con_aud_mean);
            subj{iResp,iroi,4} = cat(1, subj{iResp,iroi,4}, cmb_group(iSubj).subj(iResp).resp(iroi).low_con_vis_mean);
            subj{iResp,iroi,5} = cat(1, subj{iResp,iroi,5}, cmb_group(iSubj).subj(iResp).resp(iroi).high_con_aud_mean);
            subj{iResp,iroi,6} = cat(1, subj{iResp,iroi,6}, cmb_group(iSubj).subj(iResp).resp(iroi).high_con_vis_mean);
            
            subj{iResp,iroi,7} = cat(1, subj{iResp,iroi,7}, cmb_group(iSubj).subj(iResp).resp(iroi).pooldisp_con_aud_mean);
            subj{iResp,iroi,8} = cat(1, subj{iResp,iroi,8}, cmb_group(iSubj).subj(iResp).resp(iroi).pooldisp_con_vis_mean);
            
        end % subj
        
        roi_mean{iResp,iroi,1} = mean(subj{iResp,iroi,1});
        roi_std{iResp,iroi,1} = std(subj{iResp,iroi,1});
        roi_sem{iResp,iroi,1} = roi_std{iResp,iroi,1}/n;
        
        roi_mean{iResp,iroi,2} = mean(subj{iResp,iroi,2});
        roi_std{iResp,iroi,2} = std(subj{iResp,iroi,2});
        roi_sem{iResp,iroi,2} = roi_std{iResp,iroi,2}/n;
        
        roi_mean{iResp,iroi,3} = mean(subj{iResp,iroi,3});
        roi_std{iResp,iroi,3} = std(subj{iResp,iroi,3});
        roi_sem{iResp,iroi,3} = roi_std{iResp,iroi,3}/n;
        
        roi_mean{iResp,iroi,4} = mean(subj{iResp,iroi,4});
        roi_std{iResp,iroi,4} = std(subj{iResp,iroi,4});
        roi_sem{iResp,iroi,4} = roi_std{iResp,iroi,4}/n;
        
        roi_mean{iResp,iroi,5} = mean(subj{iResp,iroi,5});
        roi_std{iResp,iroi,5} = std(subj{iResp,iroi,5});
        roi_sem{iResp,iroi,5} = roi_std{iResp,iroi,5}/n;
        
        roi_mean{iResp,iroi,6} = mean(subj{iResp,iroi,6});
        roi_std{iResp,iroi,6} = std(subj{iResp,iroi,6});
        roi_sem{iResp,iroi,6} = roi_std{iResp,iroi,6}/n;
        
        roi_mean{iResp,iroi,7} = mean(subj{iResp,iroi,7});
        roi_std{iResp,iroi,7} = std(subj{iResp,iroi,7});
        roi_sem{iResp,iroi,7} = roi_std{iResp,iroi,7}/n;
        
        roi_mean{iResp,iroi,8} = mean(subj{iResp,iroi,8});
        roi_std{iResp,iroi,8} = std(subj{iResp,iroi,8});
        roi_sem{iResp,iroi,8} = roi_std{iResp,iroi,8}/n;        
    end % roi    
end % resp

%% Figure

% colors
cols.yelAud1  = [165 0 33]/255;
cols.yelAud2  = [255 153 153]/255;
cols.blueVis1 = [153 204 255]/255;
cols.blueVis2 = [0 51 204]/255;

indroi=1:5;

% size
positionXY = [0, 0, 450, 400];
figure('color', [1 1 1], 'Position', positionXY);
ind_x=[1.2 2.2; 4.6 5.6; 8.0 9.0; 11.4 12.4; 14.8 15.8];

for iroi=1:length(indroi)
    plot_matrix_mean2=[roi_mean{1,indroi(iroi),8} roi_mean{1,indroi(iroi),7}; roi_mean{2,indroi(iroi),8} roi_mean{2,indroi(iroi),7}];
    for i=1:2
        hold on
        plot(ind_x(iroi,:),plot_matrix_mean2(i,:),...
            'Color','k',...
            'LineStyle','-',...
            'LineWidth',1.5)
    end
end

conCols=[cols.yelAud2;cols.yelAud1;cols.blueVis2;cols.blueVis1];
conLineStyles={
    '-'
    '-'
    '-'
    '-'
    };

plot_matrix_mean = [roi_mean{1,indroi,8}; roi_mean{1,indroi,7}; roi_mean{2,indroi,8}; roi_mean{2,indroi,7}];
plot_matrix_sem = [roi_sem{1,indroi,8}; roi_sem{1,indroi,7}; roi_sem{2,indroi,8}; roi_sem{2,indroi,7}];

for iroi=1:length(indroi)
    a = [ind_x(iroi,:),ind_x(iroi,:)];
    for i=1:4
        if i<3
            k=0.05;
        else
            k=-0.05;
        end
        hold on
        line([a(i)+k a(i)+k],[plot_matrix_mean(i,iroi)+plot_matrix_sem(i,iroi) ...
            plot_matrix_mean(i,iroi)-plot_matrix_sem(i,iroi)],...
            'Color',conCols(i,:),'LineWidth',1.5);hold on;
        plot(a(i)+k,plot_matrix_mean(i,iroi),'o',...
            'Color',conCols(i,:),...
            'LineStyle',conLineStyles{i,:},...
            'LineWidth',1.5,'MarkerSize',4,...
            'MarkerEdgeColor',conCols(i,:),...
            'MarkerFaceColor',conCols(i,:));hold on;
    end    
end

% plot settings
xl = [0.5 ind_x(end)+0.5]; xlim(xl);
set(gca,'FontName', 'Helvetica');
set(gca,'FontSize', 12);
set(gca, 'XTick', [1.2 2.2 4.6 5.6 8.0 9.0 11.4 12.4 14.8 15.8]);
set(gca, 'XTickLabel', []);
set(gca,'TickLength', [0.01 0.01]);
set(gca,'LineWidth',1.2)
