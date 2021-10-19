%% Behavioural analysis: WAV

clear;
close all;
clc;

% Initialize for group analysis
WAV_con = [];

% Exp info
Resp = 'Aud'; % response modality: Aud or Vis
subjID=dir('*Exp_All_Sessions*.mat');

for iSubj = 1:length(subjID)
    
    % Get a list of .mat files to analyse
    cd(dataPath_ve);
    Files = dir('*Exp_All_Sessions*.mat');
    
    % Load current subject dataset
    load(Files.name,'tdata');
    
    %% Compute WAV
    
    WAV_AL_VC_con_aud = (tdata.Response(strcmp(tdata.ResponseModality,Resp) & ...
        strcmp(tdata.AttentionModality,'Aud') & ...
        tdata.TargetLoc==-9 & tdata.NonTargetLoc==0 & ...
        tdata.IncorResp==0 & tdata.MissedResp==0 & tdata.AntResp==0,:) ...
        -(cong_tot_mean(1)))/(cong_tot_mean(2)-(cong_tot_mean(1)));
    
    WAV_AL_VC_con_aud_mean = mean(WAV_AL_VC_con_aud);
    WAV_AL_VC_con_aud_std = std(WAV_AL_VC_con_aud);
    WAV_AL_VC_con_aud_sem = WAV_AL_VC_con_aud_std/sqrt(length(WAV_AL_VC_con_aud));
    
    WAV_AL_VC_con_vis = (tdata.Response(strcmp(tdata.ResponseModality,Resp) & ...
        strcmp(tdata.AttentionModality,'Vis') & ...
        tdata.TargetLoc==-9 & tdata.NonTargetLoc==0 & ...
        tdata.IncorResp==0 & tdata.MissedResp==0 & tdata.AntResp==0,:) ...
        -(cong_tot_mean(1)))/(cong_tot_mean(2)-(cong_tot_mean(1)));
    
    WAV_AL_VC_con_vis_mean = mean(WAV_AL_VC_con_vis);
    WAV_AL_VC_con_vis_std = std(WAV_AL_VC_con_vis);
    WAV_AL_VC_con_vis_sem = WAV_AL_VC_con_vis_std/sqrt(length(WAV_AL_VC_con_vis));
    
    WAV_AL_VR_con_aud = (tdata.Response(strcmp(tdata.ResponseModality,Resp) & ...
        strcmp(tdata.AttentionModality,'Aud') & ...
        tdata.TargetLoc==-9 & tdata.NonTargetLoc==9 & ...
        tdata.IncorResp==0 & tdata.MissedResp==0 & tdata.AntResp==0,:) ...
        -(cong_tot_mean(1)))/(cong_tot_mean(3)-(cong_tot_mean(1)));
    
    WAV_AL_VR_con_aud_mean = mean(WAV_AL_VR_con_aud);
    WAV_AL_VR_con_aud_std = std(WAV_AL_VR_con_aud);
    WAV_AL_VR_con_aud_sem = WAV_AL_VR_con_aud_std/sqrt(length(WAV_AL_VR_con_aud));
    
    WAV_AL_VR_con_vis = (tdata.Response(strcmp(tdata.ResponseModality,Resp) & ...
        strcmp(tdata.AttentionModality,'Vis') & ...
        tdata.TargetLoc==-9 & tdata.NonTargetLoc==9 & ...
        tdata.IncorResp==0 & tdata.MissedResp==0 & tdata.AntResp==0,:) ...
        -(cong_tot_mean(1)))/(cong_tot_mean(3)-(cong_tot_mean(1)));
    
    WAV_AL_VR_con_vis_mean = mean(WAV_AL_VR_con_vis);
    WAV_AL_VR_con_vis_std = std(WAV_AL_VR_con_vis);
    WAV_AL_VR_con_vis_sem = WAV_AL_VR_con_vis_std/sqrt(length(WAV_AL_VR_con_vis));
    
    WAV_AC_VL_con_aud = (tdata.Response(strcmp(tdata.ResponseModality,Resp) & ...
        strcmp(tdata.AttentionModality,'Aud') & ...
        tdata.TargetLoc==0 & tdata.NonTargetLoc==-9 & ...
        tdata.IncorResp==0 & tdata.MissedResp==0 & tdata.AntResp==0,:) ...
        -(cong_tot_mean(2)))/(cong_tot_mean(1)-(cong_tot_mean(2)));
    
    WAV_AC_VL_con_aud_mean = mean(WAV_AC_VL_con_aud);
    WAV_AC_VL_con_aud_std = std(WAV_AC_VL_con_aud);
    WAV_AC_VL_con_aud_sem = WAV_AC_VL_con_aud_std/sqrt(length(WAV_AC_VL_con_aud));
    
    WAV_AC_VL_con_vis = (tdata.Response(strcmp(tdata.ResponseModality,Resp) & ...
        strcmp(tdata.AttentionModality,'Vis') & ...
        tdata.TargetLoc==0 & tdata.NonTargetLoc==-9 & ...
        tdata.IncorResp==0 & tdata.MissedResp==0 & tdata.AntResp==0,:) ...
        -(cong_tot_mean(2)))/(cong_tot_mean(1)-(cong_tot_mean(2)));
    
    WAV_AC_VL_con_vis_mean = mean(WAV_AC_VL_con_vis);
    WAV_AC_VL_con_vis_std = std(WAV_AC_VL_con_vis);
    WAV_AC_VL_con_vis_sem = WAV_AC_VL_con_vis_std/sqrt(length(WAV_AC_VL_con_vis));
    
    WAV_AC_VR_con_aud = (tdata.Response(strcmp(tdata.ResponseModality,Resp) & ...
        strcmp(tdata.AttentionModality,'Aud') & ...
        tdata.TargetLoc==0 & tdata.NonTargetLoc==9 & ...
        tdata.IncorResp==0 & tdata.MissedResp==0 & tdata.AntResp==0,:) ...
        -(cong_tot_mean(2)))/(cong_tot_mean(3)-(cong_tot_mean(2)));
    
    WAV_AC_VR_con_aud_mean = mean(WAV_AC_VR_con_aud);
    WAV_AC_VR_con_aud_std = std(WAV_AC_VR_con_aud);
    WAV_AC_VR_con_aud_sem = WAV_AC_VR_con_aud_std/sqrt(length(WAV_AC_VR_con_aud));
    
    WAV_AC_VR_con_vis = (tdata.Response(strcmp(tdata.ResponseModality,Resp) & ...
        strcmp(tdata.AttentionModality,'Vis') & ...
        tdata.TargetLoc==0 & tdata.NonTargetLoc==9 & ...
        tdata.IncorResp==0 & tdata.MissedResp==0 & tdata.AntResp==0,:) ...
        -(cong_tot_mean(2)))/(cong_tot_mean(3)-(cong_tot_mean(2)));
    
    WAV_AC_VR_con_vis_mean = mean(WAV_AC_VR_con_vis);
    WAV_AC_VR_con_vis_std = std(WAV_AC_VR_con_vis);
    WAV_AC_VR_con_vis_sem = WAV_AC_VR_con_vis_std/sqrt(length(WAV_AC_VR_con_vis));
    
    WAV_AR_VL_con_aud = (tdata.Response(strcmp(tdata.ResponseModality,Resp) & ...
        strcmp(tdata.AttentionModality,'Aud') & ...
        tdata.TargetLoc==9 & tdata.NonTargetLoc==-9 & ...
        tdata.IncorResp==0 & tdata.MissedResp==0 & tdata.AntResp==0,:) ...
        -(cong_tot_mean(3)))/(cong_tot_mean(1)-(cong_tot_mean(3)));
    
    WAV_AR_VL_con_aud_mean = mean(WAV_AR_VL_con_aud);
    WAV_AR_VL_con_aud_std = std(WAV_AR_VL_con_aud);
    WAV_AR_VL_con_aud_sem = WAV_AR_VL_con_aud_std/sqrt(length(WAV_AR_VL_con_aud));
    
    WAV_AR_VL_con_vis = (tdata.Response(strcmp(tdata.ResponseModality,Resp) & ...
        strcmp(tdata.AttentionModality,'Vis') & ...
        tdata.TargetLoc==9 & tdata.NonTargetLoc==-9 & ...
        tdata.IncorResp==0 & tdata.MissedResp==0 & tdata.AntResp==0,:) ...
        -(cong_tot_mean(3)))/(cong_tot_mean(1)-(cong_tot_mean(3)));
    
    WAV_AR_VL_con_vis_mean = mean(WAV_AR_VL_con_vis);
    WAV_AR_VL_con_vis_std = std(WAV_AR_VL_con_vis);
    WAV_AR_VL_con_vis_sem = WAV_AR_VL_con_vis_std/sqrt(length(WAV_AR_VL_con_vis));
    
    WAV_AR_VC_con_aud = (tdata.Response(strcmp(tdata.ResponseModality,Resp) & ...
        strcmp(tdata.AttentionModality,'Aud') & ...
        tdata.TargetLoc==9 & tdata.NonTargetLoc==0 & ...
        tdata.IncorResp==0 & tdata.MissedResp==0 & tdata.AntResp==0,:) ...
        -(cong_tot_mean(3)))/(cong_tot_mean(2)-(cong_tot_mean(3)));
    
    WAV_AR_VC_con_aud_mean = mean(WAV_AR_VC_con_aud);
    WAV_AR_VC_con_aud_std = std(WAV_AR_VC_con_aud);
    WAV_AR_VC_con_aud_sem = WAV_AR_VC_con_aud_std/sqrt(length(WAV_AR_VC_con_aud));
    
    WAV_AR_VC_con_vis = (tdata.Response(strcmp(tdata.ResponseModality,Resp) & ...
        strcmp(tdata.AttentionModality,'Vis') & ...
        tdata.TargetLoc==9 & tdata.NonTargetLoc==0 & ...
        tdata.IncorResp==0 & tdata.MissedResp==0 & tdata.AntResp==0,:) ...
        -(cong_tot_mean(3)))/(cong_tot_mean(2)-(cong_tot_mean(3)));
    
    WAV_AR_VC_con_vis_mean = mean(WAV_AR_VC_con_vis);
    WAV_AR_VC_con_vis_std = std(WAV_AR_VC_con_vis);
    WAV_AR_VC_con_vis_sem = WAV_AR_VC_con_vis_std/sqrt(length(WAV_AR_VC_con_vis));
    
    % low spatial disparity
    % attended
    WAV_low_con_aud_mean = mean([WAV_AC_VL_con_aud;WAV_AC_VR_con_aud;WAV_AL_VC_con_aud;WAV_AR_VC_con_aud]);
    WAV_low_con_aud_std = std([WAV_AC_VL_con_aud;WAV_AC_VR_con_aud;WAV_AL_VC_con_aud;WAV_AR_VC_con_aud]);
    WAV_low_con_aud_sem = WAV_low_con_aud_std/sqrt(length([WAV_AC_VL_con_aud;WAV_AC_VR_con_aud;WAV_AL_VC_con_aud;WAV_AR_VC_con_aud]));
    % unattended
    WAV_low_con_vis_mean = mean([WAV_AC_VL_con_vis;WAV_AC_VR_con_vis;WAV_AL_VC_con_vis;WAV_AR_VC_con_vis]);
    WAV_low_con_vis_std = std([WAV_AC_VL_con_vis;WAV_AC_VR_con_vis;WAV_AL_VC_con_vis;WAV_AR_VC_con_vis]);
    WAV_low_con_vis_sem = WAV_low_con_vis_std/sqrt(length([WAV_AC_VL_con_vis;WAV_AC_VR_con_vis;WAV_AL_VC_con_vis;WAV_AR_VC_con_vis]));
    
    % high spatial disparity
    % attended
    WAV_high_con_aud_mean = mean([WAV_AL_VR_con_aud;WAV_AR_VL_con_aud]);
    WAV_high_con_aud_std = std([WAV_AL_VR_con_aud;WAV_AR_VL_con_aud]);
    WAV_high_con_aud_sem = WAV_high_con_aud_std/sqrt(length([WAV_AL_VR_con_aud;WAV_AR_VL_con_aud]));
    % unattended
    WAV_high_con_vis_mean = mean([WAV_AL_VR_con_vis;WAV_AR_VL_con_vis]);
    WAV_high_con_vis_std = std([WAV_AL_VR_con_vis;WAV_AR_VL_con_vis]);
    WAV_high_con_vis_sem = WAV_high_con_vis_std/sqrt(length([WAV_AL_VR_con_vis;WAV_AR_VL_con_vis]));
    
    WAV_con = cat(1, WAV_con, [WAV_low_con_aud_mean, WAV_high_con_aud_mean, WAV_low_con_vis_mean, WAV_high_con_vis_mean]);

end

if strcmp(Resp,'Vis')
    WAV_con=1-WAV_con;
end

%% Save for sign permutation testing
save([Resp '_WAV'], 'WAV_con');

WAV_con_group_mean = mean(WAV_con);
WAV_con_group_std = std(WAV_con);
WAV_con_group_sem = std(WAV_con)/sqrt(size(WAV_con,1));

%% Figure

cols.yelAud1  = [165 0 33]/255;
cols.yelAud2  = [255 153 153]/255;
cols.blueVis1 = [153 204 255]/255;
cols.blueVis2 = [0 51 204]/255;

% size
positionXY = [0, 0, 150, 400];
figure('color', [1 1 1], 'Position', positionXY);

conLineStyles={
    '-'
    '--'
    '-'
    '--'
    };

plot_matrix_mean2=[WAV_con_group_mean(3) WAV_con_group_mean(1); WAV_con_group_mean(4) WAV_con_group_mean(2);
    WAV_con_group_mean(7) WAV_con_group_mean(5); WAV_con_group_mean(8) WAV_con_group_mean(6)];
for i=1:4
    hold on
    if i<3
        k=0;%0.05;
    else
        k=0;%-0.05;
    end
    plot((1:2)+k,plot_matrix_mean2(i,:),...
        'Color','k',...
        'LineStyle',conLineStyles{i},...
        'LineWidth',1.5)
end

conCols=[cols.yelAud2;cols.yelAud1;cols.yelAud2;cols.yelAud1;cols.blueVis2;cols.blueVis1;cols.blueVis2;cols.blueVis1];

plot_matrix_mean = [WAV_con_group_mean(3); WAV_con_group_mean(1); WAV_con_group_mean(4); WAV_con_group_mean(2);
    WAV_con_group_mean(7); WAV_con_group_mean(5); WAV_con_group_mean(8); WAV_con_group_mean(6)];
plot_matrix_sem = [WAV_con_group_sem(3); WAV_con_group_sem(1); WAV_con_group_sem(4); WAV_con_group_sem(2);
    WAV_con_group_sem(7); WAV_con_group_sem(5); WAV_con_group_sem(8); WAV_con_group_sem(6)];

a=repmat([1 2],1,4);

for i=1:8
    if i<5
        k=0;%0.05;
    else
        k=0;%-0.05;
    end
    hold on
    line([a(i)+k a(i)+k],[plot_matrix_mean(i)+plot_matrix_sem(i) ...
        plot_matrix_mean(i)-plot_matrix_sem(i)],...
        'Color',conCols(i,:),'LineWidth',1.5);hold on;
    plot(a(i)+k,plot_matrix_mean(i),'o',...
        'Color',conCols(i,:),...
        'LineWidth',1.5,'MarkerSize',4,...
        'MarkerEdgeColor',conCols(i,:),...
        'MarkerFaceColor',conCols(i,:));hold on;
end

% adjust x and y ticks
xl = [0.7 2.3]; xlim(xl);
yl = [0.15 1.05]; ylim(yl);
% adjust x and y ticks
set(gca,'FontName', 'Helvetica');
set(gca,'FontSize', 12);
set(gca,'YAxisLocation','left');
set(gca,'XTick', [1 2]);
set(gca,'XTickLabel', {'attV';'attA'});
set(gca,'YTick', (0:0.1:1));
set(gca,'YTickLabel', (0:0.1:1));
set(gca,'TickLength', [0.01 0.01]);
set(gca,'LineWidth',1.2)
