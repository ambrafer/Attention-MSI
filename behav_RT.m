%% Behavioural analysis: RT

clear;
close all;
clc;

% Initialize for group analysis
RT = [];

subjID=dir('*Exp_All_Sessions*.mat');

for iSubj = 1:length(subjID)
    
    % Get a list of .mat files to put together
    Files = dir('*Exp_All_Sessions*.mat');
    
    % Initialize
    RT_subj = struct;
    
    % Load current subject dataset
    load(Files.name,'tdata');
    
    %% Statistics for auditory (avdisp)
    
    % attended
    RT_subj.aud.att.AV0 = median(tdata.ResponseTime(tdata.AVdisparity==0 & ...
        strcmp(tdata.ResponseModality,'Aud') & strcmp(tdata.AttentionValidity,'Attended') & ...
        tdata.IncorResp==0 & tdata.MissedResp==0 & tdata.AntResp==0));
    
    RT_subj.aud.att.AV9 = median(tdata.ResponseTime(tdata.AVdisparity==9 & ...
        strcmp(tdata.ResponseModality,'Aud') & strcmp(tdata.AttentionValidity,'Attended') & ...
        tdata.IncorResp==0 & tdata.MissedResp==0 & tdata.AntResp==0));
    
    RT_subj.aud.att.AV18 = median(tdata.ResponseTime(tdata.AVdisparity==18 & ...
        strcmp(tdata.ResponseModality,'Aud') & strcmp(tdata.AttentionValidity,'Attended') & ...
        tdata.IncorResp==0 & tdata.MissedResp==0 & tdata.AntResp==0));
    
    RT_subj.aud.att.AVinc = median(tdata.ResponseTime(tdata.AVdisparity~=0 & ...
        strcmp(tdata.ResponseModality,'Aud') & strcmp(tdata.AttentionValidity,'Attended') & ...
        tdata.IncorResp==0 & tdata.MissedResp==0 & tdata.AntResp==0));
    
    % unattended
    RT_subj.aud.unatt.AV0 = median(tdata.ResponseTime(tdata.AVdisparity==0 & ...
        strcmp(tdata.ResponseModality,'Aud') & strcmp(tdata.AttentionValidity,'Unattended') & ...
        tdata.IncorResp==0 & tdata.MissedResp==0 & tdata.AntResp==0));
    
    RT_subj.aud.unatt.AV9 = median(tdata.ResponseTime(tdata.AVdisparity==9 & ...
        strcmp(tdata.ResponseModality,'Aud') & strcmp(tdata.AttentionValidity,'Unattended') & ...
        tdata.IncorResp==0 & tdata.MissedResp==0 & tdata.AntResp==0));
    
    RT_subj.aud.unatt.AV18 = median(tdata.ResponseTime(tdata.AVdisparity==18 & ...
        strcmp(tdata.ResponseModality,'Aud') & strcmp(tdata.AttentionValidity,'Unattended') & ...
        tdata.IncorResp==0 & tdata.MissedResp==0 & tdata.AntResp==0));
    
    RT_subj.aud.unatt.AVinc = median(tdata.ResponseTime(tdata.AVdisparity~=0 & ...
        strcmp(tdata.ResponseModality,'Aud') & strcmp(tdata.AttentionValidity,'Unattended') & ...
        tdata.IncorResp==0 & tdata.MissedResp==0 & tdata.AntResp==0));
    
    %% Statistics for auditory (pooled over avdisp)
    
    RT_subj.aud.att.AVpool = median(tdata.ResponseTime(strcmp(tdata.ResponseModality,'Aud') & strcmp(tdata.AttentionValidity,'Attended') & ...
        tdata.IncorResp==0 & tdata.MissedResp==0 & tdata.AntResp==0));
    RT_subj.aud.unatt.AVpool = median(tdata.ResponseTime(strcmp(tdata.ResponseModality,'Aud') & strcmp(tdata.AttentionValidity,'Unattended') & ...
        tdata.IncorResp==0 & tdata.MissedResp==0 & tdata.AntResp==0));
    
    
    %% Statistics for visual (avdisp)
    
    % attended
    RT_subj.vis.att.AV0 = median(tdata.ResponseTime(tdata.AVdisparity==0 & ...
        strcmp(tdata.ResponseModality,'Vis') & strcmp(tdata.AttentionValidity,'Attended') & ...
        tdata.IncorResp==0 & tdata.MissedResp==0 & tdata.AntResp==0));
    
    RT_subj.vis.att.AV9 = median(tdata.ResponseTime(tdata.AVdisparity==9 & ...
        strcmp(tdata.ResponseModality,'Vis') & strcmp(tdata.AttentionValidity,'Attended') & ...
        tdata.IncorResp==0 & tdata.MissedResp==0 & tdata.AntResp==0));
    
    RT_subj.vis.att.AV18 = median(tdata.ResponseTime(tdata.AVdisparity==18 & ...
        strcmp(tdata.ResponseModality,'Vis') & strcmp(tdata.AttentionValidity,'Attended') & ...
        tdata.IncorResp==0 & tdata.MissedResp==0 & tdata.AntResp==0));
    
    RT_subj.vis.att.AVinc = median(tdata.ResponseTime(tdata.AVdisparity~=0 & ...
        strcmp(tdata.ResponseModality,'Vis') & strcmp(tdata.AttentionValidity,'Attended') & ...
        tdata.IncorResp==0 & tdata.MissedResp==0 & tdata.AntResp==0));
    
    % unattended
    RT_subj.vis.unatt.AV0 = median(tdata.ResponseTime(tdata.AVdisparity==0 & ...
        strcmp(tdata.ResponseModality,'Vis') & strcmp(tdata.AttentionValidity,'Unattended') & ...
        tdata.IncorResp==0 & tdata.MissedResp==0 & tdata.AntResp==0));
    
    RT_subj.vis.unatt.AV9 = median(tdata.ResponseTime(tdata.AVdisparity==9 & ...
        strcmp(tdata.ResponseModality,'Vis') & strcmp(tdata.AttentionValidity,'Unattended') & ...
        tdata.IncorResp==0 & tdata.MissedResp==0 & tdata.AntResp==0));
    
    RT_subj.vis.unatt.AV18 = median(tdata.ResponseTime(tdata.AVdisparity==18 & ...
        strcmp(tdata.ResponseModality,'Vis') & strcmp(tdata.AttentionValidity,'Unattended') & ...
        tdata.IncorResp==0 & tdata.MissedResp==0 & tdata.AntResp==0));
    
    RT_subj.vis.unatt.AVinc = median(tdata.ResponseTime(tdata.AVdisparity~=0 & ...
        strcmp(tdata.ResponseModality,'Vis') & strcmp(tdata.AttentionValidity,'Unattended') & ...
        tdata.IncorResp==0 & tdata.MissedResp==0 & tdata.AntResp==0));
    
    %% Statistics for visual (pooled over avdisp)
    
    RT_subj.vis.att.AVpool = median(tdata.ResponseTime(strcmp(tdata.ResponseModality,'Vis') & strcmp(tdata.AttentionValidity,'Attended') & ...
        tdata.IncorResp==0 & tdata.MissedResp==0 & tdata.AntResp==0));
    RT_subj.vis.unatt.AVpool = median(tdata.ResponseTime(strcmp(tdata.ResponseModality,'Vis') & strcmp(tdata.AttentionValidity,'Unattended') & ...
        tdata.IncorResp==0 & tdata.MissedResp==0 & tdata.AntResp==0));
    
    %% Group current subjects data
    
    RT = cat(1, RT, [RT_subj.aud.att.AV0, RT_subj.aud.att.AV9, RT_subj.aud.att.AV18,RT_subj.aud.unatt.AV0, RT_subj.aud.unatt.AV9, RT_subj.aud.unatt.AV18, ...
        RT_subj.vis.unatt.AV0, RT_subj.vis.unatt.AV9, RT_subj.vis.unatt.AV18, RT_subj.vis.att.AV0, RT_subj.vis.att.AV9, RT_subj.vis.att.AV18, ...
        RT_subj.aud.att.AVinc, RT_subj.aud.unatt.AVinc, RT_subj.vis.unatt.AVinc, RT_subj.vis.att.AVinc,...
        RT_subj.aud.att.AVpool, RT_subj.aud.unatt.AVpool, RT_subj.vis.unatt.AVpool, RT_subj.vis.att.AVpool]);
    
end

%% Group data
RT_group_mean = mean(RT);
RT_group_std = std(RT);
RT_group_sem = std(RT)/sqrt(size(RT,1));

%% Figure

cols.yelAud1  = [165 0 33]/255;
cols.yelAud2  = [255 153 153]/255;
cols.blueVis1 = [153 204 255]/255;
cols.blueVis2 = [0 51 204]/255;

% size
positionXY = [0, 0, 150, 200];
figure('color', [1 1 1], 'Position', positionXY);
conLineStyles={
    ':'
    '-'
    ':'
    '-'
    };
plot_matrix_mean2=[RT_group_mean(10) RT_group_mean(7);RT_group_mean(16) RT_group_mean(15);...
    RT_group_mean(4) RT_group_mean(1); RT_group_mean(14) RT_group_mean(13);];

for i=1:4
    hold on
    if i==1
        k=-0.08;
    elseif i==2
        k=-0.04;
    elseif i==3
        k=0.04;
    elseif i==4
        k=0.08;
    end
    plot((1:2)+k,plot_matrix_mean2(i,:),...
        'Color','k',...
        'LineStyle',conLineStyles{i},...
        'LineWidth',1.5)
end

conCols=[cols.blueVis2;cols.blueVis1;cols.blueVis2;cols.blueVis1;...
    cols.yelAud2;cols.yelAud1;cols.yelAud2;cols.yelAud1];

plot_matrix_mean = [RT_group_mean(10); RT_group_mean(7); RT_group_mean(16); RT_group_mean(15);...
    RT_group_mean(4); RT_group_mean(1); RT_group_mean(14); RT_group_mean(13)];
plot_matrix_sem = [RT_group_sem(10); RT_group_sem(7); RT_group_sem(16); RT_group_sem(15);...
    RT_group_sem(4); RT_group_sem(1); RT_group_sem(14); RT_group_sem(13)];

a=repmat([1 2],1,4);

for i=1:8
    if i==1 || i==2
        k=-0.08;
    elseif i==3 || i==4
        k=-0.04;
    elseif i==5 || i==6
        k=0.04;
    elseif i==7 || i==8
        k=0.08;
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
%yl = [0 1]; ylim(yl);
% adjust x and y ticks
set(gca,'FontName', 'Helvetica');
set(gca,'FontSize', 12);
set(gca,'XTick', [1 2]);
set(gca,'XTickLabel', []);
%set(gca,'YTick', (0:0.1:1));
set(gca,'YTickLabel', (40:10:80)); %% fake values
set(gca,'TickLength', [0.02 0.02]);
set(gca,'LineWidth',1.2)
