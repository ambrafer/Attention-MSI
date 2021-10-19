%% Plot BCI parameters of winning model (BCI model with separate sigmas as a function of attention)

%% Get data
clear;
close all;
clc;

cd('E:\AMBRA\UoB\Data\MAMSI_MRI\group\behav\BCI\single_fit');
file = fullfile('group_bciSimulations_attmod_best10.mat');
load(file,'group_best_results_mod');

subjN = numel(group_best_results_mod);

aud = [];
vis = [];

for iSubj = 1:subjN
    
    aud = cat(1,aud, group_best_results_mod{iSubj}.bci_details(1).parameters);
    vis = cat(1,vis, group_best_results_mod{iSubj}.bci_details(10).parameters);
    
end % End of loop over subjects

% Prepare matrix for SPSS
params_matrix_attmod = [aud, vis(:,[3 4])];
save('group_params_attmod_stats.mat','params_matrix_attmod');

%% Calculate descriptive statistics across subjects
% mean
aud_mean = mean(aud);
vis_mean = mean(vis);
% std
aud_std = std(aud);
vis_std = std(vis);
% sem
aud_sem = aud_std/sqrt(subjN);
vis_sem = vis_std/sqrt(subjN);

%% Figure
% colors
cols.grey = [166 166 166]/255;
cols.aud  = [166 166 166]/255;
cols.vis  = [166 166 166]/255;

% x axis location and shift
loc = [1 2];
loc2 = [1.5 5];
a = 0.145;

% nr columns
row = 1; % number of figure rows
col = 4; % number of figure columns

% size
positionXY = [0, 0, 800, 400];
base.left   = 60;% pixels from left
base.bottom = 60;% pixels from bottom
jump.left   = 10; % pixels between plots laterally
jump.bottom = 10; % pixels between plots vertically
figsize.x   = (positionXY(3)-base.left*2-jump.left*(col-1))/col;
figsize.y   = (positionXY(4)-base.bottom*2);

Cols=[cols.vis;cols.aud];
LineStyles={
    '-'
    '-'
    };

% locations on x axis
fig1 = figure('color', [1 1 1], 'Position', positionXY);

for c = 1:col
    % determine plot id and position
    hsp = subplot(1,col,c);
    switch c
        case 1
            plot_matrix_mean=aud_mean(1);
            plot_matrix_sem=aud_sem(1);
            for i=1
                hold on
                bar(i+0.5,plot_matrix_mean(i),0.6,...
                    'FaceColor',cols.grey,...
                    'LineStyle','-',...
                    'LineWidth',1.5)
                line([i+0.5 i+0.5],[plot_matrix_mean(i)+plot_matrix_sem(i) ...
                    plot_matrix_mean(i)-plot_matrix_sem(i)],...
                    'Color',[0 0 0],'LineWidth',1.5)
            end
        case 2
            plot_matrix_mean=aud_mean(2);
            plot_matrix_sem=aud_sem(2);
            for i=1
                hold on
                bar(i+0.5,plot_matrix_mean(i),0.6,...
                    'FaceColor',cols.grey,...
                    'LineStyle','-',...
                    'LineWidth',1.5)
                line([i+0.5 i+0.5],[plot_matrix_mean(i)+plot_matrix_sem(i) ...
                    plot_matrix_mean(i)-plot_matrix_sem(i)],...
                    'Color',[0 0 0],'LineWidth',1.5)
            end
        case 3
            plot_matrix_mean=[vis_mean(3) aud_mean(3)];
            plot_matrix_sem=[vis_sem(3) aud_sem(3)];
            for i=1:2
                hold on
                bar(i,plot_matrix_mean(i),0.6,...
                    'FaceColor',Cols(i,:),...
                    'LineStyle','-',...
                    'LineWidth',1.5)
                line([i i],[plot_matrix_mean(i)+plot_matrix_sem(i) ...
                    plot_matrix_mean(i)-plot_matrix_sem(i)],...
                    'Color',[0 0 0],'LineWidth',1.5)
            end
        case 4
            plot_matrix_mean=[vis_mean(4) aud_mean(4)];
            plot_matrix_sem=[vis_sem(4) aud_sem(4)];
            for i=1:2
                hold on
                bar(i,plot_matrix_mean(i),0.6,...
                    'FaceColor',Cols(i,:),...
                    'LineStyle','-',...
                    'LineWidth',1.5)
                line([i i],[plot_matrix_mean(i)+plot_matrix_sem(i) ...
                    plot_matrix_mean(i)-plot_matrix_sem(i)],...
                    'Color',[0 0 0],'LineWidth',1.5)
            end
    end
    
    % set x and y axes and ticks    
    if c == 1
        xl = [0.3 2.7]; xlim(xl);
        yl = [0 0.8]; ylim(yl);
        ticksY = (0:0.2:1);
        set(gca, 'YTick', ticksY);
        set(gca, 'YTickLabel', (0:0.2:1));
    elseif c == 2
        xl = [0.3 2.7]; xlim(xl);
        yl = [0 24]; ylim(yl);
        ticksY = (0:6:30);
        set(gca, 'YTick', ticksY);
        set(gca, 'YTickLabel', (0:6:30));
    elseif c == 3
        xl = [0.3 2.7]; xlim(xl);
        yl = [0 10]; ylim(yl);
        ticksY = (0:1:10);
        set(gca, 'YTick', ticksY);
        set(gca, 'YTickLabel', (0:1:10));
    elseif c == 4
        xl = [0.3 2.7]; xlim(xl);
        yl = [0 10]; ylim(yl);
        ticksY = (0:1:10);
        set(gca, 'YTick', ticksY);
        set(gca, 'YTickLabel', (0:1:10));
    end
    set(gca,'FontName', 'Helvetica');
    set(gca,'FontSize', 15);
    set(gca, 'XTickLabel', '');
    set(gca, 'XTick', []);
    set(gca,'TickLength', [0.015 0.015]);
    set(gca,'LineWidth',1.5)
end
