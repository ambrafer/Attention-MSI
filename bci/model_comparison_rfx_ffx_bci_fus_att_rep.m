%% RFX and FFX analysis on BIC
% 1. BCI attention 2. BCI no attention
% 3. fusion attention 4. fusion no attention
% 5. fusion (report) no attention 6. fusion (report) attention

clear;
close all;
clc;

addpath(genpath('E:\AMBRA\UoB\Programs\spm12'));

load('group_bci_aic_bic_attmod.mat')
bci_att_bic_group = bic_group_mod;
load('group_aic_bic.mat')
bci_noatt_bic_group = bic_group_null;
load('group_fus_aic_bic_attmod.mat')
fus_att_bic_group = bic_group_mod;
load('group_fus_aic_bic.mat')
fus_noatt_bic_group = bic_group_null;
load('group_fus_aic_bic_repmod.mat')
fus_rep_bic_group = bic_group_repmod;
load('group_fus_aic_bic_mod_switch.mat')
fus_attrep_bic_group = bic_group_mod_switch;
subjList = numel(bic_group_null);

%% Create matrix for spm rfx (cols = 6 models, row = 12 subjects)

rfx_bciA_bcinoA_fusA_fusnoA_fusR_fusRA = zeros(subjList,6);
%rfx_fusatt_fusnoatt = zeros(subjList,2);
rfx_results = cell(1,1); % 1 comparison

% Comparison 1:
% col 1: bci_att
rfx_bciA_bcinoA_fusA_fusnoA_fusR_fusRA(:,1) = bci_att_bic_group;
% col 2: bci_noatt
rfx_bciA_bcinoA_fusA_fusnoA_fusR_fusRA(:,2) = bci_noatt_bic_group;
% col 3: fus_att
rfx_bciA_bcinoA_fusA_fusnoA_fusR_fusRA(:,3) = fus_att_bic_group;
% col 4: fus_noatt
rfx_bciA_bcinoA_fusA_fusnoA_fusR_fusRA(:,4) = fus_noatt_bic_group;
% col 5: fus_rep
rfx_bciA_bcinoA_fusA_fusnoA_fusR_fusRA(:,5) = fus_rep_bic_group;
% col 6: fus_repatt
rfx_bciA_bcinoA_fusA_fusnoA_fusR_fusRA(:,6) = fus_attrep_bic_group;

BIC_group=rfx_bciA_bcinoA_fusA_fusnoA_fusR_fusRA(:,[4 3 5 6 2 1]);
save('BIC.mat','BIC_group');

[alpha,exp_r,xp,pxp,bor] = spm_BMS(rfx_bciA_bcinoA_fusA_fusnoA_fusR_fusRA);
rfx_results{1} = [alpha,exp_r,xp,pxp,bor]; % 3 x alpha + 3 x exp_r + 3 x xp + 3 x pxp + bor

fig1=figure('Position', [0 0 150 200]);
bar(pxp(1:4),'FaceColor',[166 166 166]/255,'LineWidth',1.2);
set(gca,'TickLength', [0.02 0.02]);
set(gca,'LineWidth',1.2);
xl = [0.2 3.8]; xlim(xl);
set(gca, 'XTick', []);
set(gca, 'XTickLabel', []);
box off

pxp_toplot = [pxp(4) pxp(3); pxp(5) pxp(6); pxp(2) pxp(1);];
map_vector=0:0.01:1;
map=repmat(map_vector',1,3);
positionXY = [0, 0, 240, 240];
figure('Position', positionXY);
imagesc(pxp_toplot);
colormap(map);
colorbar;
set(gca,'XTick',[1 2 3],...
    'YTick',[1 2 3],...
    'FontName', 'Helvetica', ...
    'TickLength', [0 0], ...
    'YTickLabel',[], ...
    'XTickLabel',[]);

saveas(gcf, 'pxp_bci_fus_mod_rep_plot_av09_MRI', 'emf');

% % Comparison 2:
% % col 1: fus_att
% rfx_fusatt_fusnoatt(:,1) = fus_att_bic_group;
% % col 2: fus_noatt
% rfx_fusatt_fusnoatt(:,2) = fus_noatt_bic_group;
% 
% [alpha,exp_r,xp,pxp,bor] = spm_BMS(rfx_fusatt_fusnoatt);
% rfx_results{2} = [alpha,exp_r,xp,pxp,bor]; % 3 x alpha + 3 x exp_r + 3 x xp + 3 x pxp + bor

%% Create matrix for ffx (cols = 6 models, row = 12 subjects)

ffx_comp = zeros(subjList,3);
ffx_comp(:,1) = bci_att_bic_group-fus_noatt_bic_group;
ffx_comp(:,2) = bci_att_bic_group-fus_att_bic_group;
ffx_comp(:,3) = bci_att_bic_group-fus_rep_bic_group;
ffx_comp(:,4) = bci_att_bic_group-fus_attrep_bic_group;
ffx_comp(:,5) = bci_att_bic_group-bci_noatt_bic_group;

fig2=figure;
line([0 13],[0 0], 'Color','k','LineStyle', '-');
hold on;
plot(1:subjList, ffx_comp(:,1),'g*-');
hold on;
plot(1:subjList, ffx_comp(:,2),'b*-');
hold on;
plot(1:subjList, ffx_comp(:,3),'r*-');
hold on;
plot(1:subjList, ffx_comp(:,4),'c*-');
hold on;
plot(1:subjList, ffx_comp(:,5),'k*-');
xlim([0 13]);
set(gca, 'XTick', (1:1:12));
set(gca, 'XTickLabel', (1:1:12));
set(gca,'FontName', 'Helvetica');
set(gca,'FontSize', 12);

% ffx
groupBIC = sum(ffx_comp);
