%% MVPA analysis: group sign permutation test on decoded congruent spatial locations

close all;
clear;
clc;

% General info
Tail = 'right';
scaling = 'min0max1_concat_scale';
save_notes = 'results';
savePath=['E:\AMBRA\UoB\Data\MAMSI_MRI\group\MVPA\' scaling];

% Subject ID
subjID_list = {'sub-MA01'};
% rois list
roi_list = {'V1-3';'IPS0-2';'IPS3-4_SPL1';'TE1.0-1.1';'PT'};
rho_tot_roi = nan(length(subjID_list),length(roi_list));

for iroi = 1:length(roi_list)    
    ROI = roi_list{iroi};    
    for isubj = 1:length(subjID_list)        
        subjID = subjID_list{isubj};        
        % Select the data folders of the subject
        dataPath = fullfile('E:\AMBRA\UoB\Data\MAMSI_MRI', subjID, 'derivatives\ve\MVPA', save_notes, scaling);
        cd(dataPath);
        
        Files_con = dir('res_cong_cong_pred_labels.mat');
        load(Files_con.name);
        res_con = res;        
        AV_con = AV;
        clear res
        % each line is a different subject
        % each column is a different ROI
        rho_tot_roi(isubj,iroi)=res_con(iroi).corr_rho_tot;
    end    
    clear res_con    
end

%% Plot rho
% rho = correlation coefficient between true and decoded congruent locations

rho_tot=atanh(rho_tot_roi);
rho_tot_mean_ftrans=mean(rho_tot);
rho_tot_mean=tanh(mean(rho_tot));

sigmaz=std(rho_tot)/sqrt(size(rho_tot,1));
l95=tanh(mean(rho_tot)-1.96*sigmaz);
u95=tanh(mean(rho_tot)+1.96*sigmaz);

indroi=1:5;
plot_mean=rho_tot_mean(indroi);
plot_l95=l95(indroi);
plot_u95=u95(indroi);

positionXY = [0, 0, 350, 250];
figure('color', [1 1 1], 'Position', positionXY);

plot(plot_mean,'Color','k','LineWidth',1.5);
for iroi=1:length(indroi)
    line([iroi iroi],[plot_l95(iroi) plot_u95(iroi)],...
        'Color','k','LineWidth',1.5);hold on;
    plot(iroi,plot_mean(iroi),'o',...
        'Color','k',...
        'LineWidth',1.5,'MarkerSize',4,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','k');hold on;
end
% adjust x and y ticks
set(gca,'FontName', 'Helvetica');
set(gca,'FontSize', 12);
xl = [0.5 length(indroi)+0.5]; xlim(xl);
set(gca, 'XTickLabel', []);
set(gca,'TickLength', [0.01 0.01]);
set(gca,'LineWidth',1.2)
box off
saveas(gcf,fullfile(savePath, 'rho_group.emf'));

%% Run permutation

% Gets all the possible permutations (via cartesian product)
for iSub=1:size(rho_tot,1)
    sets{iSub} = [-1 1];
end
[a, b, c, d, e, f, g, h, i, j, k, l] = ndgrid(sets{:});
ToPermute = [a(:), b(:), c(:), d(:), e(:), f(:), g(:), h(:), i(:), j(:), k(:), l(:)];

% Compute the null distributions: one for each permutation
for iPerm = 1:size(ToPermute,1)
    tmp2 = (ToPermute(iPerm,:))';
    % Each row of Perms contains the mean of a permutation; each column of Perms is a different ROI
    Perms(iPerm,:) = mean(rho_tot.*repmat(tmp2,1,size(rho_tot,2),1)); %#ok<*SAGROW>
end

for i = 1:size(rho_tot,2)
    if strcmp(Tail,'left')
        % check the proportion of permutation results that are inferior to
        % the mean of my sample
        P(i) = sum(Perms(:,i)<mean(rho_tot(:,i)))/numel(Perms(:,i));
    elseif strcmp(Tail,'right')
        % same but the other way
        P(i) = sum(Perms(:,i)>mean(rho_tot(:,i)))/numel(Perms(:,i));
    elseif strcmp(Tail,'both')
        % for the 2 tailed just compare to the distribution of absolute value of the distance
        % between the result of each permutation to the mean of all
        % permutation results
        
        % Then you check the proportion of those distances are superior to
        % the distance between the mean of your sample and the mean of all
        % permutation results
        % P(i) = sum( abs((Perms(:,i)-mean(Perms(:,i)))) > abs((mean(betas(:,i))-mean(Perms(:,i)))) ) / numel(Perms(:,i)) ;
        
        % Actually not just take the absolute values: the above assumes
        % that your null distribution is symmetric
        P(i) = sum( abs(Perms(:,i)) > abs(mean(rho_tot(:,i)) ) )  / numel(Perms(:,i)) ;        
    end
    
    ES(i)=mean(rho_tot(:,i)-mean(Perms(:,i)));
    CI(i,1)=mean(rho_tot(:,i)-mean(Perms(:,i)))-1.96*(std(rho_tot(:,i)-mean(Perms(:,i)))/sqrt(length(subjID_list)));
    CI(i,2)=mean(rho_tot(:,i)-mean(Perms(:,i)))+1.96*(std(rho_tot(:,i)-mean(Perms(:,i)))/sqrt(length(subjID_list)));    
end

for i=1:length(roi_list)
    roi_list{i,2}=rho_tot_mean(i);
    roi_list{i,3}=[l95(i) u95(i)];
    roi_list{i,4}=rho_tot_mean_ftrans(i);
    roi_list{i,5}=P(i);
    roi_list{i,6}=ES(i);
    roi_list{i,7}=CI(i,[1 2]);
end
