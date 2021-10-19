%% MVPA analysis: group sign permutation test on nWAV

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
roi_num = size(roi_list,1);

for iroi = 1:length(roi_list)
    
    load(fullfile('E:\AMBRA\UoB\Data\MAMSI_MRI\group\MVPA\',scaling,...
    [roi_list{iroi} '_nWAV.mat']));

    %% main effects and interactions
    % main effect of attention (V-A)
    mainAtt(:,iroi)=mean(results(:,[3 4 7 8]),2)-mean(results(:,[1 2 5 6]),2);
    mainAtt_mean(iroi,1)=mean(mainAtt(:,iroi));
    % main effect of report (V-A)
    mainRep(:,iroi)=mean(results(:,5:8),2)-mean(results(:,1:4),2);
    mainRep_mean(iroi,1)=mean(mainRep(:,iroi));
    % main effect of AV disparity (H-L)
    mainDisp(:,iroi)=mean(results(:,2:2:8),2)-mean(results(:,1:2:8),2);
    mainDisp_mean(iroi,1)=mean(mainDisp(:,iroi));
    % attention x report
    intAttRep(:,iroi)=mean(results(:,[7 8]),2)-mean(results(:,[5 6]),2)-...
        (mean(results(:,[3 4]),2)-mean(results(:,[1 2]),2));
    intAttRep_mean(iroi,1)=mean(intAttRep(:,iroi));
    % attention x AV disparity
    intAttDisp(:,iroi)=mean(results(:,[4 8]),2)-mean(results(:,[3 7]),2)-...
        (mean(results(:,[2 6]),2)-mean(results(:,[1 5]),2));
    intAttDisp_mean(iroi,1)=mean(intAttDisp(:,iroi));
    % report x AV disparity
    intRepDisp(:,iroi)=mean(results(:,[6 8]),2)-mean(results(:,[5 7]),2)-...
        (mean(results(:,[2 4]),2)-mean(results(:,[1 3]),2));
    intRepDisp_mean(iroi,1)=mean(intRepDisp(:,iroi));
    % attention x report x AV disparity
    intAttRepDisp(:,iroi)=(mean(results(:,8),2)-mean(results(:,7),2)-...
        (mean(results(:,4),2)-mean(results(:,3),2)))-...
        ((mean(results(:,6),2)-mean(results(:,5),2))-...
        (mean(results(:,2),2)-mean(results(:,1),2)));
    intAttRepDisp_mean(iroi,1)=mean(intAttRepDisp(:,iroi));
    
    % pooled over all conditions
    All(:,iroi)=mean(results,2);
    All_mean(iroi,1)=mean(All(:,iroi));
    
    %% follow-up simple main effects
    intAttADisp(:,iroi)=mean(results(:,[2 6]),2)-mean(results(:,[1 5]),2);
    intAttADisp_mean(iroi,1)=mean(intAttADisp(:,iroi));
    intAttVDisp(:,iroi)=mean(results(:,[4 8]),2)-mean(results(:,[3 7]),2);
    intAttVDisp_mean(iroi,1)=mean(intAttVDisp(:,iroi));
    
    intDispHAtt(:,iroi)=mean(results(:,[4 8]),2)-mean(results(:,[2 6]),2);
    intDispHAtt_mean(iroi,1)=mean(intDispHAtt(:,iroi));
    intDispLAtt(:,iroi)=mean(results(:,[3 7]),2)-mean(results(:,[1 5]),2);
    intDispLAtt_mean(iroi,1)=mean(intDispLAtt(:,iroi));
    
    intRepADisp(:,iroi)=mean(results(:,[2 4]),2)-mean(results(:,[1 3]),2);
    intRepADisp_mean(iroi,1)=mean(intRepADisp(:,iroi));
    intRepVDisp(:,iroi)=mean(results(:,[6 8]),2)-mean(results(:,[5 7]),2);
    intRepVDisp_mean(iroi,1)=mean(intRepVDisp(:,iroi));
    
    intDispHRep(:,iroi)=mean(results(:,[6 8]),2)-mean(results(:,[2 4]),2);
    intDispHRep_mean(iroi,1)=mean(intDispHRep(:,iroi));
    intDispLRep(:,iroi)=mean(results(:,[5 7]),2)-mean(results(:,[1 3]),2);
    intDispLRep_mean(iroi,1)=mean(intDispLRep(:,iroi));
    
    intRepAAtt(:,iroi)=mean(results(:,[3 4]),2)-mean(results(:,[1 2]),2);
    intRepAAtt_mean(iroi,1)=mean(intRepAAtt(:,iroi));
    intRepVAtt(:,iroi)=mean(results(:,[7 8]),2)-mean(results(:,[5 6]),2);
    intRepVAtt_mean(iroi,1)=mean(intRepVAtt(:,iroi));
    
end

% Gets all the possible permutations (via cartesian product)
for iSub=1:size(results,1)
    sets{iSub} = [-1 1];
end
[a, b, c, d, e, f, g, h, i, j, k, l] = ndgrid(sets{:});
ToPermute = [a(:), b(:), c(:), d(:), e(:), f(:), g(:), h(:), i(:), j(:), k(:), l(:)];

% Compute the null distributions: one for each permutation
for iPerm = 1:size(ToPermute,1)
    tmp = (ToPermute(iPerm,:))';
    % Each row of Perms contains the mean of a permutation; each column of Perms is a different ROI
    %% main effects and interactions
    PermsAtt(iPerm,:) = mean(mainAtt.*repmat(tmp,1,roi_num,1)); %#ok<*SAGROW>
    PermsRep(iPerm,:) = mean(mainRep.*repmat(tmp,1,roi_num,1)); %#ok<*SAGROW>
    PermsDisp(iPerm,:) = mean(mainDisp.*repmat(tmp,1,roi_num,1)); %#ok<*SAGROW>
    PermsAttRep(iPerm,:) = mean(intAttRep.*repmat(tmp,1,roi_num,1)); %#ok<*SAGROW>
    PermsAttDisp(iPerm,:) = mean(intAttDisp.*repmat(tmp,1,roi_num,1)); %#ok<*SAGROW>
    PermsRepDisp(iPerm,:) = mean(intRepDisp.*repmat(tmp,1,roi_num,1)); %#ok<*SAGROW>
    PermsAttRepDisp(iPerm,:) = mean(intAttRepDisp.*repmat(tmp,1,roi_num,1)); %#ok<*SAGROW>
    PermsAll(iPerm,:)= mean(All.*repmat(tmp,1,roi_num,1)); %#ok<*SAGROW>
    PermsAll2(iPerm,:)= 1+mean(All.*repmat(tmp,1,roi_num,1)); %#ok<*SAGROW>
    %% follow-up simple main effects
    PermsAttADisp(iPerm,:) = mean(intAttADisp.*repmat(tmp,1,roi_num,1)); %#ok<*SAGROW>
    PermsAttVDisp(iPerm,:) = mean(intAttVDisp.*repmat(tmp,1,roi_num,1)); %#ok<*SAGROW>
    PermsDispHAtt(iPerm,:) = mean(intDispHAtt.*repmat(tmp,1,roi_num,1)); %#ok<*SAGROW>
    PermsDispLAtt(iPerm,:) = mean(intDispLAtt.*repmat(tmp,1,roi_num,1)); %#ok<*SAGROW>
    PermsRepADisp(iPerm,:) = mean(intRepADisp.*repmat(tmp,1,roi_num,1)); %#ok<*SAGROW>
    PermsRepVDisp(iPerm,:) = mean(intRepVDisp.*repmat(tmp,1,roi_num,1)); %#ok<*SAGROW>
    PermsDispHRep(iPerm,:) = mean(intDispHRep.*repmat(tmp,1,roi_num,1)); %#ok<*SAGROW>
    PermsDispLRep(iPerm,:) = mean(intDispLRep.*repmat(tmp,1,roi_num,1)); %#ok<*SAGROW>
    PermsRepAAtt(iPerm,:) = mean(intRepAAtt.*repmat(tmp,1,roi_num,1)); %#ok<*SAGROW>
    PermsRepVAtt(iPerm,:) = mean(intRepVAtt.*repmat(tmp,1,roi_num,1)); %#ok<*SAGROW>
end

for i = 1:roi_num
    if strcmp(Tail,'left')
        % check the proportion of permutation results that are inferior to
        % the mean of my sample
        P(i,1) = sum(PermsAtt(:,i)<mainAtt_mean(i))/numel(PermsAtt(:,i));
        P(i,2) = sum(PermsRep(:,i)<mainRep_mean(i))/numel(PermsRep(:,i));
        P(i,3) = sum(PermsDisp(:,i)<mainDisp_mean(i))/numel(PermsDisp(:,i));
        P(i,4) = sum(PermsAttRep(:,i)<intAttRep_mean(i))/numel(PermsAttRep(:,i));
        P(i,5) = sum(PermsAttDisp(:,i)<intAttDisp_mean(i))/numel(PermsAttDisp(:,i));
        P(i,6) = sum(PermsRepDisp(:,i)<intRepDisp_mean(i))/numel(PermsRepDisp(:,i));        
        P(i,7) = sum(PermsAttADisp(:,i)<intAttADisp_mean(i))/numel(PermsAttADisp(:,i));
        P(i,8) = sum(PermsAttVDisp(:,i)<intAttVDisp_mean(i))/numel(PermsAttVDisp(:,i));
        P(i,9) = sum(PermsDispHAtt(:,i)<intDispHAtt_mean(i))/numel(PermsDispHAtt(:,i));
        P(i,10) = sum(PermsDispLAtt(:,i)<intDispLAtt_mean(i))/numel(PermsDispLAtt(:,i));
        P(i,11) = sum(PermsRepADisp(:,i)<intRepADisp_mean(i))/numel(PermsRepADisp(:,i));
        P(i,12) = sum(PermsRepVDisp(:,i)<intRepVDisp_mean(i))/numel(PermsRepVDisp(:,i));
        P(i,13) = sum(PermsDispHRep(:,i)<intDispHRep_mean(i))/numel(PermsDispHRep(:,i));
        P(i,14) = sum(PermsDispLRep(:,i)<intDispLRep_mean(i))/numel(PermsDispLRep(:,i));
        P(i,15) = sum(PermsRepAAtt(:,i)<intRepAAtt_mean(i))/numel(PermsRepAAtt(:,i));
        P(i,16) = sum(PermsRepVAtt(:,i)<intRepVAtt_mean(i))/numel(PermsRepVAtt(:,i));
        P(i,17) = sum(PermsAttRepDisp(:,i)<intAttRepDisp_mean(i))/numel(PermsAttRepDisp(:,i));
        P(i,18) = sum(PermsAll2(:,i)<All_mean(i))/numel(PermsAll2(:,i));
    elseif strcmp(Tail,'right')
        % same but the other way
        P(i,1) = sum(PermsAtt(:,i)>mainAtt_mean(i))/numel(PermsAtt(:,i));
        P(i,2) = sum(PermsRep(:,i)>mainRep_mean(i))/numel(PermsRep(:,i));
        P(i,3) = sum(PermsDisp(:,i)>mainDisp_mean(i))/numel(PermsDisp(:,i));
        P(i,4) = sum(PermsAttRep(:,i)>intAttRep_mean(i))/numel(PermsAttRep(:,i));
        P(i,5) = sum(PermsAttDisp(:,i)>intAttDisp_mean(i))/numel(PermsAttDisp(:,i));
        P(i,6) = sum(PermsRepDisp(:,i)>intRepDisp_mean(i))/numel(PermsRepDisp(:,i));
        P(i,7) = sum(PermsAttADisp(:,i)>intAttADisp_mean(i))/numel(PermsAttADisp(:,i));
        P(i,8) = sum(PermsAttVDisp(:,i)>intAttVDisp_mean(i))/numel(PermsAttVDisp(:,i));
        P(i,9) = sum(PermsDispHAtt(:,i)>intDispHAtt_mean(i))/numel(PermsDispHAtt(:,i));
        P(i,10) = sum(PermsDispLAtt(:,i)>intDispLAtt_mean(i))/numel(PermsDispLAtt(:,i));
        P(i,11) = sum(PermsRepADisp(:,i)>intRepADisp_mean(i))/numel(PermsRepADisp(:,i));
        P(i,12) = sum(PermsRepVDisp(:,i)>intRepVDisp_mean(i))/numel(PermsRepVDisp(:,i));
        P(i,13) = sum(PermsDispHRep(:,i)>intDispHRep_mean(i))/numel(PermsDispHRep(:,i));
        P(i,14) = sum(PermsDispLRep(:,i)>intDispLRep_mean(i))/numel(PermsDispLRep(:,i));
        P(i,15) = sum(PermsRepAAtt(:,i)>intRepAAtt_mean(i))/numel(PermsRepAAtt(:,i));
        P(i,16) = sum(PermsRepVAtt(:,i)>intRepVAtt_mean(i))/numel(PermsRepVAtt(:,i));
        P(i,17) = sum(PermsAttRepDisp(:,i)>intAttRepDisp_mean(i))/numel(PermsAttRepDisp(:,i));
        P(i,18) = sum(PermsAll(:,i)>All_mean(i))/numel(PermsAll(:,i));
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
        P(i,1) = sum( abs(PermsAtt(:,i)) > abs(mainAtt_mean(i) ) )  / numel(PermsAtt(:,i)) ;
        P(i,2) = sum( abs(PermsRep(:,i)) > abs(mainRep_mean(i) ) )  / numel(PermsRep(:,i)) ;
        P(i,3) = sum( abs(PermsDisp(:,i)) > abs(mainDisp_mean(i) ) )  / numel(PermsDisp(:,i)) ;
        P(i,4) = sum( abs(PermsAttRep(:,i)) > abs(intAttRep_mean(i) ) )  / numel(PermsAttRep(:,i)) ;
        P(i,5) = sum( abs(PermsAttDisp(:,i)) > abs(intAttDisp_mean(i) ) )  / numel(PermsAttDisp(:,i)) ;
        P(i,6) = sum( abs(PermsRepDisp(:,i)) > abs(intRepDisp_mean(i) ) )  / numel(PermsRepDisp(:,i)) ;        
        P(i,7) = sum( abs(PermsAttADisp(:,i)) > abs(intAttADisp_mean(i) ) )  / numel(PermsAttADisp(:,i)) ;
        P(i,8) = sum( abs(PermsAttVDisp(:,i)) > abs(intAttVDisp_mean(i) ) )  / numel(PermsAttVDisp(:,i)) ;
        P(i,9) = sum( abs(PermsDispHAtt(:,i)) > abs(intDispHAtt_mean(i) ) )  / numel(PermsDispHAtt(:,i)) ;
        P(i,10) = sum( abs(PermsDispLAtt(:,i)) > abs(intDispLAtt_mean(i) ) )  / numel(PermsDispLAtt(:,i)) ;
        P(i,11) = sum( abs(PermsRepADisp(:,i)) > abs(intRepADisp_mean(i) ) )  / numel(PermsRepADisp(:,i)) ;
        P(i,12) = sum( abs(PermsRepVDisp(:,i)) > abs(intRepVDisp_mean(i) ) )  / numel(PermsRepVDisp(:,i)) ;
        P(i,13) = sum( abs(PermsDispHRep(:,i)) > abs(intDispHRep_mean(i) ) )  / numel(PermsDispHRep(:,i)) ;
        P(i,14) = sum( abs(PermsDispLRep(:,i)) > abs(intDispLRep_mean(i) ) )  / numel(PermsDispLRep(:,i)) ;
        P(i,15) = sum( abs(PermsRepAAtt(:,i)) > abs(intRepAAtt_mean(i) ) )  / numel(PermsRepAAtt(:,i)) ;
        P(i,16) = sum( abs(PermsRepVAtt(:,i)) > abs(intRepVAtt_mean(i) ) )  / numel(PermsRepVAtt(:,i)) ;
        P(i,17) = sum( abs(PermsAttRepDisp(:,i)) > abs(intAttRepDisp_mean(i) ) )  / numel(PermsAttRepDisp(:,i)) ;
        P(i,18) = sum( abs(PermsAll(:,i)) > abs(All_mean(i)) )  / numel(PermsAll(:,i)) ;
        P(i,19) = sum( abs(PermsAll2(:,i)) < abs(All_mean(i)) )  / numel(PermsAll2(:,i)) ;
    end
    ES(i,1)=mean(mainAtt(:,i)-mean(PermsAtt(:,i)));
    ES(i,2)=mean(mainRep(:,i)-mean(PermsRep(:,i)));
    ES(i,3)=mean(mainDisp(:,i)-mean(PermsDisp(:,i)));
    ES(i,4)=mean(intAttRep(:,i)-mean(PermsAttRep(:,i)));
    ES(i,5)=mean(intAttDisp(:,i)-mean(PermsAttDisp(:,i)));
    ES(i,6)=mean(intRepDisp(:,i)-mean(PermsRepDisp(:,i)));
    ES(i,7)=mean(intAttRepDisp(:,i)-mean(PermsAttRepDisp(:,i)));
    ES(i,8)=mean(intRepAAtt(:,i)-mean(PermsRepAAtt(:,i)));
    if strcmp(Tail,'left')
        ES(i,9)=mean(All(:,i)-mean(PermsAll2(:,i)));
    else
        ES(i,9)=mean(All(:,i)-mean(PermsAll(:,i)));        
    end
    
    CI(i,1,1)=mean(mainAtt(:,i)-mean(PermsAtt(:,i)))-1.96*(std(mainAtt(:,i)-mean(PermsAtt(:,i)))/sqrt(size(results,1)));
    CI(i,1,2)=mean(mainAtt(:,i)-mean(PermsAtt(:,i)))+1.96*(std(mainAtt(:,i)-mean(PermsAtt(:,i)))/sqrt(size(results,1)));
    CI(i,2,1)=mean(mainRep(:,i)-mean(PermsRep(:,i)))-1.96*(std(mainRep(:,i)-mean(PermsRep(:,i)))/sqrt(size(results,1)));
    CI(i,2,2)=mean(mainRep(:,i)-mean(PermsRep(:,i)))+1.96*(std(mainRep(:,i)-mean(PermsRep(:,i)))/sqrt(size(results,1)));
    CI(i,3,1)=mean(mainDisp(:,i)-mean(PermsDisp(:,i)))-1.96*(std(mainDisp(:,i)-mean(PermsDisp(:,i)))/sqrt(size(results,1)));
    CI(i,3,2)=mean(mainDisp(:,i)-mean(PermsDisp(:,i)))+1.96*(std(mainDisp(:,i)-mean(PermsDisp(:,i)))/sqrt(size(results,1)));
    CI(i,4,1)=mean(intAttRep(:,i)-mean(PermsAttRep(:,i)))-1.96*(std(intAttRep(:,i)-mean(PermsAttRep(:,i)))/sqrt(size(results,1)));
    CI(i,4,2)=mean(intAttRep(:,i)-mean(PermsAttRep(:,i)))+1.96*(std(intAttRep(:,i)-mean(PermsAttRep(:,i)))/sqrt(size(results,1)));
    CI(i,5,1)=mean(intAttDisp(:,i)-mean(PermsAttDisp(:,i)))-1.96*(std(intAttDisp(:,i)-mean(PermsAttDisp(:,i)))/sqrt(size(results,1)));
    CI(i,5,2)=mean(intAttDisp(:,i)-mean(PermsAttDisp(:,i)))+1.96*(std(intAttDisp(:,i)-mean(PermsAttDisp(:,i)))/sqrt(size(results,1)));
    CI(i,6,1)=mean(intRepDisp(:,i)-mean(PermsRepDisp(:,i)))-1.96*(std(intRepDisp(:,i)-mean(PermsRepDisp(:,i)))/sqrt(size(results,1)));
    CI(i,6,2)=mean(intRepDisp(:,i)-mean(PermsRepDisp(:,i)))+1.96*(std(intRepDisp(:,i)-mean(PermsRepDisp(:,i)))/sqrt(size(results,1)));
    CI(i,7,1)=mean(intAttRepDisp(:,i)-mean(PermsAttRepDisp(:,i)))-1.96*(std(intAttRepDisp(:,i)-mean(PermsAttRepDisp(:,i)))/sqrt(size(results,1)));
    CI(i,7,2)=mean(intAttRepDisp(:,i)-mean(PermsAttRepDisp(:,i)))+1.96*(std(intAttRepDisp(:,i)-mean(PermsAttRepDisp(:,i)))/sqrt(size(results,1)));
    CI(i,8,1)=mean(intRepAAtt(:,i)-mean(PermsRepAAtt(:,i)))-1.96*(std(intRepAAtt(:,i)-mean(PermsRepAAtt(:,i)))/sqrt(size(results,1)));
    CI(i,8,2)=mean(intRepAAtt(:,i)-mean(PermsRepAAtt(:,i)))+1.96*(std(intRepAAtt(:,i)-mean(PermsRepAAtt(:,i)))/sqrt(size(results,1)));
    if strcmp(Tail,'left')
        CI(i,9,1)=mean(All(:,i)-mean(PermsAll2(:,i)))-1.96*(std(All(:,i)-mean(PermsAll2(:,i)))/sqrt(size(results,1)));
        CI(i,9,2)=mean(All(:,i)-mean(PermsAll2(:,i)))+1.96*(std(All(:,i)-mean(PermsAll2(:,i)))/sqrt(size(results,1)));
    else
        CI(i,9,1)=mean(All(:,i)-mean(PermsAll(:,i)))-1.96*(std(All(:,i)-mean(PermsAll(:,i)))/sqrt(size(results,1)));
        CI(i,9,2)=mean(All(:,i)-mean(PermsAll(:,i)))+1.96*(std(All(:,i)-mean(PermsAll(:,i)))/sqrt(size(results,1)));      
    end
end

for i=1:roi_num
    roi_list{i,2}=mainAtt_mean(i);
    roi_list{i,3}=P(i,1);
    roi_list{i,4}=mainRep_mean(i);
    roi_list{i,5}=P(i,2);
    roi_list{i,6}=mainDisp_mean(i);
    roi_list{i,7}=P(i,3);
    roi_list{i,8}=intAttRep_mean(i);
    roi_list{i,9}=P(i,4);
    roi_list{i,10}=intAttDisp_mean(i);
    roi_list{i,11}=P(i,5);
    roi_list{i,12}=intRepDisp_mean(i);
    roi_list{i,13}=P(i,6);
    roi_list{i,14}=intAttRepDisp_mean(i);
    roi_list{i,15}=P(i,17);
    roi_list{i,16}=All_mean(i);
    roi_list{i,17}=P(i,18);
    if strcmp(Tail,'both')
        roi_list{i,18}=P(i,19);
    end
    
    roi_list2{i,2}=intAttADisp_mean(i);
    roi_list2{i,3}=P(i,7);
    roi_list2{i,4}=intAttVDisp_mean(i);
    roi_list2{i,5}=P(i,8);
    roi_list2{i,6}=intDispHAtt_mean(i);
    roi_list2{i,7}=P(i,9);
    roi_list2{i,8}=intDispLAtt_mean(i);
    roi_list2{i,9}=P(i,10);
    
    roi_list2{i,10}=intRepADisp_mean(i);
    roi_list2{i,11}=P(i,11);
    roi_list2{i,12}=intRepVDisp_mean(i);
    roi_list2{i,13}=P(i,12);
    roi_list2{i,14}=intDispHRep_mean(i);
    roi_list2{i,15}=P(i,13);
    roi_list2{i,16}=intDispLRep_mean(i);
    roi_list2{i,17}=P(i,14);
    roi_list2{i,18}=intRepAAtt_mean(i);
    roi_list2{i,19}=P(i,15);
    roi_list2{i,20}=intRepVAtt_mean(i);
    roi_list2{i,21}=P(i,16);
end
