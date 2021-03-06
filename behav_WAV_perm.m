%% Behavioural analysis: WAV sign permutation test

clear;
clc;

control_cb=0;
Tail = 'both';

igroup=1;
group_num=1;

% load and concatenate group matrices generated by behav_WAV.m
load('Aud_WAV', 'WAV_con');
WAV_aud=WAV_con;
load('Vis_WAV', 'WAV_con');
WAV_vis=WAV_con;
results=[WAV_aud,WAV_vis];

%% main effects and interactions
% main effect of attention (V-A)
mainAtt(:,igroup)=mean(results(:,[3 4 7 8]),2)-mean(results(:,[1 2 5 6]),2);
mainAtt_mean(igroup,1)=mean(mainAtt(:,igroup));
% main effect of report (V-A)
mainRep(:,igroup)=mean(results(:,5:8),2)-mean(results(:,1:4),2);
mainRep_mean(igroup,1)=mean(mainRep(:,igroup));
% min effect of AV disparity (H-L)
mainDisp(:,igroup)=mean(results(:,2:2:8),2)-mean(results(:,1:2:8),2);
mainDisp_mean(igroup,1)=mean(mainDisp(:,igroup));
% attention x report
intAttRep(:,igroup)=mean(results(:,[7 8]),2)-mean(results(:,[5 6]),2)-...
    (mean(results(:,[3 4]),2)-mean(results(:,[1 2]),2));
intAttRep_mean(igroup,1)=mean(intAttRep(:,igroup));
% attention x AV disparity
intAttDisp(:,igroup)=mean(results(:,[4 8]),2)-mean(results(:,[3 7]),2)-...
    (mean(results(:,[2 6]),2)-mean(results(:,[1 5]),2));
intAttDisp_mean(igroup,1)=mean(intAttDisp(:,igroup));
% report x AV disparity
intRepDisp(:,igroup)=mean(results(:,[6 8]),2)-mean(results(:,[5 7]),2)-...
    (mean(results(:,[2 4]),2)-mean(results(:,[1 3]),2));
intRepDisp_mean(igroup,1)=mean(intRepDisp(:,igroup));
% attention x report x AV disparity
intAttRepDisp(:,igroup)=(mean(results(:,8),2)-mean(results(:,7),2)-...
    (mean(results(:,4),2)-mean(results(:,3),2)))-...
    ((mean(results(:,6),2)-mean(results(:,5),2))-...
    (mean(results(:,2),2)-mean(results(:,1),2)));
intAttRepDisp_mean(igroup,1)=mean(intAttRepDisp(:,igroup));

%% follow-up simple main effects
intAttADisp(:,igroup)=mean(results(:,[2 6]),2)-mean(results(:,[1 5]),2);
intAttADisp_mean(igroup,1)=mean(intAttADisp(:,igroup));
intAttVDisp(:,igroup)=mean(results(:,[4 8]),2)-mean(results(:,[3 7]),2);
intAttVDisp_mean(igroup,1)=mean(intAttVDisp(:,igroup));

intDispHAtt(:,igroup)=mean(results(:,[4 8]),2)-mean(results(:,[2 6]),2);
intDispHAtt_mean(igroup,1)=mean(intDispHAtt(:,igroup));
intDispLAtt(:,igroup)=mean(results(:,[3 7]),2)-mean(results(:,[1 5]),2);
intDispLAtt_mean(igroup,1)=mean(intDispLAtt(:,igroup));

intRepADisp(:,igroup)=mean(results(:,[2 4]),2)-mean(results(:,[1 3]),2);
intRepADisp_mean(igroup,1)=mean(intRepADisp(:,igroup));
intRepVDisp(:,igroup)=mean(results(:,[6 8]),2)-mean(results(:,[5 7]),2);
intRepVDisp_mean(igroup,1)=mean(intRepVDisp(:,igroup));

intDispHRep(:,igroup)=mean(results(:,[6 8]),2)-mean(results(:,[2 4]),2);
intDispHRep_mean(igroup,1)=mean(intDispHRep(:,igroup));
intDispLRep(:,igroup)=mean(results(:,[5 7]),2)-mean(results(:,[1 3]),2);
intDispLRep_mean(igroup,1)=mean(intDispLRep(:,igroup));

intRepAAtt(:,igroup)=mean(results(:,[3 4]),2)-mean(results(:,[1 2]),2);
intRepAAtt_mean(igroup,1)=mean(intRepAAtt(:,igroup));
intRepVAtt(:,igroup)=mean(results(:,[7 8]),2)-mean(results(:,[5 6]),2);
intRepVAtt_mean(igroup,1)=mean(intRepVAtt(:,igroup));

% Gets all the possible permutations (via cartesian product): might
% not be necessary if you have "a lot" of subjects then you can just randomly permutate X number of times
% but with 10 subejcts that's only 1024 permutations.

for iSub=1:nperm
    sets{iSub} = [-1 1];
end

[a, b, c, d, e, f, g, h, i, j, k, l] = ndgrid(sets{:});
ToPermute = [a(:), b(:), c(:), d(:), e(:), f(:), g(:), h(:), i(:), j(:),k(:), l(:)];    

mainAtt=mainAtt(1:nperm);
mainRep=mainRep(1:nperm);
mainDisp=mainDisp(1:nperm);
intAttRep=intAttRep(1:nperm);
intAttDisp=intAttDisp(1:nperm);
intRepDisp=intRepDisp(1:nperm);
intAttRepDisp=intAttRepDisp(1:nperm);

intAttADisp=intAttADisp(1:nperm);
intAttVDisp=intAttVDisp(1:nperm);
intDispHAtt=intDispHAtt(1:nperm);
intDispLAtt=intDispLAtt(1:nperm);
intRepADisp=intRepADisp(1:nperm);
intRepVDisp=intRepVDisp(1:nperm);
intDispHRep=intDispHRep(1:nperm);
intDispLRep=intDispLRep(1:nperm);
intRepAAtt=intRepAAtt(1:nperm);
intRepVAtt=intRepVAtt(1:nperm);

% Compute the null distributions: one for each permutation
for iPerm = 1:size(ToPermute,1)
    tmp = (ToPermute(iPerm,:))';
    % Each row of Perms contains the mean of a permutation; each column of Perms is a different ROI
    %% main effects and interactions
    PermsAtt(iPerm,:) = mean(mainAtt.*repmat(tmp,1,group_num,1)); %#ok<*SAGROW>
    PermsRep(iPerm,:) = mean(mainRep.*repmat(tmp,1,group_num,1)); %#ok<*SAGROW>
    PermsDisp(iPerm,:) = mean(mainDisp.*repmat(tmp,1,group_num,1)); %#ok<*SAGROW>
    PermsAttRep(iPerm,:) = mean(intAttRep.*repmat(tmp,1,group_num,1)); %#ok<*SAGROW>
    PermsAttDisp(iPerm,:) = mean(intAttDisp.*repmat(tmp,1,group_num,1)); %#ok<*SAGROW>
    PermsRepDisp(iPerm,:) = mean(intRepDisp.*repmat(tmp,1,group_num,1)); %#ok<*SAGROW>
    PermsAttRepDisp(iPerm,:) = mean(intAttRepDisp.*repmat(tmp,1,group_num,1)); %#ok<*SAGROW>
    %% follow-up simple main effects
    PermsAttADisp(iPerm,:) = mean(intAttADisp.*repmat(tmp,1,group_num,1)); %#ok<*SAGROW>
    PermsAttVDisp(iPerm,:) = mean(intAttVDisp.*repmat(tmp,1,group_num,1)); %#ok<*SAGROW>
    PermsDispHAtt(iPerm,:) = mean(intDispHAtt.*repmat(tmp,1,group_num,1)); %#ok<*SAGROW>
    PermsDispLAtt(iPerm,:) = mean(intDispLAtt.*repmat(tmp,1,group_num,1)); %#ok<*SAGROW>
    PermsRepADisp(iPerm,:) = mean(intRepADisp.*repmat(tmp,1,group_num,1)); %#ok<*SAGROW>
    PermsRepVDisp(iPerm,:) = mean(intRepVDisp.*repmat(tmp,1,group_num,1)); %#ok<*SAGROW>
    PermsDispHRep(iPerm,:) = mean(intDispHRep.*repmat(tmp,1,group_num,1)); %#ok<*SAGROW>
    PermsDispLRep(iPerm,:) = mean(intDispLRep.*repmat(tmp,1,group_num,1)); %#ok<*SAGROW>
    PermsRepAAtt(iPerm,:) = mean(intRepAAtt.*repmat(tmp,1,group_num,1)); %#ok<*SAGROW>
    PermsRepVAtt(iPerm,:) = mean(intRepVAtt.*repmat(tmp,1,group_num,1)); %#ok<*SAGROW>
end

ES(1)=mean(mainAtt-mean(PermsAtt));
ES(2)=mean(mainRep-mean(PermsRep));
ES(3)=mean(mainDisp-mean(PermsDisp));
ES(4)=mean(intAttRep-mean(PermsAttRep));
ES(5)=mean(intAttDisp-mean(PermsAttDisp));
ES(6)=mean(intRepDisp-mean(PermsRepDisp));
ES(7)=mean(intAttADisp-mean(PermsAttADisp));
ES(8)=mean(intAttVDisp-mean(PermsAttVDisp));
ES(9)=mean(intDispHAtt-mean(PermsDispHAtt));
ES(10)=mean(intDispLAtt-mean(PermsDispLAtt));
ES(11)=mean(intRepADisp-mean(PermsRepADisp));
ES(12)=mean(intRepVDisp-mean(PermsRepVDisp));
ES(13)=mean(intDispHRep-mean(PermsDispHRep));
ES(14)=mean(intDispLRep-mean(PermsDispLRep));
ES(15)=mean(intRepAAtt-mean(PermsRepAAtt));
ES(16)=mean(intRepVAtt-mean(PermsRepVAtt));
ES(17)=mean(intAttRepDisp-mean(PermsAttRepDisp));

CI(1,1)=mean(mainAtt-mean(PermsAtt))-1.96*(std(mainAtt-mean(PermsAtt))/sqrt(size(results,1)));
CI(1,2)=mean(mainAtt-mean(PermsAtt))+1.96*(std(mainAtt-mean(PermsAtt))/sqrt(size(results,1)));
CI(2,1)=mean(mainRep-mean(PermsRep))-1.96*(std(mainRep-mean(PermsRep))/sqrt(size(results,1)));
CI(2,2)=mean(mainRep-mean(PermsRep))+1.96*(std(mainRep-mean(PermsRep))/sqrt(size(results,1)));
CI(3,1)=mean(mainDisp-mean(PermsDisp))-1.96*(std(mainDisp-mean(PermsDisp))/sqrt(size(results,1)));
CI(3,2)=mean(mainDisp-mean(PermsDisp))+1.96*(std(mainDisp-mean(PermsDisp))/sqrt(size(results,1)));
CI(4,1)=mean(intAttRep-mean(PermsAttRep))-1.96*(std(intAttRep-mean(PermsAttRep))/sqrt(size(results,1)));
CI(4,2)=mean(intAttRep-mean(PermsAttRep))+1.96*(std(intAttRep-mean(PermsAttRep))/sqrt(size(results,1)));
CI(5,1)=mean(intAttDisp-mean(PermsAttDisp))-1.96*(std(intAttDisp-mean(PermsAttDisp))/sqrt(size(results,1)));
CI(5,2)=mean(intAttDisp-mean(PermsAttDisp))+1.96*(std(intAttDisp-mean(PermsAttDisp))/sqrt(size(results,1)));
CI(6,1)=mean(intRepDisp-mean(PermsRepDisp))-1.96*(std(intRepDisp-mean(PermsRepDisp))/sqrt(size(results,1)));
CI(6,2)=mean(intRepDisp-mean(PermsRepDisp))+1.96*(std(intRepDisp-mean(PermsRepDisp))/sqrt(size(results,1)));
CI(7,1)=mean(intAttADisp-mean(PermsAttADisp))-1.96*(std(intAttADisp-mean(PermsAttADisp))/sqrt(size(results,1)));
CI(7,2)=mean(intAttADisp-mean(PermsAttADisp))+1.96*(std(intAttADisp-mean(PermsAttADisp))/sqrt(size(results,1)));
CI(8,1)=mean(intAttVDisp-mean(PermsAttVDisp))-1.96*(std(intAttVDisp-mean(PermsAttVDisp))/sqrt(size(results,1)));
CI(8,2)=mean(intAttVDisp-mean(PermsAttVDisp))+1.96*(std(intAttVDisp-mean(PermsAttVDisp))/sqrt(size(results,1)));
CI(9,1)=mean(intDispHAtt-mean(PermsDispHAtt))-1.96*(std(intDispHAtt-mean(PermsDispHAtt))/sqrt(size(results,1)));
CI(9,2)=mean(intDispHAtt-mean(PermsDispHAtt))+1.96*(std(intDispHAtt-mean(PermsDispHAtt))/sqrt(size(results,1)));
CI(10,1)=mean(intDispLAtt-mean(PermsDispLAtt))-1.96*(std(intDispLAtt-mean(PermsDispLAtt))/sqrt(size(results,1)));
CI(10,2)=mean(intDispLAtt-mean(PermsDispLAtt))+1.96*(std(intDispLAtt-mean(PermsDispLAtt))/sqrt(size(results,1)));
CI(11,1)=mean(intRepADisp-mean(PermsRepADisp))-1.96*(std(intRepADisp-mean(PermsRepADisp))/sqrt(size(results,1)));
CI(11,2)=mean(intRepADisp-mean(PermsRepADisp))+1.96*(std(intRepADisp-mean(PermsRepADisp))/sqrt(size(results,1)));
CI(12,1)=mean(intRepVDisp-mean(PermsRepVDisp))-1.96*(std(intRepVDisp-mean(PermsRepVDisp))/sqrt(size(results,1)));
CI(12,2)=mean(intRepVDisp-mean(PermsRepVDisp))+1.96*(std(intRepVDisp-mean(PermsRepVDisp))/sqrt(size(results,1)));
CI(13,1)=mean(intDispHRep-mean(PermsDispHRep))-1.96*(std(intDispHRep-mean(PermsDispHRep))/sqrt(size(results,1)));
CI(13,2)=mean(intDispHRep-mean(PermsDispHRep))+1.96*(std(intDispHRep-mean(PermsDispHRep))/sqrt(size(results,1)));
CI(14,1)=mean(intDispLRep-mean(PermsDispLRep))-1.96*(std(intDispLRep-mean(PermsDispLRep))/sqrt(size(results,1)));
CI(14,2)=mean(intDispLRep-mean(PermsDispLRep))+1.96*(std(intDispLRep-mean(PermsDispLRep))/sqrt(size(results,1)));
CI(15,1)=mean(intRepAAtt-mean(PermsRepAAtt))-1.96*(std(intRepAAtt-mean(PermsRepAAtt))/sqrt(size(results,1)));
CI(15,2)=mean(intRepAAtt-mean(PermsRepAAtt))+1.96*(std(intRepAAtt-mean(PermsRepAAtt))/sqrt(size(results,1)));
CI(16,1)=mean(intRepVAtt-mean(PermsRepVAtt))-1.96*(std(intRepVAtt-mean(PermsRepVAtt))/sqrt(size(results,1)));
CI(16,2)=mean(intRepVAtt-mean(PermsRepVAtt))+1.96*(std(intRepVAtt-mean(PermsRepVAtt))/sqrt(size(results,1)));
CI(17,1)=mean(intAttRepDisp-mean(PermsAttRepDisp))-1.96*(std(intAttRepDisp-mean(PermsAttRepDisp))/sqrt(size(results,1)));
CI(17,2)=mean(intAttRepDisp-mean(PermsAttRepDisp))+1.96*(std(intAttRepDisp-mean(PermsAttRepDisp))/sqrt(size(results,1)));

for i = 1:group_num
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
    end
end

res_list{1}='behav';
res_list2{1}='behav';

for i=1:group_num
    res_list{i,2}=mainAtt_mean(i);
    res_list{i,3}=P(i,1);
    res_list{i,4}=mainRep_mean(i);
    res_list{i,5}=P(i,2);
    res_list{i,6}=mainDisp_mean(i);
    res_list{i,7}=P(i,3);
    res_list{i,8}=intAttRep_mean(i);
    res_list{i,9}=P(i,4);
    res_list{i,10}=intAttDisp_mean(i);
    res_list{i,11}=P(i,5);
    res_list{i,12}=intRepDisp_mean(i);
    res_list{i,13}=P(i,6);
    res_list{i,14}=intAttRepDisp_mean(i);
    res_list{i,15}=P(i,17);
    
    res_list2{i,2}=intAttADisp_mean(i);
    res_list2{i,3}=P(i,7);
    res_list2{i,4}=intAttVDisp_mean(i);
    res_list2{i,5}=P(i,8);
    res_list2{i,6}=intDispHAtt_mean(i);
    res_list2{i,7}=P(i,9);
    res_list2{i,8}=intDispLAtt_mean(i);
    res_list2{i,9}=P(i,10);
    
    res_list2{i,10}=intRepADisp_mean(i);
    res_list2{i,11}=P(i,11);
    res_list2{i,12}=intRepVDisp_mean(i);
    res_list2{i,13}=P(i,12);
    res_list2{i,14}=intDispHRep_mean(i);
    res_list2{i,15}=P(i,13);
    res_list2{i,16}=intDispLRep_mean(i);
    res_list2{i,17}=P(i,14);
    res_list2{i,18}=intRepAAtt_mean(i);
    res_list2{i,19}=P(i,15);
    res_list2{i,20}=intRepVAtt_mean(i);
    res_list2{i,21}=P(i,16);
end
