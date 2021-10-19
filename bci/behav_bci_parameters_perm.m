clear;
clc;

%% Exp info

Tail = 'left';
nperm=12;
igroup=1;
group_num=1;

cd('E:\AMBRA\UoB\Data\MAMSI_MRI\group\behav\BCI\single_fit');
load('group_params_bci_attmod_stats');
results=params_matrix_attmod(:,3:6);

%% Effect of attention (V-A)
sigA(:,igroup)=results(:,3)-results(:,1);
sigA_mean(igroup,1)=mean(sigA(:,igroup));

sigV(:,igroup)=results(:,4)-results(:,2);
sigV_mean(igroup,1)=mean(sigV(:,igroup));

% Gets all the possible permutations (via cartesian product)
for iSub=1:nperm
    sets{iSub} = [-1 1];
end
[a, b, c, d, e, f, g, h, i, j, k, l] = ndgrid(sets{:});
ToPermute = [a(:), b(:), c(:), d(:), e(:), f(:), g(:), h(:), i(:), j(:),k(:), l(:)];

sigA=sigA(1:nperm);
sigV=sigV(1:nperm);

% Compute the null distributions: one for each permutation
for iPerm = 1:size(ToPermute,1)
    tmp = (ToPermute(iPerm,:))';
    %% main effects and interactions
    sigAPerms(iPerm,:) = mean(sigA.*repmat(tmp,1,group_num,1)); %#ok<*SAGROW>
    sigVPerms(iPerm,:) = mean(sigV.*repmat(tmp,1,group_num,1)); %#ok<*SAGROW>
end

ES(1)=mean(sigA-mean(sigAPerms));
ES(2)=mean(sigV-mean(sigVPerms));

CI(1,1)=mean(sigA-mean(sigAPerms))-1.96*(std(sigA-mean(sigAPerms))/sqrt(size(results,1)));
CI(1,2)=mean(sigA-mean(sigAPerms))+1.96*(std(sigA-mean(sigAPerms))/sqrt(size(results,1)));
CI(2,1)=mean(sigV-mean(sigVPerms))-1.96*(std(sigV-mean(sigVPerms))/sqrt(size(results,1)));
CI(2,2)=mean(sigV-mean(sigVPerms))+1.96*(std(sigV-mean(sigVPerms))/sqrt(size(results,1)));

for i = 1:group_num
    if strcmp(Tail,'left')
        % check the proportion of permutation results that are inferior to
        % the mean of my sample
        P(i,1) = sum(sigAPerms(:,i)<sigA_mean(i))/numel(sigAPerms(:,i));
        P(i,2) = sum(sigVPerms(:,i)<sigV_mean(i))/numel(sigVPerms(:,i));
    elseif strcmp(Tail,'right')
        % same but the other way
        P(i,1) = sum(sigAPerms(:,i)>sigA_mean(i))/numel(sigAPerms(:,i));
        P(i,2) = sum(sigVPerms(:,i)>sigV_mean(i))/numel(sigVPerms(:,i));
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
        P(i,1) = sum( abs(sigAPerms(:,i)) > abs(sigA_mean(i) ) )  / numel(sigAPerms(:,i)) ;
        P(i,2) = sum( abs(sigVPerms(:,i)) > abs(sigV_mean(i) ) )  / numel(sigVPerms(:,i)) ;
    end
end

res_list{1}='behav';

for i=1:group_num
    res_list{i,2}=sigA_mean(i);
    res_list{i,3}=P(i,1);
    res_list{i,4}=sigV_mean(i);
    res_list{i,5}=P(i,2);
end
