%% Preparations
clc;
clear;

% Checking setup
[~,setupID] = system('hostname');
setupID = regexp(setupID,'[a-zA-Z_0-9-]*','match');
setupID = setupID{1};

if ~isempty(regexp(setupID,'^bb.+','once'))
    isServer = true;
    progrMonitor = false;
else
    isServer = false;
    progrMonitor = true;
    dbstop if error;
end

% Opening parallel pool for parfor loop during gridsearch.
% if there is no parallel pool running, open one.
currPool = gcp('nocreate');
if isempty(currPool)
    if isServer
        parpool('local',16);
    else
        parpool('local');
    end
end

% Path to BCI scripts
addpath(genpath('E:\AMBRA\UoB\Exp\misc'));

% number of minima to evaluate with fminsearch
nmin = 10;

% model
model='fus';

% initialise cell for best 10 results
group_best10_results_repmod = cell(length(subjID),1);

% initialise structure for best results
group_best_results_repmod = cell(length(subjID),1);

figure_plot=0;

for iSubj = 1:length(subjID)
    
    if strcmp(expID, 'MRI')
        cd(fullfile('E:\AMBRA\UoB\Data\MAMSI_MRI', subjID{iSubj}, 'behav\scanner'));
    elseif strcmp(expID, 'behav')
        cd(fullfile('E:\AMBRA\UoB\Data\MAMSI_MRI', 'MAMSI_MRI_behav', subjID{iSubj}));
    end
    
    origData = dir('*Exp_All_Sessions*.mat');
    load(origData.name, 'tdata');
    origData = tdata;
    
    % initialise structure for subject results
    subj_results = struct;
    
    % Reorganize dataset 'tdata' so that I get the information I need for model estimation (AttentionModality, locA, locV, respA, respV)
    % Get rid of incorrect responses (wrong keypad), missed responses (no answer provided) and anticipated responses (RT < 100ms)
    origData = origData(origData.IncorResp==0 & origData.MissedResp==0 & origData.AntResp==0,:);
    
    % Put a variable for Attention Modality
    if ~sum(strcmp(origData.Properties.VarNames(:),'AttentionModalityBCI')) % if AttentionModality is not present
        origData.AttentionModalityBCI = NaN(size(origData,1),1);
        for j = 1:size(origData,1)
            % Report Auditory = 1
            if strcmp(origData.ResponseModality(j), 'Aud')
                origData.AttentionModalityBCI(j) = 1;
                % Report Visual = 2
            elseif strcmp(origData.ResponseModality(j), 'Vis')
                origData.AttentionModalityBCI(j) = 2;
            end
        end
    end
    
    % Put locA and locV
    if ~sum(strcmp(origData.Properties.VarNames(:),'locA')) % if locA is not present (and locV as well)
        origData.locA = NaN(size(origData,1),1);
        origData.locV = NaN(size(origData,1),1);
        for j = 1:size(origData,1)
            % locA: when people have to locate A, locA = location of the target and locV = location of the distractor
            if strcmp(origData.ResponseModality(j), 'Aud')
                origData.locA(j) = origData.TargetLoc(j);
                origData.locV(j) = origData.NonTargetLoc(j);
                % locV: when people have to locate V, locV = location of the target and locA = location of the distractor
            elseif strcmp(origData.ResponseModality(j), 'Vis')
                origData.locV(j) = origData.TargetLoc(j);
                origData.locA(j) = origData.NonTargetLoc(j);
            end
        end
    end
    
    % Put respA and respV
    if ~sum(strcmp(origData.Properties.VarNames(:),'respA')) % if locA is not present (and respV as well)
        origData.respA = NaN(size(origData,1),1);
        origData.respV = NaN(size(origData,1),1);
        for j = 1:size(origData,1)
            % respA: when people have to locate A, respA = response
            if strcmp(origData.ResponseModality(j), 'Aud')
                origData.respA(j) = origData.Response(j);
                origData.respV(j) = NaN;
                % respV: when people have to locate V, respV = response
            elseif strcmp(origData.ResponseModality(j), 'Vis')
                origData.respV(j) = origData.Response(j);
                origData.respA(j) = NaN;
            end
        end
    end
    
    % Keep only useful data (AttentionModality, locaV, locA, respV, respA)
    bciData = dataset2table(origData(:,{'AttentionModalityBCI','locV','locA','respV','respA'}));
    bciData.Properties.VariableNames = {'AttentionModality' 'locV' 'locA' 'respV' 'respA'};
    
    %% Perform gridsearch on a grid of parameters
    % Response locations
    responseLoc = unique(bciData.locV);
    % Readout type: model averaging (1) vs. model selection (2)
    readout = 1;
    
    % Preparing the parameter space
    % grid resolution
    n = 5;
    if strcmp(model, 'fus')
        p_common=1;
    elseif strcmp(model, 'segA') || strcmp(model, 'segV') || strcmp(model, 'taskRel')
        p_common=0;
    else
        p_common=linspace(0.1,0.9,n);
    end
    sigP = linspace(0.1,30,n);
    sigV1 = linspace(0.1,30,n);
    sigV2 = linspace(0.1,30,n);
    sigA1 = linspace(0.1,30,n);
    sigA2 = linspace(0.1,30,n);
    
    % These parameters can be used as well, but not included in this analysis
    % kernelWidth = 1; 'kW'
    % muP = 0; % mean of spatial prior
    
    % Parameter names
    if strcmp(model, 'fus') || strcmp(model, 'taskRel')
        parameterNames = {'sigP','sigA1','sigA2','sigV1','sigV2'};
        gridVectors = {sigP,sigA1,sigA2,sigV1,sigV2};
    elseif strcmp(model, 'segA')
        parameterNames = {'sigP','sigA1','sigA2'};
        gridVectors = {sigP,sigA1,sigA2};
    elseif strcmp(model, 'segV')
        parameterNames = {'sigP','sigV1','sigV2'};
        gridVectors = {sigP,sigV1,sigV2};
    else
        parameterNames = {'p_common','sigP','sigA1','sigA2','sigV1','sigV2'};
        gridVectors = {p_common,sigP,sigA1,sigA2,sigV1,sigV2};
    end
    nParameters = numel(gridVectors);
    % Full factorial expansion of the specified parameter vectors
    coords = cell(1,nParameters);
    [coords{:}] = ndgrid(gridVectors{:});
    coords = cellfun(@(x) x(:),coords,'UniformOutput',false);
    % Matrix of all possible parameter combinations: nRows = number of
    % combinations, columns correspond to parameters
    paramCombinations = cat(2,coords{:});
    
    % Pre-allocating array for data collection
    logLike_all = NaN(size(paramCombinations,1),1);
    
    % Starting the timer and printing details
    cStart = clock;
    fprintf('Performing gridsearch... \n');
    
    % Performing grid search on the specified parameter space
    parfor in = 1:size(paramCombinations,1)
        [logLike_all(in)] = fitModelAcrossAtt(paramCombinations(in,:),parameterNames,bciData,responseLoc,model,readout,'logLike');
    end
    
    % Printing elapsed time
    cTemp = clock;
    fprintf('Gridsearch elapsed time (days hours:minutes:seconds) %s \n\n',...
        datestr(etime(cTemp,cStart)/86400,'dd HH:MM:SS'));
    
    % Saving subject specific bci simulations
    fprintf('\n\nSaving data...\n\n');
    save([subjID{iSubj} '_' model 'Simulations_repmod_gridsearch'], ...
        'subjID', 'logLike_all', 'paramCombinations', 'parameterNames', 'bciData');
    
    % Finding the n minima of the log-likelihood values and the corresponding parameter combination
    % sort gridsearch log-likelihoods from the 1st best (global minimum)
    logLike_sort = sort(logLike_all);
    
    % get id number and parameter combination for the n best log-likelihoods
    for i = 1:nmin
        idx = find(logLike_all==logLike_sort(i));
        if length(idx)>1
            idx_temp=randi(length(idx)); %randomly select one
            idx=idx(idx_temp);
        end
        bestParam = paramCombinations(idx,:);
        % save in subj-specific structure
        subj_results(i).idx = idx;
        subj_results(i).bestParam = bestParam;
    end
    % get n best log-likelihoods and save in subj-specific structure
    bestlogLike(:,1) = logLike_sort(1:nmin);
    for i=1:nmin
        subj_results(i).bestlogLike = bestlogLike(i);
    end
    
    % Plot log-likelihoods (with n best marked with a *)
    cols = {[103 0 13];[133 6 17];[165 15 21];[203 24 29];[213 34 35];...
        [239 59 44];[251 106 74];[252 146 114];[252 187 161];[254 224 210]};
    if figure_plot==1
        figure('Position', [0, 0, 1920, 1080]); plot(logLike_all,'k'); hold on;
        for i = 1:nmin
            plot(subj_results(i).idx,logLike_all(subj_results(i).idx), '*', 'Color', cols{i}/255); hold on;
        end
        saveas(gcf,[subjID{iSubj} '_gridsearch_repmod_logLike' ],'emf');
    end
    %% Perform fminsearch to refine the best parameters
    % Options for fminsearchbnd
    opts = optimset('fminsearch');
    opts.Display = 'iter';
    opts.MaxFunEvals = 4000;
    opts.TolFun = 1.e-12;
    opts.MaxIter = 1000;
    
    % set upper and lower bounds on parameters
    % 'p_common','sigP','sigA1', 'sigA2', 'sigV1', 'sigV2'
    if strcmp(model, 'fus') || strcmp(model, 'taskRel')
        LB = [0.001 0.001 0.001 0.001 0.001];
        UB = [30    30    30    30    30   ];
    elseif strcmp(model, 'segA') || strcmp(model, 'segV')
        LB = [0.001 0.001 0.001];
        UB = [30    30    30   ];
    else
        LB = [0 0.001 0.001 0.001 0.001 0.001];
        UB = [1 30    30    30    30    30   ];
    end
    
    % Creating anonymous function for input to fminsearch
    fun = @(param) fitModelAcrossAtt(param,parameterNames,bciData,responseLoc,model,readout,'logLike');
    
    % Starting the timer and printing details
    cStart = clock;
    
    % Performing fminsearch for each parameter combination
    for i = 1:nmin
        fprintf(['Performing fminsearch ' num2str(i) ' of subject ' num2str(iSubj) '... \n']);
        [fmin_parameters,errorval] = fminsearchbnd(fun,subj_results(i).bestParam,LB,UB,opts);
        subj_results(i).fmin_parameters = fmin_parameters;
        subj_results(i).errorval = errorval;
    end
    
    % Finishing timer and printing elapsed time
    fprintf('fminsearch elapsed time (days hours:minutes:seconds) %s \n\n',...
        datestr(etime(clock,cStart)/86400,'dd HH:MM:SS'));
    
    %% Evaluating model fit of best model
    
    % Save all simulation data and plot
    my_actions = 1;
    
    for i = 1:nmin
        [fmin_logLike,bic, aic, bci_details] = fitModelAcrossAtt(subj_results(i).fmin_parameters,parameterNames,bciData,responseLoc,model,readout,'full');
        
        subj_results(i).fmin_logLike = fmin_logLike;
        subj_results(i).bci_details = bci_details;
        subj_results(i).bic = bic;
        subj_results(i).aic = aic;
    end
    % save structure with current subj results into group cell
    group_best10_results_repmod{iSubj} = subj_results;
    
    % Find the best parameters combination after fminsearch
    fmin_logLike = [];
    for i = 1:nmin
        fmin_logLike  = cat(1,fmin_logLike,subj_results(i).fmin_logLike);
    end
    bestSim = find(fmin_logLike==min(fmin_logLike));
    fmin_parameters = subj_results(bestSim).fmin_parameters;
    
    % Produce bci simulations for the best parameters combination
    [fmin_bestlogLike, bic, aic, bci_details] = fitModelAcrossAtt(fmin_parameters,parameterNames,bciData,responseLoc,model,readout,'full');
    
    % Save data in structure
    group_best_results_repmod{iSubj}.fmin_bestlogLike = fmin_bestlogLike;
    group_best_results_repmod{iSubj}.bci_details = bci_details;
    group_best_results_repmod{iSubj}.bic = bic;
    group_best_results_repmod{iSubj}.aic = aic;
    
end % End of loop over subjects

%% Saving subject specific bci simulations

fprintf('\n\nSaving data...\n\n');
save(['group_' model 'Simulations_repmod_best10'], ...
    'subjID', 'parameterNames', 'group_best10_results_repmod', 'group_best_results_repmod');
