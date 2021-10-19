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

%% Loading data

expID = 'behav'; % MRI, behav

% Subjects
if strcmp(expID, 'MRI')
    subjID = {'sub-MA01';'sub-MA08';'sub-MA53';'sub-MA98';'sub-MA114';'sub-MA117';...
        'sub-MA124';'sub-MA129';'sub-MA130';'sub-MA132';'sub-MA138';'sub-MA145'};
elseif strcmp(expID, 'behav')
    subjID = {'sub-MA01';'sub-MA08';'sub-MA53';'sub-MA98';...
        'sub-MA103';'sub-MA104';'sub-MA109';'sub-MA114';...
        'sub-MA115';'sub-MA117';'sub-MA122';...
        'sub-MA124';'sub-MA125';'sub-MA126';'sub-MA129';'sub-MA130';...
        'sub-MA132';'sub-MA134';'sub-MA137';'sub-MA138';...
        'sub-MA142';'sub-MA145';'sub-MA146';'sub-MA147';...
        'sub-MA149';'sub-MA151';'sub-MA152'};
end

% number of minima to evaluate with fminsearch
nmin = 10;

% model
model='bci'; % bci fus segA segV taskRel null

figure_plot=0;

for iSubj = 1:length(subjID)
    
    if strcmp(expID, 'MRI')
        cd(fullfile('E:\Data\MAMSI_MRI', subjID{iSubj}, 'behav\scanner'));
    elseif strcmp(expID, 'behav')
        cd(fullfile('E:\Data\MAMSI_MRI', 'MAMSI_MRI_behav', subjID{iSubj}));
    end
    
    origData = dir('*Exp_All_Sessions*.mat');
    load(origData.name, 'tdata');
    origData = tdata;
    
    % Reorganize dataset 'tdata' so that I get the information I need for model estimation (locA, locV, respA, respV)
    % Get rid of incorrect responses (wrong keypad), missed responses (no answer provided) and anticipated responses (RT < 100ms)
    origData = origData(origData.IncorResp==0 & origData.MissedResp==0 & origData.AntResp==0,:);
    
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
    
    % Keep only useful data (locV, locA, respV, respA)
    bciData = origData(:,{'locV','locA','respV','respA'});
    actDataVA = dataset2table(bciData);
    
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
    sigV = linspace(0.1,30,n);
    sigA = linspace(0.1,30,n);
    
    % These parameters can be used as well, but not included in this analysis
    % kernelWidth = 1; 'kW'
    % muP = 0; % mean of spatial prior
    
    % Parameter names
    if strcmp(model, 'fus') || strcmp(model, 'taskRel')
        parameterNames = {'sigP','sigA','sigV'};
        gridVectors = {sigP,sigA,sigV};
    elseif strcmp(model, 'segA')
        parameterNames = {'sigP','sigA'};
        gridVectors = {sigP,sigA};
    elseif strcmp(model, 'segV')
        parameterNames = {'sigP','sigV'};
        gridVectors = {sigP,sigV};
    else
        parameterNames = {'p_common','sigP','sigA','sigV'};
        gridVectors = {p_common,sigP,sigA,sigV};
    end
    nParameters = numel(gridVectors);
    % Full factorial expansion of the specified parameter vectors
    coords = cell(1,nParameters);
    [coords{:}] = ndgrid(gridVectors{:});
    coords = cellfun(@(x) x(:),coords,'UniformOutput',false);
    % Matrix of all possible parameter combinations: nRows = number of
    % combinations, columns correspond to parameters
    paramCombinations = cat(2,coords{:});
    
    % This is actually relevant if you have two or more V reliabilities (it works with one nevertheless)
    sigVnames = parameterNames(~cellfun(@isempty,regexp(parameterNames,'sigV[0-9]*')));
    otherParamIdx = cellfun(@isempty,regexp(parameterNames,'sigV[0-9]*'));
    % Select visual variance parameter for this reliability level
    actSigV = paramCombinations(:,ismember(parameterNames,sigVnames));
    % generate the parameter setting for one visual reliability level
    actParameters = [paramCombinations(:,otherParamIdx),actSigV];
    actParamNames = [parameterNames(otherParamIdx),'sigV'];
    
    % Pre-allocating array for data collection
    logLike_all = NaN(size(paramCombinations,1),1);
    grid_10bestIdx = NaN(nmin,1);
    grid_10bestParameters = cell(nmin,1);
    grid_10bestlogLike = NaN(nmin,1);
    fm_10bestParameters = cell(nmin,1);
    fm_10bestErrorval = NaN(nmin,1);
    fm_10bestLogLike = NaN(nmin,1);
    bciSimulations = struct([]);
    
    % Starting the timer and printing details
    cStart = clock;
    fprintf('Performing gridsearch... \n');
    % Performing grid search on the specified parameter space
    switch model
        case 'bci'
            parfor i = 1:size(actParameters,1)
                [logLike_all(i)] = bci_fitmodel(actParameters(i,:),actParamNames,actDataVA,responseLoc,readout);
            end
        case 'fus'
            parfor i = 1:size(actParameters,1)
                [logLike_all(i)] = fitModelFus(actParameters(i,:),actParamNames,actDataVA,responseLoc);
            end
    end
    
    % Printing elapsed time
    cTemp = clock;
    fprintf('gridsearch elapsed time (days hours:minutes:seconds) %s \n\n',...
        datestr(etime(cTemp,cStart)/86400,'dd HH:MM:SS'));
    
    % Finding the n minima of the log-likelihood values and the corresponding parameter combination
    % sort gridsearch log-likelihoods from the 1st best (global minimum)
    logLike_sort = sort(logLike_all);
    
    % get id number and parameter combination for the n best log-likelihoods
    for i = 1:nmin
        grid_10bestIdx(i,1) = find(logLike_all==logLike_sort(i));
        grid_10bestParameters{i,1} = actParameters(grid_10bestIdx(i),:);
    end
    
    % get n best log-likelihoods and save in subj-specific structure
    grid_10bestlogLike(:,1) = logLike_sort(1:nmin);
    
    % Plot log-likelihoods (with n best marked with a *)
    cols = {[103 0 13];[133 6 17];[165 15 21];[203 24 29];[213 34 35];...
        [239 59 44];[251 106 74];[252 146 114];[252 187 161];[254 224 210]};
    if figure_plot==1
        figure('Position', [0, 0, 1920, 1080]); plot(logLike_all,'k'); hold on;
        for i = 1:nmin
            plot(grid_10bestIdx(i,1),logLike_all(grid_10bestIdx(i,1)), '*', 'Color', cols{i}/255); hold on;
        end
        saveas(gcf,[subjID{iSubj} '_gridsearch_logLike' ],'emf');
    end
    
    %% Perform fminsearch to refine the best parameters
    % Take the best parameters found in the gridsearch and use them as starting
    % point in the fmincon
    
    % Options for fminsearchbnd
    opts = optimset('fminsearch');
    opts.Display = 'iter';
    opts.MaxFunEvals = 4000;
    opts.TolFun = 1.e-12;
    opts.MaxIter = 1000;
    
    % set upper and lower bounds on parameters
    % 'p_common','sigP','sigA1', 'sigA2', 'sigV1', 'sigV2'
    if strcmp(model, 'fus') || strcmp(model, 'taskRel')
        LB = [0.001 0.001 0.001];
        UB = [30    30    30];
    elseif strcmp(model, 'segA') || strcmp(model, 'segV')
        LB = [0.001 0.001];
        UB = [30    30];
    else
        LB = [0 0.001 0.001 0.001 ];
        UB = [1 30    30    30    ];
    end
    
    % Creating anonymous function for input to fmincon
    switch model
        case 'bci'
            fun = @(param) bci_fitmodel(param,actParamNames,actDataVA,responseLoc,readout);
        case 'fus'
            fun = @(param) fitModelFus(param,actParamNames,actDataVA,responseLoc);
    end
    
    % Performing fminsearch for each parameter combination
    for i = 1:nmin
        fprintf(['Performing fminsearch ' num2str(i) ' of subject ' num2str(iSubj) '... \n']);
        [fm_10bestParameters{i,1},fm_10bestErrorval(i,1)] = fminsearchbnd(fun,grid_10bestParameters{i,1},LB,UB,opts);
    end
    
    % Finishing timer and printing elapsed time
    fprintf('fminsearch elapsed time (days hours:minutes:seconds) %s \n\n',...
        datestr(etime(clock,cTemp)/86400,'dd HH:MM:SS'));
    
    %% Evaluating model fit of best model
    
    % Save all simulation data and plot
    my_actions = 1;
    
    % Evaluate model fit of the n best parameters combinations
    switch model
        case 'bci'
            for i = 1:nmin
                [fm_10bestLogLike(i,1)] = bci_fitmodel(fm_10bestParameters{i,1},actParamNames,actDataVA,responseLoc,readout);
            end
        case 'fus'
            for i = 1:nmin
                [fm_10bestLogLike(i,1)] = fitModelFus(fm_10bestParameters{i,1},actParamNames,actDataVA,responseLoc);
            end
    end
    
    % Find the best parameters combination after fminsearch
    bestSim = find(fm_10bestLogLike==min(fm_10bestLogLike));
    fm_bestParameters = fm_10bestParameters{bestSim};
    
    % Produce bci simulations for the best parameters combination
    switch model
        case 'bci'
            [fm_bestlogLike,all] = bci_fitmodel(fm_bestParameters,actParamNames,actDataVA,responseLoc,readout);
        case 'fus'
            [fm_bestlogLike,all] = fitModelFus(fm_bestParameters,actParamNames,actDataVA,responseLoc);
    end
    
    bciSimulations = cat(2,bciSimulations,all);
    
    % Saving subject specific bci simulations
    fprintf('\n\nSaving data...\n\n');
    save([subjID{iSubj} '_' model 'Simulations_gridsearch'], ...
        'subjID', 'logLike_all', 'actParameters', 'actParamNames', 'actDataVA');
    
    % Computing Bayesian Information Criterion and Akaike Information Criterion
    ndata = size(bciData,1);
    k = numel(actParamNames);
    
    bic = -fm_bestlogLike-0.5*k*log(ndata);
    aic = 2*k + 2*fm_bestlogLike;
    
    %% Saving subject specific bci simulations
    fprintf('\n\nSaving data...\n\n');
    save([subjID{iSubj} '_' model 'Simulations_best10'],...
        'subjID', 'logLike_all', 'grid_10bestIdx', 'grid_10bestParameters',...
        'grid_10bestlogLike', 'fm_10bestParameters', 'fm_10bestErrorval',...
        'fm_10bestLogLike', 'fm_bestParameters', 'fm_bestlogLike', 'bciSimulations',...
        'ndata', 'k', 'bic', 'aic', '-v7.3');
    
end % End of loop over subjects
