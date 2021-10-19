function varargout = fitModelAcrossAtt(parameters,parameterNames,dataVA,responseLoc,model,varargin)

%% Checking input

% Parsing input
p = inputParser;

validModels = {'bci','fus','segA','segV','taskRel'};
validOutputs = {'logLike','full'};

addRequired(p,'parameters',@(x) validateattributes(x,{'numeric'},{'vector'}));
addRequired(p,'parameterNames',@(x) validateattributes(x,{'cell'},...
    {'numel',numel(parameters)}));
addRequired(p,'dataVA',@(x) validateattributes(x,{'table'},{'nonempty'}));
addRequired(p,'responseLoc',@(x) validateattributes(x,{'numeric'},{'vector'}));
addRequired(p,'model',@(x) any(validatestring(x,validModels)));
addOptional(p,'decisionFun',1,@(x) validateattributes(x,{'numeric'},...
    {'scalar','integer','positive','<=',3}));
addOptional(p,'output','logLike',@(x) any(validatestring(x,validOutputs)));
parse(p,parameters,parameterNames,dataVA,responseLoc,model,varargin{:});

parameters = p.Results.parameters;
parameterNames = p.Results.parameterNames;
dataVA = p.Results.dataVA;
responseLoc = p.Results.responseLoc;
model = p.Results.model;
decisionFun = p.Results.decisionFun;
output = p.Results.output;

% Making sure parameters and parameterNames are row vectors
if iscolumn(parameters), parameters = parameters'; end
if iscolumn(parameterNames), parameterNames = parameterNames'; end

%% Fit across attention levels

% get number of attention levels
attLevels = unique(dataVA.AttentionModality);
nLevelsAtt = numel(attLevels);

sigVnames = parameterNames(~cellfun(@isempty,regexp(parameterNames,'sigV[0-9]*')));
sigAnames = parameterNames(~cellfun(@isempty,regexp(parameterNames,'sigA[0-9]*')));
otherParamIdxV = cellfun(@isempty,regexp(parameterNames,'sigV'));
otherParamIdxA = cellfun(@isempty,regexp(parameterNames,'sigA'));
otherParamIdx = otherParamIdxV==1&otherParamIdxA==1;

% initialise log likelihood output
logLike_sum = 0;
% initialise structure for bci details in case they are requested
if strcmp(output,'full')
    bci_details = struct([]);
end

% loop over number of attention levels (p_common and sigP kept constant; sigA and sigV let vary with attention levels)
for iAtt = 1:nLevelsAtt
    
    % Select data for each attention level
    actDataVA = dataVA((dataVA.AttentionModality == attLevels(iAtt)),:);
    
    % Select visual variance and auditory variance parameters for this attention level
    if strcmp(model,'segA')
        actSigA = parameters(:,ismember(parameterNames,sigAnames{iAtt}));
        % generate the parameter setting for one attention level
        actParameters = [parameters(:,otherParamIdx),actSigA];
        actParamNames = [parameterNames(otherParamIdx),'sigA'];
    elseif strcmp(model,'segV')
        actSigV = parameters(:,ismember(parameterNames,sigVnames{iAtt}));
        % generate the parameter setting for one attention level
        actParameters = [parameters(:,otherParamIdx),actSigV];
        actParamNames = [parameterNames(otherParamIdx),'sigV'];
    else
        actSigV = parameters(:,ismember(parameterNames,sigVnames{iAtt}));
        actSigA = parameters(:,ismember(parameterNames,sigAnames{iAtt}));
        % generate the parameter setting for one attention level
        actParameters = [parameters(:,otherParamIdx),actSigA,actSigV];
        actParamNames = [parameterNames(otherParamIdx),'sigA','sigV'];
    end
    
    switch output
        
        case 'logLike' % one single output: log likelihood
            
            switch model
                case 'bci'
                    logLike = fitModel(actParameters,actParamNames,actDataVA,responseLoc,decisionFun);
                    % insert new cases for different models ('fus','segA','segV','taskRel')
                case 'fus'
                    logLike = fitModelFus(actParameters,actParamNames,actDataVA,responseLoc);
                case 'segA'
                    logLike = fitModelSegA(actParameters,actParamNames,actDataVA,responseLoc);
                case 'segV'
                    logLike = fitModelSegV(actParameters,actParamNames,actDataVA,responseLoc);
                case 'taskRel'
                    logLike = fitModelTaskRel(actParameters,actParamNames,actDataVA,responseLoc);
            end
            
        case 'full' % two output: log likelihood and bci model details
            
            switch model
                case 'bci'
                    [logLike,det] = fitModel(actParameters,actParamNames,actDataVA,responseLoc,decisionFun);
                    % insert new cases for different models ('fus','segA','segV','taskRel')
                case 'fus'
                    [logLike,det] = fitModelFus(actParameters,actParamNames,actDataVA,responseLoc);
                case 'segA'
                    [logLike,det] = fitModelSegA(actParameters,actParamNames,actDataVA,responseLoc);
                case 'segV'
                    [logLike,det] = fitModelSegV(actParameters,actParamNames,actDataVA,responseLoc);
                case 'taskRel'
                    [logLike,det] = fitModelTaskRel(actParameters,actParamNames,actDataVA,responseLoc);
            end
            
            temp = repmat({attLevels(iAtt)},size(det));
            [det.attMod] = temp{:};
            bci_details = cat(2,bci_details,det);
    end
    
    % sum logLikes over attention levels
    logLike_sum = logLike_sum + logLike;
    
end

varargout{1} = logLike_sum;
if strcmp(output,'full')
    % Computing Bayesian Information Criterion and Akaike Information Criterion
    n = size(dataVA,1);
    k = numel(parameterNames);
    
    bic = -logLike_sum-0.5*k*log(n);
    aic = 2*k + 2*logLike_sum;
    varargout{2} = bic;
    varargout{3} = aic;
    varargout{4} = bci_details;
end

end
