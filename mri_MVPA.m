%% MVPA SVR analysis using TDT toolbox

clear;
close all;
clc;

addpath(genpath('E:\AMBRA\UoB\Programs\spm12'));
addpath(genpath('E:\AMBRA\UoB\Programs\decoding_toolbox_v3.96'));

%% Define some settings
PID_list = {'sub-MA01'};
driveLetter = 'E:\AMBRA\UoB\Data';

sessionNr = 7; % 14 or 7 (concat2runs)

extraNotes = ['_' num2str(sessionNr) '_sessions_nosmooth_tempder_concat'];

%What kind of analysis to use
%'roi' or 'searchlight'
decoding_type = 'roi';

if strcmp(decoding_type, 'roi')
    decoding_region_list = {'V1-3';'IPS0-2';'IPS3-4_SPL1';'TE1.0-1.1';'PT'};
else
    decoding_region_list  = {'all_brain'};
end

%Select training data
%'cong_cong','cong_incong'
trainTestScheme_list = {'cong_incong'};

if sessionNr == 14
    glmFolderName = 'GLM';
elseif sessionNr == 7 && ~strcmp(trainTestScheme_list,'cong_cong_resp')
    glmFolderName = 'GLM2';
    scaling = 'min0max1_concat_scale';
elseif sessionNr == 7 && strcmp(trainTestScheme_list,'cong_cong_resp')
    glmFolderName = 'GLM3';
    scaling = 'min0max1_concat_scale';
end

for iScheme = 1:length(trainTestScheme_list)
    
    trainTestScheme = trainTestScheme_list{iScheme};
    
    if ~isempty(regexp(trainTestScheme, 'cong_cong','once'))
        condNr = 12;
    else
        condNr = 36;
    end
    
    for iSubj = 1:length(PID_list)
        
        PID = PID_list{iSubj};
        
        for iroi = 1:length(decoding_region_list)
            
            decoding_region = decoding_region_list{iroi};
            
            %% Set folders with the above settings
            
            if strcmp(decoding_type, 'roi')
                thisAnalysisName = [decoding_type '_' decoding_region ...
                    '_train_test_' trainTestScheme extraNotes];
            elseif strcmp(decoding_type, 'searchlight')
                thisAnalysisName = [decoding_type ...
                    '_train_test_' trainTestScheme extraNotes];
            end
            
            roiFolder = fullfile(driveLetter,'ROIs\MAMSI',PID);
            glmFolder = fullfile(driveLetter,'MAMSI_MRI',...
                PID,'derivatives','ve','MVPA',glmFolderName);
            workingFolder = fullfile(driveLetter,'MAMSI_MRI',...
                PID,'derivatives','ve','MVPA',scaling,thisAnalysisName);
            betasFolder = fullfile(driveLetter,'MAMSI_MRI',...
                PID,'derivatives','ve','MVPA',scaling,thisAnalysisName,'betas');
            
            if ~isfolder(roiFolder)
                mkdir(roiFolder);
            end
            
            if ~isfolder(glmFolder)
                mkdir(glmFolder);
            end
            
            if ~isfolder(workingFolder)
                mkdir(workingFolder);
            end
            
            if length(dir(glmFolder))==2
                if sessionNr == 14
                    copyfile(fullfile(driveLetter,'MAMSI_MRI',...
                        PID,'derivatives','ve', 'anl_ve_stim_native_smooth0'),glmFolder);
                elseif sessionNr == 7 && ~strcmp(trainTestScheme_list,'cong_cong_resp')
                    copyfile(fullfile(driveLetter,'MAMSI_MRI',...
                        PID,'derivatives','ve', 'anl_ve_stim_concat2runs_native_smooth0'),glmFolder);
                elseif sessionNr == 7 && strcmp(trainTestScheme_list,'cong_cong_resp')
                    copyfile(fullfile(driveLetter,'MAMSI_MRI',...
                        PID,'derivatives','ve', 'anl_ve_stim_resp_reg_concat2runs_native_smooth0'),glmFolder);
                end
            end
            
            if ~isfolder(betasFolder)
                mkdir(betasFolder);
            end
            
            cd(workingFolder);
            
            %% Set early cfg parameters
            
            clear cfg
            cfg = decoding_defaults;
            cfg.results.dir = workingFolder;
            
            %% Set ROI mask
            
            cfg.analysis = decoding_type;
            
            switch lower(decoding_type)
                
                case 'searchlight' %If searchlight, use whole-brain mask
                    cfg.files.mask = fullfile(glmFolder,'mask.nii');
                    
                case 'roi' %Otherwise use specified ROI masks
                    switch decoding_region
                        case 'V1-3'
                            cfg.files.mask = fullfile(roiFolder,'rwV1-3.nii');
                        case 'IPS0-2'
                            cfg.files.mask = fullfile(roiFolder,'rwIPS0-2.nii');
                        case 'IPS3-4_SPL1'
                            cfg.files.mask = fullfile(roiFolder,'rwIPS3-4_SPL1.nii');
                        case 'PT'
                            cfg.files.mask = fullfile(roiFolder,'rPT_aparc2009s.nii');
                        case 'TE1.0-1.1'
                            cfg.files.mask = fullfile(roiFolder,'rwTE1.0-1.1.nii');
                    end
            end
            
            %% Set up design (based on specified train/test scheme)
            
            %Get regressor names from SPM file
            regressor_names = design_from_spm(glmFolder);
            condIndex = zeros(0,length(regressor_names));
            
            switch trainTestScheme
                
                case 'cong_cong'
                    
                    %Is this cross-classification (train on different data from test)?
                    isXclass = 0;
                    
                    %Specify names of regressors we're interested in
                    inclCond = {'AL_VL';'AC_VC';'AR_VR'};
                    inclBin = '1';
                    
                    %Get an index of the location of these regressors
                    for iIncl = 1:length(inclCond)
                        tempCond = regexp(regressor_names(1,:),inclCond{iIncl}, 'once');
                        condIndex(size(condIndex,1)+1,:) = cellfun(@(x) ~isempty(x),tempCond);
                    end
                    condIndex = any(condIndex);
                    
                    % Get bin number info (1 = HRF; 2 = TEMPORAL DERIVATIVE; 3 = DISPERSION DERIVATIVE)
                    [~,bin_names] = strtok(regressor_names(1,:));
                    % Index of desired bin
                    tempBin = regexp(bin_names(1,:),inclBin, 'once');
                    binIndex = cellfun(@(x) ~isempty(x),tempBin);
                    
                    % Get beta index that takes into account conditions (specified by inclCond) and derivatives (specified by inclBin)
                    betaIndex = and(binIndex,condIndex);
                    
                    % And use the index to filter regressor_names
                    regressor_names = regressor_names(:,betaIndex);
                    
                    labels = repmat([-9 0 9],1,4);
                    
                case 'cong_incong'
                    
                    %Is this cross-classification (train on different data from test)?
                    isXclass = 1;
                    
                    %Specify names of regressors we're interested in
                    inclCond = {'RespA';'RespV'};
                    inclBin = '1';
                    
                    %Get an index of the location of these regressors
                    for iIncl = 1:length(inclCond)
                        tempCond = regexp(regressor_names(1,:),inclCond{iIncl}, 'once');
                        condIndex(size(condIndex,1)+1,:) = cellfun(@(x) ~isempty(x),tempCond);
                    end
                    condIndex = any(condIndex);
                    
                    % Get bin number info (1 = HRF; 2 = TEMPORAL DERIVATIVE; 3 = DISPERSION DERIVATIVE)
                    [~,bin_names] = strtok(regressor_names(1,:));
                    % Index of desired bin
                    tempBin = regexp(bin_names(1,:),inclBin, 'once');
                    binIndex = cellfun(@(x) ~isempty(x),tempBin);
                    
                    % Get beta index that takes into account conditions (specified by inclCond) and derivatives (specified by inclBin)
                    betaIndex = and(binIndex,condIndex);
                    
                    % And use the index to filter regressor_names
                    regressor_names = regressor_names(:,betaIndex);
                    
                    labels = repmat([-9 -9 -9 0 0 0 9 9 9],1,4); % true auditory labels (repeated 4 times to pool over AttxResp experimental conditions)
                    xclass = repmat([1 2 2 2 1 2 2 2 1],1,4); % 1 = cong; 2 = incong
            end
            
            %% Select and move the correct beta images, then generate design
            
            %Get label names from regressor list
            labelnames = regressor_names(1,cell2mat(regressor_names(2,:))==1)';
            
            %Subselect correct betas
            tempStruct = dir(fullfile(glmFolder,'beta_*.nii'));
            allBetas = {tempStruct.name}';
            relevantBetas = allBetas(betaIndex);
            
            if length(dir(betasFolder))==2
                % move these to the relevant analysis folder
                for iBeta = 1:length(relevantBetas)
                    copyfile(fullfile(glmFolder,relevantBetas{iBeta}),betasFolder);
                end
                if sessionNr == 7
                    cfg.scale.meanbetas = 0;
                elseif sessionNr == 14
                    cfg.scale.meanbetas = 1;
                end
            else
                if sessionNr == 7
                    cfg.scale.meanbetas = 0;
                elseif sessionNr == 14
                    cfg.scale.meanbetas = 1;
                end
            end
            
            %Generate the design (different depending on whether we generalise)
            if isXclass
                cfg = decoding_describe_data(cfg,{labelnames},labels,regressor_names,betasFolder,xclass);
                cfg.design = make_design_xclass_cv(cfg);
            else
                cfg = decoding_describe_data(cfg,{labelnames},labels,regressor_names,betasFolder);
                cfg.design = make_design_cv(cfg);
            end
            
            %% Set final parameters, plot design, and run
            
            cfg.scale.method = 'min0max1';
            cfg.scale.estimation = 'all'; % 'all' (recommended), 'across', 'separate', or 'none'
            
            % Analysis type (classification with accuracy_minus_chance set as default)
            cfg.decoding.method = 'regression';
            cfg.results.output = {'corr', 'predicted_labels'};
            
            % in command window
            display_design(cfg);
            % write .nii file (defauls = 1)
            if strcmp(cfg.decoding.method,'regression')
                cfg.results.write = 2;
            end
            
            cfg.results.overwrite = 1;
            
            % design summary
            events = cfg.files.descr';
            events_labels = cfg.files.label;
            events_betas = cell(size(cfg.files.name,1),1);
            if isempty(regexp(extraNotes,'mean_betas','match'))
                for i = 1:size(events_betas,1)
                    events_betas{i} = cfg.files.name{i}(end-12:end);
                end
            else
                for i = 1:size(events_betas,1)
                    events_betas{i} = cfg.files.name{i}(end-27:end);
                end
            end
            design = dataset(events,events_labels,events_betas);
            save('design_info', 'design');
            
            % run decoding
            if strcmp(cfg.scale.method,'cov')
                passed_data = [];
                [results,cfg] = decoding(cfg,passed_data,misc);
            else
                [results,cfg] = decoding(cfg);
            end
            saveas(gcf,'cv_design','png');
            
        end % roi
    end % subj
end % scheme
