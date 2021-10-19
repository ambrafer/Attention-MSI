%% GLM information (names, onsets, durations) and contrasts

clear;
clc;

% We are in the data directory
rootPath = pwd;
subjID = 'sub-MA01';

% Select subject data folder
dataPath = fullfile('E:\Data\MAMSI_MRI', subjID, 'behav\scanner');
savePath = fullfile('E:\Data\MAMSI_MRI', subjID, 'derivatives\ve\GLM_info');
if ~isfolder(savePath) % create new folder if needed
    mkdir(savePath);
end
cd(dataPath);

% Get list of .mat files for separate sessions
Files = dir('*Session_*.mat');

%% Collect data for SPM design matrix

for iFile = 1:size(Files,1)
    
    % load
    file = Files(iFile).name;
    load(file);
    
    AllCdt = unique(S.TrialList, 'rows');
    AllCdt(~all(AllCdt,2),:) = []; % remove fixation and attentional cue
    TrialList = S.TrialList;
    TrialList(~all(TrialList,2),:) = []; % remove fixation and attentional cue
    
    names = cell(1,size(AllCdt,1)+2);
    for i = 1:length(names)-2
        % Attention
        if AllCdt(i,1) == 1 % A
            Att = 'A';
        elseif AllCdt(i,1) == 2 % V
            Att = 'V';
        end
        % Response
        if AllCdt(i,2) == 1 % A
            Resp = 'A';
        elseif AllCdt(i,2) == 2 % V
            Resp = 'V';
        end
        % A location
        if AllCdt(i,3) == 1 % Left
            Aloc = 'L';
        elseif AllCdt(i,3) == 2 % Centre
            Aloc = 'C';
        elseif AllCdt(i,3) == 3 % Right
            Aloc = 'R';
        end
        % V location
        if AllCdt(i,4) == 1 % Left
            Vloc = 'L';
        elseif AllCdt(i,4) == 2 % Centre
            Vloc = 'C';
        elseif AllCdt(i,4) == 3 % Right
            Vloc = 'R';
        end
        % Name for current condition
        nameCond = ['Att' Att '_Resp' Resp '_A' Aloc '_V' Vloc];
        names{i} = nameCond;
    end
    names{size(AllCdt,1)+1} = 'AttACue';
    names{size(AllCdt,1)+2} = 'AttVCue';
    
    durations = cell(1,size(AllCdt,1)+2);
    for i = 1:length(durations)
        durations{i} = 0;
    end
    
    onsets = cell(1,size(AllCdt,1)+2);
    
    % recode onsets
    AVOnset = tdata.AVOnset(~isnan(tdata.AVOnset)) - tdata.CueOnset(1);
    % collect condition-specific onsets
    for i=1:size(AllCdt,1)
        onsets{1,i} = AVOnset(all(ismember(TrialList, AllCdt(i,:),'rows'),2));
    end
    % A Attention
    onsets{size(AllCdt,1)+1} = tdata.CueOnset(~isnan(tdata.CueOnset) & strcmp(tdata.AttentionModality,'Aud'))-tdata.CueOnset(1);
    % V Attention
    onsets{size(AllCdt,1)+2} = tdata.CueOnset(~isnan(tdata.CueOnset) & strcmp(tdata.AttentionModality,'Vis'))-tdata.CueOnset(1);
    
    cd(savePath);
    if iFile < 10
        GLM = fullfile([subjID '_task_ve_stim_GLM_0' num2str(iFile)]);
    else
        GLM = fullfile([subjID '_task_ve_stim_GLM_' num2str(iFile)]);
    end
    save(GLM, 'names', 'onsets', 'durations');
    cd(dataPath);
end

%% Create GLM constrasts for VE experiment

contrastPath = fullfile('E:\Data\MAMSI_MRI\contrasts');
con_tot = cell(size(names,2),1);

for iname = 1:size(names,2)
    con = [];
    temp = regexp(names,names{iname}, 'once');
    temp = cellfun(@(x) ~isempty(x),temp);
    for t = 1:length(temp)
        con = cat(2, con, [temp(t) 0]); % discount 1st derivative
    end
    con = cat(2, con, zeros(1,6)); % discount 6 motion realignment parameters
    con = repmat(con,1,size(Files,1)); % all runs
    con_tot{iname} = con;
end

%% Save

cd(contrastPath);
save(['ve_stim_contrasts_2lev_anova_sessions_' num2str(size(Files,1))], 'con_tot', 'names');
 
cd(rootPath);
