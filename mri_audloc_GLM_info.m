%% GLM information (names, onsets, durations) and contrasts

clear;
clc;

% We are in the data directory
rootPath = pwd;

% Select subject data folder
subjID = 'sub-MA01';
dataPath = fullfile('E:\Data\MAMSI_MRI', subjID, 'behav\audloc');
savePath = fullfile('E:\Data\MAMSI_MRI', subjID, 'derivatives\audloc\GLM_info');
if ~isfolder(savePath) % create new folder if needed
    mkdir(savePath);
end
cd(dataPath);

% Get list of .mat files
Files = dir('*Aud_Session*.mat');
runNr = size(Files,1);

names = cell(1,3);
names{1} = 'L';
names{2} = 'C';
names{3} = 'R';

onsets = cell(1,3);

durations = cell(1,3);
for i = 1:length(durations)
    durations{1,i} = 0;
end

for iFile = 1:runNr
    % load
    file = Files(iFile).name;    
    load(file);
    % recode onsets
    AVOnset = tdata.AVOnset(~isnan(tdata.AVOnset)) - S.startTask;
    % collect condition-specific onsets
    LId = find(tdata.TargetLoc(~isnan(tdata.TargetLoc)) == -9);
    onsets{:,1} = AVOnset(LId);
    CId = find(tdata.TargetLoc(~isnan(tdata.TargetLoc)) == 0);
    onsets{:,2} = AVOnset(CId);
    RId = find(tdata.TargetLoc(~isnan(tdata.TargetLoc)) == 9);
    onsets{:,3} = AVOnset(RId);
    % save
    GLM = fullfile([subjID '_task_aloc_GLM_' num2str(iFile)]);
    cd(savePath);
    save(GLM, 'names', 'onsets', 'durations');
    cd(dataPath);
end

%% Create GLM constrasts

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
save(['audloc_contrasts_2lev_anova_sessions_' num2str(size(Files,1))], 'con_tot', 'names');
cd(rootPath);
