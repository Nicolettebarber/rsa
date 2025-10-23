%% SpecifyModel
%   Editor:    Daniel Elbich
%   Updated:   7/28/23
%
%   Script designed to build multiple conditions files for later use with
%   SPM's matlabbatch system.
%
%   Assumes that the behavioral data are organized as follows:
%
%   /StudyDir/BehavDir/s001/ConcatenatedBehavioralData.csv
%   /StudyDir/BehavDir/s002/ConcatenatedBehavioralData.csv
%   /StudyDir/BehavDir/s003/ConcatenatedBehavioralData.csv
%
%   Where ConcatenatedBehavioralData.csv is a comma seperated text file
%   where all functional runs are concatenated, with a column that
%   indicates which rows (i.e., tonsetsrials) belong to which functional run.
%
%   See also createParams.m, EstimateModel.m
%
%   Updates:
%
%   5/6/19 - Number of trials, durations, and onsets saved as csv file for
%   data QA.

%   Note: runRSA and runRSASearchlight require tags to be in lowercase -
%   change accordingly prior to running this script.
%
%
% TO CHANGE BELOW:
%    Change lines 48 to indicate file name, and edit lines 113-116



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set Analysis Parameters & Paths
%Load all relevent project information
 if exist('commandFlag','var') == 0

     %Select parameter file is flag does not exist
     [file,path]=uigetfile('*.mat','Select params file');
     filename=fullfile(path,file);
     load(filename);

 end

%% Main Code
for i = 1:length(subjects)

    % DEB: Get subject-specific configuration
    current_subject_id = subjects{i};
    subject_taskInfo = subject_config(current_subject_id, taskInfo);

    % Creates path to the current subjects behavioral file
    curSubj.behavDir  = [directory.Project filesep current_subject_id filesep rawData.behavDir]; %For fMRIprep
    curSubj.behavFile = dir([curSubj.behavDir filesep '*' subject_taskInfo.Name '*.' rawData.behavFile]); %change to look for enc files or ret files

    % Creates a path to this subjects analysis directory & creates that
    % directory if it does not already exist.
    curSubj.directory = fullfile(directory.Model, current_subject_id);
    if ~isdir(curSubj.directory)
        mkdir(curSubj.directory)
    end

    fprintf('Processing subject %s with %d runs...\n', current_subject_id, subject_taskInfo.Runs);
    fprintf('Sorting Behavioral Data...\n\n')

    % Build the multiple conditions *.mat file for each run
    for curRun = 1:subject_taskInfo.Runs

        % Reads in the subjects behavioral data using the readtable command.
        fprintf('Reading in subject %s behavioral data for run %d...\n', current_subject_id, curRun);
        BehavData = readtable([curSubj.behavDir filesep curSubj.behavFile(curRun).name],'FileType','text');

        %-- Initialize the names, onsets, durations, and pmods structure arrays
		onsets = num2cell(BehavData.onset)';
		durations = num2cell(BehavData.duration)';

        %-- Filter out trials of no interest
        trialsOfInterest = ~strcmp(BehavData.trial_type, 'no_interest');
        onsets = onsets(trialsOfInterest);
        durations = durations(trialsOfInterest);
        lbin = BehavData.LBin(trialsOfInterest);

        %-- Create condition names based on LBin similarity ratings
        names = cell(length(lbin),1);
        for j = 1:length(lbin)
            names{j} = ['Similarity_' num2str(lbin(j))];
        end
        names = names';

        %-- Save the Multiple Conditions *.mat file
        % Save the names, onsets, and durations variables to .mat file
        % to be used for later model estimation in SPM. See EstimateModel.m

        matfilename = fullfile(curSubj.directory, ...
            ['Run', num2str(curRun, '%03d'), '_multiple_conditions.mat']);
        fprintf('Saving subject %s run %d multiple conditions file...\n',...
            subjects{i}, curRun)
        save(matfilename, 'names', 'onsets', 'durations');

        % Summary Trial Numbers
        if i==1
            labelHeader={['names_' num2str(curRun, '%03d')],...
                ['onsets_' num2str(curRun, '%03d')],...
                ['durations_' num2str(curRun, '%03d')]};

            finalSummary{1,1}='Subject_ID';
            index=size(finalSummary,2);
            index=[index+1 index+3];
            finalSummary(1,index(1):index(2))=labelHeader;
        end

        if curRun==1
            finalSummary{i+1,1}=subjects{i};
            finalIndex=[2:4];
        end

        finalSummary{i+1,finalIndex(1)}=num2str(length(names));
        finalSummary{i+1,finalIndex(2)}=num2str(length(onsets));
        finalSummary{i+1,finalIndex(3)}=num2str(length(durations));

        finalIndex=finalIndex+3;

    end

    %ersTagFilename = fullfile(curSubj.directory,'ersTags.mat');
    %save(ersTagFilename,'ersTags');

end


%% Write output summary file

file = fopen([directory.Model filesep 'fullTrialCountQC.csv'], 'w');

for a=1:size(finalSummary,1)
    for b=1:size(finalSummary,2)
        var = eval('finalSummary{a,b}');
        try
            fprintf(file, '%s', var);
        end
        fprintf(file, ',');
    end
    fprintf(file, '\n');
end

fclose(file);
clc;
disp('All Finished!!');
clear;
