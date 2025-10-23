%% EstimateModel
%   Editor:    Daniel Elbich
%   Updated:   2/28/19
%
%   Script designed to for estimating a GLM specified using the
%   SpecifyModel.m script.  Allows user to display the trial type and
%   onsets/durations in the SPM Batch GUI.
%
% See also:  SpecifyModel, createParams

%% Set Analysis Parameters & Paths
% Load all relevent project information
if exist('commandFlag','var') == 0

    %Select parameter file is flag does not exist
    [file,path]=uigetfile('*.mat','Select params file');
    filename=fullfile(path,file);
    load(filename);

end

% User Input Step 1: Options
% Set the following jobman_option to 'interactive' to view in SPM parameters the GUI.
% Press any key into the command window to continue to next one sample t test.
% Set the following jobman option to 'run' to skip the viewing of the
% SPM parameters in the GUI and go directly to running of the one
% sample t-tests

show          = 0; % 0 = input as multiple conditions file, 1 = input parameters for showing in GUI
jobman_option = 'run'; % interactive = show in GUI, run = run through SPM

%% Routine

spm('Defaults','FMRI')
spm_jobman('initcfg')
clc;



for curSub =1:length(subjects) %for curSub = number %1:length(Subjects)

    % DEB: Get subject-specific configuration
    current_subject_id = subjects{curSub};
    subject_taskInfo = subject_config(current_subject_id, taskInfo);

    fprintf('\n')
    fprintf('Processing subject %s with %d runs...\n', current_subject_id, subject_taskInfo.Runs);

    % Model Directory: directory containing this subject's model
    Model.directory = fullfile(directory.Model, current_subject_id);

    % If we are using a mask, create a path to the mask
    if Mask.on == 1
        Model.mask  = fullfile(Mask.dir, current_subject_id, Mask.name);
    end

    % Find the SpecModel *.mat files. These should be in the model
    % directory, defined above
    SpecModelMats   = cellstr(spm_select('List', Model.directory, '.*Run.*\.mat'));

    % If we do not find any *.mat files, give an error informing the user
    % that this has occured
    if cellfun('isempty', SpecModelMats)
        error('Could not find Specify Model *.mat files for %s', current_subject_id) %#ok<*NODEF>
    end

    % Get directory of motion files
    motionFiles = dir(fullfile(directory.Model, current_subject_id,['*.txt']));

    switch preprocPipeline
        case 'spm12'
            % This case is not being used for the current user, so it is left as-is.
            % It would need similar refactoring if it were to be used.
            procDataDir = fullfile(directory.Project, 'Func_ret_unsmoothed');
            procFuncFiles = dir(fullfile(procDataDir, subjects{curSub},'run*',['*wa*.nii']));
            for i = 1:taskInfo.Runs
                setenv('modelData',[procFuncFiles(i).folder filesep procFuncFiles(i).name]);
                !gunzip $modelData
                Model.runs{i}.scans = cellstr(spm_select('ExtFPList', procFuncFiles(i).folder, ['wa*'], Inf));
                Model.runs{i}.multicond = fullfile(Model.directory, SpecModelMats{i});
                Model.runs{i}.motion = [motionFiles(i).folder filesep motionFiles(i).name];
            end

        case 'fmriprep'
            procFuncFiles = dir(fullfile(directory.Project, 'derivatives',...
                'fmriprep', current_subject_id, 'func',...
                [Func.prefix '*' subject_taskInfo.Name '*_bold.nii']));

            % DEB: Define destDir once before the loop
            destDir = fullfile(directory.Model, current_subject_id);

            for i = 1:subject_taskInfo.Runs

                % Get the source file path
                sourceFile = fullfile(procFuncFiles(i).folder, procFuncFiles(i).name);

                % Copy functional to Model Directory (no gunzip needed)
                copyfile(sourceFile, destDir);

                % Define the path to the copied scan for the model
                Model.runs{i}.scans = {fullfile(destDir, procFuncFiles(i).name)};
                Model.runs{i}.multicond = fullfile(Model.directory, SpecModelMats{i});
                Model.runs{i}.motion = fullfile(motionFiles(i).folder, motionFiles(i).name);
            end
    end

    %% Set and Save the SPM job

    try
    onsets    = [];
    durations = [];
    names     = [];
    pmod      = [];

    % Directory
    matlabbatch{1}.spm.stats.fmri_spec.dir = {Model.directory};

    % Model Parameters
    matlabbatch{1}.spm.stats.fmri_spec.timing.units   = Model.units;
    matlabbatch{1}.spm.stats.fmri_spec.timing.RT      = Model.TR;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t  = 16;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 1;

    % Session Specific
    for curRun = 1:length(Model.runs)
        matlabbatch{1}.spm.stats.fmri_spec.sess(curRun).scans = Model.runs{curRun}.scans; %#ok<*AGROW>
        if show == 0
            matlabbatch{1}.spm.stats.fmri_spec.sess(curRun).cond  = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess(curRun).multi = {Model.runs{curRun}.multicond};
        elseif show == 1
            load(Model.runs{curRun}.multicond)
            for curTrialType = 1:length(names)
                matlabbatch{1}.spm.stats.fmri_spec.sess(curRun).cond(curTrialType).name     = names{curTrialType};
                matlabbatch{1}.spm.stats.fmri_spec.sess(curRun).cond(curTrialType).onset    = onsets{curTrialType};
                matlabbatch{1}.spm.stats.fmri_spec.sess(curRun).cond(curTrialType).duration = durations{curTrialType};
                matlabbatch{1}.spm.stats.fmri_spec.sess(curRun).cond(curTrialType).tmod     = 0;
                if isempty(pmod)
                    matlabbatch{1}.spm.stats.fmri_spec.sess(curRun).cond(curTrialType).pmod = struct('name', {}, 'param', {}, 'poly', {});
                else
                    matlabbatch{1}.spm.stats.fmri_spec.sess(curRun).cond(curTrialType).pmod.name  = pmod(curTrialType).name{1};
                    matlabbatch{1}.spm.stats.fmri_spec.sess(curRun).cond(curTrialType).pmod.param = pmod(curTrialType).param{1};
                    matlabbatch{1}.spm.stats.fmri_spec.sess(curRun).cond(curTrialType).pmod.poly  = pmod(curTrialType).poly{1};
                end
            end
        end
        matlabbatch{1}.spm.stats.fmri_spec.sess(curRun).regress   = struct('name', {}, 'val', {});
        matlabbatch{1}.spm.stats.fmri_spec.sess(curRun).multi_reg = {Model.runs{curRun}.motion};
        matlabbatch{1}.spm.stats.fmri_spec.sess(curRun).hpf       = 128;
    end

    % Misc
    matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
    matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
    matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
    matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
    if Mask.on == 1
        matlabbatch{1}.spm.stats.fmri_spec.mask = {Model.mask};
    else
        matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
    end
    matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';

    % Estimation
    matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep;
    matlabbatch{2}.spm.stats.fmri_est.spmmat(1).tname = 'Select SPM.mat';
    matlabbatch{2}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(1).name  = 'filter';
    matlabbatch{2}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(1).value = 'mat';
    matlabbatch{2}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(2).name  = 'strtype';
    matlabbatch{2}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(2).value = 'e';
    matlabbatch{2}.spm.stats.fmri_est.spmmat(1).sname = 'fMRI model specification: SPM.mat File';
    matlabbatch{2}.spm.stats.fmri_est.spmmat(1).src_exbranch = substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
    matlabbatch{2}.spm.stats.fmri_est.spmmat(1).src_output = substruct('.','spmmat');
    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
    catch
        disp('Unable to create matlabbatch!');
        fprintf('ERROR ON: %s', current_subject_id);
    end

    save(fullfile(Model.directory, 'Job.mat'), 'matlabbatch');
    clear onsets durations names pmod;

    %% Run the SPM job

    if strcmp(jobman_option,'interactive')
        fprintf('\n')
        fprintf('Displaying SPM Job...\n')
        spm_jobman(jobman_option, matlabbatch)
        pause
    elseif strcmp(jobman_option,'run')
        try
            spm_jobman(jobman_option, matlabbatch)
        catch ER %#ok<*NASGU>
            disp(ER)
            fprintf('ERROR ON: %s', current_subject_id)
        end
    end

    %% Manage resulting files
    switch preprocPipeline
        case 'spm12'
            % This case is not being used, leaving as is.
            subjFuncDir = fullfile(directory.Model,subjects{curSub});
            setenv('procDataDir',[procDataDir filesep subjects{curSub}]);
            setenv('subjFuncDir',subjFuncDir);
            !gzip $subjFuncDir/*.nii'
            !gzip $procDataDir/run*/*.nii

        case 'fmriprep'
            % DEB: Correctly identify and remove the temporary functional files
            for i = 1:subject_taskInfo.Runs
                tempFile = fullfile(destDir, procFuncFiles(i).name);
                if exist(tempFile, 'file')
                    delete(tempFile);
                end
            end

            % Gzip output model files
            modelFiles = dir(fullfile(destDir, '*.nii'));
            for i = 1:length(modelFiles)
                gzip(fullfile(destDir, modelFiles(i).name));
                delete(fullfile(destDir, modelFiles(i).name));
            end
    end

    % Create new SPM mat file for gzip nifti files
    load(fullfile(Model.directory, 'SPM.mat'));
    for i=1:length(SPM.Vbeta)
        SPM.Vbeta(i).fname = strrep(SPM.Vbeta(i).fname,'nii','nii.gz');
    end
    save(fullfile(Model.directory, 'SPM_gz.mat'),'SPM');

    clear SpecModelMats NumOfRuns curFuncDir curMotDir matlabbatch SPM;
    Model = rmfield(Model,'runs');
    Model = rmfield(Model,'directory');

end

clear;
clc;
disp('All finished!!');
