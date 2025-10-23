function run_dartel_normalization()
    % Applies DARTEL normalization to fMRIPrep outputs without smoothing.
    %
    % This script is designed to take preprocessed functional data from fMRIPrep,
    % apply a DARTEL template to normalize the images to MNI space, and save
    % the output with a 'w' prefix. It does not apply any smoothing.
    %
    % Key features:
    % - Robustly handles subjects with a variable number of runs.
    % - Skips subjects who already have the output 'w' files.
    % - Skips subjects who are missing the necessary DARTEL flowfield file.
    %
    % Before running:
    % 1. Make sure SPM is in your MATLAB path.
    % 2. Update the 'subs' list to include all the subject IDs you want to process.
    % 3. Ensure the 'base_dir' and 'template_path' are correct for your system.

    clear; clc;

    % --- USER SETTINGS --- %
    % Define subject IDs to be processed
    subs = [101:102]; % Example: [101:105, 201, 203]

    % Base directory containing the 'derivatives' folder
    base_dir = '/home/acclab/Desktop/axc/derivatives';

    % Full path to your DARTEL template file
    template_path = fullfile(base_dir, 'matlab/spm/DARTEL_templates/Full_MST/Template_AxCFull_6.nii');

    % --- END USER SETTINGS --- %

    % Initialize SPM
    spm('defaults', 'FMRI');
    spm_jobman('initcfg');

    % Process each subject in parallel
    parfor i = 1:length(subs)
        sub = subs(i);
        sub_id_str = sprintf('sub-%03d', sub);
        fprintf('\nStarting Subject: %s\n', sub_id_str);

        % Define subject-specific paths
        func_dir = fullfile(base_dir, 'fmriprep', sub_id_str, 'func');
        flowfield_path = fullfile(base_dir, 'fmriprep', sub_id_str, 'anat', ...
            sprintf('u_rc1%s_desc-preproc_T1w_Template_AxCFull.nii', sub_id_str));

        % --- Pre-flight checks --- %

        % Skip if flowfield is missing
        if ~exist(flowfield_path, 'file')
            fprintf('WARNING: Flowfield missing for %s, skipping...\n', sub_id_str);
            continue;
        end

        % Find all existing functional runs for this subject
        input_files_struct = dir(fullfile(func_dir, ...
            sprintf('%s_task-retrieval_run-*_space-MNI152NLin2009cAsym_desc-preproc_bold.nii', sub_id_str)));

        if isempty(input_files_struct)
            fprintf('WARNING: No input functional files found for %s, skipping...\n', sub_id_str);
            continue;
        end

        num_runs = length(input_files_struct);
        fprintf('Found %d runs for %s.\n', num_runs, sub_id_str);

        % Check if all output 'w' files already exist
        w_files_exist = true;
        for run = 1:num_runs
            output_file = fullfile(func_dir, sprintf('w%s_task-retrieval_run-%d_space-MNI152NLin2009cAsym_desc-preproc_bold.nii', sub_id_str, run));
            if ~exist(output_file, 'file')
                w_files_exist = false;
                break;
            end
        end

        if w_files_exist
             fprintf('Skipping %s (all w-prefixed files already exist)\n', sub_id_str);
             continue;
        end

        % --- Prepare and run batch --- %

        % Collate input file paths
        input_files = cell(num_runs, 1);
        for run = 1:num_runs
            input_files{run} = fullfile(func_dir, input_files_struct(run).name);
        end

        % Configure SPM batch for DARTEL normalization
        matlabbatch = [];
        matlabbatch{1}.spm.tools.dartel.mni_norm.template = {template_path};
        matlabbatch{1}.spm.tools.dartel.mni_norm.data.subj.flowfield = {flowfield_path};
        matlabbatch{1}.spm.tools.dartel.mni_norm.data.subj.images = input_files;

        % Apply DARTEL normalization without smoothing
        matlabbatch{1}.spm.tools.dartel.mni_norm.vox = [2.3 2.3 2.3];
        matlabbatch{1}.spm.tools.dartel.mni_norm.preserve = 0;
        matlabbatch{1}.spm.tools.dartel.mni_norm.fwhm = [0 0 0];  % No smoothing

        % Run processing
        try
            spm_jobman('run', matlabbatch);

            % Rename outputs to the desired 'w' prefix format
            for run = 1:num_runs
                % SPM appends 'w' to the *beginning* of the original filename
                spm_output_file = fullfile(func_dir, ['w' input_files_struct(run).name]);
                final_output_file = fullfile(func_dir, sprintf('w%s_task-retrieval_run-%d_space-MNI152NLin2009cAsym_desc-preproc_bold.nii', sub_id_str, run));

                if exist(spm_output_file, 'file')
                    % In this case, SPM's default output matches our desired output, so no move is needed.
                    % We just confirm it was created.
                    fprintf('Created: %s\n', spm_output_file);
                else
                    fprintf('WARNING: Expected output missing for %s, run-%d\n', sub_id_str, run);
                end
            end

            fprintf('Successfully processed %s\n', sub_id_str);

        catch ME
            fprintf('ERROR processing %s: %s\n', sub_id_str, ME.message);
        end
    end

    fprintf('\nProcessing complete!\n');
    fprintf('All new `w` files saved in their respective `fmriprep/sub-*/func/` directories.\n');
end
