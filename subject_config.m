function [subject_taskInfo] = subject_config(subject_id, default_taskInfo)
    % subject_config Returns subject-specific modifications to the taskInfo struct.
    %
    % This function checks if a given subject ID has any exceptions to the
    % default parameters (e.g., a different number of runs). If an exception
    % is found, it returns a modified taskInfo struct for that subject.
    % Otherwise, it returns the original, unmodified taskInfo struct.
    %
    % INPUTS:
    %   subject_id       - The ID of the subject (e.g., 'sub-211').
    %   default_taskInfo - The default taskInfo struct loaded from the params file.
    %
    % OUTPUT:
    %   subject_taskInfo - The potentially modified taskInfo struct for the subject.

    % By default, the subject's info is the default info
    subject_taskInfo = default_taskInfo;

    % Use a switch statement to handle specific subjects with exceptions
    switch subject_id
        case 'sub-211'
            fprintf('Found custom configuration for %s: Overriding number of runs to 3.\n', subject_id);
            subject_taskInfo.Runs = 3;

        case 'sub-241'
            fprintf('Found custom configuration for %s: Overriding number of runs to 3.\n', subject_id);
            subject_taskInfo.Runs = 3;

        case 'sub-201'
            fprintf('Found custom configuration for %s: Overriding number of runs to 5.\n', subject_id);
            subject_taskInfo.Runs = 5;

        case 'sub-242'
            fprintf('Found custom configuration for %s: Overriding number of runs to 5.\n', subject_id);
            subject_taskInfo.Runs = 5;

        % Add other subjects with exceptions here
        % case 'sub-XXX'
        %     subject_taskInfo.Runs = Y;

        otherwise
            % For any other subject, the default taskInfo is used.
            % No changes are needed.
    end
end
