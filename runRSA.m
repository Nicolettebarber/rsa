%% Run Representational Similarity Analysis (RSA)
%   Editor:    CAN Lab
%   Updated:   4/9/2025
%
% Representational similarity analysis (RSA) for a single subject for single subject.
% Flagged for either single ROI RSA or searchlight analysis, but only
% currently supports ROI RSA. For RSA Searchlight, see runRSASearchlight.m
%
% Load single-trial beta images from each subject, apply ROI mask, calculate
% correlations between trial patterns, take the mean across trial types
%
%
% Updates:
%
% 3/28/19 - Loads in parameter file created by createParams.m subscript. Mat
% file should contain paths to roi and data folders, lists of rois,
% subjects, & conditions, and analysis name. See createParams.m for full
% list.

% 4/9/25 - removed all mention of ERS -- separate script now 

%Outputs within-category similarity, between-category similarity and
%distinctiveness described by Haxby 2001

% Only works for 2 conditions

% Requires zipped betas in nifti format and co-registered ROIs 

%% Set Analysis Parameters & Paths
% The conditions to be used in the analysis are now loaded from the
% params file, e.g. taskInfo.Conditions

% Load all relevent project information
if exist('commandFlag','var') == 0
    
    %Select parameter file is flag does not exist
    [file,path]=uigetfile('*.mat','Select params file');
    filename=fullfile(path,file);
    load(filename);
    
end

% turn cosmo warnings off
cosmo_warning('off');

% Filepath for results folder
parentDir = directory.Model;

% Base output directory name
analysis = [directory.Analysis filesep 'models' filesep file(1:end-4)];

%Debug
%subjects(2)=[];

%% Main Body
for iteration=1:length(subjects)
    %% Subject-Specific Directories
    
    % Current subject data paths:
    %  dataPath = fullpath to this subject's Single Trial Model directory
    %  spmFile  = fullpath to this subject
    %%%% For Just model folders w/Masks in's SPM.mat file. Note: the
    %                 :beta appended to the end tells cosmo to pull the beta
    %                 information from the SPM.mat file.
    
    dataPath   = fullfile(directory.Model, subjects{iteration});
    outputPath = fullfile(analysis, subjects{iteration});   
    spmFile = [dataPath '/SPM_gz.mat']; %For fMRIprep   
    
    % create the output path if it doesn't already exist
    if ~exist(outputPath, 'dir')
        mkdir(outputPath)
    end
    
    %% Load Mask Data
    % DEB: Path directly to the user-specified mask file.
    mask_path = '/home/acclab/Desktop/axc/multivariate/whole_aal_mask_resliced.nii';
    [mask_dir, mask_name, mask_ext] = fileparts(mask_path);
    masks = struct('folder', mask_dir, 'name', [mask_name, mask_ext]);
          
    for curMask = 1:length(masks)
        
        switch classType
            case 'RSA'
                % Path to current region mask
                curROI = fullfile(masks(curMask).folder, masks(curMask).name);
                
                % Current region name
                regionName = erase(masks(curMask).name, '.nii');
                
                %%% Loading ROI data into CosmoMVPA - use SPM betas
                % Note that loading data through the SPM.mat file will automatically
                % "chunk" by runs, which is what we want
                fprintf('Current subject: %s\n',subjects{iteration});
                fprintf('Loading data from ROI: %s\n',masks(curMask).name);
                currDataset=cosmo_fmri_dataset([spmFile ':beta'],'mask',curROI);
                

                % Clear errant Not a Number (NaN) values
                % Remove constant features
                currDataset=cosmo_remove_useless_data(currDataset);
                
                   
                % Trial by Trial Correlation Matrix - All Conditions --
                % Script Condition: Working (CAT CARPENTER 5/9/2025)
                %Flip Dataset
                    currDataset_flipped = currDataset.samples.';

                %Run Correlations
                    EntireMatrix=(atanh(corrcoef(currDataset_flipped)));
                        
                    %%%% You can save this out if you'd like too! %%%%
                        %Save out correlation matrices by subject and ROI in csv:
                        %csvwrite([analysis filesep subjects{iteration} filesep subjects{iteration} '_' regionName '_EntireCorrMatrix.csv'],
                        %EntireMatrix);

                %Plot 
                     % figure;
                     % imagesc(atanh(corrcoef(currDataset_flipped))); 
                     % colorbar; 
                     % colormap(jet); 
                     % caxis([-1 1]);
                     % numVars = size(currDataset_flipped, 2);  
                     % set(gca, 'XTick', 1:numVars, 'YTick', 1:numVars);
                     % set(gca, 'XTickLabel', currDataset.sa.labels, 'YTickLabel', currDataset.sa.labels);


                switch regressRT.flag
                    case 'Yes'
                        files = dir([study_path filesep subject filesep 'Run*']);
                        for i=1:length(files)
                            curMat(i) = load([files(i).folder filesep files(i).name]);
                            if i==length(files)
                                rtCell = [curMat.RT];
                                
                                % Convert from cell to double for regression
                                for ii=1:length(rtCell)
                                    
                                    % Flag outlier RT greater than 4 seconds
                                    if double(rtCell{ii}) >= regressRT.trialSec
                                        rtDouble(ii,1) = regressRT.trialSec;
                                    else
                                        rtDouble(ii,1) = double(rtCell{ii});
                                    end
                                end
                                
                                % Replace with trial duration (set in params)
                                rtDouble(isnan(rtDouble))=regressRT.trialSec;
                            end
                        end
                        
                        for i=1:length(currDataset.samples)
                            model = LinearModel.fit(rtDouble,currDataset.samples(:,i));
                            if i==1
                                allResiduals = model.Residuals.Raw;
                            else
                                allResiduals = [allResiduals model.Residuals.Raw];
                            end
                        end
                        
                        zscoreResid = zscore(allResiduals);
                        currDataset.samples = zscoreResid;
                        
                        clear files curMat rtCell rtDouble model allResiduals zscoreResid;
                end

        end
        
        %% Define trial information
        switch classType
            case 'RSA'
                try
                    % Find all trials that match any of the conditions to compare
                    is_target_condition = false(size(currDataset.sa.labels));
                    for k = 1:numel(taskInfo.Conditions)
                        is_target_condition = is_target_condition | ~cellfun(@isempty, strfind(currDataset.sa.labels, taskInfo.Conditions{k}));
                    end

                    % Keep only the samples that match the target conditions
                    currDataset = cosmo_slice(currDataset, is_target_condition);

                    % Rebuild the CondList and Cond structures based on the filtered dataset
                    CondList = zeros(size(currDataset.samples,1),1);
                    for ii=1:length(taskInfo.Conditions)
                        Cond(ii).labels = ~cellfun(@isempty, strfind...
                            (currDataset.sa.labels, taskInfo.Conditions{ii}));
                        Cond(ii).idx = find(Cond(ii).labels == 1);
                        CondList(Cond(ii).idx) = ii;
                    end
                    
                    currDataset.sa.targets = CondList;
                    
                    % Codes trials/conditions of no interest as 0 (see SpecifyModel script
                    % for trial tag information)
                    Zeroidx = find(CondList == 0);
                    
                    % Removes all trials of no interest from analysis
                    if isempty(Zeroidx)==0
                        currDataset.samples(Zeroidx,:)=[];
                        currDataset.sa.targets(Zeroidx)=[];
                        currDataset.sa.beta_index(Zeroidx)=[];
                        currDataset.sa.chunks(Zeroidx)=[];
                        currDataset.sa.fname(Zeroidx)=[];
                        currDataset.sa.labels(Zeroidx)=[];
                    end
                    
                    fprintf('Number of possible targets: %i\n',length...
                        (unique(currDataset.sa.targets)));
                    
                    % Print dataset
                    fprintf('Dataset input:\n');
                    
                    % Print dataset
                    fprintf('Number of samples: %i\n',...
                        size(currDataset.samples,1));
                    fprintf('Number of features (voxels): %i\n',...
                        size(currDataset.samples,2));
                    fprintf('Number of chunks (runs): %i\n',...
                        length(unique(currDataset.sa.chunks)));
                    
                    % RSA ROI analysis
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    switch analysisType
                        %case 'Searchlight'  
                        case 'ROI'
                            
                            % Calculate separate matrices for each condition & compare
                            % those
                            
                            % Split dataset into separate condition
                            % variables
                            for i=1:length(Cond)
                                index = find(currDataset.sa.targets == i);
                                Conditions(i) = currDataset;
                                Conditions(i).samples = Conditions(i).samples(index,:);
                                Conditions(i).sa.beta_index = Conditions(i).sa.beta_index(index);
                                Conditions(i).sa.chunks = Conditions(i).sa.chunks(index);
                                Conditions(i).sa.fname = Conditions(i).sa.fname(index);
                                Conditions(i).sa.labels = Conditions(i).sa.labels(index);
                                
                                % Re-index trials within condition into
                                % separate targets
                                Conditions(i).sa.targets = [1:length(index)]';
                                Conditions(i).samples = Conditions(i).samples';
                               
                            end
                    end
                    
                    % Compute RDM
                    rdm = zeros(length(Cond));
                    for i = 1:length(Cond)
                        for j = 1:length(Cond)
                            % Correlate the mean pattern for each condition
                            pattern_i = mean(Conditions(i).samples, 2);
                            pattern_j = mean(Conditions(j).samples, 2);
                            rdm(i,j) = 1 - corr(pattern_i, pattern_j);
                        end
                    end
                    
                    % Save the RDM
                    rdm_filename = fullfile(outputPath, sprintf('%s_rdm.csv', regionName));
                    csvwrite(rdm_filename, rdm);

                    clear targetDSM Conditions;
                end
             
        end
        
