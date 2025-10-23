%% RSA Searchlight

% This script was created by Cat Carpenter - contact cmc84@psu.edu for any
% questions
% Last updated 7/25/23

% This script implements the cosmoMVPA toolbox and uses a dissimilarity
% matrix within the searchlight

%This is multi-FRAME compatible as of 7/26/23

% This script is currently set up to run a binary similarity searchlight, 
% wherein the similarity between conditions are weighted to 1, 
% and similarity within conditions is weighted to  0
% from here, similarity is computed on each voxel within the specified searchlight.

%%  demo from http://www.cosmomvpa.org/_static/publish/demo_fmri_searchlight_rsm.html 

%path to cosmo-mvpa

%% Set data paths
% DEB: Define which conditions to compare for this specific analysis.
% The names must exactly match the labels in your behavioral files.
conditions_to_compare = {'cr_new', 'hit_same', 'fa_sim', 'cr_sim'};

% The function cosmo_config() returns a struct containing paths to tutorial
% data. (Alternatively the paths can be set manually without using
% cosmo_config.)
config=cosmo_config();

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
study_path = directory.Model;

% Base output directory name
analysis = [directory.Analysis filesep 'models' filesep file(1:end-4)];

% reset citation list
cosmo_check_external('-tic');

%% Load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iteration = 1:length(subjects)

  % define data filenames & load data   
  data_path=fullfile(study_path, subjects{iteration}); 
  output_path = fullfile(analysis, subjects{iteration});  

  % DEB: Path directly to the user-specified mask file.
  mask_path = '/home/acclab/Desktop/axc/multivariate/whole_aal_mask_resliced.nii';
  [mask_dir, mask_name, mask_ext] = fileparts(mask_path);
  masks = struct('folder', mask_dir, 'name', [mask_name, mask_ext]);

  data_fn=[data_path '/SPM_gz.mat']; %For fMRIprep  

  % create the output path if it doesn't already exist
    if ~exist(output_path, 'dir')
        mkdir(output_path)
    end

for curMask = 1:length(masks)

curROI = fullfile(masks(curMask).folder, masks(curMask).name); 
regionName=erase(masks(curMask).name,'.nii'); 

% Set conditions
Conds = taskInfo.Conditions;

ds=cosmo_fmri_dataset([data_fn ':beta'],'mask',curROI);

%ds=cosmo_fmri_dataset([data_fn ':beta'],'mask',curROI,'targets', 1:144, 'chunks', 1);

% Clear errant Not a Number (NaN) values
    % Remove constant features
   ds=cosmo_remove_useless_data(ds);  

%% Simple RSM searchlight
% define 'simple' target structure where trials from the same condition are
% equal to 0. (distance = 0)
% all other pairs are equally dissimilar (distance=1).
curROI

% DEB: Use the user-defined conditions_to_compare to filter the dataset
is_target_condition = false(size(ds.sa.labels));
for k = 1:numel(conditions_to_compare)
    is_target_condition = is_target_condition | ~cellfun(@isempty, strfind(ds.sa.labels, conditions_to_compare{k}));
end

% Keep only the samples that match the target conditions
ds = cosmo_slice(ds, is_target_condition);

% Rebuild the CondList based on the filtered dataset
CondList = zeros(size(ds.samples,1),1);
for ii=1:length(conditions_to_compare)
    Condition(ii).labels = ~cellfun(@isempty, strfind...
        (ds.sa.labels, conditions_to_compare{ii}));
    Condition(ii).idx = find(Condition(ii).labels == 1);
    CondList(Condition(ii).idx) = ii;
end

ds.sa.targets = CondList;
   
nsamples=size(ds.samples,1);
target_dsm=zeros(nsamples);
target_class=ds.sa.targets;

%% Linear RSA Searchlight

% simple sanity check to ensure all attributes are set properly
cosmo_check_dataset(ds);

% print dataset
fprintf('Dataset input:\n');
cosmo_disp(ds);

%% Define feature neighorhoods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Important note! 
       %%%Your searchlight size will be
       %%%influenced by your voxel dimensions. So if your voxel is
       %%%3mm isotropic and your searchlight radius is 2 voxels, your
       %%%searchlight radius will be 6mm. But if you have 2mm isotropic
       %%%voxels, then your radius would be 4mm. Searchlight
       %%%size will need to be study specific

%   'radius', r         } either use a radius of r, or select
%   'count', c          } approximately c voxels per searchlight
%                       Notes:
%                       - These two options are mutually exclusive
%                       - When using this option for an fmri dataset, the
%                         radius r is expressed in voxel units; for an meeg
%                         source dataset, the radius r is in whatever units
%                         the source dataset uses for the positions
                    
            %can determine searchlight size by count or radius in voxels,
            %radius is more common
            
            %searchlightSize=1.5;  %Average 17.5 voxels per light
            %searchlightSize=2.5;  %Average 70.9 voxels per light
            %searchlightSize=3;  %Average 120.3 voxels per light
            
  searchlightSize = searchlight.Size;
  searchlightMetric = searchlight.Metric; 
  
% The neighborhood defined here is used three times (one for each target
% similarity matrix), so it is not recomputed for every searchlight call.
fprintf('Defining neighborhood for each feature\n');
%nbrhood=cosmo_spherical_neighborhood(ds,'radius',searchlightSize);

nbrhood=cosmo_spherical_neighborhood(ds,searchlightMetric,searchlightSize);

% print neighborhood

fprintf('Searchlight neighborhood definition:\n');
cosmo_disp(nbrhood);

%% For linear RSA Searchlight
%computes absolute difference between each pair of samples
%target_dsm=abs(bsxfun(@minus,target_class,target_class'));

%% For binary RSA Searchlight 
for row=1:nsamples
    for col=1:nsamples
        same_target_class=target_class(row)==target_class(col);

      if same_target_class
            target_dsm(row,col)= 0; %% weighted as dissimilar
        else
            target_dsm(row,col)= 1; %% weighted as similar
        end
    end
end

fprintf('Using the following target dsm\n');
disp(target_dsm);
imagesc(target_dsm)
set(gca,'XTick',1:nsamples,'XTickLabel',ds.sa.labels,...
        'YTick',1:nsamples,'YTickLabel',ds.sa.labels)

% set measure
measure=@cosmo_target_dsm_corr_measure; 
measure_args=struct();
measure_args.target_dsm=target_dsm;

% print measure and arguments
fprintf('Searchlight measure:\n');
cosmo_disp(measure);
fprintf('Searchlight measure arguments:\n');
cosmo_disp(measure_args);     

% run searchlight
ds_rsm=cosmo_searchlight(ds,nbrhood,measure,measure_args);

% Note: when these results are used for group analysis across multiple
% participants, it may be good to Fisher-transform the correlation values,
% so that they are more normally distributed. This can be done by:
ds_rsm.samples=atanh(ds_rsm.samples);

% show results
cosmo_plot_slices(ds_rsm);

% store results
output_fn=strcat(output_path,'/', 'rsm','_',...
                                    num2str(searchlightSize),'_',...
                                    regionName,'.nii');
cosmo_map2fmri(ds_rsm, output_fn)

save([sprintf([output_path '/', 'rsm','_',...
                                    num2str(searchlightSize),'_',...
                                    regionName,'.mat'])]);
 end

end