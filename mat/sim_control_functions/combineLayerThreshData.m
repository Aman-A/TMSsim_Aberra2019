function combineLayerThreshData(model_prefix,cell_ids)
% COMBINELAYERTHRESHDATA Combines threshold data for each cell
%(<cell_model_name>_thresh.mat) into one data file in
%data/<model_prefix>/<nrn_model_ver>/<model_prefix>.mat
%
% AUTHOR    : Aman Aberra
if nargin == 0
   model_prefix = 'utms_maxH_w1';    
   cell_ids = 1:25;
end
par_on = 1; % use parfor loop
mat_dir = addPaths;
data_dir = fullfile(mat_dir,'nrn_sim_data',model_prefix); 
num_cells = length(cell_ids);
cell_model_names = cellModelNames(cell_ids); 
if num_cells == 1
   cell_model_names = {cell_model_names}; 
end
threshEs = cell(1,num_cells); init_inds = cell(1,num_cells); 
E_mags = cell(1,num_cells); 
incomplete_cellids = []; % cell ids of incomplete cells
incomplete_inds = []; % indices of incomplete cells
saved_params = 0; 
%% set up parpool
numCPUs = feature('numCores'); 
fprintf('Number of CPUs requested = %g\n',numCPUs);  
if numCPUs > 1 && par_on % parallelized
    %% Set up parpool
    if numCPUs > 4 % on cluster
        pc_storage_dir = fullfile(mat_dir,'pc_storage',getenv('SLURM_JOB_ID'));        
        mkdir(pc_storage_dir); 
        pc = parcluster('local');    
        pc.JobStorageLocation =  pc_storage_dir;    
    else
       pc = parcluster('local'); % use default 
    end      
    poolobj = parpool(pc,numCPUs);
    % initialize arrays    
    params = cell(1,num_cells);
    parfor i = 1:num_cells
        % load data
        data_filei = fullfile(data_dir,cell_model_names{i},sprintf('%s_thresh.mat',cell_model_names{i}));
        if exist(data_filei,'file') == 2
            datai = load(data_filei);                         
            params{i} = datai.params;
            threshEs{i} =  datai.threshEs; 
            init_inds{i} = datai.init_inds; 
            E_mags{i} = datai.E_mag; 
        else
           fprintf('cell %g: %s data not found\n',cell_ids(i),cell_model_names{i}); 
           incomplete_cellids = [cell_ids(i),incomplete_cellids];            
           incomplete_inds = [incomplete_inds,i]; 
        end
    end
    delete(poolobj); 
    rmdir(pc_storage_dir,'s');     
    % get params from first completed simulation
    not_empty_inds = find(~cellfun(@isempty,params));    
    params = params{not_empty_inds(1)};     
else % serial
    for i = 1:num_cells
       % load data
        data_filei = fullfile(data_dir,cell_model_names{i},sprintf('%s_thresh.mat',cell_model_names{i}));
        if exist(data_filei,'file') == 2
            datai = load(data_filei); 
            if saved_params == 0 % save once for completed simulation                
                params = datai.params;
                saved_params = 1; 
            end
            threshEs{i} =  datai.threshEs; 
            init_inds{i} = datai.init_inds;  
            E_mags{i} = datai.E_mag; 
        else            
            fprintf('cell %g: %s data not found\n',cell_ids(i),cell_model_names{i});
            incomplete_cellids = [cell_ids(i),incomplete_cellids];
            incomplete_inds = [incomplete_inds, i];
        end
    end
end
% fill in unfinished runs with nans
for i = 1:length(incomplete_inds)
   threshEs{incomplete_inds(i)} = nan;  
   init_inds{incomplete_inds(i)} = nan; 
end
% print incomplete runs
if ~isempty(incomplete_cellids)
    incomplete_cellids = sort(incomplete_cellids); 
    fprintf('Incomplete cells:\n'); 
    fprintf('%g ',incomplete_cellids);
    fprintf('\n'); 
end
% Save data
save_filename = [model_prefix '.mat'];
save(fullfile(mat_dir,'nrn_sim_data',save_filename),'cell_model_names','params','threshEs','init_inds','E_mags'); 
fprintf('Saved data to %s\n',save_filename); 