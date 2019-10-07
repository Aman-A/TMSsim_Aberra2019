% Creates directory for data storage of given run in Data/ of current
% directory
% Used for MATLAB and NEURON data folders, assumes model directory is in a
% folder called NEURON
% MATLAB data organization
    % matlab_direc -> Data -> cell_name -> model_prefix -> model_name ->
    % .mat files
% NEURON data organization
    % nrn_direc -> tmp -> current .txt files
    % nrn_direc -> Data -> cell_name -> model_prefix -> model_name -> .txt
    % files
% NEW NEURON parallelized data organization
    % nrn_direc -> tmp -> tmp_<cell_name>_<model_name> -> current .txt files
    % no Data folder anymore    
function [data_folder_full,run_tmp_fold,run_params_fold] = createDataFolder(parent_direc,simulation,cell_model)    
cell_model_name = cell_model.cell_model_name;
model_prefix = simulation.model_prefix; % model name for name of tmp folder and MATLAB data folder
nrn_model_ver = cell_model.nrn_model_ver;    
model_name = simulation.model_name;   
run_name = simulation.run_name;     
run_tmp_fold = 'blah'; run_params_fold = 'blah'; data_folder_full = 'blah'; % overwrite if creating NEURON folders
% if parent_direc is nrn check if tmp directory is made already
% and make this run's tmp folder and Params folder
if ~isempty(regexp(parent_direc,'nrn','ONCE'))
% *** Make nrn params/data/cell_data folders ***
    % Create general tmp directory if necessary
    % tmp/<run_name>/
    if exist(fullfile(parent_direc,'tmp'),'dir') == 0  
        mkdir(fullfile(parent_direc,'tmp')); 
        fprintf('Created tmp directory in %s\n',parent_direc);
    else
%             fprintf('tmp directory already made in %s\n',parent_direc);
    end
    % Create current simulation's specific tmp directory in the general tmp        
    run_tmp_fold = fullfile(parent_direc,'tmp',run_name);
    if exist(run_tmp_fold,'dir') == 0
        mkdir(run_tmp_fold);
%             fprintf('Created run tmp directory\n');
    else
%             fprintf('Overwriting run tmp directory\n');            
%             rmdir(run_tmp_fold,'s');
%             mkdir(run_tmp_fold);            
    end
    % Create general Params directory if necessary
    % params/<run_name>/
    if exist(fullfile(parent_direc,'params'),'dir') == 0  
        mkdir(fullfile(parent_direc,'params')); 
%             fprintf('Created params directory in %s\n',parent_direc);                    
    end
    % Create current simulation's specific Params directory in the
    % general Params folder
    run_params_fold = fullfile(parent_direc,'params',run_name);
    if exist(run_params_fold,'dir') == 0
        mkdir(run_params_fold);
        fprintf('Created run params directory\n');
    else
%             fprintf('Overwriting run params directory\n');
%             rmdir(run_params_fold,'s');
%             mkdir(run_params_fold);            
    end         
%         fprintf('Created nrn folders\n');
else 
% *** Make mat data folders *** 
% nrn_sim_data/<model_prefix>/<nrn_model_ver>/<cell_model_name>/
    data_direc = fullfile(parent_direc,'nrn_sim_data');
    if exist(data_direc,'dir') == 0
        mkdir(data_direc);
        fprintf('Created data directory in %s\n',parent_direc);
    else
        %fprintf('Data directory already made in %s\n',parent_direc);
    end
    model_prefix_path = fullfile(data_direc,model_prefix); %directory containing all results of run
    if exist(model_prefix_path,'dir') == 0         
        mkdir(model_prefix_path);
        fprintf('Created %s directory\n',model_prefix_path);
    else % Data directory already made
        %fprintf('Writing results to %s\n',model_prefix);         
    end    
    cell_fold = fullfile(model_prefix_path,cell_model_name);
    if exist(cell_fold,'dir') == 0
        mkdir(cell_fold)
%             fprintf('Created cell model directory %s\n',cell_fold);
    else
        %fprintf('Cell model directory %s already made\n',cell_fold);
    end    
    data_folder_full = fullfile(cell_fold,model_name); % full path name of data folder, containing either  .mat files     
    if exist(data_folder_full,'dir') == 0
        mkdir(data_folder_full)
%             fprintf('Created data folder %s\n',data_folder_full);
    else
       %fprintf('Data folder %s already made\n',data_folder_full); 
    end
    % Create cell_data folder
    if exist(fullfile(parent_direc,'cell_data'),'dir') == 0  
        mkdir(fullfile(parent_direc,'cell_data'));
        fprintf('Created cell_data directory in %s\n',parent_direc);                    
    end
    nrn_model_ver_fold = fullfile(parent_direc,'cell_data',nrn_model_ver);
    if exist(nrn_model_ver_fold,'dir') == 0
       mkdir(nrn_model_ver_fold);
       mkdir([nrn_model_ver_fold filesep 'coordinates']); 
       mkdir([nrn_model_ver_fold filesep 'comp_types']); 
       mkdir([nrn_model_ver_fold filesep 'diams']); 
       mkdir([nrn_model_ver_fold filesep 'parent_inds']); 
       mkdir([nrn_model_ver_fold filesep 'secnames']); 
       mkdir([nrn_model_ver_fold filesep 'sectypes']); 
       mkdir([nrn_model_ver_fold filesep 'areas']); 
       mkdir([nrn_model_ver_fold filesep 'branchorders']); 
       fprintf('Created nrn_model_ver directories in cell_data/%s\n',nrn_model_ver);                    
    end 
%         fprintf('Created mat folders\n');
end            
end