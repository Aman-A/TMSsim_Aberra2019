function ctmsTimeConstants(cTMS_durs,data_files,data_output_file,varargin)
%GETCTMSTIMECONSTANTS estimate time constants for all cells  
%
%   Inputs 
%   ------ 
%   cTMS_durs : vector
%               pulse-widths in ms, same order as data_files
%   data_files : cell array of strings
%               file names of simulated thresholds to use to estimate time
%               constants 
%   data_output_file : string
%                      name of output file to save time constants to
%   Optional Inputs 
%   --------------- 
%   data_folder : string
%                 name of folder to load data from and save output file to

% AUTHOR    : Aman Aberra 
if nargin == 0
   cTMS_durs = [30,60,120];
   data_files = {'tms_maxH_w4_30us_ls_1_E_M1_PA_Magstim70mm_P_nrn_pop1-nrn_pop6_all';
       'tms_maxH_w4_60us_ls_1_E_M1_PA_Magstim70mm_P_nrn_pop1-nrn_pop6_all';
       'tms_maxH_w4_120us_ls_1_E_M1_PA_Magstim70mm_P_nrn_pop1-nrn_pop6_all'};
   data_output_file = 'tms_maxH_w4_30-120us_ls_1_E_M1_PA_Magstim70mm_P_nrn_pop1-nrn_pop6_all_timeconstants';
end
par_on = 1; % 1 to use parfor loop
mat_dir = addPaths; 
in.data_folder = fullfile(mat_dir,'nrn_sim_data');
in = sl.in.processVarargin(in,varargin); 
%% Load data
num_durs = length(cTMS_durs);
% Load data
threshEs_all = cell(num_durs,1);
cell_model_names_all = cell(num_durs,1); 
for i = 1:num_durs
    data = load(fullfile(in.data_folder,[data_files{i} '.mat']));
    fprintf('Loaded %s\n',data_files{i}); 
    cell_model_names_all{i} = data.cell_model_names;
    threshEs_all{i} = data.threshEs;    
end
cell_model_names = cell_model_names_all{1}; % all are same, extract single array of names
num_cells = length(cell_model_names_all{1});
% array with num_cell elements, each with num_pos vectors with dimension
% num_durs x 1, i.e. str-dur curve for each position
threshEs_cell = cell(1,num_cells); 
for n = 1:num_cells
    num_positions = size(threshEs_all{1}{n},1);
    num_rotations = size(threshEs_all{1}{n},2);
    threshEs_cell{n} = repmat({zeros(num_durs,num_rotations)},num_positions,1);    
    for j = 1:num_durs
        for i = 1:num_positions
            % nth cell, jth duration, ith position (all rotations)
            threshEs_cell{n}{i}(j,:) = threshEs_all{j}{n}(i,:); 
        end        
    end    
end
% Load waveforms
% Load cTMS waveforms
cTMS_data = load(fullfile(mat_dir,'run_nrn_functions/cTMSwaves.mat')); 
t = cTMS_data.t; % time in sec
pw = cTMS_data.pw;
[~,pw_inds,~] = intersect(pw,cTMS_durs*1e-3);
W = cTMS_data.W(:,pw_inds); 
%% Calculate timeconstants
timeconstants = cell(1,num_cells);
rheobases = cell(1,num_cells);
residuals = cell(1,num_cells);
%% Set up parpool
numCPUs = feature('numCores'); 
fprintf('Number of CPUs requested = %g\n',numCPUs);    
if numCPUs > 1 && par_on
   if numCPUs > 4 % on cluster
       pc_storage_dir = fullfile(mat_dir,'pc_storage',getenv('SLURM_JOB_ID'));
       mkdir(pc_storage_dir);
       pc = parcluster('local');
       pc.JobStorageLocation =  pc_storage_dir;
   else
      pc = parcluster('local');  
   end
   poolobj = parpool(pc,numCPUs); 
   parfor n = 1:num_cells
       fprintf('Calculating time constants for %s\n',cell_model_names{n});
       num_positions = length(threshEs_cell{n}); % # of positions for this cell
       num_rotations = size(threshEs_cell{n}{1},2); 
       timeconstants{n} = zeros(num_positions,num_rotations);
       rheobases{n} = zeros(num_positions,num_rotations);
       residuals{n} = zeros(num_positions,num_rotations);
       for i = 1:num_positions 
           for j =  1:num_rotations
               mt = threshEs_cell{n}{i}(:,j);
               [tau,rb,resnorm] = estTimeConstant(cTMS_durs,mt,t,W);
               timeconstants{n}(i,j) = tau; % time constant (in us)
               rheobases{n}(i,j) = rb;
               residuals{n}(i,j) = resnorm; % r^2 of ith orientation
           end
%            fprintf('pos %g: tau = %.3f us. Rh = %.3f A/us\n',i,timeconstants{n}(i),rheobases{n}(i))
       end
   end
   delete(poolobj);
else % serial loop
    for n = 1:num_cells
        fprintf('Calculating time constants for %s\n',cell_model_names{n});
        num_positions = length(threshEs_cell{n}); % # of positions for this cell
        num_rotations = size(threshEs_cell{n}{1},2);
        timeconstants{n} = zeros(num_positions,num_rotations);
        rheobases{n} = zeros(num_positions,num_rotations);
        residuals{n} = zeros(num_positions,num_rotations);
        for i = 1:num_positions
            for j =  1:num_rotations
                mt = threshEs_cell{n}{i}(:,j);
                [tau,rb,resnorm] = estTimeConstant(cTMS_durs,mt,t,W);
                timeconstants{n}(i,j) = tau; % time constant (in us)
                rheobases{n}(i,j) = rb;
                residuals{n}(i,j) = resnorm; % r^2 of ith orientation
            end
%             fprintf('pos %g: tau = %.3f us. Rh = %.3f A/us\n',i,timeconstants{n}(i),rheobases{n}(i))
        end
    end
end
%% Save data
save_file = fullfile(in.data_folder,[data_output_file '.mat']); 
save(save_file,'cell_model_names','timeconstants','rheobases','residuals'); 
fprintf('Saved time constants to %s\n',data_output_file); 
