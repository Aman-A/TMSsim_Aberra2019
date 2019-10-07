function combineNrnPopData(model_prefix_pre,nrn_pop_names,varargin)
%COMBINENRNPOPDATA Combine neuron simulation data from multiple populations
%(i.e. azimuthal rotations) into single data file
%Assumes last part of model_prefix is name of neuron population, e.g.
% if model_prefix_pre = 'tms_sim', full model_prefix names would be:
% 'tms_sim_nrn_pop1', 'tms_sim_nrn_pop2', etc. 
%
%   Inputs 
%   ------ 
%   Optional Inputs 
%   --------------- 

%   Examples 
%   --------------- 

% AUTHOR    : Aman Aberra 
in.data_fields = {'threshEs','init_inds'}; % default is for threshold simulations
in = sl.in.processVarargin(in,varargin); 
%% Load data
mat_dir = addPaths; 
data_fold = fullfile(mat_dir,'nrn_sim_data'); 
num_pops = length(nrn_pop_names); 
data_all = cell(num_pops,1); 
for i = 1:num_pops
    data_filei = [model_prefix_pre '_P_' nrn_pop_names{i} '.mat'];
    data_all{i} = load(fullfile(data_fold,data_filei));
end
% convert cell array of structs to struct array
data_all = cell2mat(data_all); 
data_out.cell_model_names = data_all(1).cell_model_names; % same for all
data_out.E_mags = data_all(1).E_mags; % same for all
num_data_fields = length(in.data_fields); 
num_cells = length(data_out.cell_model_names); 
% combine data from separate simulation into arrays, grouped by
% cell_model_name
for i = 1:num_data_fields   
   data_alli = [data_all.(in.data_fields{i})]; % concatenate all data for this field
   % separate by cell model
   datai = cell(1,num_cells); 
   for j = 1:num_cells
       datai{j} = cell2mat(data_alli(j:num_cells:end));
   end
   data_out.(in.data_fields{i}) = datai; % 
end
%% Save data
data_out_file = fullfile(data_fold,[model_prefix_pre '_P_' nrn_pop_names{1} '-' nrn_pop_names{end} '_all.mat']);
save(data_out_file,'-struct','data_out'); 
