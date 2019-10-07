% outputs raw data from .mat files in data folder for specific model run,
% i.e. mat/data/<model_prefix>/<nrn_model_ver>/<cell_model_name>/<run_name>
function [t,ap_times,ap_netids,vm_matrix,params,threshE] = loadRawData(data_folder)
ap_file = fullfile(data_folder,'APtimes.mat'); % net id of each event in ap_file (ap_netids)    
vm_file = fullfile(data_folder,'Vm.mat'); % membrane potentials at each time point    
params_file = fullfile(data_folder,'params.mat'); % sim parameters
threshE_file = fullfile(data_folder,'threshE.mat'); 
% Load AP data
ap_data = load(ap_file); 
ap_times = ap_data.ap_times; 
ap_netids = ap_data.ap_netids; 
% Load Vm data
vm_data = load(vm_file); 
t = vm_data.tvec;
vm_matrix = vm_data.vm_matrix; 
% Load params data
params_data = load(params_file); 
params = params_data.params; 
% Load threshE or deltaVm
if exist(threshE_file,'file') == 2
    threshE_data = load(threshE_file); 
    threshE = threshE_data.threshE; 
else
    threshE = nan; 
end

end