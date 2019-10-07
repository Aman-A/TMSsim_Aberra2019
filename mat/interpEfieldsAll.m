function interpEfieldsAll()
addPaths;
% Runs interpEfield for both E-field solutions and all 6 neuron populations
% (6 azimuthal rotations for each neuron)
layer_set_num = 1; 
nrn_model_ver = 'maxH';
pop_inds = 1:6; % set to 1 to only interpolate E-fields on single population (1 rotation)
nrn_pops = arrayfun(@(x) sprintf('nrn_pop%g',x),pop_inds,'UniformOutput',0)';
over_write = 0; 
%% Interpolate E-fields - Efield for conventional TMS simulations 
% Warning takes ~2 hours per neural population if run serially 
Efield_solution='M1_PA_MCB70'; 
for i = 1:length(nrn_pops)
    interpEfield(layer_set_num,nrn_pops{i},nrn_model_ver,Efield_solution,over_write); 
end
% Reverse E-fields to get E-fields for AP stimulation
for i = 1:length(nrn_pops)
    reverseEfield(layer_set_num,nrn_pops{i},nrn_model_ver,Efield_solution); 
end
% Do same for L2/3 PC axon models
nrn_model_ver_linax = 'maxHlin'; 
interpEfield(layer_set_num,nrn_pops{1},nrn_model_ver_linax,Efield_solution,over_write); 
reverseEfield(layer_set_num,nrn_pops{1},nrn_model_ver,Efield_solution); 
%% Efield for cTMS simulations
Efield_solution='M1_PA_Magstim70mm';
for i = 1:length(nrn_pops)
    interpEfield(layer_set_num,nrn_pops{i},nrn_model_ver,Efield_solution,over_write); 
end
end