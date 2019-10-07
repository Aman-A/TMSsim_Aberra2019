function writeParamsHoc(run_params_fold,simulation,cell_model,Efield,waveform)
if nargin == 0
   mat_dir = addPaths;
   run_params_fold = fullfile(mat_dir,'../nrn/params/test');  
end
file_name = 'defineParams.hoc';
fid = fopen(fullfile(run_params_fold,file_name),'w'); 
fprintf(fid,'NSTACK_size = 100000 // always call with NSTACK at least 100000\n'); 
%% Simulation parameters
fprintf(fid,'//*****Simulation Parameters****\n'); 
fprintf(fid,'steps_per_ms = 1/%f\n',simulation.dt);
sim_params = fieldnames(simulation);
% remove parameters not used by NEURON
sim_params = sim_params(~strcmp(sim_params,'model_prefix') & ~strcmp(sim_params,'model_name') & ~strcmp(sim_params,'model_type')); 
for i = 1:length(sim_params)    
    if isnumeric(simulation.(sim_params{i}))
        if abs(simulation.(sim_params{i})) > 1e6
            fprintf(fid,'%s = %g\n',sim_params{i},simulation.(sim_params{i}));
        else
            fprintf(fid,'%s = %f\n',sim_params{i},simulation.(sim_params{i}));
        end
    elseif ischar(simulation.(sim_params{i}))
        fprintf(fid,'strdef %s\n%s = "%s"\n',sim_params{i},sim_params{i},simulation.(sim_params{i}));
    else
        error('Parameter %s is not string or numeric',sim_params{i}); 
    end
end
%% Cell model parameters
fprintf(fid,'//*****Cell model Parameters****\n'); 
cell_params = fieldnames(cell_model); 
for i = 1:length(cell_params)    
    if isnumeric(cell_model.(cell_params{i}))
        fprintf(fid,'%s = %f\n',cell_params{i},cell_model.(cell_params{i}));     
    elseif ischar(cell_model.(cell_params{i}))
        fprintf(fid,'strdef %s\n%s = "%s"\n',cell_params{i},cell_params{i},cell_model.(cell_params{i}));
    else
        error('Parameter %s is not string or numeric',sim_params{i}); 
    end
end
%% E-field parameters
fprintf(fid,'//*****E-field Parameters****\n'); 
Efield_params = fieldnames(Efield); 
% remove parameters not used by NEURON (if load_potentials = 1)
Efield_params = Efield_params(~strcmp(Efield_params,'E_position') & ~strcmp(Efield_params,'E_file')); 
Efield_params = Efield_params(~strcmp(Efield_params,'layer_set_num') & ~strcmp(Efield_params,'nrn_pop_name')); 
for i = 1:length(Efield_params)
    fprintf(fid,'%s = %f\n',Efield_params{i},Efield.(Efield_params{i}));
end
%% Waveform parameters
% Set initial E-field stimulus amplitude (from waveform struct)
fprintf(fid,'AMP = %f\n',waveform.amp);
%%
fclose(fid); 
fprintf('Wrote defineParams.hoc to %s\n',run_params_fold);
end