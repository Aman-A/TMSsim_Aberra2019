%% Run single neuron simulation
% set neuron model parameters in nrn/params/test/defineParams.hoc
% Stimulation waveform parameters (tvec.txt, Evec.txt)
mat_dir = addPaths; 
dt = 0.005; 
tstop = 1;
waveform.mode = 1; % 1. monophasic ,2. biphasic,3. half sine MagProX100  
                   % 4 - cTMS1 waveforms (30 µs - 160 µs) 
waveform.del = dt; 
waveform.dur = 0.03; 
plot_wave = 1; 
[tvec,Evec] = TMSwave(dt,tstop,waveform,plot_wave); 
%% Write Er.txt
% E-field parameters (for load_potentials = 1)
cell_model.cell_id = 1; % L23_PC_cADpyr229_1
cell_model.nrn_model_ver = 'maxH';
Efield.E_position = 1133; % 1500
Efield.nrn_pop_name = 'nrn_pop1';
Efield.layer_set_num = 1;
Efield.E_file = 'M1_PA_MCB70';
nrn_dir = [mat_dir '/../nrn']; 
run_params_fold = fullfile(nrn_dir,'params','test'); 
Er = writeEr(mat_dir,run_params_fold,cell_model,Efield); 
%% Plot Er on morphology
cell_data = loadCellData(cell_model.cell_id,cell_model.nrn_model_ver);
figure; 
plotCellLines(cell_model.cell_id,cell_model.nrn_model_ver,'vals',vmag(Er),'lw',1);
axis equal; axis tight; 
view([0 0]); 
colorbar; 
title(sprintf('E magnitude on cell: %s at position: %g',...
    cellModelNames(cell_model.cell_id),Efield.E_position),'Interpreter','none'); 