%% Setup - run this function once to set up data files necessary for simulations
% uses defaults from Aberra et al. 2019
function initializeCellData()
mat_dir = addPaths; 

%% Check that nrnmech.dll file exists, if not display a warning message
nrnmech_file = fullfile(mat_dir, '../nrn', 'nrnmech.dll');
if ~isfile(nrnmech_file)
    msg = 'nrniv.exe does not exist. Compile mechanisms and move nrniv.exe to nrn folder';
    f = warndlg(msg, 'modal');
    return
end

%% Generate cell data (based on neuron models placed in '../nrn/cells/'
writeCellModelTable(); % creates cell_data/cell_data.mat, called by cellModelNames
cell_model_names = cellModelNames; % outputs all cell_model_names
num_cells = length(cell_model_names); 
nrn_model_ver = 'maxH'; % adult, human, myelinated (see outputMorphParams.m)
for cell_id = 1:num_cells
    saveCellData(cell_id,nrn_model_ver);
end
% Generate cell data for L2/3 PC with linear axon (maxHlin)
cell_id_l23pcs = 6:10; 
nrn_model_ver_linax = 'maxHlin'; 
for cell_id = cell_id_l23pcs
    saveCellData(cell_id,nrn_model_ver_linax); 
end
%% Generate neuron populations (in output_data/layer_data/<nrn_model_ver>/)
layer_set_num = 1; 
nrn_pop_name = 'nrn_pop1';
nrn_model_ver = 'maxH'; 
cell_ids = {1:5;6:10;11:15;16:20;21:25}; 
% Load first neuron population 
NeuronPop_file = fullfile(mat_dir,'output_data','layer_data',nrn_model_ver,...
                            [nrn_pop_name '_' nrn_model_ver '.mat']); 
if exist(NeuronPop_file,'file')
    NeuronPop_data = load(NeuronPop_file); 
    NeuronPop1 = NeuronPop_data.NeuronPop; 
else
   % create first neuron population with random azimuthal orientations
    NeuronPop1 = generateNeuronPop(layer_set_num,nrn_pop_name,nrn_model_ver,cell_ids); 
end
phis1 = NeuronPop1.phis; 
%% Generate additional populations by rotating each neuron in azimuthal direction
% by fixed amount
phi_rot = 60:60:300; 
nrn_pops = cell(length(phi_rot)+1,1); 
nrn_pops{1} = nrn_pop_name; 
for i = 1:length(phi_rot)
    rot_i = phi_rot(i); 
    phisi = rotateNeuronPopPhis(phis1,rot_i);
    nrn_pop_name = sprintf('nrn_pop%g',i+1); 
    nrn_pops{i+1} = nrn_pop_name; 
    generateNeuronPop(layer_set_num,nrn_pop_name,nrn_model_ver,cell_ids,'phis',phisi);         
end
%% Generate population of L2/3 PCs with linear axon (maxHlin)
layer_set_num = 1; 
nrn_pop_name = 'nrn_pop1';
nrn_model_ver = 'maxHlin'; 
cell_ids = {[];cell_id_l23pcs;[];[];[]}; 
% Saves neuron population for linear axon L2/3 PCs
generateNeuronPop(layer_set_num,nrn_pop_name,nrn_model_ver,cell_ids);                         
end