function NeuronPop = generateNeuronPop(layer_set_num,nrn_pop_name,nrn_model_ver,cell_ids,varargin)
%generateNeuronPop Places and orients neuron population in layers   
%
%   Inputs
%   ------
%   layer_set_num : integer 
%                   specify layer set to populate with neurons                       
%   nrn_pop_name : string
%                  name of neuron population (set of neuron models in
%                  layer and their azimuthal orientations)
%   nrn_model_ver : string
%                   name corresponding to set of neuron model parameters,
%                   specified in outputMorphParams
%   cell_ids : cell array
%              each element contains vector of cell ids within that layer
%
%   Optional Inputs
%   ---------------
%   phis : cell array (default empty)
%          cell array containing vectors of azimuthal rotations to apply to
%          each neuron in the population (num_layers x 1)
%   over_write : 1 or 0 (default 0)
%                determines whether existing population is overwritten
%   Outputs
%   -------
%   Examples
%   ---------- 
% 
% AUTHOR    : Aman Aberra
in.phis = {};
in.over_write = 0; 
in = sl.in.processVarargin(in,varargin); 
%% Load files
% load MeshROI
mat_dir = addPaths;
layer_data_folder = fullfile('output_data','layer_data');
MeshROI_file = load(fullfile(mat_dir,layer_data_folder,'MeshROI.mat'));
MeshROI = MeshROI_file.MeshROI;
% load layers
layer_set_name = sprintf('layer_set_%g',layer_set_num);
layers_file = load(fullfile(mat_dir,layer_data_folder,sprintf('%s.mat',layer_set_name)));
layers = layers_file.layers;
%% Create phis
num_layers = length(layers);     
if length(cell_ids) ~= num_layers
    error('Number of cell_id entries does not match number of layers');
end
if isempty(in.phis)
    phis = cell(num_layers,1);
    new_phis = 1;
else
    new_phis = 0; 
end
if new_phis % if placing cells with random phis        
    for i = 1:num_layers   
        layer_cells = cell_ids{i}; 
        num_layer_cells = length(layer_cells);
        % Extract cell positions according to cells_per_layer
        num_cell_positionsi = layers(i).num_elem;            
        for j = 1:num_layer_cells                        
            phis{i}{j} = rand(num_cell_positionsi,1)*360; % random angle between 0 and 360
        end    
    end
else
    phis = in.phis; 
end
%% Create NeuronPop struct
NeuronPop.layers = layers;
NeuronPop.MeshROI = MeshROI;
NeuronPop.cell_ids = cell_ids;
NeuronPop.phis = phis;
NeuronPop.nrn_model_ver = nrn_model_ver;
NeuronPop.name = nrn_pop_name;
nrn_pop_folder = fullfile(layer_data_folder,nrn_model_ver);
if exist(fullfile(mat_dir,nrn_pop_folder),'dir') == 0
   mkdir(fullfile(mat_dir,nrn_pop_folder));  
end
nrn_pop_file = fullfile(nrn_pop_folder,[nrn_pop_name '_' nrn_model_ver '.mat']);
if ~exist(fullfile(mat_dir,nrn_pop_file),'file') || in.over_write
   save(fullfile(mat_dir,nrn_pop_file),'NeuronPop');
   fprintf('Saved %s NeuronPop to %s\n',nrn_pop_name,nrn_pop_folder);
else
   fprintf('%s already exists, over_write = %g\n',nrn_pop_file,in.over_write);  
end
end