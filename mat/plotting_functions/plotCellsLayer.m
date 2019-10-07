function plotCellsLayer(cell_id,l_num,NeuronPop,inds,vals)
% Plots one neuron model population in its corresponding layer
if nargin < 5
   vals = 0;  % num_compartment x 1 vector
end
nrn_model_ver = NeuronPop.nrn_model_ver;
layers = NeuronPop.layers; 
cell_ids = NeuronPop.cell_ids;
% check that cell_id was placed in this layer
if ~ismember(cell_id,cell_ids{l_num})
   error('Cell %g was not placed in layer %g (%s)',cell_id,l_num,NeuronPop.name) 
end
cell_origins = layers(l_num).cell_origins;
cell_normals = layers(l_num).cell_normals;
phis = NeuronPop.phis{l_num}{cell_id==cell_ids{l_num}};
cell_data = loadCellData(cell_id,nrn_model_ver); 
cell_data.C = cell_data.C*1e-3; % convert to mm
if iscolumn(inds)
   inds = inds'; 
end
celli_opts.cell_data = cell_data; % same for every neuron
for i = inds
    celli_opts.cell_origin = cell_origins(i,:); 
    celli_opts.cell_normal = cell_normals(i,:); 
    celli_opts.phi = phis(i); 
    celli_opts.vals = vals; 
    plotCellLines(cell_id,nrn_model_ver,celli_opts);
end
drawnow;