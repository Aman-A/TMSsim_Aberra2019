function [Call,cell_data,axon_inds] = plotCellsLayerAxon(cell_id,l_num,NeuronPop,inds,vals)
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
cell_data = loadCellData(cell_id,nrn_model_ver,'coordinates','comp_types','parent_inds'); 
cell_data.C = cell_data.C*1e-3; % convert to mm
if iscolumn(inds)
   inds = inds'; 
end
% extract just axon
axon_inds = find(cell_data.comp_types == 1 | cell_data.comp_types == 2 | cell_data.comp_types == 3 | cell_data.comp_types == 4);
cell_data.C = cell_data.C([1;axon_inds],:); % soma (1) and axon indices
cell_data.comp_types = [0;cell_data.comp_types(axon_inds)];
parent_inds_ax0 = cell_data.parent_inds(axon_inds-1); 
parent_inds_ax = ones(length(axon_inds),1);
for j = 2:length(parent_inds_ax)
    parent_inds_ax(j) = find(parent_inds_ax0(j) == axon_inds)+1;
end
cell_data.parent_inds = parent_inds_ax; 
celli_opts.cell_data = cell_data; % same for every neuron
celli_opts.lw = 0.5; 
Call = cell(length(inds),1); 
for i = 1:length(inds)
    ii = inds(i); 
    celli_opts.cell_origin = cell_origins(ii,:); 
    celli_opts.cell_normal = cell_normals(ii,:); 
    celli_opts.phi = phis(ii); 
    celli_opts.vals = vals; 
    Call{i} = plotCellLines(cell_id,nrn_model_ver,celli_opts);
end
drawnow;