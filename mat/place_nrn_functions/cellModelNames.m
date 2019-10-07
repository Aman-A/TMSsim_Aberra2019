function cell_model_names = cellModelNames(cell_ids,mat_dir)
% CELLMODELNAMES opens table with names of neuron models
%input scalar or vector of cell_ids, outputs string or cell array,
%respectively, of cell_model_names
if nargin < 2
    mat_dir = addPaths; 
end
cell_model_data_file = fullfile(mat_dir,'cell_data','cell_models.mat'); 
cell_model_data = load(cell_model_data_file); 
T = cell_model_data.T; 
if nargin == 0
   cell_ids = 1:size(T,1); % output all cell ids  
end
[~,inds] = ismember(cell_ids,T.cell_id);
if ~isempty(inds) && ~any(inds == 0)
    cell_model_names = T.Properties.RowNames(inds);
    if length(cell_model_names) == 1
        cell_model_names = cell_model_names{1}; 
    end
else
   error('All cell_ids should be between %g and %g',T.cell_id(1),T.cell_id(end))
end
end