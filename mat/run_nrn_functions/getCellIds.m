% outputs cell_ids and cell_model_names for given class of cell types
% type_list - cell array of cell classifications to output
% query_type can be 'layer','mtype', 'layer_mtype','etype'
% cell_models_or_data_file is optional input, defaults to 25 cells included in
% Aberra et al 2019. Input name of file or cell table itself
% ex: type_list = {'L1','L4','L6'}; query_type = {'layer'}; 
function [cell_ids,cell_model_names] = getCellIds(type_list,query_type,cell_models_or_data_file)
if nargin < 3
   cell_models_or_data_file = 'cell_models';  % load
end
if ischar(cell_models_or_data_file)
    mat_dir = addPaths; 
    cell_model_data = load(fullfile(mat_dir,'cell_data',[cell_models_or_data_file '.mat'])); 
    T = cell_model_data.T;
else
    T = cell_models_or_data_file; % input cell table
end
var_names = T.Properties.VariableNames;
if nargin < 2
   query_type = 'Row'; % default input is 'cell_model_name'   
end
if ~ismember(query_type,var_names) && ~strcmp(query_type,'Row')
   error('%s is not a valid query_type\n',query_type);  
end
% initialize vectors
cell_ids = []; 
cell_model_names = {}; 
if ischar(type_list)
   type_list = {type_list};  
end
for i = 1:length(type_list)
    [~,inds] = ismember(T.(query_type),type_list{i});
    inds = logical(inds); 
    cell_ids = [cell_ids;T.cell_id(inds)]; % append to list
    cell_model_names = [cell_model_names;T.Row(inds)];
end
end