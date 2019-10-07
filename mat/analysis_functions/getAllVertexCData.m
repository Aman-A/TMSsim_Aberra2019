function data_vertices = getAllVertexCData(data,layers,cell_ids,cell_model_names)      
% Input data calculated for cells placed at center of faces, layers
% structure, and cell_ids being plotted from plotThreshLayers, outputs new
% threshEs cell array with data converted to thresholds at each vertex in
% surface using weighted averages of face data
mat_dir = addPaths; 
data_vertices = cell(size(data)); 
num_layers = length(layers); 
for i = 1:num_layers
  num_cells_in_layer = length(cell_ids{i}); 
  for j = 1:num_cells_in_layer
       cell_model_name = cellModelNames(cell_ids{i}(j),mat_dir);            
       data_ind = strcmp(cell_model_name,cell_model_names);
       data_ij = data{data_ind}; % extract threshold for this cell
       % convert threshold face data to vertex data using weighted
       % averages
       data_ij = getVertexCData(data_ij,layers(i).surface,layers(i).cell_origins); 
       data_vertices{data_ind} = data_ij; 
  end
end
fprintf('Converted data from face data to vertex data\n'); 
end