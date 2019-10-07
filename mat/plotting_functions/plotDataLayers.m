function p_all = plotDataLayers(layers,data,cell_model_names,cell_ids,varargin)
% Plots median value of simulation data (thresholds,time constants, etc.)
% on layer surfaces, shifts each layer by shift_dir   
in.shift_dir = [0,-35,0]; % default shift
in.plot_vertices = 1; % interpolate data onto vertices (cells at element centers)
in.norm_mode = 'none'; % no normalization
in = sl.in.processVarargin(in,varargin);
%% Plot single figure with all layer surfaces
num_layers = length(layers);        
data_layer = cell(num_layers,1);
for i = 1:num_layers
    if ~isempty(cell_ids{i})
        cell_model_names_i = cellModelNames(cell_ids{i}); % cell names in layer
        [~,~,thresh_inds] = intersect(cell_model_names_i,cell_model_names); % get indices of layer cells in threshEs    
        data_layer{i} = median(cell2mat(data(thresh_inds)),2);            
    else
        data_layer{i} = 1e6; % placeholder large number to allow UniformOutput in cellfun below, no cells in this layer    
    end
end
fprintf('Computed layer medians\n')
%% Normalize
if strcmp(in.norm_mode,'min') % divide by minimum threshold (min = 1)
    global_min = min(cellfun(@min,data_layer));
    data_layer = cellfun(@(x) x/global_min,data_layer,'UniformOutput',0);
    fprintf('Applied normalization: %s\n',in.norm_mode)
elseif strcmp(in.norm_mode,'min_layer') % divide by minimum threshold in each layer (min = 1)
    data_layer = cellfun(@(x) x/min(x),data_layer,'UniformOutput',0); 
    fprintf('Applied normalization: %s\n',in.norm_mode)
end
%% Get data points on vertices using weighted average
if in.plot_vertices    
    for i = 1:num_layers
        if length(cell_ids{i}) >= 1
            data_layer{i} = getVertexCData(data_layer{i},layers(i).surface,layers(i).cell_origins); 
        end
    end
    fprintf('Computed data on vertices\n'); 
end
%% Convert to vertex data
% if in.plot_vertices
%     data = getAllVertexCData(data,layers,cell_ids,cell_model_names);
% end
%%
figure;     
for i = 1:num_layers
    num_cells_in_layer = length(cell_ids{i});    
    if num_cells_in_layer >= 1                                    
        layer_surf = layers(i).surface;
        layer_surf.vertices = layer_surf.vertices + (i-1)*repmat(in.shift_dir,size(layer_surf.vertices,1),1);
        p = patch(layer_surf);
        p.FaceVertexCData = data_layer{i};
        if in.plot_vertices
            p.FaceColor = 'interp';
        else
            p.FaceColor = 'flat';
        end
        p.CDataMapping = 'scaled';
        p.EdgeColor = 'none';
        hold on;          
    else
        fprintf('No cells in layer\n');
    end
end            
camlight; lighting gouraud;
axis equal; axis tight; 
axis off; 
end