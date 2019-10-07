function layers_or_layersE = loadLayers(layer_set_num,opt,Efield_name)
% Load layers struct or layersE struct
%
mat_dir = addPaths;
if nargin < 1
   layer_set_num = 1;   
end
if nargin < 2
   opt = 0; % just cell placement layers (between layer boundaries)
end
if opt == 1 % layersP includes with both cell placement layer and boundaries
    layer_name = sprintf('layer_set_%gp.mat',layer_set_num);
elseif opt == 2 % layersE includes E-field data for Efield_name within layers
    layer_name = sprintf('layer_set_%g_E_%s.mat',layer_set_num,Efield_name);    
else
    layer_name = sprintf('layer_set_%g.mat',layer_set_num);
end
fprintf('Loading %s\n',layer_name);  
layer_folder = fullfile(mat_dir,'output_data','layer_data'); 
layer_data = load(fullfile(layer_folder,layer_name));
if opt < 2
    layers_or_layersE = layer_data.layers; 
else
    layers_or_layersE = layer_data.layersE; 
end