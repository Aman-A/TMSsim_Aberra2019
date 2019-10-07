function plotFig6a_bottom(save_fig)
% PLOTFIG2B Plots straight axon L2/3 population between GM/WM layer meshes
% for Fig. 6
%  
% AUTHOR    : Aman Aberra
if nargin == 0
   save_fig = 0;  
end
mat_dir = addPaths; 
layer_set_num = 1; 
nrn_pop_name = 'nrn_pop1';
nrn_model_ver = 'maxHlin';
cell_id = 6; % first cell in each layer
plot_layer = 2;
% Plot options
axlims = [-42 -30 -23 -10 33 46.8]; 
% Load data
layers = loadLayers(layer_set_num,0); % load layers 
% layersP = loadLayers(layer_set_num,1); % layersP
MeshROI_data = load(fullfile(mat_dir,'output_data','layer_data','MeshROI.mat')); 
MeshROI = MeshROI_data.MeshROI; 
gm = MeshROI.surfaces.GrayMatter; % gray matter mesh
wm = MeshROI.surfaces.WhiteMatter; % white matter mesh
NeuronPop = loadNeuronPop(nrn_pop_name,nrn_model_ver); 
%% Plot neuron populations
fig = figure('Units','Normalized','Position',[0.2639 0.0689 0.7354 0.8244]); 
patch(gm,'FaceColor',[0.8 0.8 0.8],'EdgeColor','none','FaceAlpha',0.5,...
            'FaceLighting','gouraud'); hold on; 
patch(wm,'FaceColor',[0.8 0.8 0.8],'EdgeColor','none','FaceAlpha',1,...
            'FaceLighting','gouraud'); 
ax = gca;
axis(ax,'equal'); axis(ax,'off');
ax.View = [-82.27 53.2]; 
camlight('head');
axis(ax,axlims); 
% Plot cell populations
[~,inds] = clipPoints3d(layers(plot_layer).cell_origins,axlims);
plotCellsLayer(cell_id,plot_layer,NeuronPop,inds,[])
if save_fig
   fig_name = 'Fig6a_bottom'; 
   fig_fold = fullfile(mat_dir,'figures');    
   print(fig,fullfile(fig_fold,[fig_name '.tiff']),'-dtiff','-cmyk','-r300');    
end
end