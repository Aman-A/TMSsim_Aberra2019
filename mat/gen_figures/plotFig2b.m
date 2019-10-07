function plotFig2b(save_fig)
% PLOTFIG2B Plots neuron populations between GM/WM layer meshes for Fig2b
%  
% REQUIRES 'CURVE FITTING' MATLAB TOOLBOX FOR MAKESLICES
% AUTHOR    : Aman Aberra
if nargin == 0
   save_fig = 0;  
end
mat_dir = addPaths; 
layer_set_num = 1; 
nrn_pop_name = 'nrn_pop1';
nrn_model_ver = 'maxH';
cell_ids = 1:5:25; % first cell in each layer
% slice options
plane = [-35.3816  -17.0009 44.3816 -0.5 0.5 0]; % Slice analysis plane
uvROI = [0 0.6 0 1]; % bounds for 45� plane
normROI = [-40 -35.4 0 0 33.93 45]; % posterior part of sulcus plus little lip
start_dim = 2; % + y direction
% Plot options
axlims = [-42 -30 -23 -10 33 46.8]; 
cmap = [0 0.5 1;0 1 1; 0.5 1 0.5; 0.5961 0.3059 0.6392; 1 0.5 0];
caxlims =[0.4 5.5]; 
% Load data
layers = loadLayers(layer_set_num,0); % load layers 
layersP = loadLayers(layer_set_num,1); % layersP
MeshROI_data = load(fullfile(mat_dir,'output_data','layer_data','MeshROI.mat')); 
MeshROI = MeshROI_data.MeshROI; 
gm = MeshROI.surfaces.GrayMatter; % gray matter mesh
wm = MeshROI.surfaces.WhiteMatter; % white matter mesh
NeuronPop = loadNeuronPop(nrn_pop_name,nrn_model_ver); 
phis = NeuronPop.phis;
%% Plot neuron populations
fig = figure('Units','Normalized','Position',[0.2639 0.0689 0.7354 0.8244]); 
pGM = patch(gm,'FaceColor',[0.8 0.8 0.8],'EdgeColor','none','FaceAlpha',0.5,...
            'FaceLighting','gouraud'); hold on; 
pWM = patch(wm,'FaceColor',[0.8 0.8 0.8],'EdgeColor','none','FaceAlpha',1,...
            'FaceLighting','gouraud'); 
ax = gca;
axis(ax,'equal'); axis(ax,'off');
ax.View = [-82.27 53.2]; 
camlight head; 
axis(ax,axlims); 
% Plot cell populations
num_layers = length(layers); % plot one cell pop per layer
for i = 1:num_layers
    celli = cell_ids(i);
    [~,inds] = clipPoints3d(layers(i).cell_origins,axlims);
    val = ((i-1)+0.8); % same color for all compartments
     plotCellsLayer(celli,i,NeuronPop,inds,val)
end
colormap(ax,cmap); 
caxis(ax,caxlims); 
% plot plane and layer boundaries in plane
p3d = drawPlane3d(createPlane(plane(1:3),plane(4:6))); 
p3d.FaceColor = 'r'; p3d.FaceAlpha = 0.5; p3d.EdgeColor = 'k'; p3d.LineWidth = 1; 
slices = makeSlices(layersP,MeshROI,plane,start_dim,0,uvROI); 
lgm = slices(1).slice_line;
lwm = slices(end).slice_line;
% plot line denoting gm boundary
plot3(ax,lgm(:,1),lgm(:,2),lgm(:,3),'k','LineWidth',3);
plot3(ax,lwm(:,1),lwm(:,2),lwm(:,3),'k','LineWidth',6);
if save_fig
   fig_name = 'Fig2b'; 
   fig_fold = fullfile(mat_dir,'figures');    
   print(fig,fullfile(fig_fold,[fig_name '.tiff']),'-dtiff','-cmyk','-r300'); 
end
end