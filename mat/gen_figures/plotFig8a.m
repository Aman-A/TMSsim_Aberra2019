function plotFig8a(save_fig)
%PLOTFIG8A Plot median timeconstant across clones/rotations on L2/3,L4,L5 
%for Fig. 8a
%  
% AUTHOR    : Aman Aberra 
if nargin==0
   save_fig = 0;
end
mat_dir = addPaths; 
% Simulation settings
cTMS_durs = [30,60,120]; % us
cell_ids = {6:10;11:15;16:20}; 
layer_set_num = 1;
Efield_name = 'M1_PA_Magstim70mm';
nrn_model_ver = 'maxH';
mode = 4; % cTMS1
nrn_pop='nrn_pop1-nrn_pop6_all'; 
model_prefix = sprintf('tms_%s_w%g_%g-%gus_ls_%g_E_%s_P_%s_timeconstants',nrn_model_ver,mode,...
                            cTMS_durs(1),cTMS_durs(end),layer_set_num,Efield_name,nrn_pop); 
% Plot settings
plot_region_name = 'FDI_rep_inds';
plot_layers = 2:4; 
z_lims = [22 52.4057]; % or []
ax_view = [-89.2 45]; % [-89.2 70.8]
lt_pos = [-411 -807 836]; % [170.3353 52.0246 1.2195e3] 
shift_dir = [0 -35 0]; 
cax_lims = [100 400]; % set by min, max
cmap_name = 'jet'; 
cmap = jet(1000); 
lw = 1; 
%% Load data
layers = loadLayers(layer_set_num); 
data_fold = fullfile(mat_dir,'nrn_sim_data');
data_struct = load(fullfile(data_fold,model_prefix)); 
timeconstants = data_struct.timeconstants; 
cell_model_names = data_struct.cell_model_names;
% FDI representation
ROI_filename = fullfile(mat_dir,'gen_figures',[plot_region_name '.mat']); 
ROI_data = load(ROI_filename);
ROIi = ROI_data.ROIi(plot_layers); 
layersROIs = [layers.surface];
layersROIs = layersROIs(plot_layers); 
%% Plot
plotDataLayers(layers(plot_layers),timeconstants,cell_model_names,cell_ids); 
fig = gcf;
fig.Color = 'w';
% add FDI representation outlines
for i = 1:length(plot_layers)
    layersROIs(i).faces = layersROIs(i).faces(ROIi{i},:);
    [~,vb] = meshboundary(layersROIs(i));
    plot3(vb(1,:)+shift_dir(1)*(i-1),vb(2,:)+shift_dir(2)*(i-1),vb(3,:)+shift_dir(3)*(i-1),'k','LineWidth',lw);
end
ax = gca;
axis(ax,'off'); 
lt_ind = length(plot_layers)+1; 
ax.Children(lt_ind).Position = lt_pos;
ax.Children(lt_ind).Style = 'local';
ax.View = ax_view; 
colormap(ax,cmap); 
zlim(ax,z_lims);  
caxis(ax,cax_lims); 
fig.Units = 'centimeters';
fig.Position(3:4) = [16 7.3]*1.3;
if save_fig
   fig_name = 'Fig8a';
   fig_fold = fullfile(mat_dir,'figures');
   savefig(fig,fullfile(fig_fold,[fig_name '.fig']));
   print(fig,fullfile(fig_fold,[fig_name '.tiff']),'-dtiff','-cmyk','-r300','-painters');         
end
