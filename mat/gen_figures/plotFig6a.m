function plotFig6a(save_fig,fig_name,ax_view)
% Plots median threshold across clones and rotations at each location on
% layers 1-5, adjacent to each other, as in Fig. 3c
if nargin==0
   save_fig = 0;   
   fig_name = 'Fig6a';
   ax_view = [-89.2 45]; % [-94.78 68.66] or [-89.2 70.8]
end
mat_dir = addPaths; 
% Simulation settings
nrn_model_ver = 'maxHlin';
mode = 1; % monophasic MagProX100 pulse
layer_set_num = 1;
Efield_name = 'M1_PA_MCB70'; 
cell_ids = {[];6:10;[];[];[]}; 
nrn_pop = 'nrn_pop1'; 
model_prefix1 = sprintf('tms_%s_w%g_ls_%g_E_%s_P_%s',nrn_model_ver,mode,...
                            layer_set_num,Efield_name,nrn_pop); 
model_prefix2 = sprintf('tms_%s_w%g_ls_%g_E_%s_P_%s',nrn_model_ver,mode,...
                            layer_set_num,[Efield_name '_r'],nrn_pop); 
%% Plot settings
cmap = [flipud(fake_parula(1000));0.8 0.8 0.8]; % add gray for values above cutoff
% clims = [90 500]; % A/us
% clims = [90 1010]; % A/us
clims = [70 230]; 
z_lims = [22 52.4057]; % or []
lt_pos = [-411 -807 836]; % [170.3353 52.0246 1.2195e3] 
%% Load data
layers = loadLayers(layer_set_num); 
data_fold = fullfile(mat_dir,'nrn_sim_data');
data_struct = load(fullfile(data_fold,model_prefix1)); 
threshEs = data_struct.threshEs; 
cell_model_names = data_struct.cell_model_names;
data_struct = load(fullfile(data_fold,model_prefix2)); 
threshEs2 = data_struct.threshEs; 
cell_model_names2 = data_struct.cell_model_names;
%% Plot - PA
plotDataLayers(layers,threshEs,cell_model_names,cell_ids); 
fig = gcf;
fig.Units = 'normalized';
fig.OuterPosition = [0 0 1 1];
fig.Color = 'w';
colormap(cmap); caxis(clims); 
view(ax_view); 
% change light
ax = gca;
ax.Children(1).Position = lt_pos;
ax.Children(1).Style = 'local';
% cut off sulcus
if ~isempty(z_lims)
   ax.ZLim = z_lims; 
end
%% Plot - AP
plotDataLayers(layers,threshEs2,cell_model_names2,cell_ids); 
fig2 = gcf;
fig2.Units = 'normalized';
fig2.OuterPosition = [0 0 1 1];
fig2.Color = 'w';
colormap(cmap); caxis(clims); 
view(ax_view); 
% change light
ax2 = gca;
ax2.Children(1).Position = lt_pos;
ax2.Children(1).Style = 'local';
% cut off sulcus
if ~isempty(z_lims)
   ax2.ZLim = z_lims; 
end
%% Save figs
if save_fig   
   fig_fold = fullfile(mat_dir,'figures');    
   savefig(fig,fullfile(fig_fold,[fig_name '-1.fig']));
   print(fig,fullfile(fig_fold,[fig_name '-1.tif']),'-dtiff','-cmyk','-r250');
   savefig(fig2,fullfile(fig_fold,[fig_name '-2.fig']));
   print(fig,fullfile(fig_fold,[fig_name '-2.tif']),'-dtiff','-cmyk','-r250');
end
