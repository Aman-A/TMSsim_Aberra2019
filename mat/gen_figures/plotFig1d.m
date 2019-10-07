function plotFig1d(save_fig)
% PLOTFIG1D Plots grid of threshold-direction maps (Mollweide projection) 
% Data should already be generated and saved in
% mat/nrn_sim_data/<model_prefix>.mat
if nargin == 0
   save_fig = 0;  
end
mat_dir = addPaths; 
nrn_model_ver = 'maxH';
tms_mode = 1; % monophasic
cell_ids = 1:25;
model_prefix = sprintf('utms_%s_w%g',nrn_model_ver,tms_mode); 
% Plot settings
cmap = flipud(fake_parula(1000)); 
mw_settings.display_scale = 'log';
mw_settings.colorbar_auto=0;
mw_settings.colorbar_norm=1;
mw_settings.colorbar_min=1;
mw_settings.colorbar_max=4;
mw_settings.colorbar_on=0; 
mw_settings.title_on=0; 
mw_settings.display_grid=0; 
%% Load data
data_fold = fullfile(mat_dir,'nrn_sim_data'); % save .mat file here
data_file = sprintf('%s_%g-%g.mat',model_prefix,cell_ids(1),cell_ids(end));
% Load composite data file (rax=0 for all cells, except rax=2 for cells 32-36,52-56)
data = load(fullfile(data_fold,data_file));
thetas = data.thetas;
phis = data.phis; 
threshEs = data.threshEs;
init_inds = data.init_inds; 
%% reorder
min_threshEs = reshape(cellfun(@min,threshEs),5,5); 
reorder = reshape(1:25,5,5); 
for i = 1:5
   [~,indi] = sort(min_threshEs(:,i),'descend');
   reorder(:,i) = reorder(indi,i); 
end
% reorder = [1:5:25,2:5:25,3:5:25,4:5:25,5:5:25];
reorder = reorder'; reorder = reorder(:); 
cell_ids = cell_ids(reorder); 
% cell_model_names = cell_model_names(reorder); 
threshEs = threshEs(reorder); 
init_inds = init_inds(reorder); 
%% Plot
fig = plotMWmorphGrid(threshEs,thetas,phis,init_inds,cell_ids,nrn_model_ver,mw_settings);
for i = (length(fig.Children)/2 + 1):length(fig.Children) % just MW maps
   colormap(fig.Children(i),cmap);
   if strcmp(mw_settings.display_scale,'log')
        caxis(fig.Children(i),log10([mw_settings.colorbar_min mw_settings.colorbar_max])) % log scale
   else
       caxis(fig.Children(i),[mw_settings.colorbar_min mw_settings.colorbar_max]) % linear scaling
   end
end
% adjust cell 36 plot
fig.Children(2).XLim = [-1000 1400];
% save fig
if save_fig
   fig_name = 'Fig1d';
   fig_fold = fullfile(mat_dir,'figures'); 
   print(fig,fullfile(fig_fold,[fig_name '.tiff']),'-dtiff','-r350','-cmyk'); 
end