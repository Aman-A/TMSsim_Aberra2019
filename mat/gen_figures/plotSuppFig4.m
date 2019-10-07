function plotSuppFig4(save_fig)
%PLOTSUPPFIG4 Plot threshold-direction maps for L2/3 PCs with and without
%straight axons for Supp. Fig 4
% 
%   Inputs 
%   ------ 
%   save_fig : 1 or 0
%             Set to 1 to save .fig and .tiff figure files
%   fig_mode : 1 or 2
%             Set to 1 to plot Supp. Fig. 5, set to 2 to plot Supp. Fig. 6
%
% AUTHOR    : Aman Aberra 
if nargin==0
   save_fig = 0;   
end
mat_dir = addPaths; 
nrn_model_ver1 = 'maxHlin'; 
nrn_model_ver2 = 'maxH'; 
tms_mode = 1; % monophasic
cell_ids = 6:10; 
model_prefix1 = sprintf('utms_%s_w%g_rax5',nrn_model_ver1,tms_mode); % original
data_cell_ids1 = 6:10; % cell_ids for model_prefix1
model_prefix2 = sprintf('utms_%s_w%g',nrn_model_ver2,tms_mode); % axon disabled
data_cell_ids2 = 1:25; % cell_ids for model_prefix2
% Plot settings
cmap1 = flipud(fake_parula(1000)); % threshold maps
clims_lin = [200 1000];
clims_real = [175 500]; 
clims2 = [-40 900];
mw_settings.display_scale = 'log';
mw_settings.colorbar_auto=1;
mw_settings.colorbar_norm=0;
mw_settings.colorbar_min=clims_lin(1);
mw_settings.colorbar_max=clims_lin(2);
mw_settings.colorbar_on=0; 
mw_settings.title_on=0; 
mw_settings.display_grid=0; 
%% Load data
data_fold = fullfile(mat_dir,'nrn_sim_data'); % save .mat file here
data_file1 = sprintf('%s_%g-%g.mat',model_prefix1,data_cell_ids1(1),data_cell_ids1(end));
% Load rax = 0 data for L5 PCs
data = load(fullfile(data_fold,data_file1));
thetas = data.thetas;
phis = data.phis;
cell_model_names1 = data.cell_model_names; 
[~,~,data_inds] = intersect(cellModelNames(cell_ids),cell_model_names1); % get indices of layer cells in threshEs
threshEs1 = data.threshEs(data_inds); 
% Load composite data file (rax=0 for all cells, except rax=2 for cells 32-36,52-56)
data_file2 = sprintf('%s_%g-%g.mat',model_prefix2,data_cell_ids2(1),data_cell_ids2(end));
data = load(fullfile(data_fold,data_file2));
cell_model_names2 = data.cell_model_names; 
[~,~,data_inds] = intersect(cellModelNames(cell_ids),cell_model_names2); % get indices of layer cells in threshEs
threshEs2 = data.threshEs(data_inds); 
% Get percent difference (linear axon - full axon)/full axon
threshEs_diff = cellfun(@(x,y) 100*(x-y)./y,threshEs1,threshEs2,'UniformOutput',0);
% Combine all
threshEs_all = [threshEs1,threshEs2,threshEs_diff]; % [original;modified;difference]
%% Plot 
num_rows = 3; num_cols = 5; 
% Generate grid and axes handles on gcf
ax_handles = cell(num_rows,num_cols);
fprintf('Printing with %g rows and %g cols, %g entries\n',num_rows,num_cols,num_rows*num_cols);
fig = figure('units','normalized','outerposition',[0 0 1 1],'Color',[1 1 1]); 
for r = 1:num_rows
    for c = 1:num_cols
        ax_handles{r,c} = axes('Position',[(0.9/num_cols)*(c-1),(0.95/num_rows)*(num_rows-r),0.9/num_cols,0.9/num_rows ]);
        axis(ax_handles{r,c},'off');
        hold(ax_handles{r,c},'all');
    end
end
plotMWGrid(ax_handles(1,:),threshEs_all(1:num_cols),thetas,phis,1,num_cols,mw_settings);
mw_settings.colorbar_min=clims_real(1);
mw_settings.colorbar_max=clims_real(2);
plotMWGrid(ax_handles(2,:),threshEs_all((num_cols+1):(num_cols*2)),thetas,phis,1,num_cols,mw_settings);
% ratio
mw_settings.display_scale = 'lin';
mw_settings.colorbar_min=clims2(1);
mw_settings.colorbar_max=clims2(2);
mw_settings.plot_min_point=0; 
plotMWGrid(ax_handles(3,:),threshEs_all((num_cols*2+1):end),thetas,phis,1,num_cols,mw_settings);
for r = 1:num_rows
   for c = 1:num_cols
       if r <= 2
           colormap(ax_handles{r,c},cmap1);
       else
           caxis(ax_handles{r,c},clims2); 
           colormap(ax_handles{r,c},bluewhitered(1000));           
       end
   end
end
%%
if save_fig
   fig_name = 'SuppFig4';   
   fig_fold = fullfile(mat_dir,'figures');
   savefig(fig,fullfile(fig_fold,[fig_name '.fig']));
   print(fig,fullfile(fig_fold,[fig_name '.tiff']),'-dtiff','-r250','-cmyk');
end
