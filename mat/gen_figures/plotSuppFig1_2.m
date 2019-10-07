function plotSuppFig1_2(save_fig,fig_mode)
%PLOTSUPPFIG1_2 Plot threshold-direction maps for L5/6 PCs with and without
%axon disabled for Supp. Fig. 1/2
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
   fig_mode = 1; % 1 - Supp Fig 5 - L5 PCs. 2 - Supp Fig 6 - L6 PCs  
end
mat_dir = addPaths; 
nrn_model_ver = 'maxH'; 
tms_mode = 1; % monophasic
cell_ids = {16:20;21:25}; 
cell_ids = cell_ids{fig_mode}; 
model_prefix1 = sprintf('utms_%s_w%g_rax0',nrn_model_ver,tms_mode); % original
data_cell_ids1 = 16:25; % cell_ids for model_prefix1
model_prefix2 = sprintf('utms_%s_w%g',nrn_model_ver,tms_mode); % axon disabled
data_cell_ids2 = 1:25; % cell_ids for model_prefix2
% Plot settings
cmap1 = flipud(fake_parula(1000)); % threshold maps
clims = {[120 500];[200 800]};
clims = clims{fig_mode}; 
cmap2 = jet(1000); % difference maps
clims2 = {[-5 120];[-20 200]};
clims2 = clims2{fig_mode}; 
mw_settings.display_scale = 'lin';
mw_settings.colorbar_auto=0;
mw_settings.colorbar_norm=0;
mw_settings.colorbar_min=clims(1);
mw_settings.colorbar_max=clims(2);
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
% Get percent difference  - (full axon - disabled term)/(disabled term)
threshEs_diff = cellfun(@(x,y) 100*(x-y)./y,threshEs2,threshEs1,'UniformOutput',0);
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
plotMWGrid(ax_handles(1:2,:),threshEs_all(1:(num_cols*2)),thetas,phis,num_rows-1,num_cols,mw_settings);
% ratio
% mw_settings.display_scale = 'lin';
mw_settings.colorbar_min=clims2(1);
mw_settings.colorbar_max=clims2(2);
mw_settings.plot_min_point=0; 
plotMWGrid(ax_handles(3,:),threshEs_all((num_cols*2+1):end),thetas,phis,1,num_cols,mw_settings);
for r = 1:num_rows
   for c = 1:num_cols
       if r <= 2
           colormap(ax_handles{r,c},cmap1);
       else
           colormap(ax_handles{r,c},cmap2);
       end
   end
end
%%
if save_fig
   fig_names = {'SuppFig1','SuppFig2'};
   fig_name = fig_names{fig_mode}; 
   fig_fold = fullfile(mat_dir,'figures');
   savefig(fig,fullfile(fig_fold,[fig_name '.fig']));
   print(fig,fullfile(fig_fold,[fig_name '.tiff']),'-dtiff','-r250','-cmyk');
end
