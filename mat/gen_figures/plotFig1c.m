function plotFig1c(save_fig)
%PLOTFIG1C Plots single threshold-direction map for example L2/3 PC,
%monophasic TMS pulse
%  
% AUTHOR    : Aman Aberra
if nargin == 0
   save_fig = 0;  
end
mat_dir = addPaths; 
nrn_model_ver = 'maxH';
tms_mode = 1; % monophasic
cell_id = 6; % 'L23_PC_cADpyr229_1'
model_prefix = sprintf('utms_%s_w%g',nrn_model_ver,tms_mode); 
data_cell_ids = 1:25; % cell_ids in data file
%% Plot settings
cmap = flipud(fake_parula(1000)); 
display_scale = 'lin';
plot_min_point = 1; 
title_on = 0; 
%% Load data
fig_fold = fullfile(mat_dir,'figures'); 
data_fold = fullfile(mat_dir,'nrn_sim_data'); % save .mat file here
data_file = sprintf('%s_%g-%g.mat',model_prefix,data_cell_ids(1),data_cell_ids(end));
data = load(fullfile(data_fold,data_file)); 
theta = data.thetas;
phi = data.phis; 
if iscell(theta)
    theta = theta{1};
end
if iscell(phi)
   phi = phi{1};  
end
threshEs = data.threshEs{cell_id};
%% Plot
fig = figure('units','normalized','outerposition',[0 0 1 1],'Color','w');
plotThreshMapMW(threshEs,theta,phi,'plot_min_point',plot_min_point,...
                'display_scale',display_scale,'title_on',title_on);
ax = gca;
colormap(ax,cmap); 
if save_fig
    fig_name = fullfile(fig_fold,'Fig1c.tiff');    
    print(fig,fig_name,'-dtiff','-r300','-cmyk');
    delete(fig.Children(1)); % remove colorbar    
end
end
