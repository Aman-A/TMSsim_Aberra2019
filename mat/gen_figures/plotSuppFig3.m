function plotSuppFig3(save_fig)
%PLOTSUPPFIG3 Plots threshold differences between different E-field
%directions for Supp. Fig. 3
%  
% AUTHOR    : Aman Aberra 
if nargin==0
   save_fig = 0;
end
mat_dir = addPaths; 
% Data settings
nrn_model_ver = 'maxH';
tms_mode = 1; % monophasic
cell_ids = 1:25;
model_prefix = sprintf('utms_%s_w%g',nrn_model_ver,tms_mode); 
% Plot settings
theta_up = 60; theta_mid = 120; % cutoffs for E-field directions (thirds)
fig_size = [9 9];
markers = {'+','^','x','o','s'}; 
m_size = 50; % 8 for plot, 50 for scatter 
m_lw = 0.5;
font_size = 10; 
x_lims = [0.7 5.3];
y_lims = [-40 40]; 
font_name = 'Times';
cell_labels = {'L1 NGC','L2/3 PC','L4 LBC','L5 PC','L6 PC'}; 
num_cells = length(cell_labels); 
num_clones = length(cell_ids)/num_cells; 
colors = jet(num_cells); colors(4,:) = [0.5961 0.3059 0.6392]; 
rand_jit_max = 0.1; % range of rand distribution for point jitter (for marker_mode = 2)
rng(99);     % seed rn generator
%% Load data
data_fold = fullfile(mat_dir,'nrn_sim_data'); 
data_file = sprintf('%s_%g-%g.mat',model_prefix,cell_ids(1),cell_ids(end));
% Load composite data file (rax=0 for all cells, except rax=2 for cells 16-25 (L5/L6 PC))
data = load(fullfile(data_fold,data_file)); 
threshEs = data.threshEs;
phis = data.phis; 
thetas = data.thetas;
num_phi = length(unique(phis)); 
num_theta = length(unique(thetas)); 
thetasu = unique(thetas); 
%% Analysis
% Compute mean threshold in top, middle, and bottom third of threshold maps
% (upward/outward, transverse, and downward/inward E-field directions, respectively)
meanthr = zeros(num_theta,length(threshEs)); % threshold averaged across azimuthal rotations at each polar angle
thr_anis_all =  zeros(length(threshEs),1); % anisotropy across all E-field directions
% mean threshold within 3rd
meanthr_up = zeros(length(threshEs),1); % <= theta_up
meanthr_mid = zeros(length(threshEs),1); % > theta_up & <= theta_mid
meanthr_down = zeros(length(threshEs),1); % >theta_mid (up to 180)
minthr = zeros(length(threshEs),1); 
for i = 1:length(threshEs)
    threshEsi = threshEs{i}; 
    thr = [repmat(threshEsi(1),1,num_phi);
        reshape(threshEsi(2:end-1),num_phi,num_theta-2)';
        repmat(threshEsi(end),1,num_phi)];     
    meanthr(:,i) = mean(thr,2); 
    thr_anis_all(i) = max(max(thr))/min(min(thr));     
    meanthr_up(i) = mean(meanthr(thetasu <= theta_up,i)); 
    meanthr_mid(i) = mean(meanthr(thetasu > theta_up & thetasu <= theta_mid,i)); 
    meanthr_down(i) = mean(meanthr(thetasu > theta_mid,i));         
    minthr(i) = min(threshEsi); 
end
% calc directional anisotropy
meanthr_anis_mid_down = 100*(meanthr_mid-meanthr_down)./meanthr_down;
meanthr_anis_mid_up = 100*(meanthr_mid-meanthr_up)./meanthr_up;
%% Plot transverse - inward E-field: (mid - down)/down
meanthr_anis_mid_down2 = reshape(meanthr_anis_mid_down,num_cells,num_clones); 
fig1 = figure('Color','w'); 
fig1.Units = 'centimeters';% down/mid
fig1.Position(3:4) = fig_size; 
for i = 1:num_cells
   for j = 1:num_clones
       scatter(i,meanthr_anis_mid_down2(j,i),m_size,'MarkerEdgeColor',colors(i,:),...
           'Marker',markers{j},'MarkerFaceColor','none','LineWidth',m_lw,...
           'jitter','on','jitterAmount',rand_jit_max); hold on; 
   end  
end
ax = gca;
ax.FontSize = font_size; 
ax.FontName = font_name;
ax.YColor = 'k'; ax.XColor = 'k';
ax.XTick = 1:num_cells;
ax.XTickLabel = cell_labels; 
ax.XLim = x_lims; 
ax.YLim = y_lims;
ax.YGrid = 'on';
box off;
plot(x_lims,[0 0],'k','LineWidth',0.5)
ylabel('Percent threshold difference (%)');
%% Plot transverse - outward E-field: (mid-up)/up
meanthr_anis_mid_up2 = reshape(meanthr_anis_mid_up,num_cells,num_clones); 
fig2 = figure('Color','w'); 
fig2.Units = 'centimeters';% down/mid
fig2.Position(3:4) = fig_size; 
for i = 1:num_cells
   for j = 1:num_clones
       scatter(i,meanthr_anis_mid_up2(j,i),m_size,'MarkerEdgeColor',colors(i,:),...
           'Marker',markers{j},'MarkerFaceColor','none','LineWidth',m_lw,...
           'jitter','on','jitterAmount',rand_jit_max); hold on; 
   end      
end
ax = gca;
ax.FontSize = font_size; 
ax.FontName = font_name;
ax.YColor = 'k'; ax.XColor = 'k';
ax.XTick = 1:num_cells;
ax.XTickLabel = cell_labels; 
ax.XLim = x_lims; 
ax.YLim = y_lims; 
ax.YGrid = 'on';
box off;
plot(x_lims,[0 0],'k','LineWidth',0.5)
ylabel('Percent threshold difference (%)');
%% Save
if save_fig
   fig_name1 = 'SuppFig3a'; fig_name2 = 'SuppFig3b';
   fig_fold = fullfile(mat_dir,'figures');
   savefig(fig1,fullfile(fig_fold,[fig_name1 '.fig']));
   print(fig1,fullfile(fig_fold,[fig_name1 '.eps']),'-depsc2','-cmyk');   
   savefig(fig2,fullfile(fig_fold,[fig_name2 '.fig']));
   print(fig2,fullfile(fig_fold,[fig_name2 '.eps']),'-depsc2','-cmyk');   
end
