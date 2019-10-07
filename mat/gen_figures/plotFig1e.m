function plotFig1e(save_fig)
% PLOTFIG1E Plots horizontal box plot of minimum thresholds for monophasic
% TMS for Fig. 1e
if nargin == 0
   save_fig = 0;  
end
% Data settings
mat_dir = addPaths; 
nrn_model_ver = 'maxH';
tms_mode = 1; % monophasic
cell_ids = 1:25;
model_prefix = sprintf('utms_%s_w%g',nrn_model_ver,tms_mode); 
fig_fold = fullfile(mat_dir,'figures'); 
data_fold = fullfile(mat_dir,'nrn_sim_data'); % save .mat file here
% plot settings
lw = 1; 
% marker_mode = 1 - all same marker (1st in markers), 2 - all diff marker
marker_mode = 2; 
markers = {'+','^','x','o','s'}; 
m_size = 50;  % marker size
m_lw = 0.5; % marker line width
font_size = 8; 
font_name = 'Times';
font_weight = 'normal'; % or 'bold'
fig_dims = [20.8 3.8]; 
cell_labels = {'L1 NGC','L2/3 PC','L4 LBC','L5 PC','L6 PC'}; 
numCells = 5; numClones = 5; 
x_lim = [0.7 5.3];
include_ave = 0; % 0 - none, 1 - mean, 2 - median line
xoff = 0.2; % offset from tick for mean line
rand_jit_max = 0.2; % range of rand distribution for point jitter (for marker_mode = 2)
rng(50);     % seed rn generator 
colors = jet(numCells); colors(4,:) = [0.5961 0.3059 0.6392]; 
%% Load Data
data_file = sprintf('%s_%g-%g.mat',model_prefix,cell_ids(1),cell_ids(end));
% Load composite data file (rax=0 for all cells, except rax=2 for cells 16-25 (L5/L6 PC))
data = load(fullfile(data_fold,data_file)); 
threshEs = data.threshEs;
%%
min_threshEs = cellfun(@min,threshEs); 
min_threshEs = reshape(min_threshEs,numCells,numClones);
% plot_cell_data(cell_model_names,thetas,phis,threshEs,'Threshold','V/m',default_celltypes,which_plot); 
fig = figure('Color','w');
for i = 1:numCells
    if marker_mode == 1
       plot(i*ones(numCells,1),min_threshEs(:,i),'LineStyle','none',...
           'Marker',markers{1},'Color',colors(i,:),'MarkerSize',m_size,...
           'LineWidth',m_lw);  hold on;
    elseif marker_mode == 2               
       for j = 1:numClones          
          scatter(i,min_threshEs(j,i),m_size,'MarkerEdgeColor',colors(i,:),...
           'Marker',markers{j},'MarkerFaceColor','none',...
           'LineWidth',m_lw,'jitter','on','jitterAmount',rand_jit_max);  hold on; 
       end
    end
    if include_ave == 1 % mean
        plot([i-xoff,i+xoff],mean(min_threshEs(:,i))*ones(1,2),'Color',colors(i,:),...
            'LineWidth',lw);
    elseif include_ave == 2 % median
        plot([i-xoff,i+xoff],median(min_threshEs(:,i))*ones(1,2),'Color',colors(i,:),...
            'LineWidth',lw);
    end
end
ax = gca;
ax.FontSize=font_size;
ax.LineWidth=0.75;
ax.FontName=font_name;
ax.FontWeight=font_weight;
ax.XGrid='off';
ax.YGrid = 'on';
ax.YMinorGrid = 'off';
ax.YScale = 'linear';
ax.XLim = x_lim; 
ax.YLim = [150 400];
ax.YColor = 'k'; ax.XColor = 'k';
ylabel('Minimum threshold (V/m)'); 
% ax.TickLabelInterpreter = 'tex';
ax.XTick = 1:numCells; 
ax.XTickLabel = cell_labels; 
ax.Box = 'off';
fig.Units = 'centimeters';
fig.Position(3:4) = fig_dims; 
if save_fig
   fig_name = 'Fig1e'; 
   savefig(fig,fullfile(fig_fold,fig_name));       
   print(fig,fullfile(fig_fold,[fig_name '.eps']),'-depsc2','-cmyk'); 
end