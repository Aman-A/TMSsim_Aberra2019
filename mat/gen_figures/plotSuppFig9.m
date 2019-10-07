function plotSuppFig9(save_fig)
%PLOTSUPPFIG9 plot boxplots of threshold distributions within FDI 
% subROIs
%  
% AUTHOR    : Aman Aberra 
if nargin==0
   save_fig = 1;
end
mat_dir = addPaths; 
% Simulation settings
nrn_model_ver = 'maxH';
modes = [1,3,2]; % 1 monophasic, 3 half sine, 2 biphasic MagProX100 pulses
layer_set_num = 1;
Efield_name = 'M1_PA_MCB70'; 
plot_layer = 4; % l4 of 5 (L5 PCs)
cell_ids = 16:20;
nrn_pop = 'nrn_pop1-nrn_pop6_all'; 
model_prefixes = [arrayfun(@(x) sprintf('tms_%s_w%g_ls_%g_E_%s_P_%s',nrn_model_ver,x,...
                            layer_set_num,Efield_name,nrn_pop),modes,'UniformOutput',0);
                   arrayfun(@(x) sprintf('tms_%s_w%g_ls_%g_E_%s_r_P_%s',nrn_model_ver,x,...
                            layer_set_num,Efield_name,nrn_pop),modes,'UniformOutput',0)]; 
% reverse biphasic wave, b/c 2nd (dominant) phase is in opposite direction
model_prefixes([5,6]) = model_prefixes([6,5]);
model_prefixes = model_prefixes(:);  
model_names = repmat({'P-A','A-P'},1,3); 
model_names([5,6]) = model_names([6,5]); 
% Plot settings
plot_region_name = 'FDI_rep_inds_subROIs';
lw = 1; 
exp_color = 0.8*ones(1,3); 
lw_exp = 2;
font_size = 12;
font_name = 'Times';
yaxlim = [30 1000]; % A/us
colors = (1/255)*[27,158,119;217,95,2;117,112,179;231,41,138];
% colors = jet(length(cell_ids)); 
% colors(4,:) = [0.5961 0.3059 0.6392]; 
%% Load data
% FDI representation
ROI_filename = fullfile(mat_dir,'gen_figures',[plot_region_name '.mat']); 
ROI_data = load(ROI_filename);
roi_names = fieldnames(ROI_data); 
ROIi_all = struct2cell(ROI_data);
% reorder to left/right sup, left/right inf
roi_names = roi_names([2,4,1,3]);
ROIi_all = ROIi_all([2,4,1,3]); 
ROIi_all = cellfun(@(x) x{plot_layer},ROIi_all,'UniformOutput',0); % extract indices in single layer for each ROI
num_rois = length(ROIi_all);
% Load threshold data
num_models = length(model_prefixes); 
threshEs_all_models = cell(num_models,1); 
% loads data
for i = 1:num_models
    datai = load(fullfile(mat_dir,'nrn_sim_data',[model_prefixes{i} '.mat']));     
    threshEs_all_models{i} = datai.threshEs;     
end
cell_model_names = datai.cell_model_names;
%% Plot boxplots
fig = figure('Color','w'); 
plotDistsModelROIs(threshEs_all_models,ROIi_all,cell_ids,cell_model_names,colors);
ax = gca;
ax.YScale = 'log';
ax.YLim = yaxlim;
ax.YTick = [50 100 200 500 1000]; 
box(ax,'off'); 
delete(fig.Children(1)); % remove legend
% add fill boxes
addBoxPlotFillAll(colors,num_rois,num_rois*num_models);
ax.FontName = font_name; 
ax.FontSize = font_size;
ax.YLim = yaxlim; 
ax.XTickLabel = model_names; 
fig.Units = 'centimeters';
fig.Position(3:4) = [14.9,11.2];
ax.YGrid = 'off'; ax.YMinorGrid = 'off'; 
ax.XColor = 'k'; ax.YColor = 'k';
ylabel('Threshold coil current (A/\mu s)');
if save_fig
   fig_name = 'SuppFig9';
   fig_fold = fullfile(mat_dir,'figures');
   savefig(fig,fullfile(fig_fold,[fig_name '.fig']));
   print(fig,fullfile(fig_fold,[fig_name '.eps']),'-depsc2','-cmyk','-painters');
end
