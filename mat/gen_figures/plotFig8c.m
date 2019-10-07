function plotFig8c(save_fig)
%PLOTFIG8C Plots time constants for range of threshold percentiles of L2/3,
%L4, and L5 in FDI representation for Fig. 8c
%  
% AUTHOR    : Aman Aberra 
if nargin==0
   save_fig = 0;
end
mat_dir = addPaths; 
% Simulation settings
cTMS_durs = [30,60,120]; % us
layer_set_num = 1; 
Efield_name = 'M1_PA_Magstim70mm';
nrn_model_ver = 'maxH';
mode = 4; % cTMS1
nrn_pop='nrn_pop1-nrn_pop6_all'; 
model_prefixes = arrayfun(@(x) sprintf('tms_%s_w%g_%gus_ls_%g_E_%s_P_%s',nrn_model_ver,mode,...
                            x,layer_set_num,Efield_name,nrn_pop),cTMS_durs,'UniformOutput',0)';
cell_ids = {6:10;11:15;16:20}; 
% Plot settings
plot_layers = 2:4;
plot_region_name = 'FDI_rep_inds';
cutoffs = [.025,0.05,0.1,0.25,0.5];
exp_lw = 1; % experimental data line width
exp_msize = 8; % experimental data marker size ('o')
font_size = 10; 
font_name = 'Times'; 
colors = [0 0.5 1;0 1 1; 0.5 1 0.5; 0.5961 0.3059 0.6392; 1 0.5 0];
bar_col = colors(plot_layers,:);
exp_color = 'w';
bar_width = 0.8;
%% Load data
% Load cTMS waveforms
cTMS_data = load(fullfile(mat_dir,'run_nrn_functions/cTMSwaves.mat')); 
t = cTMS_data.t; % time in sec
pw = cTMS_data.pw;
[~,pw_inds,~] = intersect(pw,cTMS_durs*1e-3);
W = cTMS_data.W(:,pw_inds); 
% FDI representation
ROI_filename = fullfile(mat_dir,'gen_figures',[plot_region_name '.mat']); 
ROI_data = load(ROI_filename);
ROIi = ROI_data.ROIi(plot_layers); 
% Load threshEs and extract thresholds within ROI
num_pws = length(cTMS_durs); 
num_layers = length(plot_layers); 
num_cutoffs = length(cutoffs);
quants = zeros(num_pws,num_layers,num_cutoffs); 
threshEsROI = cell(num_pws,1); 
for i = 1:num_pws   
   threshEsROI{i} = cell(1,num_layers); 
   data_fold = fullfile(mat_dir,'nrn_sim_data');
   data_struct = load(fullfile(data_fold,model_prefixes{i}));
   threshEs = data_struct.threshEs;
   cell_model_names = data_struct.cell_model_names;
   for j = 1:num_layers       
       [~,cell_inds] = intersect(cell_model_names,cellModelNames(cell_ids{j}));
        threshEsj = cell2mat(threshEs(cell_inds));                  
        threshEsj =  threshEsj(ROIi{j},:);
        threshEsROI{i}{j} =  threshEsj;
        quants(i,j,:) = quantile(threshEsj(:),cutoffs); 
   end
end
% Peterchev 2013 data 
exp_range = [200 - 33,200,200 + 33]; % 200 +/- 33 individual mean +/- SD
%% Get time constants
taus = zeros(num_layers,num_cutoffs); % time constants - us
rbs = zeros(num_layers,num_cutoffs); % rheobase - A/us
for i = 1:num_layers
    for j = 1:length(cutoffs)
       mt = quants(:,i,j);
       [taus(i,j),rbs(i,j)] = estTimeConstant(cTMS_durs,mt,t,W);
    end
end
%% Plot
ingroup_space = 0.1;
btgroup_space = bar_width*2;
x = 0; % initialize
alphas = fliplr(1/(num_cutoffs):1/(num_cutoffs):1);
fig = figure('Color','w');
bar(x,exp_range(2),'BarWidth',bar_width,'FaceColor',exp_color); hold on; 
errorbar(x,exp_range(2),exp_range(2)-exp_range(1),exp_range(3)-exp_range(2),'-k','Marker','none',...
    'MarkerFaceColor','none','LineWidth',exp_lw,'MarkerSize',exp_msize);
for i = 1:length(plot_layers)   
    x = x(end) + btgroup_space + (0:(ingroup_space+bar_width):((ingroup_space+bar_width)*(num_cutoffs-1)));
    for j = 1:length(x)
        bar(x(j),taus(i,j),'BarWidth',bar_width,...
            'FaceColor',bar_col(i,:),'FaceAlpha',alphas(j),'EdgeColor','k');
    end
end
ax = gca;
box(ax,'off');
ax.YGrid = 'on';
ax.FontSize = font_size;
ax.FontName = font_name;
ax.YColor = 'k'; ax.XColor =  'k';
fig.Units = 'centimeters';
fig.Position(3:4) = [10 6.5];
ax.XTickLabel = {}; 
xlim([-bar_width/2-0.5,ax.Children(1).XData + bar_width/2+0.5]);  
ax.XTick = [0,(btgroup_space+(bar_width+ingroup_space)*floor(length(cutoffs)/2)):(btgroup_space+(bar_width+ingroup_space)*2*floor(length(cutoffs)/2)):ax.XLim(2)];
if save_fig
   fig_name = 'Fig8c';
   fig_fold = fullfile(mat_dir,'figures');
   savefig(fig,fullfile(fig_fold,[fig_name '.fig']));
   print(fig,fullfile(fig_fold,[fig_name '.png']),'-dpng','-cmyk');   
end
