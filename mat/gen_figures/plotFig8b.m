function plotFig8b(save_fig)
%PLOTFIG8B Plots strength-duration curves for L2/3, L4, and L5 with 50% threshold 
%percentile for Fig. 8b
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
exp_data_file = 'Peterchev2013_data';
cell_ids = {6:10;11:15;16:20}; 
% Plot settings
plot_layers = 2:4;
plot_region_name = 'FDI_rep_inds';
norm_ind = 2; % normalize to 60
cutoff = 0.5; 
lw = 1; % model line width
msize = 16; % model marker size ('.')
exp_lw = 1; % experimental data line width
exp_msize = 8; % experimental data marker size ('o')
font_size = 10; 
font_name = 'Times'; 
colors = [0 0.5 1;0 1 1; 0.5 1 0.5; 0.5961 0.3059 0.6392; 1 0.5 0];
line_col = colors(plot_layers,:);
layer_labels = {'L2/3 PC','L4 LBC','L5 PC'}; 
%% Load data
% FDI representation
ROI_filename = fullfile(mat_dir,'gen_figures',[plot_region_name '.mat']); 
ROI_data = load(ROI_filename);
ROIi = ROI_data.ROIi(plot_layers); 
% Load threshEs and extract thresholds within ROI
threshEsROI = cell(length(model_prefixes),1); 
num_pws = length(cTMS_durs); 
num_layers = length(plot_layers); 
quants = zeros(num_pws,num_layers); 
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
        quants(i,j) = quantile(threshEsj(:),cutoff); 
   end
end
% Load Peterchev 2013 data motor threshold (MT) data
exp_data = load(fullfile(mat_dir,'gen_figures',[exp_data_file,'.mat'])); 
MT = exp_data.MT; 
MTs = []; 
for i =1:length(MT)
    if (size(MT(i).MTs,1) == 3)
        MTs = [MTs;MT(i).MTs(:,2)']; 
    end
end
MTs(31,:) = []; 
meanMT = mean(MTs,1); 
stdMT = std(MTs,1); 
%% Plot
fig = figure('Color','w'); 
models = cell(length(plot_layers),1);
for i = 1:length(plot_layers)        
    models{i} = plot(cTMS_durs,quants(:,i)/quants(norm_ind,i),'Marker','.','LineWidth',lw,'Color',line_col(i,:));
    [models{i}.MarkerSize] = deal(msize);
    hold on;
end
exp = errorbar(cTMS_durs,meanMT/meanMT(norm_ind),stdMT/meanMT(norm_ind),'ok',...
        'Marker','o','MarkerFaceColor','none','LineWidth',exp_lw,'MarkerSize',exp_msize); 
ax = gca;
box(ax,'off');         
ax.YGrid = 'on';
ax.FontSize = font_size;
ax.FontName = font_name;
ax.YColor = 'k'; ax.XColor = 'k';
leg = legend([[models{:}],exp],[cellfun(@(x) sprintf('Model - %s',x),layer_labels,'UniformOutput',0),...
                    'Experiment - Mean �SD'],'FontName',font_name,'FontSize',font_size,'Box','off'); 
ax.XLim = [20 120]; 
ax.YLim = [0.5 2]; 
fig.Units = 'centimeters';
fig.Position(3:4) = [9 6.5]; 
leg.Position(1:2) = [0.4 0.65]; 
leg.Visible = 'off'; 
xlabel('Pulse Duration (\mu s)');  
ylabel(num2str(cTMS_durs(norm_ind),'Threshold (normalized to %g \\mu s)')) 
if save_fig
   fig_name = 'Fig8b';
   fig_fold = fullfile(mat_dir,'figures');
   savefig(fig,fullfile(fig_fold,[fig_name '.fig']));
   print(fig,fullfile(fig_fold,[fig_name '.eps']),'-depsc2','-cmyk');    
end
