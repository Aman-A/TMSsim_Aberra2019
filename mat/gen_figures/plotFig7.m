function plotFig7(save_fig)
%PLOTFIG7 Plot boxplots of threshold distributions within FDI
%representation with experimental data overlaid for Fig. 7
%  
% AUTHOR    : Aman Aberra 
if nargin==0
   save_fig = 0;
end
mat_dir = addPaths; 
% Simulation settings
nrn_model_ver = 'maxH';
modes = [1,3,2]; % 1 monophasic, 3 half sine, 2 biphasic MagProX100 pulses
layer_set_num = 1;
Efield_name = 'M1_PA_MCB70'; 
cell_ids = {1:5;6:10;11:15;16:20;21:25};  
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
exp_data_file = 'Sommer2006_data';
plot_region_name = 'FDI_rep_inds';
lw = 1; 
exp_color = 0.8*ones(1,3); 
lw_exp = 2;
font_size = 12;
font_name = 'Times';
yaxlim = [30 1000]; % A/us
colors = jet(length(cell_ids)); colors(4,:) = [0.5961 0.3059 0.6392]; 
%% Load data
% FDI representation
ROI_filename = fullfile(mat_dir,'gen_figures',[plot_region_name '.mat']); 
ROI_data = load(ROI_filename);
ROIi = ROI_data.ROIi; 
% experimental data from Sommer et al., 2006
som_mt_data = load(fullfile(mat_dir,'gen_figures',[exp_data_file '.mat'])); 
som_mt = som_mt_data.som_mt; 
quants = [0.01,0.25,0.5,0.75,0.99]; % [min,25th,median,75th,max]
som_mt = quantile(som_mt,quants); % [APm PAm APh PAh APb PAb]
som_mt = [som_mt(:,2),som_mt(:,1),som_mt(:,4),som_mt(:,3),som_mt(:,5:6)];
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
plotDistsModelLayer(threshEs_all_models,ROIi,cell_ids,cell_model_names,colors);
ax = gca;
ax.YScale = 'log';
ax.YLim = yaxlim;
ax.YTick = [50 100 200 500 1000]; 
box(ax,'off'); 
delete(fig.Children(1)); % remove legend
xvals = 3:5:28;
% add fill boxes
addBoxPlotFillAll(colors);
% Plot experimental data
xoff=0.5;
medoff = 0.2;
for i = 1:size(som_mt,2)   
    % plot median
    plot([0.5+5*(i-1)+medoff,5.5+5*(i-1)-medoff],[som_mt(3,i),som_mt(3,i)],'w','LineWidth',lw_exp); 
    % plot 75th - max
    plot([xvals(i),xvals(i),nan,0.5+5*(i-1)+xoff,5.5+5*(i-1)-xoff],...
        [som_mt(4,i),som_mt(5,i),nan,som_mt(5,i),som_mt(5,i)],'Color',exp_color,'LineWidth',lw_exp);
    % plot 25th - min
    plot([xvals(i),xvals(i),nan,0.5+5*(i-1)+xoff,5.5+5*(i-1)-xoff],...
        [som_mt(2,i),som_mt(1,i),nan,som_mt(1,i),som_mt(1,i)],'Color',exp_color,'LineWidth',lw_exp)
end
ax.Children = [ax.Children(3*size(som_mt,1):end);ax.Children(1:size(som_mt,2)*3)]; 
box_off = 0.1;
for i = 1:size(som_mt,2)    
%     pi = patch([0.5+5*(i-1),5.5+5*(i-1),5.5+5*(i-1),0.5+5*(i-1)],...
%         [som_mt(2,i),som_mt(2,i),som_mt(4,i),som_mt(4,i)],exp_color);
%     pi.EdgeColor = 'none';
    rectangle('Position',[0.5+5*(i-1)+box_off,som_mt(2,i),5-box_off*2,som_mt(4,i)-som_mt(2,i)],...        
        'FaceColor',exp_color,'EdgeColor','none','Curvature',[0.08,0.08]);    
end 
ax.Children = [ax.Children(size(som_mt,1):end);ax.Children(1:size(som_mt,2))]; 
set(ax,'layer','top')
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
   fig_name = 'Fig7';
   fig_fold = fullfile(mat_dir,'figures');
   savefig(fig,fullfile(fig_fold,[fig_name '.fig']));
   print(fig,fullfile(fig_fold,[fig_name '.eps']),'-depsc2','-cmyk','-painters');
end
