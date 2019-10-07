function plotSuppFig8(save_fig)
%PLOTSUPPFIG8 Plot model thresholds vs. experimental MT from Sommer et al
%2006 for Supp Fig. 8
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
% Plot settings
exp_data_file = 'Sommer2006_data';
cutoff = 0.5; % 50th threshold percentile
plot_region_name = 'FDI_rep_inds';
som_ind = 5; % 1 - exp mean, 5 - exp median
x_lim = [45 90]; % A/us - exp
y_lim = []; 
fig_dims = [9.0 9.5];% mean - [9.2 9.5] 
fig_dims2 = [9.0 10.5]; 
font_size = 12;
font_name = 'Times';
colors = jet(length(cell_ids)); colors(4,:) = [0.5961 0.3059 0.6392]; 
markers = {'^','^','o','o','d','d'}; % mPA,mAP,hPA,hAP,bAP,bPA
m_size = 6; 
units = '(A / \mu s)';
%% Load data
% FDI representation
ROI_filename = fullfile(mat_dir,'gen_figures',[plot_region_name '.mat']); 
ROI_data = load(ROI_filename);
ROIi = ROI_data.ROIi; 
% experimental data from Sommer et al., 2006
som_mt_data = load(fullfile(mat_dir,'gen_figures',[exp_data_file '.mat'])); 
som_mt = som_mt_data.som_mt; 
som_mt = quantile(som_mt,cutoff); % [APm PAm APh PAh APb PAb]
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
% Reshape data
num_layers = length(cell_ids); 
threshEs_all = cell(num_layers,num_models);
cell_model_names_i_all = cellfun(@(x) cellModelNames(x),cell_ids,'UniformOutput',0); 
for j = 1:num_models
   for i = 1:num_layers
       cell_model_names_i = cell_model_names_i_all{i}; % cell names in layer
       [~,~,data_inds] = intersect(cell_model_names_i,cell_model_names); % get indices of layer cells in threshEs
       if ~isempty(data_inds)
           data_layer =  cell2mat(threshEs_all_models{j}(data_inds));
           threshEs_all{i,j} = data_layer(ROIi{i},:); % combine all thresholds in layer
           threshEs_all{i,j} = threshEs_all{i,j}(:); % collapse into single vector
       else
           threshEs_all{i,j} = 1e6; % placeholder large number to allow UniformOutput in cellfun below, no cells in this layer
       end
   end
end
%% Analysis
% Extract cutoff value from ROI thresholds 
model_thresh = cellfun(@(x) quantile(x,cutoff),threshEs_all,'UniformOutput',1); % numLayers x numModels
model_thresh(model_thresh==1e6) = []; % clear empty elements
%% Plot a
fig = figure('Color','w'); 
fig.Units = 'centimeters';
fig.Position(3:4) = fig_dims; 
for i = 1:size(model_thresh,1)
    for j = 1:num_models
        if mod(j,2) == 0 % even: AP, leave open
            plot(som_mt(j),model_thresh(i,j),'Marker',markers{j},'MarkerFaceColor','none',...
            'MarkerEdgeColor',colors(i,:),'LineStyle','none','MarkerSize',m_size);
        else % odd: PA, fill
            plot(som_mt(j),model_thresh(i,j),'Marker',markers{j},'MarkerFaceColor',colors(i,:),...
                'MarkerEdgeColor',colors(i,:),'LineStyle','none','MarkerSize',m_size);
        end
        hold on; 
    end   
end
% get regression lines
layer_fits = zeros(size(model_thresh,1),2); 
rsq = zeros(size(model_thresh,1),1); 
ps = zeros(size(model_thresh,1),1); 
xexp = som_mt';
for i = 1:size(model_thresh,1)    
    ymodel = model_thresh(i,:)';
    layer_fits(i,:) = polyfit(xexp,ymodel,1);
    % calc R^2
    yfit = polyval(layer_fits(i,:),xexp);
    [Ri,pi] = corrcoef(ymodel,yfit);
    rsq(i) = Ri(1,2)^2; ps(i) = pi(1,2);    
    plot(xexp,yfit,'Color',colors(i,:));    
end
ax = gca;
if ~isempty(x_lim)
   ax.XLim = x_lim; 
end
if ~isempty(y_lim)
    ax.YLim = y_lim;
end
ax.FontSize = font_size;
ax.FontName = font_name; 
ax.Box = 'off';
ax.YGrid = 'on';
ax.YColor = 'k';
ax.XColor = 'k';
if som_ind == 1 % mean
    xlabel(sprintf('Mean experimental MT %s',units)); 
elseif som_ind == 5 % median
    xlabel(sprintf('Median experimental MT %s',units)); 
end
if som_ind == 1
    ylabel(sprintf('Mean model population threshold %s',units)); 
else
    if cutoff == 0.5
        ylabel(sprintf('Median model population threshold %s',units)); 
    else
        ylabel(sprintf('Model population threshold - %g %% %s',cutoff*100,units)); 
    end
end
layer_names = {'L1','L2/3','L4','L5','L6'}; 
fprintf('Cutoff = %.3f\n',cutoff); 
for i = 1:size(model_thresh,1)
   fprintf('%s - Rsq = %.4f, p = %.4f\n',layer_names{i},rsq(i),ps(i));  
end
%% get r^2 for all cutoffs for b
cutoffs_all = [0.025,0.5,0.1,0.25,0.5];
num_cutoffs = length(cutoffs_all);
rsq_all = zeros(num_cutoffs,size(model_thresh,1));
ps_all = zeros(num_cutoffs,size(model_thresh,1));
for i = 1:num_cutoffs
    model_threshi = cellfun(@(x) quantile(x,cutoffs_all(i)),threshEs_all,'UniformOutput',1); % numLayers x numModels    
    model_threshi(model_threshi==1e6) = []; % clear empty elements
    som_mti = quantile(som_mt_data.som_mt,cutoff);
    som_mti = [som_mti(:,2),som_mti(:,1),som_mti(:,4),som_mti(:,3),som_mti(:,5:6)];
    xexp = som_mti'; 
    for j = 1:size(model_threshi,1)    
        ymodel = model_threshi(j,:)';
        layer_fits(j,:) = polyfit(xexp,ymodel,1);
        % calc R^2
        yfit = polyval(layer_fits(j,:),xexp);
        [Ri,pi] = corrcoef(ymodel,yfit);
        rsq_all(i,j) = Ri(1,2)^2; ps_all(i,j) = pi(1,2);           
    end
end
%% Plot b
bar_width = 0.8;
fig2 = figure('Color','w'); 
fig2.Units = 'centimeters';
fig2.Position(3:4) = fig_dims2; 
ingroup_space = 0.1;
btgroup_space = bar_width*2;
bar_col = colors; 
x = 0; % initialize
alphas = fliplr(1/(num_cutoffs):1/(num_cutoffs):1);
for i = 1:size(rsq_all,2)
    x = x(end) + btgroup_space + (0:(ingroup_space+bar_width):((ingroup_space+bar_width)*(num_cutoffs-1)));
    for j = 1:length(x)
        b = bar(x(j),rsq_all(j,i),'BarWidth',bar_width,...
            'FaceColor',bar_col(i,:),'FaceAlpha',alphas(j),'EdgeColor','k');
        hold on;
        if ps_all(j,i) > 0.05
            h = hatchfill2(b,'single'); 
        end
    end
end
drawnow;
ax = fig2.Children(1);
box(ax,'off');
ax.YGrid = 'on';
ax.FontSize = font_size;
ax.FontName = font_name;
ax.YColor = 'k'; ax.XColor =  'k';
ax.XTickLabel = {'L1','L2/3','L4','L5','L6'}; 
% ax.XTickLabel = {'L1 NGC','L2/3 PC','L4 LBC','L5 PC','L6 PC'}; 
ax.XLim = [ax.Children(end).XData - bar_width/2-0.5,ax.Children(1).XData + bar_width/2+0.5];  
ax.XTick = (btgroup_space+(bar_width+ingroup_space)*floor(num_cutoffs/2)):(btgroup_space+(bar_width+ingroup_space)*2*floor(num_cutoffs/2)):ax.XLim(2);
ax.XLabel.String = 'Model Layer';
ax.YLabel.String = 'R^{2}'; 
ax.XColor = 'k';
ax.YColor = 'k';
ax.YGrid = 'on';
ax.YLim = [0 1]; 
%% Save
if save_fig
    fig_fold = fullfile(mat_dir,'figures');
    fig_name = 'SuppFig8a';
    savefig(fig,fullfile(fig_fold,[fig_name '.fig']));
    print(fig,fullfile(fig_fold,[fig_name '.eps']),'-depsc2','-cmyk');
    fig_name2 = 'SuppFig8b';
    savefig(fig2,fullfile(fig_fold,[fig_name2 '.fig']));
    print(fig2,fullfile(fig_fold,[fig_name2 '.png']),'-dpng','-cmyk','-r600');    
end
