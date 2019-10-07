function plotDistsModelROIs(data_all,ROIi_all,cellids,cell_model_names,colors)
%PLOTBOXPLOTS Plots boxplots for models across model, subgrouped by layer
%  
%   Inputs 
%   ------ 
%   Optional Inputs 
%   --------------- 
%   Outputs 
%   ------- 
%   Examples 
%   --------------- 

% AUTHOR    : Aman Aberra 
if nargin < 4
   colors = jet(length(cellids)); 
end
num_models = length(data_all); 
num_rois = length(ROIi_all); 
%%
% combine data within each layer
data = cell(num_rois,num_models); % rows - layer, columns - model
cell_model_names_i = cellModelNames(cellids); % cell names 
[~,~,data_inds] = intersect(cell_model_names_i,cell_model_names); % get indices of layer cells in threshEs
for j = 1:num_models    
    for i = 1:num_rois
        data_roi =  cell2mat(data_all{j}(data_inds));
        data{i,j} = data_roi; % combine all thresholds in layer        
    end
end
% Extract data within desired region
roi_names = strcat('R',strsplit(num2str(1:num_rois)))';
layer_ids_all = {}; 
for j = 1:num_models 
    for i = 1:num_rois  % layer specific ROI indices       
        data{i,j} = data{i,j}(ROIi_all{i},:);
        data{i,j} = data{i,j}(:); % turn into 1D array
        layer_ids_all = [layer_ids_all; repmat(strcat(roi_names(i), '_m',num2str(j)),length(data{i}),1)];           
    end
end
% threshEs_all = threshEs_all'; 
data = cell2mat(data(:)); % collapse into single vector
%% Plot
boxplot(data,layer_ids_all,'Symbol','','Width',0.5); 
hold on; 
ax = gca;
bp = ax.Children(1); 
% yaxlim = ax.YLim;
yaxlim = [0 1000]; 
if num_models > 1 % plot dividing lines
    ax.XTick = ((num_rois+1)/2:num_rois:(num_models*num_rois-1));    
    ax.XTickLabel = strcat('M',strsplit(num2str(1:num_models)))';
    plot(repmat(num_rois:num_rois:num_rois*(num_models-1),2,1)+0.5,...
            [1e-3*min(data)*ones(1,num_models-1);yaxlim(2)*ones(1,num_models-1)],...
            'Color','k','LineWidth',ax.LineWidth)
    leg_objs = bp.Children(num_models*num_rois+1:2*num_models*num_rois);
    leg_objs = leg_objs(num_rois:-1:1); 
    legend(leg_objs,roi_names);
end
colors = flipud(colors); 
for i = 1:num_rois:(7*num_models*num_rois)
   for j = 0:num_rois-1
       bp.Children(i+j).Color = colors(j+1,:);
       bp.Children(i+j).LineWidth = 1;
       bp.Children(i+j).LineStyle = '-';
   end
end
% remove outliers
delete(ax.Children(end).Children(1:num_rois*num_models)) % delete outliers
