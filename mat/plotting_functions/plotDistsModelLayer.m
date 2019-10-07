function plotDistsModelLayer(data_all,ROIi,cell_ids,cell_model_names,colors)
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
   colors = jet(length(cell_ids)); 
end
num_models = length(data_all); 
num_layers = length(cell_ids); 
%%
% combine data within each layer
data = cell(num_layers,num_models); % rows - layer, columns - model
for j = 1:num_models
    for i = 1:num_layers
        cell_model_names_i = cellModelNames(cell_ids{i}); % cell names in layer
        [~,~,data_inds] = intersect(cell_model_names_i,cell_model_names); % get indices of layer cells in threshEs
        if ~isempty(data_inds)
            data_layer =  cell2mat(data_all{j}(data_inds));
            data{i,j} = data_layer; % combine all thresholds in layer
        else
            data{i,j} = 1e6; % placeholder large number to allow UniformOutput in cellfun below, no cells in this layer
        end
    end
end
% Extract data within desired region
layer_names = strcat('L',strsplit(num2str(1:num_layers)))';
layer_ids_all = {}; 
for j = 1:num_models 
    for i = 1:num_layers  % layer specific ROI indices       
        data{i,j} = data{i,j}(ROIi{i},:);
        data{i,j} = data{i,j}(:); % turn into 1D array
        layer_ids_all = [layer_ids_all; repmat(strcat(layer_names(i), '_m',num2str(j)),length(data{i}),1)];           
    end
end
% threshEs_all = threshEs_all'; 
data = cell2mat(data(:)); % collapse into single vector
%% Plot
% boxplot(data,layer_ids_all,'Symbol','','Width',0.5); 
boxplot(data,layer_ids_all,'Width',0.5); 
hold on; 
ax = gca;
bp = ax.Children(1); 
yaxlim = ax.YLim;
if num_models > 1 % plot dividing lines
    ax.XTick = ((num_layers+1)/2:num_layers:(num_models*num_layers-1));    
    ax.XTickLabel = strcat('M',strsplit(num2str(1:num_models)))';
    plot(repmat(num_layers:num_layers:num_layers*(num_models-1),2,1)+0.5,...
            [1e-3*min(data)*ones(1,num_models-1);yaxlim(2)*ones(1,num_models-1)],...
            'Color','k','LineWidth',ax.LineWidth)
    leg_objs = bp.Children(num_models*num_layers+1:2*num_models*num_layers);
    leg_objs = leg_objs(num_layers:-1:1); 
    legend(leg_objs,layer_names);
end
colors = flipud(colors); 
for i = 1:num_layers:(7*num_models*num_layers)
   for j = 0:num_layers-1
       bp.Children(i+j).Color = colors(j+1,:);
       bp.Children(i+j).LineWidth = 1;
       bp.Children(i+j).LineStyle = '-';
   end
end
% remove outliers
delete(ax.Children(end).Children(1:length(cell_ids)*length(data_all))) % delete outliers
