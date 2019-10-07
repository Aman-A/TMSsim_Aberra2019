function data_med = plotDataSlice(data,cell_model_names,cell_ids,layers,slices,plane,varargin)
%PLOTDATASLICE ... 
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
in.cutoff = [];
in.norm_to_max = 0;
in = sl.in.processVarargin(in,varargin);
num_layers = length(layers); 
data_med = cell(num_layers,1);
plot_inds = 1:num_layers; 
% Get median in each layer
for i = 1:num_layers
    cell_model_names_i = cellModelNames(cell_ids{i}); % cell names in layer
    [~,~,data_inds] = intersect(cell_model_names_i,cell_model_names); % get indices of layer cells in threshEs
    if ~isempty(data_inds)
        data_med{i} = median(cell2mat(data(data_inds)),2); % average each position across cell models
    else
        data_med{i} = 1e6; % placeholder large number to allow UniformOutput in cellfun below, no cells in this layer
        plot_inds(plot_inds==i) = []; % remove this index
    end
end
% Extract data in slice
data_med = extractDataSlice(slices,data_med,layers,plane,plot_inds); 
if in.norm_to_max % set max to 1
   max_data = max(cellfun(@max,data_med,'UniformOutput',1)); 
   data_med = cellfun(@(x) x/max_data,data_med,'UniformOutput',0); 
end
% Plot on current axis
plotSlices(slices,data_med); 
if ~isempty(in.cutoff)
    for j = 1:5
        ax = gca;
        pi = ax.Children(1+length(ax.Children)-j); % l1 to l5
        dataj = data_med{j};
        dataj = [dataj; flipud(dataj)];
        pi.FaceVertexCData = dataj; 
        above_thresh = dataj > in.cutoff; 
        pi.FaceVertexCData(above_thresh) = 1000;
    end
end