% plot 2d slice data
% input:
%   slices: struct array with numLayers*2+1 elements, includes boundaries
%   between layers and outer boundaries (GM and WM)
%   data: numLayers x 1 cell array with Nvertices x 1 arrays in each
%   element
function plotSlices(slices,color_data)
maxpts = length(slices(1).slice_line); 
num_layers = (length(slices)-1)/2; % don't count GM and WM boundaries
if nargin < 2
    color_data = repmat({linspace(0,1,maxpts)},num_layers,1);             
end
start_ind = 1;
for j = 1:num_layers
    i = j*2;     
%     plot3(slices(i).slice_line(:,1),slices(i).slice_line(:,2),slices(i).slice_line(:,3),'--k'); 
    l1 = slices(i-1).slice_line'; l2 = slices(i+1).slice_line';
    %        color_dat = i*ones(length(l1)*2,1);    
    color_dataj = color_data{j};
    if isempty(color_dataj)
       face_color = 'none';
       color_dataj = zeros(1,maxpts);
    else
        face_color = 'interp';
    end
    if isrow(color_dataj)
        color_dataj = [color_dataj(start_ind:end), fliplr(color_dataj(start_ind:end))];        
    else
        color_dataj = [color_dataj(start_ind:end); flipud(color_dataj(start_ind:end))];
    end
    patch([l1(1,start_ind:end),fliplr(l2(1,start_ind:end))],[l1(2,start_ind:end),fliplr(l2(2,start_ind:end))],[l1(3,start_ind:end),fliplr(l2(3,start_ind:end))],...
        color_dataj,'FaceColor',face_color,'LineWidth',1);  hold on; 
end
end