function [data,ROI] = extractDataSlice(slices,data,layers,plane,plot_inds,ROI)
%EXTRACTTHRESHSLICE ... 
%  
%   Inputs 
%   ------ 
%   Optional Inputs 
%   --------------- 
%   Outputs 
%   ------- 

% AUTHOR    : Aman Aberra 
num_layers = length(layers); 
if nargin < 6
    fig = figure('units','normalized','outerposition',[0 0 1 1],'paperunits','in','paperposition',[0 0 16 9],'Color',[1 1 1],'PaperPositionMode','auto');                                  
    plotLayers(layers,1:num_layers,[0,0,0]);    
    plane3d = createPlane(plane(1:3),plane(4:6)); 
    p3d = drawPlane3d(plane3d); p3d.FaceAlpha = 0.5;
    ROI = [min(p3d.Vertices(:,1)),max(p3d.Vertices(:,1)),min(p3d.Vertices(:,2)),max(p3d.Vertices(:,2)),min(p3d.Vertices(:,3)),max(p3d.Vertices(:,3))]; 
    close(fig); 
end
cnt = 1;
% data should have one vector per cell element (single layer)
for i = 1:num_layers
    [c,inds] = clipPoints3d(layers(i).cell_origins,ROI);
    if length(data{i})> 1
        data{i} = data{i}(inds); 
        [cu,cia,~] = unique(c,'rows','stable'); % avoid duplicate points
        data{i} = data{i}(cia); % use same threshold points
        Fi = scatteredInterpolant(cu(:,1),cu(:,2),cu(:,3),data{i});
        % interpolate on actual layer surface, (even elements)
        data{i} = Fi(slices(2*plot_inds(cnt)).slice_line(:,1),slices(2*plot_inds(cnt)).slice_line(:,2),slices(2*plot_inds(cnt)).slice_line(:,3)); 
        cnt = cnt + 1; 
    else
       data{i} = [];        
    end    
end
if num_layers > 1
    fprintf('Interpolated data on %g slice lines\n',num_layers); 
end

end