function slices = makeSlices(layers,MeshROI,plane,start_dim,plot_fig,uvROI)
% MAKESLICES extract 2d slice and relevant threshold data from 3D layers and data 
% input:
%   layers: numLayers x 1 structure array
%   data: numLayers x 1 cell array with Nvertices x 1 arrays in each
%   element
%   plane: [x0 y0 z0 Nx Ny Nz] plane defined by point and normal
%   layer_nums: vector with indices of layers to plot
% MeshROI: struct with surfaces field - MUST have GrayMatter/WhiteMatter mesh surfaces 
%   plot_fig: 1 or 0 to indicate whether to plot the slices on 1 figure
if nargin < 6
   plane_filt = 0;
elseif isempty(uvROI)
    plane_filt = 0;
else % uvROI given
    plane_filt = 1;    
    %uvROI = [0 0.6 0 1]; % [umin umax, vmin vmax]
end
if nargin < 5
   plot_fig = 1; 
end
if nargin < 4
   start_dim = 2;  
end
if nargin < 3    
    plane = [-35.44 -16.46 42.82 -0.3 0.7 0]; % [x0, y0, z0, nx, ny, nz]
end
numLayers = length(layers); 
slices(numLayers+2).slice_line = []; % lines from intersection of plane with layer
% Add GM/WM slice lines
slices(1).slice_line = intersectPlaneMeshSurf(plane,MeshROI.surfaces.GrayMatter,start_dim);
slices(numLayers+2).slice_line = intersectPlaneMeshSurf(plane,MeshROI.surfaces.WhiteMatter,start_dim);
for i = 2:numLayers+1
   layer_surf = layers(i-1).surface; %1 - numLayers   
   slice_linei = intersectPlaneMeshSurf(plane,layer_surf,start_dim); 
   slices(i).layer = layers(i-1);   
    slices(i).slice_line = slice_linei; 
end
slices(1).layer = MeshROI.surfaces.GrayMatter;
slices(numLayers+2).layer = MeshROI.surfaces.WhiteMatter;
%% Filter within plane
if plane_filt    
    fig = figure; 
    plot3(slices(1).slice_line(:,1),slices(1).slice_line(:,2),slices(1).slice_line(:,3));
    hold on; 
    plot3(slices(numLayers+2).slice_line(:,1),slices(numLayers+2).slice_line(:,2),slices(numLayers+2).slice_line(:,3));
    view(plane(4:6)); axis equal;
    p3d = drawPlane3d(createPlane(plane(1:3),plane(4:6)));
    if length(p3d.Vertices)==6
        p3d_verts = p3d.Vertices([1,2,4,5],:);
    else
       p3d_verts = p3d.Vertices; 
    end
    close(fig); 
    Op = [p3d_verts(4,:)-p3d_verts(1,:);p3d_verts(2,:)-p3d_verts(1,:)];
    Opn = bsxfun(@rdivide,Op,sqrt(sum(Op.^2,2))); %orthonormal basis for plane ROI
    Co = p3d_verts(1,:); % position of orthonormal basis origin (corner)    

    oROIp = [uvROI(1)*norm(Op(1,:)) uvROI(2)*norm(Op(1,:)) uvROI(3)*norm(Op(2,:)) uvROI(4)*norm(Op(2,:))];
%     fprintf('**Filtering within plane**\n')
    for i = 1:length(slices)
       slice_linei = slices(i).slice_line;
       slice_lineip = (Opn*bsxfun(@minus,slice_linei,Co)')'; % convert to plane basis (u,v)
       slice_lineip = clipPoints3d([slice_lineip,zeros(length(slice_lineip),1)],[oROIp 0 0]); % filter in plane
       slice_lineip = slice_lineip(:,1:2); % remove z 
       slice_linei = bsxfun(@plus,Co',Opn'*slice_lineip')'; % convert back to original basis
       slices(i).slice_line = slice_linei;
       slices(i).oROIp = oROIp; 
    end
    %% Make sure points in curves are ordered
    for i = 2:numLayers+2
       slice_linei = slices(i).slice_line;
       slice_lineip = (Opn*bsxfun(@minus,slice_linei,Co)')'; % convert to plane basis (u,v)
       [u,v] = points2contour(slice_lineip(:,1),slice_lineip(:,2),1,'cw');
       slices(i).slice_line = bsxfun(@plus,Co',Opn'*[u;v])'; % convert back to original basis
    end
end
%% Resample lines using spline curves to use max number of points for all curves
maxpts = max(arrayfun(@(x) size(x.slice_line,1),slices,'UniformOutput',1));
% maxpts = 1000; 
for i = 1:numLayers+2
   splinec = cscvn(slices(i).slice_line'); % create spline curve
   interv = fnbrk(splinec,'interval'); % get min and max
   t = linspace(interv(1),interv(2),maxpts); % evenly spaced points 0-1
   slices(i).slice_line = fnval(splinec,t)'; % resample spline curve evenly with max # of pts
end
%%
if plot_fig
   figure('units','normalized','outerposition',[0 0 1 1],'paperunits','in','paperposition',[0 0 16 9],'Color',[1 1 1],'PaperPositionMode','auto');                              
   plotSlices(slices); 
   view(plane(4:6)); axis equal; 
end

end