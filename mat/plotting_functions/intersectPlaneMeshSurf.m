function slice_line = intersectPlaneMeshSurf(plane,surface,start_dim)
% Returns ordered points defining line of intersection between plane and
% mesh surface
% Input:
% plane - [x0 y0 z0 Nx Ny Nz] - point on plane and normal vector
% surface structure with fields:
%   surface.vertices
%   surface.faces
% start_dim = -1 xmin, +1 xmax, -2 ymin, +2 ymax, -3 zmin, +3 zmax
% defines dimension to start line order from 
% Adapted from Mehmet Ozturk post on MATLAB file exchange:
% https://www.mathworks.com/matlabcentral/fileexchange/22377-intersectplanesurf
% requires planeNormal2Imp3d from 
% http://people.sc.fsu.edu/~jburkardt/m_src/geometry/geometry.html
% (plane_normal2imp_3d.m)
surface.faces = surface.faces';
surface.vertices = surface.vertices';
[A,B,C,D] = planeNormal2Imp3d(plane(1:3),plane(4:6));
segment_start=nan(3,round(size(surface.faces,2)/2)); 
segment_finish=segment_start; 
count=1; 
%%
% figure;
for s=1:size(surface.faces,2)
    t(:,1)=surface.vertices(:,surface.faces(1,s));
    t(:,2)=surface.vertices(:,surface.faces(2,s));
    t(:,3)=surface.vertices(:,surface.faces(3,s));
    [ num_int, pi ] = planeImpTriangleInt3d ( A, B, C, D, t );
    if num_int==2
        segment_start(:,count)=pi(:,1);
        segment_finish(:,count)=pi(:,2);
        count=count+1;
%         hold on,plot3(pi(1,:),pi(2,:),pi(3,:))
    end
end
segment_start(:,all(isnan(segment_start),1))=[]; % remove unused poritons 
segment_finish(:,all(isnan(segment_finish),1))=[]; % remove unused poritons
%%
thr=1; nol=1; lin = {}; 
while ~isempty(segment_start)
    lin{nol}=[segment_start(:,1) segment_finish(:,1)];
    segment_start(:,1)=[];
    segment_finish(:,1)=[];
    while 1
        testDistStart1=sum((lin{nol}(1,end)-segment_start(1,:)).^2 + ...
            (lin{nol}(2,end)-segment_start(2,:)).^2 + ...
            (lin{nol}(3,end)-segment_start(3,:)).^2,1); % distance from newest segment, all remaining segment starts
        testDistStart2=sum((lin{nol}(1,end)-segment_finish(1,:)).^2 + ...
            (lin{nol}(2,end)-segment_finish(2,:)).^2 + ...
            (lin{nol}(3,end)-segment_finish(3,:)).^2,1);% distance from newest segment, all remaining segment ends
        [minDist1, best_ind1]=min(testDistStart1);
        [minDist2, best_ind2]=min(testDistStart2);   
        if isempty(minDist1)
           minDist1 = thr+1; 
        end
        if isempty(minDist2)
            minDist2 = thr+1;
        end
        if minDist1<thr || minDist2<thr
            if minDist1<minDist2
                lin{nol}=[lin{nol} segment_finish(:,best_ind1)];
                segment_start(:,best_ind1)=[];
                segment_finish(:,best_ind1)=[];
            else
                lin{nol}=[lin{nol} segment_start(:,best_ind2)];
                segment_start(:,best_ind2)=[];
                segment_finish(:,best_ind2)=[];
            end
        else
            nol=nol+1;
            break
        end
    end
end
% Combine disconnected line segments   
% get points on segments closest to desired starting dimension
if start_dim > 0  
%     start_pt = [0;0;35]; % try ordering by distance from pont
%     [max_dim_pt,max_dim_ind] = cellfun(@(x) min(sqrt(sum(bsxfun(@minus,x,start_pt).^2,1))),lin);
    [max_dim_pt,max_dim_ind] = cellfun(@(x) max(x(start_dim,:)),lin); 
    lin_lens = cellfun(@length,lin);
    flip_cells = max_dim_ind ~= 1;   
    remove_cells = flip_cells & max_dim_ind ~= lin_lens; 
    [~,cell_order] = sort(max_dim_pt,'descend');    
%     [~,cell_order] = sort(max_dim_pt,'ascend');    
else
    [min_dim_pt,min_dim_ind] = cellfun(@(x) min(x(abs(start_dim),:)),lin); 
    lin_lens = cellfun(@length,lin);
     flip_cells = min_dim_ind ~= 1;   
     remove_cells = flip_cells & min_dim_ind ~= lin_lens; 
     [~,cell_order] = sort(min_dim_pt,'ascend');
end
slice_line = [];
for i = cell_order
    if ~remove_cells(i)
        if flip_cells(i) 
%             add_line = fliplr(lin{i});
%             add_line_r = sqrt(sum(bsxfun(@minus,add_line,add_line(:,1)).^2)); % get dist of each point from starting point
%             [~,pts_order] = sort(add_line_r,'descend');            
%             slice_line = [slice_line add_line(:,pts_order)];
            slice_line = [slice_line fliplr(lin{i})];            
        else
           slice_line = [slice_line lin{i}]; 
        end
    end
end
slice_line = slice_line'; % column vectors
slice_line = unique(slice_line,'rows','stable'); 
end