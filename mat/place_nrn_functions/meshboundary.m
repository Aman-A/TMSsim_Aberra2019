% input struct p with fields: 'vertices' and 'faces' 
% outputs:
%   b: all boundaries in Nb x 2 array
%   vb: boundary points in (Nb-1)*3 array for plotting, e.g.
% 
% where Nb is number of boundary edges
function [b,vb] = meshboundary(p)
v = p.vertices;
f = p.faces; 

% [~,ed] = computeMeshEdges(f); % function in geom3d/mesh3d
ed = sort([f(:,1) f(:,2); f(:,2) f(:,3); f(:,3) f(:,1)]')';
% determine uniquess of edges
[e,~,jx] = unique(ed,'rows'); % edges in only 1 face
% determine counts for each unique edge
counts = accumarray(jx(:),1); 
b = e(counts==1,:); 

vb1 = v(b(:,1),:); % Nb x 3 - coordinates of 1st boundary nodes
vb2 = v(b(:,2),:); % Nb x 3 - coordinates of 2nd boundary nodes
vb = zeros(3,(length(vb1-1)*3)); 
vb(:,1:3:end) = vb1';
vb(:,2:3:end) = vb2';
vb(:,3:3:end) = nan(3,length(vb1-1)); 

% figure; plot3(vb(1,:),vb(2,:),vb(3,:)); % simpler method
% figure; surface([vb(1,:);vb(1,:)],[vb(2,:);vb(2,:)],[vb(3,:);vb(3,:)],...
%     'FaceColor','flat','EdgeColor','k','LineWidth',1);