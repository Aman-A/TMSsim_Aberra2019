% give vertices and faces vector and inds corresponding to face indices
% Outputs surface areas of these faces
function areas = getPatchAreas(vertices,faces,inds)   
    if nargin < 3
       inds = 1:length(faces);  
    end
    % extract x,y,z coordinates associated with element
    vi1 = vertices(faces(inds,1),:);
    vi2 = vertices(faces(inds,2),:);
    vi3 = vertices(faces(inds,3),:);
    areas = abs(vi1(:,1).*(vi2(:,2) - vi3(:,2)) +...
        vi2(:,1).*(vi3(:,2) - vi1(:,2)) +...
        vi3(:,1).*(vi1(:,2) - vi2(:,2)))/2;
    
end