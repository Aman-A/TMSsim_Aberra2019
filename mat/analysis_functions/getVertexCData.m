% Input face_data with dimension N_faces x j and surface structure with
% fields surface.faces and surface.vertices
% outputs vertex_data which contains data points for every vertex
function vertex_data = getVertexCData(face_data,surface,cell_origins)
    numVertices = size(surface.vertices,1); 
    vertices = surface.vertices;
    faces = surface.faces;
    vertex_data = zeros(numVertices,size(face_data,2)); 
    if isrow(face_data)
       face_data = face_data';  % convert to column vector
    end
    for i = 1:numVertices        
        % find all element faces adjacent to vertex_i
        r_i = vertices(i,:); % position of vertex_i
        face_inds = find(sum(faces == i,2)); % collapse into vector, get element indices        
        face_data_i = face_data(face_inds,:); % face data adjacent to vertex
        cell_origins_i = cell_origins(face_inds,:); 
        dr = cell_origins_i - repmat(r_i,length(face_inds),1); % vectors from cell_normal to vertex_i
        distances = sqrt(dr(:,1).^2+dr(:,2).^2+dr(:,3).^2); 
        inv_distances = abs(1./distances); %1/ri - weights for weighted average
        vertex_data(i,:) = inv_distances'*face_data_i/sum(inv_distances); % vertex is weighted average of thresholds of adjacent face centers      
    end
end