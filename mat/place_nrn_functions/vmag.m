function vector_mag = vmag(vector)    
% input N x 3  or N x 2 column vector, returns the
% magnitudes in N x 1 column vector
% Ex:
% vector_mag = vmag([1 2 3; 4 5 6]);
    vector_mag = sqrt(vector(:,1).^2 + vector(:,2).^2 + vector(:,3).^2);    
end