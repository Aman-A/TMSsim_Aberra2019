% Enter two vectors, A, B, returns rotation matrix R that rotates A to B
function R = getRmatrix(A,B)
    if (isrow(A))
        A = A';
    end
    if (isrow(B))
       B = B'; 
    end
    A = A./norm(A); B = B./norm(B); % convert to unit vectors 
    d = dot(A,B);
    err = 1e-5; % minimum allowable deviation in dot product
    if abs(d - 1) < err % vectors are parallel, rotation is identity
        R = eye(3); 
    elseif abs(d - -1) < err % vectors are anti-parallel, flip vector
        R = -eye(3); 
    else
        c = norm(cross(A,B));
        G = [ d, -c, 0;
            c, d, 0;
            0, 0, 1];
        F = inv([A , ( (B - d*A)/norm(B - d*A) ), cross(B,A) ]);    
        R = F\G*F;
    end
end