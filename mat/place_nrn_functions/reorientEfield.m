% Input Efield vectors calculated for translated/rotated cell, rotates back
% to z-oriented coordinates from column normal, outputs just E-field
% components in same order as NEURON (no coordinates necessary)
function E = reorientEfield(cell_normal,phi,Ecell)        
    z = [0,0,1]'; % Default alignment of somato-dendritic axis
    Rn = getRmatrix(cell_normal,z); % get rotation matrix between z-axis and cell_normal
    E = (Rn*Ecell')'; % just rotate to align cell_normal to z-axis
    phi = -phi*pi/180; % convert to radians, rotate opposite direction
    Rz = [cos(phi) -sin(phi) 0; % rotate about z CW phi degrees(not CCW)
           sin(phi) cos(phi) 0;
           0 0 1];
    E = (Rz*E')';        
end