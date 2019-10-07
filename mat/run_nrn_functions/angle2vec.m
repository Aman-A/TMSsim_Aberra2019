% generates components of normal vector given a theta and phi
% where theta is polar angle, phi is azimuthal angle in degrees
function [x, y, z] = angle2vec(theta,phi,mag)  
    if nargin < 3
       mag = 1;  
    end        
    theta = theta*pi/180; phi = phi*pi/180; % convert to radians
    % Calculate components of vector based on input angles
    x = mag.*sin(theta).*cos(phi); 
    y = mag.*sin(theta).*sin(phi); 
    z = mag.*cos(theta);    
end