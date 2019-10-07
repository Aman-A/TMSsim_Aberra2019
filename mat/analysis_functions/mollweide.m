% This function converts latitude and longitude infomation to cooridates of 
% Mollweide projection.
% 
% [X,Y]=mollweide(LAT,LONG,Center_Meridian,Radius)
% 
% Inputs:
% LAT:              latitudes to be converted in degrees. Can be matrix or
%                   vector.
% LONG:             longitudes to be converted in degrees. Must be same 
%                   size as LAT.
% CenterMeridian:   where projection is centered, in degrees. Default at 0.
% Radius:           to specify size of map, so that map area equals surface  
%                   area of sphere with given radius
% 
% Outputs:
% X,Y: map cooridates. Same size as LAT, and LONG.

function [X,Y]=mollweide(LAT,LONG,Center_Meridian,Radius)

% move projection center to center meridian
if nargin<2
    error('Not enough input arguments.');
end
if ~isequal(size(LAT), size(LONG))
    error('Dimensions of the latitude and longitude inputs don''t match.');
end
ind_lat=LAT(:)>180;
while any(ind_lat(:))
    LAT(ind_lat)=LAT(ind_lat)-360;
    ind_lat=LAT(:)>180;
end
ind_lat=LAT(:)<-180;
while any(ind_lat(:))
    LAT(ind_lat)=LAT(ind_lat)+360;
    ind_lat=LAT(:)<-180;
end
ind_lat=LAT(:)>90;
if any(ind_lat(:))
    LAT(ind_lat)=LAT(ind_lat)-180;
    LONG(ind_lat)=LONG(ind_lat)+180;
end
ind_lat=LAT(:)<-90;
if any(ind_lat(:))
    LAT(ind_lat)=LAT(ind_lat)+180;
    LONG(ind_lat)=LONG(ind_lat)+180;
end
ind_long=LONG(:)>180;
while any(ind_long)
    LONG(ind_long)=LONG(ind_long)-360;
    ind_long=LONG(:)>180;
end
ind_long=LONG(:)<-180;
while any(ind_long)
    LONG(ind_long)=LONG(ind_long)+360;
    ind_long=LONG(:)<-180;
end

if nargin==2
    NewLONG=LONG;
elseif nargin>=3            
    NewLONG=(LONG-Center_Meridian);
end
NewLONG(NewLONG>180)=NewLONG(NewLONG>180)-360;
NewLONG=NewLONG*pi/180;

% can specify size of map, so that area equals sphere of givenradius
if nargin<=3                
    Radius=1;
end

% first step of iteration
THETA_N=LAT*pi/180;
THETA_N1=THETA_N- (2*THETA_N+sin(2*THETA_N)-pi*sin(LAT*pi/180)) ./ (2+2*cos(2*THETA_N));
% adjust points where division of zero occured
ind=isnan(THETA_N1);
THETA_N1(ind)= LAT(ind)*pi/180;

% iteration
while max(max(abs(THETA_N-THETA_N1)))>0.001
    THETA_N=THETA_N1;
    THETA_N1=THETA_N- (2*THETA_N+sin(2*THETA_N)-pi*sin(LAT*pi/180)) ./ (2+2*cos(2*THETA_N));
    ind=isnan(THETA_N1);
    THETA_N1(ind)= LAT(ind)*pi/180;
end

ind=isnan(THETA_N1);
THETA_N1(ind)= LAT(ind)*pi/180;

% conversion to map cooridinates in X and Y
X=2*sqrt(2)/pi.*NewLONG.*cos(THETA_N1)*Radius;
Y=sqrt(2).*sin(THETA_N1)*Radius;
end