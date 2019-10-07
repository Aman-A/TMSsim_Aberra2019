% This function plots Mollweide projection of threshold orientation map.
% Necessary inputs:
%   threshEs:           vector of threshold values in V/m; size: 1 x num_thresholds.
%                       num_thresholds= 2+ (num_theta-2)* num_phi
%   theta:              vector of theta angles in degrees; size: 1 x num_thresholds
%   phi:                vector of phi angles in degrees; size: 1 x num_thresholds

% 
% Optional inputs
%   cell_model_name:    string containing name of cell.
%   display_scale:      string, either 'lin' or 'log'. Default: 'lin'
%   interp_interval:    interval of interpolation in degrees. Default: 1 degree intervals.
%   center_at_min:      if 1, projection centered at lowest threshold direction 
%   center_theta:       if center_at_min=0, theta of projection center. Default: 90 degree.
%   center_phi:         if center_at_min=0, phi of projection center. Default: 0 degree.
%   colorbar_auto:      if 1, colorbar is automatically set to max and min of thresholds   
%   colorbar_max:       if colorbar_auto=0, sets maximum of colorbar
%   colorbar_min:       if colorbar_auto=0, sets minimum of colorbar
%   colorbar_norm:      if colorbar_norm=1, normalizes minimum of data to 1     
%   flat:               if 1, plots 2D, otherwise 3D map (in 2D view but viewable via rotation)
%   display_grid:       if 1, display grid lines for Mollweide map in 15 degrees intervals
%   display_old_grid:   if 1, display grid lines of original sampling points
%   plot_min_point:      if 1, adds star at minimum data point
%   title_on:           if 1, display title
%   colorbar_on:        if 1, display colorbar
%   variable:           string, name of variable. Default: 'Threshold E'
%   units:              string, name of units. Default: 'V/m'
%   plot_in_grid:       if 1, changes settings for grid display
%   ax:                 axis handle to plot to, otherwise plots to gca
%   
function plotThreshMapMW(threshEs,theta,phi,varargin)
%% check input arguments and set defaults
if ~isvector(threshEs) || ~isvector(theta) || ~isvector(phi)
    error('Threshold, theta, and phi should be vectors. Plot aborted');
end
if length(threshEs)~=length(theta) || length(theta)~=length(phi) || length(threshEs)~=length(phi)
    error('Threshold, theta, and phi should be vectors of same length. Plot aborted');
end
if isrow(threshEs)
    threshEs=threshEs';
end
if isrow(theta)
    theta=theta';
end
if isrow(phi)
    phi=phi';
end
num_thresholds=length(threshEs);
num_theta = length(unique(theta)); 
num_phi = length(unique(phi)); 
if num_thresholds~= (2+(num_theta-2)* num_phi)
    error('Number of samples inconsistent with size of specified grid. Plot aborted');
end
num_lat=num_theta;
num_long=num_phi+1;

in.cell_model_name='test_cell';
in.display_scale='lin';
in.interp_interval = 0.5;
in.center_at_min = 0; 
in.center_theta = 90;
in.center_phi = -90; 
in.colorbar_auto = 1; % default auto colorbar
in.colorbar_max = []; 
in.colorbar_min = [];
in.colorbar_norm = 0; % default leave in absolute units
in.flat = 0; % default leave 3D
in.display_grid = 1; 
in.display_old_grid = 0;
in.plot_min_point = 1; 
in.title_on = 1; 
in.colorbar_on = 1;
in.variable = 'Threshold E';
in.units = 'V/m'; 
in.plot_in_grid = 0; 
in.ax = []; 
in = sl.in.processVarargin(in,varargin); 

display_cell_name=in.cell_model_name;
ind=(display_cell_name=='_');
display_cell_name(ind)=' '; %remove _ for proper title display with latex
if ~strcmp(in.display_scale,'lin') && ~strcmp(in.display_scale,'log')
    disp('Display scale invalid. Default scale: linear.');
    in.display_scale='lin';
end

if in.interp_interval<0
    in.interp_interval=0.5;
    disp('Negative interpolation interval provided. Default interval: 0.5.');
elseif in.interp_interval>15
    in.interp_interval=15;
    disp('Interpolation interval too large. Interval: 15.');
end
if ~in.center_at_min   
    flag_phi=0;flag_theta=0;
    while in.center_theta>360
        in.center_theta=in.center_theta-360;
        flag_theta=1;
    end
    while in.center_theta<0
        in.center_theta=in.center_theta+360;
        flag_theta=1;
    end
    if in.center_theta>180
        in.center_theta=360-in.center_theta;
        flag_theta=1;
        in.center_phi=in.center_phi+180;
        flag_phi=1;
    end
    while in.center_phi>180
        in.center_phi=in.center_phi-360;
        flag_phi=1;
    end
    while in.center_phi<-180
        in.center_phi=in.center_phi+360;
        flag_phi=1;
    end
    if flag_phi || flag_theta
        disp(['Projection center coordinates out of range. Coordinates adjusted set to (theta,phi)=(',num2str(in.center_theta),',',num2str(in.center_phi),').']);
    end
end
if ~in.colorbar_auto
    if isempty(in.colorbar_max)
        in.colorbar_auto=1;
        disp('No maximum value provided for manual colorbar mode. Default mode: auto.');
    elseif isempty(in.colorbar_min)
        in.colorbar_auto=1;
        disp('No minimum value provided for manual colorbar mode. Default mode: auto.');
    end
    if in.colorbar_max<in.colorbar_min 
        in.colorbar_auto=1;
        disp('Invalid range for manual colorbar mode. Default mode: auto.');
    end
end
if isempty(in.ax)
   in.ax = gca; 
end

eps=10^-4;
%% Plot threshold map
min_th = min(threshEs);
min_th_abs = min(abs(threshEs));
if in.colorbar_norm % normalize threshold values to minimum
   threshEs = threshEs/min_th_abs;
end
% ALL CAPS variables are grids
THETA=zeros(num_lat,num_long);
PHI=zeros(num_lat,num_long);
THRESHOLD=zeros(num_lat,num_long);

THETA(2:num_lat-1,1:num_long-1)=reshape(theta(2:end-1),num_long-1,num_lat-2)';
THETA(2:num_lat-1,num_long)=THETA(2:num_lat-1,1);
THETA(1,:)=theta(1);
THETA(num_lat,:)=theta(end);

PHI(2:num_lat-1,1:num_long-1)=reshape(phi(2:end-1),num_long-1,num_lat-2)';
PHI(1,:)=PHI(2,:);
PHI(num_lat,:)=PHI(num_lat-1,:);
ind=(PHI(:,1:num_long-1)>=180);
PHI(ind)=PHI(ind)-360;
[~,sort_ind]=sort(PHI(1,1:num_long-1));
PHI=[PHI(:,sort_ind),PHI(:,sort_ind(1))+360-eps];

THRESHOLD(2:num_lat-1,1:num_long-1)=reshape(threshEs(2:end-1),num_long-1,num_lat-2)';
THRESHOLD(1,:)=threshEs(1);
THRESHOLD(num_lat,:)=threshEs(end);
THRESHOLD=[THRESHOLD(:,sort_ind),THRESHOLD(:,sort_ind(1))];

THETA=[THETA(:,sort_ind),THETA(:,sort_ind(1))];
%% find center of projection if centered at lowest threshold

% max_thresh=max(THRESHOLD(:));
% get min
[min_temp,ind_temp]=min(THRESHOLD);
[min_thresh,ind_2]=min(min_temp);
ind_1=ind_temp(ind_2);
minLong=PHI(ind_1,ind_2);     %degrees
minTheta=THETA(ind_1,ind_2);
% get max
[max_temp,ind_max_temp] = max(THRESHOLD);
[max_thresh,ind_max_2] = max(max_temp);
ind_max_1 = ind_max_temp(ind_max_2);
maxLong = PHI(ind_max_1,ind_max_2);
maxTheta = THETA(ind_max_1,ind_max_2);
if in.center_at_min % only cetner azimuthal direction on min,leave center_theta=90 
    in.center_phi=minLong;     %degrees
    in.center_theta=90;
end   

%% Set coordinate grid for projection and display, origin at projection center
long=linspace(-180+eps,180-eps,360/in.interp_interval+1);
lat=linspace(-90+eps,90-eps,180/in.interp_interval+1);
[LONG,LAT]=meshgrid(long,lat);
% Map Projection
[X_map,Y_map]=mollweide(LAT,LONG);
front_bounds = [value2ind(-90,LONG(1,:)'),value2ind(90,LONG(1,:)')];
trans_bounds = [value2ind(-60,LAT(:,1)'),value2ind(60,LAT(:,1)')];
X_line = X_map(:,front_bounds); Y_line = Y_map(:,front_bounds);
X_linet = X_map(trans_bounds,:); Y_linet = Y_map(trans_bounds,:); 
%% Rotate new coordinates to original coordination of measurement
alpha=-(90-in.center_theta);   %rotation angle of Z axis

Z=cos((90-LAT)*pi/180);
X=sin((90-LAT)*pi/180).*cos(LONG*pi/180);
Y=sin((90-LAT)*pi/180).*sin(LONG*pi/180);
oldX=X*cos(alpha*pi/180)+Z*sin(alpha*pi/180);
oldZ=Z*cos(alpha*pi/180)-X*sin(alpha*pi/180);

%  _int means interpolation; applies to grid and threshold
PHI_int=atan2(Y,oldX)*180/pi+in.center_phi;
THETA_int=atan2(sqrt(oldX.^2+Y.^2),oldZ)*180/pi;
if in.center_phi>0
    ind=(PHI_int>180);
    PHI_int(ind)=PHI_int(ind)-360;
else
    ind=(PHI_int<-180);
    PHI_int(ind)=PHI_int(ind)+360;
end

%% Interpolation 
THRESHOLD_F=griddedInterpolant(THETA,PHI,THRESHOLD,'linear');
THRESHOLD_int=THRESHOLD_F(THETA_int,PHI_int);


%% Rotate original grid to new coordination, if display is needed
if in.display_old_grid || ~in.center_at_min || plot_sec_types
    alpha=(90-in.center_theta);
    
    temp_PHI=(PHI(:,1:end-1)-in.center_phi);
    if any(abs(temp_PHI(:))<eps)
        flag_on_grid=1;
        if in.center_phi<0
            ind=(temp_PHI>180);
            temp_PHI(ind)=temp_PHI(ind)-360;
            [~,sort_ind]=sort(temp_PHI(1,1:num_long-1));
            temp_PHI=[temp_PHI(:,sort_ind),temp_PHI(:,sort_ind(1))-360+eps];
        else
            ind=(temp_PHI)<-180;
            temp_PHI(ind)=temp_PHI(ind)+360;
            [~,sort_ind]=sort(temp_PHI(1,1:num_long-1));
            temp_PHI=[temp_PHI(:,sort_ind(1))+eps,temp_PHI(:,sort_ind(2:end)),temp_PHI(:,sort_ind(1))+360-eps];
        end
        temp_PHI=[temp_PHI(:,1:(num_long-1)/2),temp_PHI(:,(num_long+1)/2)-eps,...
            temp_PHI(:,(num_long+1)/2)+eps,temp_PHI(:,(num_long+3)/2:end)];
        temp_THETA=[THETA(:,sort_ind),THETA(:,sort_ind(1))];
        temp_THETA=[temp_THETA(:,1:(num_long-1)/2),temp_THETA(:,(num_long+1)/2),...
            temp_THETA(:,(num_long+1)/2),temp_THETA(:,(num_long+3)/2:end)];
        temp_THRESHOLD=[THRESHOLD(:,sort_ind),THRESHOLD(:,sort_ind(1))];
        temp_THRESHOLD=[temp_THRESHOLD(:,1:(num_long-1)/2),temp_THRESHOLD(:,(num_long+1)/2),...
            temp_THRESHOLD(:,(num_long+1)/2),temp_THRESHOLD(:,(num_long+3)/2:end)];        
    else
        flag_on_grid=0;
        if in.center_phi<0
            ind=(temp_PHI>180);
            temp_PHI(ind)=temp_PHI(ind)-360;
            [~,sort_ind]=sort(temp_PHI(1,1:num_long-1));
            temp_PHI=[-180*ones(num_lat,1)+eps,temp_PHI(:,sort_ind),180*ones(num_lat,1)-eps];
        else
            ind=(temp_PHI)<-180;
            temp_PHI(ind)=temp_PHI(ind)+360;
            [~,sort_ind]=sort(temp_PHI(1,1:num_long-1));
            temp_PHI=[-180*ones(num_lat,1)+eps,temp_PHI(:,sort_ind),180*ones(num_lat,1)-eps];
        end
        temp_PHI=[temp_PHI(:,1:(num_long+1)/2),zeros(num_lat,1)-eps,...
            zeros(num_lat,1)+eps,temp_PHI(:,(num_long+3)/2:end)];
        temp_THETA=[THETA(:,sort_ind(1)),THETA(:,sort_ind),THETA(:,sort_ind(1))];
        temp_THETA=[temp_THETA(:,1:(num_long+1)/2),kron([1,1],(temp_THETA(:,(num_long+1)/2)+...
            temp_THETA(:,(num_long+3)/2))/2),temp_THETA(:,(num_long+3)/2:end)];
        temp_THRESHOLD=[THRESHOLD(:,sort_ind(1)),THRESHOLD(:,sort_ind),THRESHOLD(:,sort_ind(1))];
        temp_THRESHOLD=[temp_THRESHOLD(:,1:(num_long+1)/2),kron([1,1],(temp_THRESHOLD(:,(num_long+1)/2)+...
            temp_THRESHOLD(:,(num_long+3)/2))/2),temp_THRESHOLD(:,(num_long+3)/2:end)];
        if plot_sec_types
           temp_SECTYPES=[SECTYPES(:,sort_ind(1)),SECTYPES(:,sort_ind),SECTYPES(:,sort_ind(1))];
           temp_SECTYPES=[temp_SECTYPES(:,1:(num_long+1)/2),kron([1,1],(temp_SECTYPES(:,(num_long+1)/2)+...
                temp_SECTYPES(:,(num_long+3)/2))/2),temp_SECTYPES(:,(num_long+3)/2:end)]; 
        end
    end
    Z=cos(temp_THETA*pi/180);
    X=sin(temp_THETA*pi/180).*cos(temp_PHI*pi/180);
    Y=sin(temp_THETA*pi/180).*sin(temp_PHI*pi/180);
    newX=X*cos(alpha*pi/180)+Z*sin(alpha*pi/180);
    newZ=Z*cos(alpha*pi/180)-X*sin(alpha*pi/180);
    
    PHI_new=atan2(Y,newX)*180/pi;
    THETA_new=atan2(sqrt(newX.^2+Y.^2),newZ)*180/pi;
    
    % Map Projection of original sampling grid
    [X_grid,Y_grid]=mollweide((90-THETA_new),PHI_new);
    
    % Grid correction
    if flag_on_grid
        if in.center_theta<=90
            kk=num_lat;
            while Y_grid(kk,(num_long+1)/2)>Y_grid(kk-1,(num_long+1)/2)
                temp_THRESHOLD(kk,(num_long+1)/2:(num_long+3)/2)=NaN;
                kk=kk-1;
            end
        else
            kk=1;
            while Y_grid(kk,(num_long+1)/2)<Y_grid(kk+1,(num_long+1)/2)
                temp_THRESHOLD(kk,(num_long+1)/2:(num_long+3)/2)=NaN;
                kk=kk+1;
            end
        end
    else
        if in.center_theta<=90
            kk=num_lat;
            while Y_grid(kk,(num_long+3)/2)>Y_grid(kk-1,(num_long+3)/2)
                temp_THRESHOLD(kk,(num_long+3)/2:(num_long+5)/2)=NaN;
                kk=kk-1;
            end
        else
            kk=1;
            while Y_grid(kk,(num_long+3)/2)<Y_grid(kk+1,(num_long+5)/2)
                temp_THRESHOLD(kk,(num_long+3)/2:(num_long+5)/2)=NaN;
                kk=kk+1;
            end
        end
    end
    % Grid correction end
end
%% Plotting
%%%% figure format
unflat=~in.flat;

if strcmp(in.display_scale,'lin')
    thresh_surf = surf(in.ax,X_map,Y_map,THRESHOLD_int*unflat,THRESHOLD_int,'LineStyle','none','Marker','none');
else
    thresh_surf = surf(in.ax,X_map,Y_map,log10(THRESHOLD_int)*unflat,log10(THRESHOLD_int),'LineStyle','none','Marker','none');
end
hold(in.ax,'on'); 
set(get(get(thresh_surf,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
if in.colorbar_auto
    if strcmp(in.display_scale,'lin')        
        caxis(in.ax,[min_thresh max_thresh]);        
    else
        if min_thresh >= 0
            caxis(in.ax,[floor(log10(min_thresh)*10)/10,ceil(log10(max_thresh)*10)/10]);
        else
            caxis(in.ax,[floor(-log10(abs(min_thresh))*10)/10,ceil(log10(max_thresh)*10)/10]);
        end
    end
else
    if strcmp(in.display_scale,'lin')
        caxis(in.ax,[in.colorbar_min,in.colorbar_max]);        
    else
        caxis(in.ax,[log10(in.colorbar_min),log10(in.colorbar_max)]);
    end
end
hold on;
xlabel(in.ax,'$$X$$ (a.u.)','Interpreter','latex');
ylabel(in.ax,'$$Y$$ (a.u.)','Interpreter','latex');
zlabel(in.ax,sprintf('%s $$\rm{(%s)}$$',in.variable,in.units),'Interpreter','latex');
if in.title_on
    if in.colorbar_norm
        if in.center_at_min
            title(in.ax,{ sprintf('%s of %s. Projection centered at minimum threshold $$\\phi$$.', in.variable,display_cell_name),...
                sprintf('Normalized to Minimum %s: $$%3.2f \\:\\rm{%s}$$ at $$\\phi=%3.1f ^{\\circ}$$, $$\\theta=%3.1f^{\\circ}$$',in.variable,min_th,in.units,mod(minLong,360),minTheta)...
                },'Interpreter','latex');
        else
            title(in.ax,{ sprintf('%s of %s. Projection centered at $$\\phi=%3.1f ^{\\circ}$$, $$\\theta=%3.1f^{\\circ}$$', in.variable,display_cell_name,in.center_phi,in.center_theta),...
                sprintf('Normalized to Minimum %s: $$%3.2f \\:\\rm{%s}$$ at $$\\phi=%3.1f ^{\\circ}$$, $$\\theta=%3.1f^{\\circ}$$',in.variable,min_th_abs,in.units,mod(minLong,360),minTheta)...
                },'Interpreter','latex');
        end
    else
        if in.center_at_min
            title(in.ax,{ sprintf('%s of %s. Projection centered at minimum threshold $$\\phi$$.', in.variable,display_cell_name),...
                sprintf('Minimum %s: $$%3.2f \\:\\rm{%s}$$ at $$\\phi=%3.1f ^{\\circ}$$, $$\\theta=%3.1f^{\\circ}$$',in.variable,min_th,in.units,mod(minLong,360),minTheta)...
                },'Interpreter','latex');
        else
            if min_th > 0
                title(in.ax,{ sprintf('%s of %s. Projection centered at $$\\phi=%3.1f ^{\\circ}$$, $$\\theta=%3.1f^{\\circ}$$', in.variable,display_cell_name,in.center_phi,in.center_theta),...
                    sprintf('Minimum %s: $$%3.2f \\:\\rm{%s}$$ at $$\\phi=%3.1f ^{\\circ}$$, $$\\theta=%3.1f^{\\circ}$$',in.variable,min_th,in.units,mod(minLong,360),minTheta)...
                    },'Interpreter','latex');
            else % negative value means deltaVm is being plotted
                title(in.ax,{ sprintf('%s of %s. Projection centered at $$\\phi=%3.1f ^{\\circ}$$, $$\\theta=%3.1f^{\\circ}$$', in.variable,display_cell_name,in.center_phi,in.center_theta),...
                    sprintf('Maximum %s: $$%1.3f \\:\\rm{%s}$$ at $$\\phi=%3.1f ^{\\circ}$$, $$\\theta=%3.1f^{\\circ}$$',in.variable,max_thresh,in.units,mod(maxLong,360),maxTheta)...
                    },'Interpreter','latex');
            end
        end
    end    
end
if in.colorbar_on
    hcb=colorbar(in.ax);
    set(hcb,'TickDirection','out');
    set(hcb,'TickLabelInterpreter','latex');
    set(hcb,'FontSize',20);
    if in.colorbar_norm
        ticklabels = (get(hcb,'Ticklabels'));    
    else    
        ticklabels = get(hcb,'TickLabels');
    end
    for ii=1:length(ticklabels)
        tickvalue=str2double(ticklabels{ii});
        if in.colorbar_norm==0
            if strcmp(in.display_scale,'log')
                ticklabels{ii}=sprintf('$$%4.1f  \\:\\rm{%s}$$',10^tickvalue,in.units);
            else % linear
                if tickvalue<0.01 % use for small deltaVms
                    ticklabels{ii}=sprintf('$$%4.3f  \\:\\rm{%s}$$',tickvalue,in.units);
                else
                    ticklabels{ii}=sprintf('$$%4.2f  \\:\\rm{%s}$$',tickvalue,in.units);
                end
            end
        else %normalized colorbar to minimum
            if strcmp(in.display_scale,'log')
                ticklabels{ii}=sprintf('$$%4.1f $$',10^tickvalue);
            else
                ticklabels{ii}=sprintf('$$%4.2f $$',tickvalue);
            end
        end
    end
    if max_thresh > 1e-2
        set(hcb,'TickLabels',ticklabels);
    else
        ylabel(hcb,'mV','FontSize',16);
    end
end
if in.display_grid == 1            % if 1, display grid lines for Mollweide map          
    downF=15/in.interp_interval;        
    if strcmp(in.display_scale,'lin')       
        grid_lw = 1;
        plot3(in.ax,X_map(:,[1:downF:front_bounds(1)-downF,front_bounds(2)+downF:downF:end]),Y_map(:,[1:downF:front_bounds(1)-downF,front_bounds(2)+downF:downF:end]),...
            THRESHOLD_int(:,[1:downF:front_bounds(1)-downF,front_bounds(2)+downF:downF:end])*unflat,...
            'LineStyle','-','LineWidth',grid_lw,'Color',[1,1,1]*0.4,'Marker','none'); hold on; % longitude lines
        plot3(in.ax,X_map(1:downF:end,1:front_bounds(1))',Y_map(1:downF:end,1:front_bounds(1))',...
            THRESHOLD_int(1:downF:end,1:front_bounds(1))'*unflat,...
            'LineStyle','-','LineWidth',grid_lw,'Color',[1,1,1]*0.4,'Marker','none'); % left latitude lines
        plot3(in.ax,X_map(1:downF:end,front_bounds(2):end)',Y_map(1:downF:end,front_bounds(2):end)',...
            THRESHOLD_int(1:downF:end,front_bounds(2):end)'*unflat,...
            'LineStyle','-','LineWidth',grid_lw,'Color',[1,1,1]*0.4,'Marker','none'); % right latitude lines
        
        % front grid lines (-90 to +90)        
        plot3(in.ax,X_map(:,front_bounds(1):downF:front_bounds(2)),Y_map(:,front_bounds(1):downF:front_bounds(2)),...
            THRESHOLD_int(:,front_bounds(1):downF:front_bounds(2))*unflat,...
            'LineStyle','-','LineWidth',grid_lw,'Color',[1,1,1]*0.9,'Marker','none'); hold on; % longitude
        plot3(in.ax,X_map(1:downF:end,front_bounds(1):front_bounds(2))',Y_map(1:downF:end,front_bounds(1):front_bounds(2))',...
            THRESHOLD_int(1:downF:end,front_bounds(1):front_bounds(2))'*unflat,...
            'LineStyle','-','LineWidth',grid_lw,'Color',[1,1,1]*0.9,'Marker','none'); % latitude
        
        % Plot -90 and 90 deg Longitude lines
        plot3(X_line,Y_line,THRESHOLD_int(:,front_bounds)'*unflat,...
            'LineStyle','-','LineWidth',grid_lw+2,'Color',[1,1,1]*0.6,'Marker','none');
    else        
        % back
        grid_lw = 1;
        plot3(in.ax,X_map(:,[1:downF:front_bounds(1)-downF,front_bounds(2)+downF:downF:end]),Y_map(:,[1:downF:front_bounds(1)-downF,front_bounds(2)+downF:downF:end]),...
            log10(THRESHOLD_int(:,[1:downF:front_bounds(1)-downF,front_bounds(2)+downF:downF:end]))*unflat,...
            'LineStyle','-','LineWidth',grid_lw,'Color',[1,1,1]*0.3,'Marker','none'); % longitude lines
        plot3(in.ax,X_map(1:downF:end,1:front_bounds(1))',Y_map(1:downF:end,1:front_bounds(1))',...
            log10(THRESHOLD_int(1:downF:end,1:front_bounds(1))')*unflat,...
            'LineStyle','-','LineWidth',grid_lw,'Color',[1,1,1]*0.3,'Marker','none'); % left latitude lines
        plot3(in.ax,X_map(1:downF:end,front_bounds(2):end)',Y_map(1:downF:end,front_bounds(2):end)',...
            log10(THRESHOLD_int(1:downF:end,front_bounds(2):end)')*unflat,...
            'LineStyle','-','LineWidth',grid_lw,'Color',[1,1,1]*0.3,'Marker','none'); % right latitude lines
        
        % front
        plot3(in.ax,X_map(:,front_bounds(1):downF:front_bounds(2)),Y_map(:,front_bounds(1):downF:front_bounds(2)),...
            log10(THRESHOLD_int(:,front_bounds(1):downF:front_bounds(2)))*unflat,...
            'LineStyle','-','LineWidth',grid_lw,'Color',[1,1,1]*0.9,'Marker','none'); % longitude
        plot3(in.ax,X_map(1:downF:end,front_bounds(1):front_bounds(2))',Y_map(1:downF:end,front_bounds(1):front_bounds(2))',...
            log10(THRESHOLD_int(1:downF:end,front_bounds(1):front_bounds(2))')*unflat,...
            'LineStyle','-','LineWidth',grid_lw,'Color',[1,1,1]*0.9,'Marker','none'); % latitude
        %Plot 180 and 0 deg Longitude lines
        plot3(in.ax,X_line,Y_line,log10(THRESHOLD_int(:,front_bounds)')*unflat,...
            'LineStyle','-','LineWidth',grid_lw+2,'Color',[1,1,1]*0.6,'Marker','none');
    end
    
elseif in.display_grid == 2 % plots latitude lines for 60 and 120 deg
    if strcmp(in.display_scale,'lin')
        plot3(in.ax,X_linet(1,:)',Y_linet(1,:)',THRESHOLD_int(trans_bounds(1),:),'--k')
        plot3(in.ax,X_linet(2,:)',Y_linet(2,:)',THRESHOLD_int(trans_bounds(2),:),'--k')
        plot3(in.ax,X_map(:,1),Y_map(:,1),THRESHOLD_int(:,1),'k')
        plot3(in.ax,X_map(:,end),Y_map(:,end),THRESHOLD_int(:,end),'k')
    else % log scale
        plot3(in.ax,X_linet(1,:)',Y_linet(1,:)',log10(THRESHOLD_int(trans_bounds(1),:)),'--k')
        plot3(in.ax,X_linet(2,:)',Y_linet(2,:)',log10(THRESHOLD_int(trans_bounds(2),:)),'--k')
        plot3(in.ax,X_map(:,1),Y_map(:,1),log10(THRESHOLD_int(:,1)),'k')
        plot3(in.ax,X_map(:,end),Y_map(:,end),log10(THRESHOLD_int(:,end)),'k')
    end
end

if in.display_old_grid         % if 1, display grid lines of original sampling points
    if strcmp(in.display_scale,'lin')
        plot3(in.ax,X_grid,Y_grid,(temp_THRESHOLD+max_thresh)*unflat,...
            'LineStyle','--','LineWidth',1.5,'Color',[1,1,1]*0.3,'Marker','none','MarkerSize',20);
        plot3(in.ax,X_grid',Y_grid',(temp_THRESHOLD'+max_thresh)*unflat,...
            'LineStyle','-.','LineWidth',1,'Color',[1,1,1]*0.1,'Marker','.','MarkerSize',20);
    else
        plot3(in.ax,X_grid,Y_grid,(log10(temp_THRESHOLD)+1)*unflat,...
            'LineStyle','--','LineWidth',1.5,'Color',[1,1,1]*0.3,'Marker','none','MarkerSize',20);
        plot3(in.ax,X_grid',Y_grid',(log10(temp_THRESHOLD')+1)*unflat,...
            'LineStyle','-.','LineWidth',1,'Color',[1,1,1]*0.1,'Marker','.','MarkerSize',20);
    end
end
view(in.ax,0,90);
axis(in.ax,'off'); 

set(in.ax,'DataAspectRatio',[repmat(max(diff(get(in.ax,'XLim')),diff(get(in.ax,'YLim'))),[1,2]),diff(get(in.ax,'ZLim'))]);

if ~in.center_at_min && in.plot_min_point
    min_marker = 'p';
    [min_temp,ind_temp]=min(temp_THRESHOLD);
    [~,ind_2]=min(min_temp);
    if min(min_temp) > 0 % indicates this is polarization data
        ind_1=ind_temp(ind_2);
    else
        [max_temp,ind_temp] = max(temp_THRESHOLD);
        [~,ind_2] = max(max_temp);
        ind_1 = ind_temp(ind_2);
    end        
    
    if strcmp(in.display_scale,'lin')
        plot3(in.ax,X_grid(ind_1,ind_2),Y_grid(ind_1,ind_2),(temp_THRESHOLD(ind_1,ind_2)+max_thresh)*unflat,...
            'Marker',min_marker,'MarkerSize',15,'MarkerEdgeColor',[0,0,0],'MarkerFaceColor',[1,1,1],'LineWidth',1.5);
    else
        plot3(in.ax,X_grid(ind_1,ind_2),Y_grid(ind_1,ind_2),(log10(temp_THRESHOLD(ind_1,ind_2))+1)*unflat,...
            'Marker',min_marker,'MarkerSize',15,'MarkerEdgeColor',[0,0,0],'MarkerFaceColor',[1,1,1],'LineWidth',1.5);
    end
    
end
drawnow
end