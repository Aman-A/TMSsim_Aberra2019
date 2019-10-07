% Enter cell_id, array of thetas and phis from full threshold simulation,
% and cell array of initial spike secnames corresponding to the
% orientations in thetas, phis. Plots AP initiation point on cell
% morphology using compartment coordinates

function plotInitPoint(ax1,cell_id,threshEs,thetas,phis,init_inds,...
         arrow_size,expand_ax,xlimits,zlimits,mark_size,nrn_model_ver)
    if isempty(init_inds)
       plot_sec = 0;
    else
        plot_sec = 1; 
    end        
    [~,min_ind] = min(threshEs);
    theta_plot = thetas(min_ind); phi_plot = phis(min_ind);         
    % Get cell coordinates    
    sec_ind = init_inds(min_ind); 
    cell_data = loadCellData(cell_id,nrn_model_ver);
    C = cell_data.C; 
    % Plot     
    [Ex,Ey,Ez] = angle2vec(theta_plot,phi_plot,1);
    axes(ax1); % set ax1 as gca
    plotCellLines(cell_id,nrn_model_ver,'cell_data',cell_data);    
    view(ax1,[0 -1 0]);                
    ax1.DataAspectRatio = [1 1 1];    
    clims = caxis(ax1);      
    caxis(ax1,'manual'); % fix colormap scaling     
    hold(ax1,'on');
    grid_mult = 1;
    [X,Y,~] = meshgrid(grid_mult*linspace(min(C(:,1)),max(C(:,1)),3),grid_mult*linspace(min(C(:,2)),max(C(:,2)),3),grid_mult*linspace(min(C(:,3)),max(C(:,3)),3));
    %[X,Y,Z] = meshgrid(grid_mult*linspace(-1365,1169,3),grid_mult*linspace(min(C(:,2)),max(C(:,2)),3),grid_mult*linspace(-1170,390.2,3));
    U = Ex*ones(size(X)); V = Ey*ones(size(X)); W = Ez*ones(size(X));
     
    if plot_sec
        % Plot AP initiation point
        if isempty(sec_ind)
           disp('Min threshold not found'); 
        end        
        plot3(ax1,C(sec_ind,1),min(C(:,2)),C(sec_ind,3),'Marker','p','MarkerSize',mark_size,'MarkerFaceColor','white','MarkerEdgeColor','black',...
            'LineWidth',0.5);
    end
    sm = 1;

    if ~strcmp(xlimits,'auto')
        ax1.XLim = xlimits; 
    end
    if ~strcmp(zlimits,'auto')
        ax1.ZLim = zlimits;
    end    
%     [~,max_dim] = max([Ex,Ey,Ez]); 
%     arrow_size = 0.6*range(C(:,max_dim));
    if Ex < 0
        x_shift = -Ex*arrow_size;
    else
       x_shift = 0; 
    end
    if Ez > 0
        z_shift = -Ez*arrow_size;
    else
        z_shift = 0;
    end
    q = quiver3D([ax1.XLim(1)*0.8+x_shift,0,0.8*ax1.ZLim(2)+z_shift],arrow_size*[Ex Ey Ez]); 
    camlight head; 
    % Plot scale bar in bottom left of figure
    plot3(ax1,[sm*ax1.XLim(1),sm*ax1.XLim(1),sm*ax1.XLim(1)+250],[min(min(min(Y))),min(min(min(Y))),min(min(min(Y)))],[ax1.ZLim(1)+250,ax1.ZLim(1),ax1.ZLim(1)],...
        'Color','k','LineWidth',3);
    hold off;
    % reset color limits to original limits
    ax1.CLim = clims;    
    if expand_ax   
        ax1.ZLim = [ax1.ZLim(1)*1.1 ax1.ZLim(2)]; 
    end       
end