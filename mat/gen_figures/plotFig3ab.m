function plotFig3ab(save_fig,plot_figs,fig_names,ax_view)
%PLOTFIG3AB Plot E-magnitude and normal component on layers for Fig3a and 
%Fig3b, respectively
%  
% AUTHOR    : Aman Aberra
if nargin==0
   save_fig = 0; 
   plot_figs = [1 2];
   fig_names = {'Fig3a','Fig3b'}; 
   ax_view = [-89.2 45]; % for Fig 3 plots
end
mat_dir = addPaths; 
layer_set_num = 1;
Efield_name = 'M1_PA_MCB70';
plot_layers = 1:5;
% Plot settings
n_cmap = 1000; % number of colormap levels
cmap_name = 'jet'; % flipped parula
normE = 1; % normalize E-magnitude to max
cmap2_name = 'bwr'; 
z_lims = [22 52.4057]; % or []

lt_pos = [-411 -807 836]; % [170.3353 52.0246 1.2195e3] 
shift_dir = [0 -35 0]; 
vecs_on = 0; % no vectors
%% Load data
% layersE
layersE = loadLayers(layer_set_num,2,Efield_name); 
%% Plot E magnitude
if any(plot_figs==1)
    fig1 = figure('units','normalized','outerposition',[0 0 1 1],'Color',[1 1 1]); 
    E_mag = plotLayersE(layersE,plot_layers,'mag',normE,'',ax_view,shift_dir,vecs_on);
    ax1 = gca;
    axis(ax1,'off'); 
    title('');     
    ax1.Children(1).Position = lt_pos;
    ax1.Children(1).Style = 'local';    
    % cut off sulcus
    if ~isempty(z_lims)
       ax1.ZLim = z_lims; 
       ROI = axis; 
       inds_all = cell(length(plot_layers),1);    
       for i = 1:length(plot_layers)
          [~,inds_all{i} ] = clipPoints3d(layersE(plot_layers(i)).surface.vertices,ROI);       
       end
       % min E in each layer within clipped ROI
       minEs = cellfun(@(x,y) min(x(y)),E_mag,inds_all,'UniformOutput',1);
       cax_lims = [min(minEs) 1];       
    else
        caxis([min(cellfun(@min,E_mag)) 1]);
    end
    if strcmp(cmap_name,'jet')
        cmap = colormap(jet(n_cmap));
    end    
end
%% Plot E-normal
if any(plot_figs==2)
    fig2 = figure('units','normalized','outerposition',[0 0 1 1],'Color',[1 1 1]); 
    E_norm = plotLayersE(layersE,plot_layers,'norm',normE,'',ax_view,shift_dir,vecs_on);
    ax2 = gca;
    axis(ax2,'off'); 
    title('');     
    ax2.Children(1).Position = lt_pos;
    ax2.Children(1).Style = 'local';
    % cut off sulcus
    if ~isempty(z_lims)
        ax2.ZLim = z_lims;
        cax2_lims = [-1 1];
        caxis(cax2_lims)
    else
        cax2_lims = [min(cellfun(@min,E_norm)) max(cellfun(@max,E_norm))];
        caxis(cax2_lims);   
    end
    if strcmp(cmap2_name,'bwr')
        colormap(bluewhitered(n_cmap));
    end    
end
%% Save figures
if save_fig
   fig_fold = fullfile(mat_dir,'figures');
   if any(plot_figs == 1)
       fig_name1 = fig_names{1};
       savefig(fig1, fullfile(fig_fold,[fig_name1 '.fig']));
       print(fig1, fullfile(fig_fold,[fig_name1 '.tiff']),'-dtiff','-r250','-cmyk');
   end   
   if any(plot_figs == 2)
       fig_name2 = fig_names{2};
       savefig(fig2, fullfile(fig_fold,[fig_name2 '.fig']));
       print(fig2, fullfile(fig_fold,[fig_name2 '.tiff']),'-dtiff','-r250','-cmyk');
   end               
end
end

