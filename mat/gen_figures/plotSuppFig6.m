function plotSuppFig6(save_figs)
%PLOTSUPPFIG6 Plots threshold distributions and analysis plots comparing
%simulations with homogeneous intracranial conductivities in FEM to those
%from heterogeneous case
%
% AUTHOR    : Aman Aberra 
if nargin==0
   save_figs = 0;
   plot_figs = {'a','b','c','d'}; 
end
mat_dir = addPaths; 
% Simulation settings
layer_set_num = 1;
Efield_namePA_het = 'M1_PA_MCB70';
Efield_namePA_hom = 'M1_PA_MCB70_hom';
Efield_nameAP_het = 'M1_PA_MCB70_r';
Efield_nameAP_hom = 'M1_PA_MCB70_hom_r';
nrn_model_ver = 'maxH';
mode = 1; % monophasic MagProX100 pulse
cell_ids = {[];[];[];16:20;[]}; 
plot_layers = 4; 
% cell_ids = {1:5;6:10;11:15;16:20;21:25}; 
nrn_pop = 'nrn_pop1-nrn_pop6_all';  
plot_region_name = 'FDI_rep_inds'; 
% Analysis settings
cutoffs = 1.1:0.1:1.5; % 10% up to 100% above minimum threshold
%% Plot settings
cmap = [flipud(fake_parula(1000));0.8 0.8 0.8]; % add gray for values above cutoff
clims = [1 1.5]; % A/us
z_lims = [22 52.4057]; % or []
lt_pos = [-411 -807 836]; % [170.3353 52.0246 1.2195e3] 
ax_view = [-89.2 45]; % [-89.2 70.8]
plot_opts.shift_dir = [0,0,0]; % vector of direction of shift
plot_opts.norm_mode = 'min_layer'; % normalize to min threshold in layer
fig_size = [9 5.4];
font_size = 10; 
font_name = 'Times';
%% Load data
layers = loadLayers(layer_set_num); 
data_fold = fullfile(mat_dir,'nrn_sim_data');
%% Load ROI
ROI_filename = fullfile(mat_dir,'gen_figures',[plot_region_name '.mat']); 
ROI_data = load(ROI_filename);
ROIi = ROI_data.ROIi; 
layersROIs = [layers.surface];
%% Load data
model_prefixa1 = sprintf('tms_%s_w%g_ls_%g_E_%s_P_%s',nrn_model_ver,mode,...
                            layer_set_num,Efield_namePA_het,nrn_pop); 
data_struct = load(fullfile(data_fold,model_prefixa1)); 
threshEsPAhet = data_struct.threshEs; 
cell_model_names = data_struct.cell_model_names;
model_prefixa2 = sprintf('tms_%s_w%g_ls_%g_E_%s_P_%s',nrn_model_ver,mode,...
                            layer_set_num,Efield_namePA_hom,nrn_pop); 
data_struct = load(fullfile(data_fold,model_prefixa2)); 
threshEsPAhom = data_struct.threshEs; 
model_prefixc1 = sprintf('tms_%s_w%g_ls_%g_E_%s_P_%s',nrn_model_ver,mode,...
                            layer_set_num,Efield_nameAP_het,nrn_pop); 
data_struct = load(fullfile(data_fold,model_prefixc1)); 
threshEsAPhet = data_struct.threshEs; 
model_prefixc2 = sprintf('tms_%s_w%g_ls_%g_E_%s_P_%s',nrn_model_ver,mode,...
                            layer_set_num,Efield_nameAP_hom,nrn_pop); 
data_struct = load(fullfile(data_fold,model_prefixc2)); 
threshEsAPhom = data_struct.threshEs; 
%% Plot a (normalized threshold plots on surface - PA)
if any(strcmp(plot_figs,'a'))
% Plot heterogeneous    
    plotDataLayers(layers,threshEsPAhet,cell_model_names,cell_ids,plot_opts); 
    figa1 = gcf;
    figa1.Units = 'normalized';
    figa1.OuterPosition = [0 0 1 1];
    figa1.Color = 'w';
    colormap(cmap); caxis(clims); 
    view(ax_view); 
    % change light
    ax = gca;
    ax.Children(1).Position = lt_pos;
    ax.Children(1).Style = 'local';
    % cut off sulcus
    if ~isempty(z_lims)
       ax.ZLim = z_lims; 
    end
    % Add hand muscle representation outline
    layersROIs_plot = layersROIs;
    for i = 1:length(plot_layers)  
        ii = plot_layers(i); 
        layersROIs_plot(ii).faces = layersROIs_plot(ii).faces(ROIi{ii},:);
        [~,vb] = meshboundary(layersROIs_plot(ii));
        plot3(vb(1,:)+plot_opts.shift_dir(1)*(i-1),vb(2,:)+plot_opts.shift_dir(2)*(i-1),vb(3,:)+plot_opts.shift_dir(3)*(i-1),'k','LineWidth',1);
    end
    title(ax,'Heterogeneous conductivities - monophasic P-A TMS');
% Plot homogeneous    
    plotDataLayers(layers,threshEsPAhom,cell_model_names,cell_ids,plot_opts); 
    figa2 = gcf;
    figa2.Units = 'normalized';
    figa2.OuterPosition = [0 0 1 1];
    figa2.Color = 'w';
    colormap(cmap); caxis(clims); 
    view(ax_view); 
    % change light
    ax = gca;
    ax.Children(1).Position = lt_pos;
    ax.Children(1).Style = 'local';
    % cut off sulcus
    if ~isempty(z_lims)
       ax.ZLim = z_lims; 
    end
    % Add hand muscle representation outline    
    layersROIs_plot = layersROIs;
    for i = 1:length(plot_layers)  
        ii = plot_layers(i); 
        layersROIs_plot(ii).faces = layersROIs_plot(ii).faces(ROIi{ii},:);
        [~,vb] = meshboundary(layersROIs_plot(ii));
        plot3(vb(1,:)+plot_opts.shift_dir(1)*(i-1),vb(2,:)+plot_opts.shift_dir(2)*(i-1),vb(3,:)+plot_opts.shift_dir(3)*(i-1),'k','LineWidth',1);
    end
    title(ax,'Homogeneous conductivities - monophasic P-A TMS');
    % Save
    if save_figs
        fig_fold = fullfile(mat_dir,'figures');
        fig_namea1 = 'SuppFig6a-1';         
        savefig(figa1,fullfile(fig_fold,[fig_namea1 '.fig']));
        print(figa1,fullfile(fig_fold,[fig_namea1 '.tif']),'-dtiff','-cmyk','-r250');
        fig_namea2 = 'SuppFig6a-2';
        savefig(figa2,fullfile(fig_fold,[fig_namea2 '.fig']));
        print(figa2,fullfile(fig_fold,[fig_namea2 '.tif']),'-dtiff','-cmyk','-r250');
    end
end
%% Plot b (activated layer area het vs. hom)
if any(strcmp(plot_figs,'b'))
    med_threshEsL5 = median(cell2mat(threshEsPAhet(16:20)),2);
    med_threshEsL5 = med_threshEsL5/min(med_threshEsL5); 
    med_threshEsL5_hom = median(cell2mat(threshEsPAhom(16:20)),2);
    med_threshEsL5_hom =med_threshEsL5_hom/min(med_threshEsL5_hom); 
    area_below = zeros(length(cutoffs),1); 
    area_below_hom = zeros(length(cutoffs),1); 
    for i = 1:length(cutoffs)
       area_below(i) = sum(getPatchAreas(layers(4).surface.vertices,layers(4).surface.faces,...
                            find(med_threshEsL5 <= cutoffs(i)*min(med_threshEsL5)))); 
       area_below_hom(i) = sum(getPatchAreas(layers(4).surface.vertices,layers(4).surface.faces,...
                            find(med_threshEsL5_hom <= cutoffs(i)*min(med_threshEsL5_hom))));
    end
    figb = figure('Color','w'); 
    figb.Units = 'centimeters';
    figb.Position(3:4) = fig_size; 
    b = bar(cutoffs,[area_below,area_below_hom]); 
    b(1).FaceColor = 'k';
    b(2).FaceColor = 'w';
    leg = legend('Heterogeneous','Homogeneous');     
    leg.FontSize = font_size; 
    leg.FontName = font_name;
    leg.Position = [0.1971 0.7386 0.4108 0.1797]; 
    leg.Box = 'off';
    box off;
    ax = gca;
    ax.XLim = [cutoffs(1)-0.1 cutoffs(end)+0.1]; 
    ax.XColor = 'k';
    ax.YColor = 'k';
    ax.YGrid = 'on';
    xlabel('Activation cutoff relative to minimum threshold'); 
    % xlabel({'Activation cutoff relative', 'to minimum threshold'}); 
    ylabel('Activated layer area (mm^{2})'); 
    ax.FontSize = font_size; 
    ax.FontName = font_name;
    ax.XTick = cutoffs;   
    if save_figs
       fig_fold = fullfile(mat_dir,'figures');
       fig_name = 'SuppFig6b';
       savefig(figb,fullfile(fig_fold,[fig_name '.fig']));
       print(figb,fullfile(fig_fold,[fig_name '.eps']),'-depsc2','-cmyk');
    end
end
%% Plot c (normalized threshold plots on surface - AP)
if any(strcmp(plot_figs,'c'))
% Plot heterogeneous    
    plotDataLayers(layers,threshEsAPhet,cell_model_names,cell_ids,plot_opts); 
    figc1 = gcf;
    figc1.Units = 'normalized';
    figc1.OuterPosition = [0 0 1 1];
    figc1.Color = 'w';
    colormap(cmap); caxis(clims); 
    view(ax_view); 
    % change light
    ax = gca;
    ax.Children(1).Position = lt_pos;
    ax.Children(1).Style = 'local';
    % cut off sulcus
    if ~isempty(z_lims)
       ax.ZLim = z_lims; 
    end
    % Add hand muscle representation outline
    layersROIs_plot = layersROIs;
    for i = 1:length(plot_layers)  
        ii = plot_layers(i); 
        layersROIs_plot(ii).faces = layersROIs_plot(ii).faces(ROIi{ii},:);
        [~,vb] = meshboundary(layersROIs_plot(ii));
        plot3(vb(1,:)+plot_opts.shift_dir(1)*(i-1),vb(2,:)+plot_opts.shift_dir(2)*(i-1),vb(3,:)+plot_opts.shift_dir(3)*(i-1),'k','LineWidth',1);
    end
    title(ax,'Heterogeneous conductivities - monophasic A-P TMS');
% Plot homogeneous    
    plotDataLayers(layers,threshEsAPhom,cell_model_names,cell_ids,plot_opts); 
    figc2 = gcf;
    figc2.Units = 'normalized';
    figc2.OuterPosition = [0 0 1 1];
    figc2.Color = 'w';
    colormap(cmap); caxis(clims); 
    view(ax_view); 
    % change light
    ax = gca;
    ax.Children(1).Position = lt_pos;
    ax.Children(1).Style = 'local';
    % cut off sulcus
    if ~isempty(z_lims)
       ax.ZLim = z_lims; 
    end
    % Add hand muscle representation outline  
    layersROIs_plot = layersROIs;
    for i = 1:length(plot_layers)
        ii = plot_layers(i); 
        layersROIs_plot(ii).faces = layersROIs_plot(ii).faces(ROIi{ii},:);
        [~,vb] = meshboundary(layersROIs_plot(ii));
        plot3(vb(1,:)+plot_opts.shift_dir(1)*(i-1),vb(2,:)+plot_opts.shift_dir(2)*(i-1),vb(3,:)+plot_opts.shift_dir(3)*(i-1),'k','LineWidth',1);
    end
    title(ax,'Heterogeneous conductivities - monophasic A-P TMS');
    % Save
    if save_figs
        fig_fold = fullfile(mat_dir,'figures');
        fig_namea1 = 'SuppFig6c-1';         
        savefig(figa1,fullfile(fig_fold,[fig_namea1 '.fig']));
        print(figa1,fullfile(fig_fold,[fig_namea1 '.tif']),'-dtiff','-cmyk','-r250');
        fig_namea2 = 'SuppFig6c-2';
        savefig(figa2,fullfile(fig_fold,[fig_namea2 '.fig']));
        print(figa2,fullfile(fig_fold,[fig_namea2 '.tif']),'-dtiff','-cmyk','-r250');
    end
end
%% Plot d (change in median threshold AP - PA)
if any(strcmp(plot_figs,'b'))
    all_data = {threshEsPAhet; % PA
            threshEsAPhet;
            threshEsPAhom; % PA
            threshEsAPhom}; 
    all_med = cell(length(all_data),1); 
    num_layers = length(layers); 
    for k = 1:length(all_data)
        threshEsk = all_data{k};
        all_med{k} = cell(num_layers,1);
        for i = 1:num_layers
            inds = (1+num_layers*(i-1)):i*num_layers;
            all_med{k}{i} = cell2mat(threshEsk(inds));
            all_med{k}{i} = all_med{k}{i}(ROIi{i},:); % get ROI
        end
    end
    % (AP-PA)/PA
    med_diff = cellfun(@(x,y) 100*(median(x(:))-median(y(:)))/median(y(:)),all_med{2},all_med{1},'UniformOutput',1);
    med_diff_hom = cellfun(@(x,y) 100*(median(x(:))-median(y(:)))/median(y(:)),all_med{4},all_med{3},'UniformOutput',1);    
    figd = figure; 
    figd.Units = 'centimeters';
    figd.Position(3:4) = fig_size;
    b2 = bar([med_diff,med_diff_hom]); 
    b2(1).FaceColor = 'k';
    b2(2).FaceColor = 'w';
    ax = gca;
    % ax.XLim = [1 5]; 
    ax.XLim = [0.5 5.5]; 
    box off; 
    ax.XTickLabels = {'L1 NGC','L2/3 PC','L4 LBC','L5 PC','L6 PC'}; 
    ylabel({'Change in median threshold (%)',' A?P - P?A'}); 
    ax.YGrid = 'on';
    ax.FontSize = font_size;  
    ax.FontName = font_name;
    ax.XColor = 'k'; ax.YColor = 'k'; 
    if save_figs
       fig_fold = fullfile(mat_dir,'figures');
       fig_name = 'SuppFig6d';
       savefig(figd,fullfile(fig_fold,[fig_name '.fig']));
       print(figd,fullfile(fig_fold,[fig_name '.eps']),'-depsc2','-cmyk');
    end
end
