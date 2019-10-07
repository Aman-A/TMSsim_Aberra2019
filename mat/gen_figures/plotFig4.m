function plotFig4(save_fig)
%PLOTFIG4 Plots activated axons overlaid on morphologies for monophasic PA
%and AP TMS in single L2/3 PC, L4 LBC, and L5 PC populations
%   
% REQUIRES 'CURVE FITTING' MATLAB TOOLBOX FOR MAKESLICES
% AUTHOR    : Aman Aberra 
%%
if nargin==0
   save_fig = 0;
end
mat_dir = addPaths; 
% Simulation settings
nrn_model_ver = 'maxH';
mode = 1; % monophasic MagProX100 pulse
layer_set_num = 1;
Efield_name = 'M1_PA_MCB70'; 
cell_ids = [6,11,16]; % L23 PC 1, L4 LBC 1, L5 PC 1
cell_layers = [2,3,4]; 
nrn_pop = 'nrn_pop1'; 
model_prefixes = {sprintf('tms_%s_w%g_ls_%g_E_%s_P_%s',nrn_model_ver,mode,...
                            layer_set_num,Efield_name,nrn_pop);
                   sprintf('tms_%s_w%g_ls_%g_E_%s_r_P_%s',nrn_model_ver,mode,...
                            layer_set_num,Efield_name,nrn_pop)}; 
% Plot settings
opts.cell_inds_file = sprintf('gen_figures/cellplot_inds_ls_%g_%s',layer_set_num,nrn_pop);
opts.morph_alpha = 0.5; % axon morphology transparency 
opts.morph_color = 'k'; % axon morphology color
opts.term_act_col = [0.1 0.8 0]; % green 
opts.term_act_col2 = [1 0 1]; % magenta
opts.term_lw = 3.5; 
bcol = 'k';% boundary color
blw = 1; % boundary linewidth
normROI = [-40 -35.4 0 0 33.93 45]; % 45� slice - postsulcusandlip2 - posterior part of sulcus plus little lip
ax_zoom = [-38.6 -32.0057 -19.8599 -13.6250 normROI(5)-1 44.6475]; % for plane 45� - deeper
title_str = {'L2/3 PC','L4 LBC','L5 PC'}; 
%% Load data
layersP = loadLayers(layer_set_num,1); 
MeshROI_file = fullfile(mat_dir,'output_data','layer_data','MeshROI.mat'); 
MeshROI_data = load(MeshROI_file); 
MeshROI = MeshROI_data.MeshROI; 
uvROI = [0.271 0.52 0.42 1]; % FDI ROI for 45 deg plane (newer)
plane = [-35.3816  -17.0009 44.3816 -0.5 0.5 0]; % **NEW 45 deg PLANE**
start_dim = 2;
slices = makeSlices(layersP,MeshROI,plane,start_dim,0,uvROI);
%% Plot
model_prefixPA = model_prefixes{1};
model_prefixAP = model_prefixes{2};
for c = 1:length(cell_ids)            
    cell_id = cell_ids(c); cell_layer = cell_layers(c);         
    fig = figure('Units','normalized','Position',[0.4236 0.05 0.572 0.843],'Color','w'); 
    plot3(slices(1).slice_line(:,1),slices(1).slice_line(:,2),slices(1).slice_line(:,3),'Color',bcol,'LineWidth',blw);
    hold on; 
    plot3(slices(end).slice_line(:,1),slices(end).slice_line(:,2),slices(end).slice_line(:,3),'Color',bcol,'LineWidth',blw);       
    view(plane(4:6)); axis equal; axis off; axis tight; 
    axis(ax_zoom); 
    ax = gca;
    opts.ax = ax; 
    plotActTermsOnMorphs2(cell_id,cell_layer,nrn_model_ver,nrn_pop,...
        model_prefixPA,model_prefixAP,opts);    
    ax.ClippingStyle = 'rectangle';
%     title(title_str{c},'FontSize',18); 
    if save_fig
        fig_name = sprintf('Fig4_%g',c);
        fig_fold = fullfile(mat_dir,'figures');
        savefig(fig,fullfile(fig_fold,[fig_name '.fig']));
        print(fig,fullfile(fig_fold,[fig_name '.tiff']),'-dtiff','-r600');       
    end
end    
