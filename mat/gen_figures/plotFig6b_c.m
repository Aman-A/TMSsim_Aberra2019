function plotFig6b_c(save_fig)
%PLOTFIG5 Plots thresholds and threshold differences on 2d analysis plane for 
%Fig. 5a-c and Fig5d
%
% REQUIRES 'CURVE FITTING' MATLAB TOOLBOX FOR MAKESLICES
% AUTHOR    : Aman Aberra 
if nargin==0
   save_fig = 0;
end
mat_dir = addPaths; 
% Simulation settings
nrn_model_ver = 'maxHlin';
modes = [1]; % 1 monophasic, 3 half sine, 2 biphasic MagProX100 pulses
layer_set_num = 1;
Efield_name = 'M1_PA_MCB70'; 
cell_ids = {[];6:10;[];[];[]}; 
nrn_pop = 'nrn_pop1'; 
model_prefixes = [arrayfun(@(x) sprintf('tms_%s_w%g_ls_%g_E_%s_P_%s',nrn_model_ver,x,...
                            layer_set_num,Efield_name,nrn_pop),modes,'UniformOutput',0)';
                   arrayfun(@(x) sprintf('tms_%s_w%g_ls_%g_E_%s_r_P_%s',nrn_model_ver,x,...
                            layer_set_num,Efield_name,nrn_pop),modes,'UniformOutput',0)']; 
% Plot settings
% Threshold plot settings
cmap = [flipud(fake_parula(1000));0.8 0.8 0.8]; % add gray for values above cutoff
clims = [70 230]; 
cutoff = clims(2)*ones(1,3);
% Threshold diff plot settings
clims2 = [-70 260];
% Sice settings
uvROI = [0.271 0.52 0.42 1]; % FDI ROI for 45deg plane
plane = [-35.3816  -17.0009 44.3816 -0.5 0.5 0]; % 
normROI = [-40 -35.4 0 0 33.93 45]; % 45deg slice - posterior part of sulcus plus little lip
start_dim = 2;
%% Load data
layers = loadLayers(layer_set_num,0);
layersP = loadLayers(layer_set_num,1);
MeshROI_file = fullfile(mat_dir,'output_data','layer_data','MeshROI.mat'); 
MeshROI_data = load(MeshROI_file); 
MeshROI = MeshROI_data.MeshROI; 
% Make slices
slices = makeSlices(layersP,MeshROI,plane,start_dim,0,uvROI);
%% Plot thresholds
% fig1 = figure('units','normalized','outerposition',[0 0 1 1],'paperunits','in','paperposition',[0 0 16 9],'Color',[1 1 1],'PaperPositionMode','auto');
fig1 = figure('units','centimeters');
fig1.Position(3:4) = [41.0986 17.0039]; 
num_rows = 1; num_cols = 2; splots = cell(num_rows,num_cols);
fprintf('Printing with %g rows and %g cols, %g entries\n',num_rows,num_cols,num_rows*num_cols);
for r = 1:num_rows
    for c = 1:num_cols
        splots{r,c} = axes('Position',[(1/num_cols)*(c-1),(1/num_rows)*(num_rows-r),1/num_cols,0.9/num_rows ]);
        axis(splots{r,c},'off');
        hold(splots{r,c},'all');        
    end
end  
threshEs_med_all = cell(length(model_prefixes),1); 
ax_zoom = [-38.6 -32.0057 -19.8599 -13.6250 normROI(5)-1 44.6475]; % for plane 45 deg - deeper
for i = 1:length(model_prefixes)
    [c,r] = ind2sub([num_cols,num_rows],i);
    ax = splots{r,c}; 
    axes(ax); % set gca
    model_prefix = model_prefixes{i};
    datai = load(fullfile(mat_dir,'nrn_sim_data',[model_prefix '.mat'])); 
    threshEs = datai.threshEs; 
    cell_model_names = datai.cell_model_names; 
    threshEs_med_all{i} = plotDataSlice(threshEs,cell_model_names,...                                        
                                            cell_ids,layers,slices,plane,'cutoff',cutoff(1),'above_cutoff_val',cutoff(1));
   view(ax,plane(4:6));
   axis(ax,'off','equal','tight',ax_zoom);       
   colormap(ax,cmap); 
   caxis(ax,clims);            
   drawnow; 
end
%% Plot threshold differences
% fig2 = figure('units','normalized','outerposition',[0 0 1 1],'paperunits','in','paperposition',[0 0 16 9],'Color',[1 1 1],'PaperPositionMode','auto');
fig2 = figure('units','centimeters');
fig2.Position(3:4) = [34 20]; 
num_rows = 1; num_cols = 1; 
splots2 = cell(num_rows,num_cols);
fprintf('Printing with %g rows and %g cols, %g entries\n',num_rows,num_cols,num_rows*num_cols);
for r = 1:num_rows
    for c = 1:num_cols
        splots2{r,c} = axes('Position',[(1/num_cols)*(c-1),(1/num_rows)*(num_rows-r),1/num_cols,0.9/num_rows ]);
        axis(splots2{r,c},'off');
        hold(splots2{r,c},'all');        
    end
end
% Plot
threshEs_pdiff_all = cell(length(model_prefixes)/2,1); 
for i = 1:length(model_prefixes)/2
    threshEsPAi = threshEs_med_all{i};
    threshEsAPi = threshEs_med_all{i+1};
    threshEs_pdiff_all{i} = cellfun(@(x,y) 100*(x-y)./y,threshEsAPi,threshEsPAi,'UniformOutput',0); 
    [c,r] = ind2sub([num_cols,num_rows],i);
    ax = splots2{r,c};
    axes(ax); % set current axis
    plotSlices(slices,threshEs_pdiff_all{i});     
    view(ax,plane(4:6));
    axis(ax,'off','equal','tight',ax_zoom);   
    caxis(ax,clims2); 
    cmap2 = [redwhiteblue(1000);0.8 0.8 0.8]; 
    colormap(ax,cmap2);        
    for j = 1:5       
       pi = ax.Children(1+length(ax.Children)-j); % l1 to l5
       layerjPA = threshEsPAi{j}; layerjAP = threshEsAPi{j}; 
       layerjPA = [layerjPA; flipud(layerjPA)];
       layerjAP = [layerjAP; flipud(layerjAP)];
       pi.FaceVertexCData = [threshEs_pdiff_all{i}{j};flipud(threshEs_pdiff_all{i}{j})]; % reset to original percent diff
       above_thresh = layerjPA > cutoff(i) & layerjAP > cutoff(i); % both PA/AP above threshold              
       pi.FaceVertexCData(above_thresh) = 1e4; 
       pi.LineWidth = 1;        
   end
    drawnow; 
end
splots{1,1}.Units = 'centimeters';
splots2{1,1}.Units = 'centimeters';
splots2{1,1}.Position(3:4) = splots{1,1}.Position(3:4); 
%%
if save_fig
   fig_fold = fullfile(mat_dir,'figures');
   fig_name1 = 'Fig6b';
   fig_name2 = 'Fig6c';   
   savefig(fig1,fullfile(fig_fold,[fig_name1 '.fig']));
   print(fig1,fullfile(fig_fold,[fig_name1 '.tiff']),'-dtiff','-r300');   
   savefig(fig2,fullfile(fig_fold,[fig_name2 '.fig']));
   print(fig2,fullfile(fig_fold,[fig_name2 '.tiff']),'-dtiff','-r280');    
end
