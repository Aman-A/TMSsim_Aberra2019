function plotActTermsOnMorphs2(cell_id,cell_layer,...
    nrn_model_ver,nrn_pop_name,model_prefix,model_prefix2,varargin)

mat_dir = addPaths; 
in.morph_alpha = 0.5;
in.morph_color = 'k';
in.term_act_col = 'r';
in.term_act_col2 = 'b';
in.term_lw = 1.5; 
in.cell_inds_file = 'gen_figures/cellplot_inds_ls_1_nrn_pop1';
in = sl.in.processVarargin(in,varargin); 
% Threshold data
if nargin == 0
    cell_id = 6; 
    cell_layer = 2; 
    layer_set_num = 1; 
    nrn_model_ver = 'maxH'; 
    nrn_pop_name = 'nrn_pop1';    
    mode = 1; % w1 - monophasic
    Efield_name = 'M1_PA_MCB70'; 
    model_prefix = sprintf('tms_%s_w%g_ls_%g_E_%s_P_%s',nrn_model_ver,mode,...
                            layer_set_num,Efield_name,nrn_pop);
    model_prefix2 = sprintf('tms_%s_w%g_ls_%g_E_%s_r_P_%s',nrn_model_ver,mode,...
                            layer_set_num,Efield_name,nrn_pop); 
end
%% Load data
indsi_data = load(fullfile(mat_dir,[in.cell_inds_file '.mat']));
indsi = indsi_data.inds{indsi_data.cellids==cell_id};   
if isempty(indsi)
   error('Cell %g indices not in %s',cell_id,inds_file)
end
% threshold data
threshEs_data = load(fullfile(mat_dir,'nrn_sim_data',[model_prefix '.mat'])); 
cell_model_names = threshEs_data.cell_model_names;
[~,~,celli] = intersect(cellModelNames(cell_id),cell_model_names);
init_indsi = threshEs_data.init_inds{celli}(indsi); % init_inds for cells to plot 
threshEs_data2 = load(fullfile(mat_dir,'nrn_sim_data',[model_prefix2 '.mat'])); 
cell_model_names2 = threshEs_data2.cell_model_names;
[~,~,celli2] = intersect(cellModelNames(cell_id),cell_model_names2);
init_indsi2 = threshEs_data2.init_inds{celli2}(indsi); 
% cell data
NeuronPop = loadNeuronPop(nrn_pop_name,nrn_model_ver);
%% Plot
% Plot cell morphologies
color_data = in.morph_color;
[Calli,cell_data,axon_inds] = plotCellsLayerAxon(cell_id,cell_layer,NeuronPop,indsi,color_data); % outputs coordinates of cells plotted
% make graph
edge_dists = vmag(cell_data.C([1;axon_inds]) - cell_data.C([1;cell_data.parent_inds],:)); % 1 for soma
gi = graph(cell_data.parent_inds,2:length(cell_data.C),edge_dists(2:end)); 
ax = gca;
[ax.Children(1:length(indsi)).EdgeAlpha] = deal(in.morph_alpha);
% Plot paths from soma to activated axon terminal
for j = 1:length(indsi)
    tj = Calli{j}(shortestpath(gi,1,1+find(axon_inds==init_indsi(j))),:);
    tj2 = Calli{j}(shortestpath(gi,1,1+find(axon_inds==init_indsi2(j))),:);
    % remove overlapping points to plot from bifurcation
    [~,ia,ib] = intersect(tj,tj2,'rows');
    inds = true(length(tj),1); inds(ia) = false; 
    inds2 = true(length(tj2),1); inds2(ib) = false; 
    tj = tj(inds,:); tj2 = tj2(inds2,:);         
    plot3(tj(:,1),tj(:,2),tj(:,3),'Color',in.term_act_col,'LineWidth',in.term_lw);
    plot3(tj2(:,1),tj2(:,2),tj2(:,3),'Color',in.term_act_col2,'LineWidth',in.term_lw);    
    % soma
    plot3(Calli{j}(1,1)-0.05,Calli{j}(1,2)+0.05,Calli{j}(1,3),'k.',...
       'MarkerFaceColor','k','MarkerSize',24);   
end
end
