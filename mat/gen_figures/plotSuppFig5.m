function plotSuppFig5(save_fig)
%PLOTSUPPFIG5 Plots minimum thresholds of full Blue Brain library after
%applying adult, human modifications described in Aberra et al. 2018 for
%Supp. Fig. 4
%
% AUTHOR    : Aman Aberra 
if nargin == 0
   save_fig = 0;  
end
mat_dir = addPaths;
model_prefix = 'utms_maxH_w1_ss1uv_vr';
data_cell_ids = 1:1035; % used all 1035 publically available Blue Brain models
cell_model_file1 = 'cell_models'; % cells included in main simulations
cell_model_file2 = 'cell_models_all'; % all cells from Blue Brain library
plot_layer_mtypes = {{'L23_LBC','L23_NBC','L23_SBC','L23_PC'};%L2/3
                      {'L4_LBC','L4_NBC','L4_SBC','L4_PC','L4_SP','L4_SS'};%L4
                      {'L5_LBC','L5_NBC','L5_SBC','L5_TTPC1','L5_TTPC2','L5_UTPC',};%L5
                      {'L6_LBC','L6_NBC','L6_SBC','L6_BPC','L6_IPC','L6_TPC_L1',... %L6
                                'L6_TPC_L4','L6_UTPC'}}; 
num_layers = length(plot_layer_mtypes);
%% Plot options
exc_marker = 'p'; % excitatory marker
m_size_big_exc = 105; % marker size for included, excitatory cells
m_size_small_exc = 55;  % marker size for not included, excitatory cells
inh_marker = 'o'; % inhibitory marker
m_size_big_inh = 70;  % marker size for included, inhibitory cells
m_size_small_inh = 25; % marker size for not included, inhibitory cells
incl_edge_col = 'k'; % edge color for included cells
incl_edge_alpha = 1;
not_incl_edge_col = 0.3*ones(1,3); % edge color for non included cells
not_incl_edge_alpha = 1;
m_lw_big = 0.5; % edge line width for included cells
m_lw_small = 0.1; % edge line width for not included cells (if on)
rand_jit_max = 0.2; % range of rand distribution for point jitter (for marker_mode = 2)
rng(50);     % seed rn generator 
colors = jet(num_layers+1); colors(4,:) = [0.5961 0.3059 0.6392]; 
colors = colors(2:end,:); % remove L1
y_font_size = 10; 
x_font_size = 8; 
font_name = 'Times';
y_lims = [100 500]; 
ax_lw = 1; 
fig_size = [21.25,11.5]; % cm
%% load cell data table
cell_model_data1 = load(fullfile(mat_dir,'cell_data',[cell_model_file1 '.mat'])); 
T = cell_model_data1.T;
T.exc_or_inh = [repmat('Inhibitory',5,1);repmat('Excitatory',5,1);
                repmat('Inhibitory',5,1);repmat('Excitatory',10,1)];
incl_cell_model_names = T.Row;
cell_model_data2 = load(fullfile(mat_dir,'cell_data',[cell_model_file2 '.mat'])); 
T_all = cell_model_data2.T;

cell_ids = cell(num_layers,1);
exc_or_inh = cell(num_layers,1); % 1 exc, 0 inh
incl_cells = cell(num_layers,1); 
for i = 1:num_layers
    num_cells_in_layer = length(plot_layer_mtypes{i});
    cell_ids{i} = cell(1,num_cells_in_layer); 
    exc_or_inh{i} = zeros(num_cells_in_layer,1);
    incl_cells{i} = cell(num_cells_in_layer,1); 
   for j = 1:num_cells_in_layer       
       [cell_ids{i}{j},cell_model_names] = getCellIds(plot_layer_mtypes{i}{j},'layer_mtype',T_all); % get cell ids                     
       incl_ind = ismember(cell_model_names,incl_cell_model_names);
       if all(incl_ind) % all cells of this layer_mtype included          
           incl_cells{i}{j} = 1; % else leave 0
       elseif any(incl_ind) % some cells but not all included
           incl_cells{i}{j} = incl_ind; % plot these separately
       else
           incl_cells{i}{j} = 0; 
       end
       if strcmp(T_all.exc_or_inh(cell_ids{i}{j}(1)),'Excitatory')
           exc_or_inh{i}(j) = 1; % else leave 0
       end
   end
end
%% Load data
data_fold = fullfile(mat_dir,'nrn_sim_data'); % save .mat file here
data_file = sprintf('%s_%g-%g.mat',model_prefix,data_cell_ids(1),data_cell_ids(end));
data = load(fullfile(data_fold,data_file)); 
threshEs = data.threshEs;
minthreshEs = min(threshEs,[],1)'; % cell_id from T_all indexes threshEs (goes from 1 to 1035)
%% Plot data
fig = figure('Color','w'); 
fig.Units = 'centimeters';
fig.Position(3:4) = fig_size; 
cnt = 1; % plot index
for i = 1:num_layers
    num_cells_in_layer = length(cell_ids{i});
    col = colors(i,:);
    for j = 1:num_cells_in_layer                     
        if sum(incl_cells{i}{j}) % included, plot bigger 
            m_lw = m_lw_big; 
            edge_alpha = incl_edge_alpha;
            if exc_or_inh{i}(j)
               m_size = m_size_big_exc;               
               edge_col = incl_edge_col;
               marker = exc_marker; 
            else
                m_size = m_size_big_inh;
                edge_col = incl_edge_col;
                marker = inh_marker; 
            end
            if length(incl_cells{i}{j}) > 1 % two separate scatter plots for this cell                
                % plot not included cells in this layer_mtype (assume inh)
                scatter(cnt*ones(size(minthreshEs(cell_ids{i}{j}(~incl_cells{i}{j})))),minthreshEs(cell_ids{i}{j}(~incl_cells{i}{j})),...
                    m_size_small_inh,'MarkerEdgeColor',not_incl_edge_col,'Marker',inh_marker,'MarkerFaceColor',col,...
                    'LineWidth',m_lw_small,'jitter','on','jitterAmount',rand_jit_max,'MarkerEdgeAlpha',not_incl_edge_alpha); hold on;
                % plot included cells (on top)
                scatter(cnt*ones(size(minthreshEs(cell_ids{i}{j}(incl_cells{i}{j})))),minthreshEs(cell_ids{i}{j}(incl_cells{i}{j})),...
                    m_size,'MarkerEdgeColor',edge_col,'Marker',marker,'MarkerFaceColor',col,...
                    'LineWidth',m_lw,'jitter','on','jitterAmount',rand_jit_max,'MarkerEdgeAlpha',incl_edge_alpha); hold on; 
            else
                scatter(cnt*ones(size(minthreshEs(cell_ids{i}{j}))),minthreshEs(cell_ids{i}{j}),...
                    m_size,'MarkerEdgeColor',edge_col,'Marker',marker,'MarkerFaceColor',col,...
                    'LineWidth',m_lw,'jitter','on','jitterAmount',rand_jit_max,'MarkerEdgeAlpha',edge_alpha); hold on; 
            end
        else % not included plot smaller
            m_lw = m_lw_small; 
            edge_alpha = not_incl_edge_alpha;
            if exc_or_inh{i}(j)
               m_size = m_size_small_exc;               
               edge_col = not_incl_edge_col;
               marker = exc_marker; 
            else
                m_size = m_size_small_inh;
                edge_col = not_incl_edge_col;
                marker = inh_marker; 
            end
            scatter(cnt*ones(size(minthreshEs(cell_ids{i}{j}))),minthreshEs(cell_ids{i}{j}),...
            m_size,'MarkerEdgeColor',edge_col,'Marker',marker,'MarkerFaceColor',col,...
            'LineWidth',m_lw,'jitter','on','jitterAmount',rand_jit_max,'MarkerEdgeAlpha',edge_alpha); hold on; 
        end        
        cnt = cnt+1; % next layer_mtype
    end
    if i < num_layers
        plot(cnt-0.5*ones(1,2),y_lims,'k','LineWidth',ax_lw);
    end
end
ylabel('Minimum threshold (V/m)'); 
ax = gca;
ax.YLabel.FontSize = y_font_size;
ax.XLabel.FontSize = x_font_size;
ax.FontName = font_name; 
ax.XTick = 1:(cnt-1); 
cell_names = [plot_layer_mtypes{:}];
ax.XTickLabel = cell_names; 
ax.XTickLabelRotation = 45; 
ax.TickLabelInterpreter = 'none';
ax.TickLength(2) = 0.025; 
ax.LineWidth = ax_lw; 
ax.YColor = 'k'; ax.XColor = 'k';
ax.YLim = y_lims;
ax.YGrid = 'on';
ax.XGrid = 'off';
if save_fig
   fig_fold = fullfile(mat_dir,'figures');
   fig_name = 'SuppFig5';   
   savefig(fig,fullfile(fig_fold,[fig_name '.fig']));    
   print(fig,fullfile(fig_fold,[fig_name '.eps']),'-depsc2','-cmyk','-painters');
end
end
