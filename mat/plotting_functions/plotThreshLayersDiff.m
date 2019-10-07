% (1 - 2)/2
function plotThreshLayersDiff(layers,cell_ids,threshEs,plot_exc,norm_data,plot_percent,plot_cell_layers,plot_layers_single_fig,save_figs)
    current_dir = '/Users/amanaberra/Documents/MATLAB/PlaceNeurons3D';  
    % (model_prefix - model_prefix2)/model_prefix2
%     model_prefix = 'thresh_TMS_maxH_rax2PCs_w9_120us_LS_7_E_M1_DostB_PA_P8_2mm_P_nrn_pop1';
%     model_prefix2 = 'thresh_TMS_maxH_rax2PCs_w9_120us_LS_7_E_M1_DostB_PA_P8_2mm_r_P_nrn_pop1';
    model_prefix = 'thresh_TMS_maxH_rax2PCs_w3_LS_9_E_M1_SomB_PA_P1_2mm_P_nrn_pop1';
    model_prefix2 = 'thresh_TMS_maxH_rax2PCs_w3_LS_9_E_M1_SomB_PA_P1_2mm_s_P_nrn_pop1';
    % takes difference of model_prefix - model_prefix2
    if nargin == 0
        cell_ids = {[];[];[];[32];[]}; 
%         cell_ids = { [42]; [17 18]; [27 41]; [32 57];  [52]};  
%         cell_ids = {42; 17; 27; 32; 52}; 
%         cell_ids = {42:46;17:21;27:31;32:36;52:56};
%         cell_ids = {[];17:21;27:31;32:36;[]};
        plot_cell_layers = 1;
        plot_layers_single_fig = 0;        
        save_figs = 0; 
       layer_set_num = 9;
       plot_exc= 0; % set to 1 to convert to 1/threshold = excitability
       norm_data  = 0; % set to 1 to normalize to min threshold/max excitability (within model_prefix)
       plot_percent = 1; % set to 1 to divide difference by reference (model_prefix2) and convert to percent 
%        shift_dir = [0 0 0]; 
       shift_dir = [0 -35 0]; 
       % load layers
       layerTable = readtable(fullfile(current_dir,'Layer_sets.csv'));
        Mesh_File = layerTable.Mesh_File{layer_set_num};
        ROI_name = layerTable.ROI_name{layer_set_num};
        layer_set = sprintf('Layer_set_%g',layer_set_num);
        if ispc
            output_folder = 'E:\PlaceNeurons3D\OutputData';
        else
            output_folder = fullfile(current_dir,'OutputData');
        end
        MeshROI_folder = fullfile(output_folder,Mesh_File,[Mesh_File '_' ROI_name]); 
        layer_set_folder = fullfile(MeshROI_folder,layer_set);
        layer_file = fullfile(layer_set_folder,[layer_set '.mat']);
        fprintf('Loading layer file %s\n',layer_file); 
        layer_data = load(layer_file);
        layers = layer_data.layers; 
        % load threshE data        
        if ispc
            input_folder = 'E:\PlaceNeurons3D\InputData';
        else
            input_folder = fullfile(current_dir,'InputData');
        end
        % load first data set
        nrn_sim_data_file = fullfile(input_folder,'NeuronFEMsimData',[model_prefix '_sectypes.mat']);
        fprintf('Loading %s NEURON data\n',nrn_sim_data_file);
        nrn_sim_data = load(nrn_sim_data_file); 
        threshEs = nrn_sim_data.threshEs;   
        cell_model_names = nrn_sim_data.cell_model_names;
        % load 2nd data set
        nrn_sim_data_file = fullfile(input_folder,'NeuronFEMsimData',[model_prefix2 '_sectypes.mat']);
        fprintf('Loading %s NEURON data\n',nrn_sim_data_file);
        nrn_sim_data = load(nrn_sim_data_file); 
        threshEs2 = nrn_sim_data.threshEs;   
        %cell_model_names2 = nrn_sim_data.cell_model_names;
    else
        layer_set_num = layers(1).layer_set_num; 
    end    
    fprintf('Plotting difference map on layer set %g\n',layer_set_num);     
    %% Convert to vertex data
    threshEs = getAllVertexCData(threshEs,layers,cell_ids,cell_model_names); 
    threshEs2 = getAllVertexCData(threshEs2,layers,cell_ids,cell_model_names); 
    fprintf('Plotting thresholds on layer set %g\n',layer_set_num);   
    %% Plot separate figures for each cell on layer surface
    num_layers = length(layers);
    if plot_cell_layers 
        if plot_exc
           threshEs = cellfun(@(x) 1./x,threshEs,'UniformOutput',0); 
           threshEs2 = cellfun(@(x) 1./x,threshEs2,'UniformOutput',0); 
           fprintf('Converting to excitability\n');
           if norm_data
               max_exc_cells = cellfun(@max,threshEs);
               max_exc_cells2 = cellfun(@max,threshEs2);
               global_max = max(max_exc_cells);
               global_max2 = max(max_exc_cells2);
               threshEs = cellfun(@(x) x/global_max,threshEs,'UniformOutput',0); 
               threshEs2 = cellfun(@(x) x/global_max2,threshEs2,'UniformOutput',0); 
               fprintf('Normalized to max excitability\n');
               title_str = sprintf('Normalized Excitability');
               units = 'normalized to max';
           else
               title_str = sprintf('Excitability');
               units = '1/threshold';
           end
        else
            if norm_data
                min_thresh_cells = cell2mat(cellfun(@min, threshEs,'UniformOutput',0));
                min_thresh_cells2 = cell2mat(cellfun(@min, threshEs2,'UniformOutput',0));
                global_min = min(min_thresh_cells); % min threshold of all cells
                global_min2 = min(min_thresh_cells2); % min threshold of all cells
                threshEs = cellfun(@(x) x/global_min,threshEs,'UniformOutput',0);
                threshEs2 = cellfun(@(x) x/global_min2,threshEs2,'UniformOutput',0);
                fprintf('Normalized to min threshold\n');
                title_str = sprintf('Threshold stimulator output'); 
                units = 'normalized to min';
            else
                title_str = sprintf('Threshold stimulator output'); 
                units = 'A/us';
            end
        end    
        % Plot layers
        for i = 1:num_layers        
           num_cells_in_layer = length(cell_ids{i}); 
            for j = 1:num_cells_in_layer               
               cell_model_name = cellmodelnames(cell_ids{i}(j));            
               thresh_ind = strcmp(cell_model_name,cell_model_names); 
               threshEs_ij = threshEs{thresh_ind} - threshEs2{thresh_ind}; % take difference of data
               if plot_percent
                   threshEs_ij = 100*threshEs_ij./threshEs2{thresh_ind}; % convert to percent
                   units = '%';
               end
               display_cell_name = cell_model_name; 
               ind=(display_cell_name=='_');
               display_cell_name(ind)=' '; % remove underscore for figure title
%                figure('units','normalized','outerposition',[0 0 1 1],'paperunits','in','paperposition',[0 0 16 9],'Color',[1 1 1],'PaperPositionMode','auto');                      
                figure; 
               layer_surf = layers(i).surface;           
               p = patch(layer_surf);
               p.FaceVertexCData = threshEs_ij;
               p.FaceColor = 'interp';
               p.CDataMapping = 'scaled'; 
               p.EdgeColor = 'none';
               camlight; lighting gouraud;
               axis equal; axis tight; axis off; 
%                view([ -137.5 28]);
               view([-98.5 32]); 
               hold on;
    %            if normals_on
    %                cell_origins = layers(i).cell_origins;           
    %                cell_normals = layers(i).cell_normals;
    %                quiver3(cell_origins(:,1),cell_origins(:,2),cell_origins(:,3),cell_normals(:,1),cell_normals(:,2),cell_normals(:,3),'k')
    %            end
               xlabel('x'); ylabel('y'); zlabel('z');
               title(sprintf('%s difference (%s):\n %s - %s\nLayer %g of %g - Cell: %s',title_str,units,model_prefix,model_prefix2,i,num_layers,display_cell_name),'Interpreter','none');
               if plot_exc                                      
                   colormap(bluewhitered); 
               else
%                    colormap(flipud(hot))                   
%                    colormap(flipud(bluewhitered))
%                     caxis([-10 40]); 
                    colormap(bluewhitered);                     
               end           
               colorbar;
               if save_figs
                   model_prefix_fig_fold = fullfile(current_dir,'Figures',model_prefix);
                   if exist(model_prefix_fig_fold,'dir') == 0
                      mkdir(model_prefix_fig_fold)
                      fprintf('Made figure folder %s\n',model_prefix_fig_fold);
                   else
                       fprintf('Saving %s of %s\n',cell_model_name,model_prefix); 
                   end
                   savefig(fullfile(model_prefix_fig_fold,sprintf('%s_norm%g.fig',cell_model_name,norm_data)));
                   saveas(gcf,fullfile(model_prefix_fig_fold,sprintf('%s_norm%g.png',cell_model_name,norm_data)));                
               end
           end                
        end
    end
    %% Plot single figure with all layer surfaces
    num_layers = length(layers);
    if plot_layers_single_fig         
%         shift_dir = [0,-14,0]; % vector of direction of shift 
        threshEs_ave = cell(num_layers,1); threshEs_ave2 = cell(num_layers,1); 
        for i = 1:num_layers
            cell_model_names_i = cellmodelnames(cell_ids{i}); % cell names in layer
            [~,~,thresh_inds] = intersect(cell_model_names_i,cell_model_names); % get indices of layer cells in threshEs
            if ~isempty(thresh_inds)
                threshEs_ave{i} = median(cell2mat(threshEs(thresh_inds)),2); % average each position across cell models                     
                threshEs_ave2{i} = median(cell2mat(threshEs2(thresh_inds)),2); % same for 2nd data set 
            else
                threshEs_ave{i} = 1e6; % placeholder large number to allow UniformOutput in cellfun below, no cells in this layer
                threshEs_ave2{i} = 1e6; 
            end
        end
        fprintf('Computed layer averages\n')
        if plot_exc
           threshEs_ave = cellfun(@(x) 1./x,threshEs_ave,'UniformOutput',0); 
           threshEs_ave2 = cellfun(@(x) 1./x,threshEs_ave2,'UniformOutput',0); 
           fprintf('Converting to excitability\n');
           if norm_data
               max_exc_cells = cellfun(@max,threshEs_ave);
               max_exc_cells2 = cellfun(@max,threshEs_ave2);
               global_max = max(max_exc_cells);
               global_max2 = max(max_exc_cells2);
               threshEs_ave = cellfun(@(x) x/global_max,threshEs_ave,'UniformOutput',0); 
               threshEs_ave2 = cellfun(@(x) x/global_max2,threshEs_ave2,'UniformOutput',0); 
               fprintf('Normalized to max excitability\n');
               title_str = sprintf('Normalized Excitability');
               units = 'normalized to max';
           else
               title_str = sprintf('Excitability');
               units = '1/threshold';
           end
        else
            if norm_data
                min_thresh_cells = cellfun(@min, threshEs_ave);
                min_thresh_cells2 = cellfun(@min, threshEs_ave2);
                global_min = min(min_thresh_cells); % min threshold of all cells
                global_min2 = min(min_thresh_cells2); % min threshold of all cells
                threshEs_ave = cellfun(@(x) x/global_min,threshEs_ave,'UniformOutput',0);
                threshEs_ave2 = cellfun(@(x) x/global_min2,threshEs_ave2,'UniformOutput',0);
                fprintf('Normalized to min threshold\n');
                title_str = sprintf('Threshold stimulator output'); 
                units = 'normalized to min';
            else
                title_str = sprintf('Threshold stimulator output'); 
                units = 'A/us';
            end
        end    
        
        figure('units','normalized','outerposition',[0 0 1 1],'paperunits','in','paperposition',[0 0 16 9],'Color',[1 1 1],'PaperPositionMode','auto');                      
        for i = 1:num_layers        
           num_cells_in_layer = length(cell_ids{i});
           fprintf('Layer %g of %g:\n',i,num_layers);
            if num_cells_in_layer >= 1                       
                cell_model_name1 = cellmodelnames(cell_ids{i}(1));
                cell_model_name2 = cellmodelnames(cell_ids{i}(num_cells_in_layer));
                fprintf('Plotting difference maps for %s to %s\n',cell_model_name1,cell_model_name2);               
                % format display names
                %             display_cell_name1 = cell_model_name1;
                %             display_cell_name2 = cell_model_name2;
                %             ind1=(display_cell_name1=='_');
                %             ind2=(display_cell_name2=='_');
                %             display_cell_name1(ind1)=' '; % remove underscore for figure title
                %             display_cell_name2(ind2)=' ';
                % plot surface
                layer_surf = layers(i).surface;
                layer_surf.vertices = layer_surf.vertices + (i-1)*repmat(shift_dir,size(layer_surf.vertices,1),1);
                p = patch(layer_surf);                
                threshEs_ave_i = threshEs_ave{i} - threshEs_ave2{i};
                if plot_percent
                   threshEs_ave_i = 100*threshEs_ave_i./threshEs_ave2{i}; 
                   units = '%';
                end
                p.FaceVertexCData = threshEs_ave_i;
                p.FaceColor = 'interp';
                p.CDataMapping = 'scaled';
                p.EdgeColor = 'none';
                hold on;
                %            if normals_on
                %                cell_origins = layers(i).cell_origins;
                %                cell_normals = layers(i).cell_normals;
                %                quiver3(cell_origins(:,1),cell_origins(:,2),cell_origins(:,3),cell_normals(:,1),cell_normals(:,2),cell_normals(:,3),'k')
                %            end
            else
                fprintf('No cells in layer\n');
            end                                     
        end
        title(sprintf('Mean Layer %s difference (%s):\n%s - %s. Layer set %g',title_str,units,model_prefix,model_prefix2,layer_set_num),'Interpreter','none');
        if plot_exc 
%             caxis([-30 30]) 
%             caxis([-20 13])
            colormap(bluewhitered(1000));
%               colormap(flipud(bluewhitered));
        else
            colormap(redwhiteblue(1000));
%             caxis([1 20])            
        end
        ax = gca; ax.FontSize = 14;
        camlight; lighting gouraud;
        axis equal; axis tight;
%         view([ -137.5 28]);
        view([-111.9 34.4]); 
        xlabel('x'); ylabel('y'); zlabel('z');                                                        
        cb = colorbar; cb.FontSize = 20; 
        if save_figs
            model_prefix_fig_fold = fullfile(current_dir,'Figures',model_prefix);
            if exist(model_prefix_fig_fold,'dir') == 0
                mkdir(model_prefix_fig_fold)
                fprintf('Made figure folder %s\n',model_prefix_fig_fold);
            else
                fprintf('Saving %s of %s\n',cell_model_name,model_prefix);
            end
            savefig(fullfile(model_prefix_fig_fold,sprintf('allLayers_norm%g.fig',norm_data)));
            saveas(gcf,fullfile(model_prefix_fig_fold,sprintf('allLayers_norm%g.png',norm_data)));
        end
    end    


end