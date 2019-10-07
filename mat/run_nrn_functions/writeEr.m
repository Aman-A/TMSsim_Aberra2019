% Input matlab directory, layer set number, name of E_file, and cell_id,
% writes Er.txt to 
function Er = writeEr(mat_dir,run_params_fold,cell_model,Efield)    
    % Open E-field data for cell    
    Edata_dir = fullfile(mat_dir,'output_data','nrn_efields',...
        sprintf('layer_set_%g',Efield.layer_set_num),Efield.E_file);    
    nrn_pop_name_full = [Efield.nrn_pop_name '_' cell_model.nrn_model_ver];
    cellfile_name = sprintf('cell%g.mat',cell_model.cell_id);         
    Ecell_data = load(fullfile(Edata_dir,nrn_pop_name_full,cellfile_name)); 
    Ecell = Ecell_data.Ecell;   
    % get phis  
    layer_set_folder = fullfile(mat_dir,'output_data','layer_data');
    nrn_pop_folder = fullfile(layer_set_folder,cell_model.nrn_model_ver);    
    NeuronPop_data = load(fullfile(nrn_pop_folder,[nrn_pop_name_full '.mat'])); 
    NeuronPop = NeuronPop_data.NeuronPop; 
    for i = 1:length(NeuronPop.cell_ids)
        if any(NeuronPop.cell_ids{i}==cell_model.cell_id)
            cell_layer = i; 
            cell_ind = NeuronPop.cell_ids{i}==cell_model.cell_id; 
            fprintf('phis from cell_id = %g, cell %g in layer %g of %g\n',cell_model.cell_id,find(cell_ind),cell_layer,length(NeuronPop.cell_ids)); 
            break; 
        end
    end    
    phi = NeuronPop.phis{cell_layer}{cell_ind}(Efield.E_position); % CCW rotation from x = 0 about z 
    % Extract E-field vectors for this position
    cell_normal = Ecell{Efield.E_position}(1,:); % first row is cell_normal
    Er = Ecell{Efield.E_position}(2:end,:); % rest of array is compartment E-vecs     
    Er = reorientEfield(cell_normal,phi,Er); % rotate E-field into local coordinate system 
    % Save to Er.txt in run_params_fold
    fid = fopen(fullfile(run_params_fold,'Er.txt'),'w');
    fprintf(fid,'%f %f %f\n',Er');
    fclose(fid); 
    fprintf('Saved E-field to %s\n',run_params_fold);    
end