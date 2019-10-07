function [tvec, ap_times,ap_netids,vm_matrix,params,threshE] = runNrn(simulation,cell_model,waveform,Efield,over_write)
%RUNNRN Calls NEURON to run simulation using parameters given in simulation,
%waveform, current_inj, and Efield structures
% 
% AUTHOR    : Aman Aberra
if nargin < 5
    over_write = 0; 
end
% set up directories            
mat_dir = addPaths();     
data_dir = mat_dir; % can change this to different directory for data storage    
nrn_direc = fullfile(mat_dir,'../nrn');          
% Create data folder in MATLAB and NEURON directories       
[data_folder_path_mat,~,~] = createDataFolder(data_dir,simulation,cell_model); % make and get matlab data directory
[~,run_tmp_fold,run_params_fold] = createDataFolder(nrn_direc,simulation,cell_model); % make and get neuron data directory         
run_name = simulation.run_name; 
nrn_file = fullfile(nrn_direc,'init_field_stim.hoc'); 
%% Check if simulation was already run, skip if yes
run_mat_data_files = dir(data_folder_path_mat); % .mat data output files
run_mat_data_files = run_mat_data_files(arrayfun(@(x) ~strcmp(x.name(1),'.'),run_mat_data_files)); % remove '.' hidden files    
num_data_files = 4; % 4 output files, includes threshEs/deltaVms output
if length(run_mat_data_files) >= num_data_files && ~over_write     
    if Efield.load_potentials == 0
        fprintf('%g output files found for theta = %g, phi = %g, skipping\n',length(run_mat_data_files),Efield.theta,Efield.phi);        
    elseif Efield.load_potentials == 1
        fprintf('%g output files found for position %g, skipping\n',length(run_mat_data_files),Efield.E_position);        
    end
    [tvec,ap_times,ap_netids,vm_matrix,params,threshE] = loadRawData(data_folder_path_mat);            
else
    if length(run_mat_data_files) >= num_data_files && over_write
        fprintf('Overwriting existing data\n');  
    end
    % Define Waveform Parameters and write t and E vectors to tvec.txt and
    % Evec.txt
    if strcmp(waveform.type,'tms') % TMS waveform
        [tvec,Evec] = TMSwave(simulation.dt,simulation.tstop,waveform,0); % create time and E field vectors
        writeVectorBin(run_params_fold,tvec,Evec);
%             fprintf('Saved TMS waveform to model directory\n');
    else
        error('Unrecognized waveform type\n')
    end            
    % Write model parameters to defineParams.hoc
    writeParamsHoc(run_params_fold,simulation,cell_model,Efield,waveform)                                
    % Write E-field/V data (if necessary)
    if Efield.load_potentials
        % write E-field vecs to params/<run_name>/Er.txt
        writeEr(mat_dir,run_params_fold,cell_model,Efield);         
    end
    cd(nrn_direc); % launch init_field_stim.hoc from here      
    % NEURON command %        
    commands = sprintf(' -c "run_name=\\"%s"\\" -c "loadFiles()" -c "find_thresh()" -c "quit()"',run_name);                    
    tic
    if ispc
        options = sprintf('-NSTACK %g -NFRAME %g ',100000,20000);
        system(['C:\nrn\bin\nrniv.exe -nobanner ' options nrn_file commands]);
    elseif ismac
        options = sprintf('-nopython -NSTACK %g -NFRAME %g ',100000,20000);
        system(['./special ' options nrn_file commands]);
    else
        options = sprintf('-NSTACK %g -NFRAME %g ',100000,20000);
        system(['./special ' options nrn_file commands]);
    end
    time_sim = toc;
    fprintf('Simulation time was %.3f sec\n',time_sim); 
    cd(mat_dir)  
    % ------ %        
    % saves data to .mat files in MATLAB directory
    [tvec, ap_times,ap_netids,vm_matrix,params,threshE] = saveDataThresh(run_tmp_fold,data_folder_path_mat,...
                    simulation,cell_model,waveform,Efield); 
end    
rmdir(run_tmp_fold,'s') % delete old tmp folder
rmdir(run_params_fold,'s'); % delete current run's params folder
end