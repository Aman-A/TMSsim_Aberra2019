% Opens tvec.txt, ap_times_all.txt,ap_netids_all.txt,vm_matrix.txt,
% secnames.txt, and threshE_angle.txt
% Reads vm vs. t vector for each section as vector in 1 x numSect cell
% array - vm_matrix
% Now reads ap times of every section in a single vector and corresponding
% netids in a second vector
function [tvec,ap_times,ap_netids,vm_matrix,params,threshE] = saveDataThresh(run_tmp_folder,data_folder_path_mat,...
                                                    simulation,cell_model,waveform,Efield)    
    time_file = fullfile(run_tmp_folder,'tvec.txt');  % time vector (tvec)
    ap_file = fullfile(run_tmp_folder,'ap_times_all.txt'); % ap times in all sections (ap_times)
    ap_netids_file = fullfile(run_tmp_folder,'ap_netids_all.txt'); % net id of each event in ap_file (ap_netids)    
    vm_bin_file = fullfile(run_tmp_folder,'vm_data_bin.txt'); % membrane potentials at each time point    
    threshE_file = fullfile(run_tmp_folder,'threshE.txt');
    % Read time vector
    tvec = nrnVread(time_file); % Read binary time vector
    % Read ap_times vector
    ap_times = nrnVread(ap_file); % Read binary spike time data 
    ap_times(ap_times<1e-5) = []; % if no spikes, leave empty
    ap_times = ap_times'; % Change to row vector
    % Read netids vector
    ap_netids = nrnVread(ap_netids_file);
    ap_netids(ap_netids<1e-5) = []; % same as above 
    ap_netids = ap_netids'; %change to row vector to match ap_times
    % open membrane potential over time for each section
    vmraw = nrnVread(vm_bin_file);
    num_timepoints = vmraw(1);    
    num_sections = vmraw(2);
    % all arrays but first have extraneous number at beginning
    % first array has num_timepoints and num_sections, remove all numbers
    vmraw = vmraw(3:end);
    %vmraw(abs(vmraw)<=1e-20) = [];
    vmraw(num_timepoints+1:num_timepoints+1:(num_timepoints+1)*(num_sections-1))=[]; % remove 8.48e-314 junk
    vmraw = reshape(vmraw,num_timepoints,num_sections);    
    if length(tvec) == size(vmraw,1)
        vm_matrix = mat2cell(vmraw,num_timepoints,ones(1,num_sections));    
        vm_matrix = cellfun(@(x) assign1(x,x(2)),vm_matrix,'UniformOutput',false); % set first value to first point
    elseif length(tvec) == (size(vmraw,1)-1) % in case (on cluster) vm vectors have extra element
        vmraw = vmraw(2:end,:); % remove first duplicate value
        vm_matrix = mat2cell(vmraw,num_timepoints-1,ones(1,num_sections));            
    else
        error('Mismatch between Vm time points and tvec time points');
    end
    fprintf('Read %g Vm time points from %g sections\n',num_timepoints,num_sections);        
    %% Read Threshold  
    fprintf('threshE_file = %s\n',threshE_file); 
    fid3 = fopen(threshE_file,'r');
    threshE_angle = fscanf(fid3,'%f');
    fclose(fid3);
    threshE = threshE_angle(1);
    threshE = threshE; %#ok variable used in save call below
    fprintf('Threshold E at theta = %g, phi = %g deg: %.3f\n',threshE_angle(2), threshE_angle(3),threshE);
    % save to .mat file
    save(fullfile(data_folder_path_mat,'threshE.mat'),'threshE');        
    % save to .mat files
    save(fullfile(data_folder_path_mat,'Vm.mat'),'tvec','vm_matrix');
    save(fullfile(data_folder_path_mat,'APtimes.mat'),'ap_times','ap_netids');    
    params = struct('simulation',simulation,'waveform',waveform,'cell_model',...
        cell_model,'Efield',Efield);       
    save(fullfile(data_folder_path_mat,'params.mat'),'params');
    fprintf('Saved threshE, Vm, APtimes, and params to %s\n',data_folder_path_mat);
end

function x = assign1(x,v)
    x(1) = v;
end
