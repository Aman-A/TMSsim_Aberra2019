function simTMSlayerThresh(cell_id,model_prefix,varargin)
% SIMTMSLAYERTHRESH Threshold loop for TMS waveforms using FEM E-field that
%spans full range of cell positions within layer for given cell_id
%   Inputs 
%   ------ 

%   Optional Inputs 
%   --------------- 
%   Outputs 
%   ------- 
%   Examples 
%   --------------- 
% AUTHOR    : Aman Aberra
par_on = 1;  % set to 1 to use parfor loop
if nargin == 0
   cell_id = 1;
   model_prefix = 'tms_maxH_w1_LS_1_E_M1_PA_MCB70_P_nrn_pop1';   
end
in.layer_set_num = 1; % default layer set
in.Efield_solution = 'M1_PA_MCB70';
in.nrn_model_ver = 'maxH';
in.nrn_pop_name = 'nrn_pop1'; 
in.ss_init = 1; % default steady state initialization
in.dt = 0.005; 
in.tstop = 1;
in.mode = 1; % TMS waveform (default monophasic)
in.del = 0.005; 
in.dur = 0.03; % doesn't matter for conventional TMS pulses
in.vm_rmode = 1; % record vm from soma
in.spike_rmode = 4; % record spikes from all compartments
in.over_write = 0; 
in.replace_axon = 0; % default leave axon intact 
in = sl.in.processVarargin(in,varargin);
fprintf('Running model: %s\n',model_prefix); 
%% Set up simulation
mat_dir = addPaths; 
Edata_dir = fullfile(mat_dir,'output_data','nrn_efields',sprintf('layer_set_%g',in.layer_set_num),in.Efield_solution);
nrn_pop_name_full = [in.nrn_pop_name '_' in.nrn_model_ver];
cellfile_name = sprintf('cell%g.mat',cell_id); 
Ecell_data = load(fullfile(Edata_dir,nrn_pop_name_full,cellfile_name)); 
Ecell = Ecell_data.Ecell;    
num_runs = length(Ecell); % or length(thetas) - every phi orientation tested at num_thetas-2, one orientation tested at 0 and 180�
Esomas = cell2mat(cellfun(@(x) x(2,:), Ecell,'UniformOutput',0)); % E-field vectors at soma for every position (1st row is cell_normal)   
[phis,lambdas,mags] = cart2sph(Esomas(:,1),Esomas(:,2),Esomas(:,3)); % outputs azimuth and elevation of vectors in rad, and magnitude
% Specify model parameters in m structure   
% Static Parameters
   m.cell_id = cell_id; % 1 - 25
   m.cell_model_name = cellModelNames(cell_id);   
   % Simulation Parameters
   m.ss_init = in.ss_init;
   m.temp = -100; % set by loadParams
   m = outputMorphParams(in.nrn_model_ver,m); % adds morph params to m
   m.replace_axon = in.replace_axon; % override replace_axon from outputMorphParams
   m.vm_record_mode = in.vm_rmode; 
   m.spike_record_mode = in.spike_rmode; 
   m.load_potentials = 1; % get pseudo-Ve using Er.txt      
   m.v_init = -75; % mV (initialization before steady state init())
   m.dt = in.dt; % ms
   m.tstop = in.tstop; % ms
   m.mode = in.mode;
   m.dur = in.dur;
   m.del = in.del;
   m.amp = 100; % amplification of E-field - starting amplitude         
   m.E_file = in.Efield_solution;
   m.layer_set_num = in.layer_set_num;  
   m.nrn_pop_name = in.nrn_pop_name;   
% Variable Parameters  
   E_position = 1:num_runs;
   E_theta = 90-lambdas*180/pi; % convert to polar angle (angle from z) in deg
   E_phi = phis*180/pi; % convert to deg
   E_mag = mags; 
   m_array = mArrayBuilder(m,E_position,E_theta,E_phi,E_mag); % generates m structure array for each loop iteration
%% Set up parpool
numCPUs = feature('numCores'); 
fprintf('Number of CPUs requested = %g\n',numCPUs);    
over_write = in.over_write; % avoid broadcasting in struct
if numCPUs > 1 && par_on
    if numCPUs > 4 % on cluster
        pc_storage_dir = fullfile(mat_dir,'pc_storage',[getenv('SLURM_JOB_ID') '_' num2str(cell_id)]);        
        mkdir(pc_storage_dir); 
        pc = parcluster('local');    
        pc.JobStorageLocation =  pc_storage_dir;    
    else
       pc = parcluster('local'); % use default 
    end
    poolobj = parpool(pc,numCPUs);
    iter = 1:num_runs; % iterator
    succeeded = false(size(iter)); % which runs have succeeded
    threshEs = zeros(num_runs,1);
    init_inds = zeros(num_runs,1);
    % initialize structs to save after parfor 
    i = 1; 
    model_name = sprintf('%s_pos_%g',model_prefix,m_array(i).E_position);
    [simulation,cell_model,waveform,Efield] = loadParams(model_prefix,model_name,m_array(i)); % load parameter structures
    while ~all(succeeded)
        todo = iter(~succeeded);
        part_threshEs = zeros(length(todo),1); % remaining results
        part_init_inds = zeros(length(todo),1);
        partsucceeded = false(size(todo)); % flags for iterations which succeeded
        parfor (falsei = 1:sum(~succeeded),poolobj.NumWorkers)
            i = todo(falsei); % get original index of this run
            try
                model_namei = sprintf('%s_pos_%g',model_prefix,m_array(i).E_position);
                % load parameter structures
                [simulationi,cell_modeli,waveformi,Efieldi] = ...
                    loadParams(model_prefix,model_namei,m_array(i)); 
                [~, ~,ap_netids,~,~,threshE] = runNrn(simulationi,cell_modeli,...
                                            waveformi,Efieldi,over_write);
                if threshE ~= 0
                    part_threshEs(falsei) = threshE;
                    part_init_inds(falsei) = ap_netids(1);
                    partsucceeded(falsei)=true;
                else % simulation exited because taking too long, set threshE to 0
                    part_threshEs(falsei) = nan;
                    part_init_inds(falsei) = nan;
                    partsucceeded(falsei)=true;
                end
                fprintf('Finished run %g\n',i); 
            catch
                fprintf('Run %g failed\n',i);
                fprintf('falsei = %g. todo(falsei) = %g\n',falsei,todo(falsei)); 
            end
        end        
        threshEs(~succeeded) = part_threshEs;
        init_inds(~succeeded) = part_init_inds;
        succeeded(~succeeded) = partsucceeded;
    end    
    delete(poolobj); 
    if numCPUs > 4 % on cluster
        rmdir(pc_storage_dir,'s'); % delete parfor temporary files           
    end
else % SERIAL LOOP - 1 cpu or nan (if env variable doesn't exist)
    threshEs = zeros(num_runs,1);
    init_inds = zeros(num_runs,1);
    model_name = sprintf('%s_pos_%g',model_prefix,m_array(1).E_position);
    [simulation,cell_model,waveform,Efield] = loadParams(model_prefix,model_name,m_array(1)); % load parameter structures
    for i = 1:num_runs
        model_namei = sprintf('%s_pos_%g',model_prefix,m_array(i).E_position);
        % load parameter structures
        [simulationi,cell_modeli,waveformi,Efieldi] = ...
            loadParams(model_prefix,model_namei,m_array(i));
        [~, ~,ap_netids,~,~,threshE] = runNrn(simulationi,cell_modeli,...
            waveformi,Efieldi,over_write);
        if threshE ~= 0 %simulation successfully completed and gave positive threshold
            threshEs(i) = threshE;
            init_inds(i) = ap_netids(1); % get init index
        else % simulation exited because taking too long, set threshE to 0
            threshEs(i) = nan;
            init_inds(i) = nan; % get init index
        end
    end
end
cell_model_name = m.cell_model_name; 
params = struct('simulation',simulation,'cell_model',cell_model,'waveform',...
                waveform,'Efield',Efield); 
cell_data_file = fullfile(mat_dir,'nrn_sim_data',model_prefix,cell_model_name,sprintf('%s_thresh.mat',cell_model_name)); 
save(cell_data_file,'cell_model_name','params','threshEs','init_inds','E_mag'); 
end
