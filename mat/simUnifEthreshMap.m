function simUnifEthreshMap(cell_id,model_prefix,varargin)
par_on = 1;  % set to 1 to use parfor loop
if nargin == 0
   cell_id = 1;
   model_prefix = 'utms_maxH_w1';   
end
in.nrn_model_ver = 'maxH'; %  cell model settings
in.dt = 0.005; 
in.tstop = 1; 
in.ss_init = 1; % steady state initialization
in.mode = 1; 
in.del = 0.005; 
in.dur = 0.03;                
in.vm_rmode = 1; % record vm from soma
in.spike_rmode = 4; % record spikes from all compartments       
in.dtheta = 15; % spacing for spherical threshold map
in.dphi = 10;         
in.over_write = 0; 
in.replace_axon = 0; % default leave axon intact
in = sl.in.processVarargin(in,varargin);
%% Set up simulation
mat_dir = addPaths;       
% generate full vectors of theta and phi values
%dtheta = 15;
thetas1 = 0:in.dtheta:180; % polar angle - note at 0 and 180, only 0deg is run, b/c of symmetry   
%dphi = 10;
phis1 = 0:in.dphi:360-in.dphi;
num_thetas = length(thetas1);
num_phis = length(phis1);   
thetas = []; 
for i = 2:num_thetas-1
    thetas = [thetas,thetas1(i)*ones(1,num_phis)]; 
end    
thetas = [0, thetas, 180]; % add in pole
phis = repmat(phis1,1,num_thetas-2);
phis = [0, phis, 0];
num_runs = length(phis); % or length(thetas) - every phi orientation tested at num_thetas-2, one orientation tested at 0 and 180�
% Specify model parameters in m structure   
% Static Parameters
   m.cell_id = cell_id; % 1 - 25
   m.cell_model_name = cellModelNames(cell_id);
   m.nrn_model_ver = in.nrn_model_ver;      
   % Simulation Parameters                     
   m.ss_init = in.ss_init;
   m.v_init = -75;    
   m.temp = -100;  % sets according to replace_ax   
   m = outputMorphParams(in.nrn_model_ver,m); % adds morph params to m
   m.replace_axon = in.replace_axon; % override replace_axon from outputMorphParams
   m.vm_record_mode = in.vm_rmode; 
   m.spike_record_mode = in.spike_rmode; 
   m.load_potentials = 0; % calc potentials in NEURON       
   m.dt = in.dt; % ms
   m.tstop = in.tstop; % ms
   m.mode = in.mode;
   m.dur = in.dur;
   m.del = in.del;
   m.amp = 100; % V/m - starting amplitude         
% Variable Parameters  
   E_theta = thetas;
   E_phi = phis;
   m_array = mArrayBuilder(m,E_theta,E_phi); % generates m structure array for each loop iteration   
% numCPUs = str2double(getenv('SLURM_CPUS_PER_TASK')); % if running on cluster   
numCPUs = feature('numCores'); 
fprintf('Number of CPUs requested = %g\n',numCPUs);    
if numCPUs > 1 && par_on % - PARALLELIZE LOOP
    %% Set up parpool
    if numCPUs > 4 % on cluster
        pc_storage_dir = fullfile(mat_dir,'pc_storage',[getenv('SLURM_JOB_ID') '_' num2str(cell_id)]); 
        mkdir(pc_storage_dir); 
        pc = parcluster('local');
        pc.JobStorageLocation =  pc_storage_dir;
    else
       pc = parcluster('local');  
    end
    poolobj = parpool(pc,numCPUs);
    iter = 1:num_runs; % iterator
    succeeded = false(size(iter)); % which runs have succeeded
    threshEs = zeros(num_runs,1);
    init_inds = zeros(num_runs,1);
    % initialize structs to save after parfor 
    i = 1; 
    model_name = sprintf('%s_%s_th%g_ph%g',model_prefix,m_array(i).nrn_model_ver,m_array(i).E_theta,m_array(i).E_phi);                    
    [simulation,cell_model,waveform,Efield] = loadParams(model_prefix,model_name,m_array(i)); % load parameter structures
    while ~all(succeeded)
        todo = iter(~succeeded); 
        part_threshEs = zeros(length(todo),1); % remaining results
        part_init_inds = zeros(length(todo),1);
        partsucceeded = false(size(todo)); % flags for iterations which succeeded                         
        parfor (falsei = 1:sum(~succeeded),poolobj.NumWorkers)
            i = todo(falsei); % get original index of this run
%             try
                model_namei = sprintf('%s_%s_th%g_ph%g',model_prefix,m_array(i).nrn_model_ver,...
                    m_array(i).E_theta,m_array(i).E_phi);                                          
                [simulationi,cell_modeli,waveformi,Efieldi] = ...
                    loadParams(model_prefix,model_namei,m_array(i)); % load parameter structures
                [~, ~,ap_netids,~,~,threshE] = runNrn(simulationi,cell_modeli,...
                                             waveformi,Efieldi,in.over_write);                    
                if threshE ~= 0                        
                    part_threshEs(falsei) = threshE;
                    part_init_inds(falsei) = ap_netids(1);
                    partsucceeded(falsei)=true;
                else % simulation exited because taking too long, set threshE to 0                        
                    part_threshEs(falsei) = nan;
                    part_init_inds(falsei) = nan;
                    partsucceeded(falsei)=true;
                end
%             catch
%                 fprintf('Run %g failed\n',i);                    
%             end
        end
        threshEs(~succeeded) = part_threshEs;
        init_inds(~succeeded) = part_init_inds;
        succeeded(~succeeded) = partsucceeded;
    end
    delete(poolobj); 
    if exist(pc_storage_dir,'dir')
        rmdir(pc_storage_dir,'s'); % delete parfor temporary files if necessary
    end
else % SERIAL LOOP - 1 cpu or nan (if env variable doesn't exist)
    threshEs = zeros(num_runs,1);
    init_inds = zeros(num_runs,1);
    for i = 1:num_runs
        model_namei = sprintf('%s_%s_th%g_ph%g',model_prefix,...
            in.nrn_model_ver,m_array(i).E_theta,m_array(i).E_phi);
        %model_names{i} = model_namei;
        [simulation,cell_model,waveform,Efield] = ...
            loadParams(model_prefix,model_namei,m_array(i)); % load parameter structures
        [~, ~,ap_netids,~,~,threshE] = runNrn(simulation,cell_model,...
                                        waveform,Efield,in.over_write);
        if threshE ~= 0 %simulation successfully completed and gave positive threshold
            threshEs(i) = threshE;
            init_inds(i) = ap_netids(1); % get init index
        else % simulation exited because taking too long, set threshE to 0
            threshEs(i) = nan;
            init_inds(i) = nan; % get init index
        end
    end        
end    
cell_model_name = m.cell_model_name; % save string
% Note: Efield, model uses last run's theta,phi
params = struct('simulation',simulation,'cell_model',cell_model,'waveform',waveform,...
                  'Efield',Efield); %#ok
cell_data_file = fullfile(mat_dir,'nrn_sim_data',model_prefix,cell_model_name,sprintf('%s_thresh.mat',cell_model_name));    
save(cell_data_file,'cell_model_name','thetas','phis','threshEs','init_inds','params');    
fprintf('Saved cell_data_file: %s\n',sprintf('%s_thresh.mat',cell_model_name)); 
end