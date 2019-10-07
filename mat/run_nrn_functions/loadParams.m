% Takes model_prefix and individual model run name (model_name) as well as
% variable number of parameter arguments and outputs the simulation,
% waveform, and Efield data structures used by run_nrn.m to
% specify and run a simulation. Define inputs using ('<var name>',
% <value>) pairs. Uses inputParser to assign default values to unspecified
% variables. 
% Aman Aberra
function [simulation,cell_model,waveform,Efield] = loadParams(model_prefix,model_name,varargin)        
p = inputParser;
% Simulation parameters
    p.addRequired('model_prefix',@(x) ~isnumeric(x)); % must input prefix name as string
    p.addRequired('model_name',@(x) ~isnumeric(x)); % must input individual model name as string            
    p.addParameter('dt',0.005,@isnumeric); % default dt = 0.005 ms for sampling TMS wave
    p.addParameter('tstop',1,@isnumeric); % default tstop for curr loop
    p.addParameter('v_init',-70,@isnumeric); % v_init for L3 and L5 Pyr Main cells
    p.addParameter('temp',-100,@isnumeric); % temperature, leave unset, determine default below
    p.addParameter('ss_init',1,@isnumeric); % default is steady state initialization
    p.addParameter('init_tstart',-1e11,@isnumeric); % if ss_init is 1, t start for initialization
    p.addParameter('init_dt',1e9,@isnumeric); % if ss_init is 1, time step for initialization        
    p.addParameter('vm_record_mode',1,@isnumeric); % default is to record Vm from soma
    p.addParameter('spike_record_mode',4,@isnumeric); % default is to record spikes from all comps
    p.addParameter('record_dt',0.025,@isnumeric); % default is 25 us
% cell_model parameters
    p.addParameter('cell_id',6,@isnumeric); % default is L2/3 PC 1
    p.addParameter('cell_model_name','L23_PC_cADpyr229_1',@ischar); % default cell name        
    p.addParameter('nrn_model_ver','maxH',@ischar); % default adult human myelinated (outputMorphParams overwrites defaults below)
    p.addParameter('replace_axon',0,@isnumeric); % default leave axon intact 
    p.addParameter('myelinate_axon',0,@isnumeric); % default don't myelinate axon
    p.addParameter('scale_axon_diam',1,@isnumeric); % default keep axon diameter scaling = 1
    p.addParameter('scale_mainaxon_diam',1,@isnumeric); % default keep main axon diameter scaling = 1
    p.addParameter('prune_axon',0,@isnumeric); % default don't prune axon 
    p.addParameter('scale_apic_diam',1,@isnumeric); % default don't scale apical dendrite diameters
    p.addParameter('scale_basal_diam',1,@isnumeric); % default don't scale basal dendrite diameters
    p.addParameter('scale_soma_area',1,@isnumeric); % default don't scale soma area
    p.addParameter('scale_basal_L',1,@isnumeric); % default don't scale basal dendritic length
    p.addParameter('min_myelinD',0.2,@isnumeric); % default 0.2 um min diameter to myelinate
    p.addParameter('min_myelinL',20,@isnumeric); % default 20 um min branch length to myelinate
    p.addParameter('max_myelin_order',0,@isnumeric); % default 0 below max order, myelinate any branch (above min_myelinD and     
    p.addParameter('total_axonL',1000,@isnumeric); % default 1000 um axon length when using replace_axon=5
% Waveform parameters
    p.addParameter('waveform_type','tms',@ischar); % tms only
    p.addParameter('del',1,@isnumeric); % default for thresh loop        
    p.addParameter('amp',1,@isnumeric); % default 1 V/m         
    p.addParameter('dur',0.3,@isnumeric); % default 0.3 (1/msec or kHz) (biphasic)                           
    p.addParameter('mode',1,@isnumeric); % default monophasic               
% E field parameters
    p.addParameter('E_theta',180,@isnumeric); % default polar angle 180 deg
    p.addParameter('E_phi',0,@isnumeric); % default azimuthal angle (rotation about z)
    p.addParameter('load_potentials',0,@isnumeric); % default is 0, calculate potentials for uniformE in NEURON     
    % if load_potentials is 1 (using quasipotentials from FEM E-field solution)
    p.addParameter('E_mag',1,@isnumeric); % magnitude of E at soma
    p.addParameter('E_position',1,@isnumeric); % index of cell position in layer        
    p.addParameter('E_file','none',@ischar); % name of E-field solution file  
    p.addParameter('layer_set_num',1,@isnumeric); % name of layer set, source of Efield vecs
    p.addParameter('nrn_pop_name','none',@ischar); % name of neuron population, source of Efield vecs
p.parse(model_prefix,model_name,varargin{:});    
% Create simulation, waveform, and Efield data structures
% Simulation
    % Define simulation name/type
    simulation.model_prefix = model_prefix; % general name of model, name of folder with all model runs
    simulation.model_name = model_name; % name of specific data folder of this run, e.g. 1Vpm, 2Vpm, etc., for .mat files    
    simulation.run_name = [p.Results.cell_model_name '_' simulation.model_name]; % for tmp data folder
    % Simulation parameters
    simulation.dt = p.Results.dt; % ms 
    simulation.tstop = p.Results.tstop; % ms - simulation time (for 2 pulse trains, use 4000 ms)        
    simulation.v_init = p.Results.v_init; % mV - if explicitly chosen, set v_init to input
    simulation.ss_init = p.Results.ss_init; % 1 - initialize to steady state, just use v_init
    simulation.init_tstart = p.Results.init_tstart; % 1 - ss initialization by starting at init_tsart ms (should be negative)
    simulation.init_dt = p.Results.init_dt; % 1 - take init_dt time steps from init_tstart to 0 (should be positive)
    % Recording settings
    simulation.vm_record_mode = p.Results.vm_record_mode; % 1 - record from centers, 2 record from x = 0.5 and 1
    simulation.spike_record_mode = p.Results.spike_record_mode; % 1 - record from centers, 2 record from x = 0.5 and 1
    if p.Results.record_dt >= simulation.dt && mod(p.Results.record_dt,simulation.dt) == 0
        % Check that recording time step is equal to or a multiple of
        % the simulation time step
        simulation.record_dt = p.Results.record_dt;
    else % either record_dt was unset, leaving it as 0, or incorrectly set to below time step
        simulation.record_dt = simulation.dt;
        fprintf('Recording time step set lower than simulation time step. Reset to %.3f ms\n',simulation.record_dt);
    end
    % set temperature
    if p.Results.temp == -100 % temperature unset
        if p.Results.replace_axon == 1
           simulation.temp = 34; % BB default temp 
        else 
            simulation.temp = 37; % myelinated axon default temp                
        end
        fprintf('Temp set to %g C\n',simulation.temp); 
    else
        simulation.temp = p.Results.temp; % set custom temp from input            
    end                        
    % Define cell model
    cell_model.cell_id = p.Results.cell_id;               
    cell_model.cell_model_name = p.Results.cell_model_name; % manually input cell_model_name                                
    cell_model.nrn_model_ver = p.Results.nrn_model_ver; % version name, referenced in output_morph_params.m
    cell_model.replace_axon = p.Results.replace_axon; % if 1, replace Blue-Brain axon with initial segment (axon[0] and axon[1])        
    cell_model.myelinate_axon = p.Results.myelinate_axon; % if 1, replace default Blue-brain axon with myelinated aoxon
    cell_model.scale_axon_diam = p.Results.scale_axon_diam; % if !=1, scale all axonal comp diameters by this factor
    cell_model.scale_mainaxon_diam = p.Results.scale_mainaxon_diam; % if !=1, scale all main axonal comp diameters by this factor
    cell_model.prune_axon = p.Results.prune_axon; % if >=1, prune collaterals of order max_order+1 - prune_axon
    if (cell_model.myelinate_axon == 1 && cell_model.replace_axon == 1)
       error('Cannot myelinate axon when replacing with initial segment'); % ensures conflict is not reached in NEURON
    end
    if (cell_model.prune_axon >=1 && cell_model.replace_axon == 1)
       error('Cannot prune axon when replacing with initial segment'); % ensures conflict is not reached in NEURON
    end                            
    cell_model.scale_apic_diam = p.Results.scale_apic_diam; % apical dendrite diameter scaling
    cell_model.scale_basal_diam = p.Results.scale_basal_diam; % basal dendrite diameter scaling
    cell_model.scale_soma_area = p.Results.scale_soma_area; % soma surface area scaling - via scaling diam(0.5)
    cell_model.scale_basal_L = p.Results.scale_basal_L; % basal dendritic length scaling 
    cell_model.min_myelinD = p.Results.min_myelinD; % min diameter of axon branch to myelinate
    cell_model.min_myelinL = p.Results.min_myelinL; % min length of axon branch to myelinate
    cell_model.max_myelin_order = p.Results.max_myelin_order; % if myelinate_ax, myelinates branches with order up to max_order - max_myelin_order         
    if cell_model.replace_axon == 5
       cell_model.total_axonL = p.Results.total_axonL; % um total length for artificial, linear myelinated axon (replace_axon = 5)
    end
    fprintf('Running %s cell model\n',cell_model.cell_model_name);        
% Waveform Parameters
    waveform.type = p.Results.waveform_type; 
    waveform.del = p.Results.del; % V/m   
    waveform.amp = p.Results.amp; % ms period of single pulse          
    waveform.mode = p.Results.mode; % ms - 
    TMS_waveforms = {'MagProX100 Monophasic','MagProX100 Biphasic','MagproX100 Half-sine',...        
                    sprintf('cTMS1 %g us',p.Results.dur*1000)};        
    waveform.desc = TMS_waveforms{waveform.mode}; % description of waveform wource                       
    waveform.dur =  p.Results.dur;                     
% Define E-field parameters
    Efield.theta = p.Results.E_theta; % polar angle of uniform E-field (default = 90 => in x-y plane)
    Efield.phi = p.Results.E_phi; % azimuthal angle of uniform E-field
    Efield.load_potentials = p.Results.load_potentials; % 1 - calculate potentials in MATLAB, load in NEURON, 0 - use theta,phi to calculate potentials in NEURON
    if Efield.load_potentials
       Efield.E_mag = p.Results.E_mag; 
       Efield.E_position = p.Results.E_position; 
       Efield.E_file = p.Results.E_file;   
       Efield.layer_set_num = p.Results.layer_set_num; 
       Efield.nrn_pop_name = p.Results.nrn_pop_name;        
    end            
end