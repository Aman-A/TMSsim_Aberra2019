% input nrn_model_ver and current params_struct m, or no params struct
function mout = outputMorphParams(nrn_model_ver,min)
if nargin < 2
   min = struct; % initialize to empty struct
end
mout = min;
if strcmp(nrn_model_ver,'umaxr')
    % P14 rat, unmyelinated, no scaling
    mout.nrn_model_ver = nrn_model_ver;   
    mout.myelinate_axon=0;
    mout.prune_axon=0;
    mout.replace_axon=0;
    mout.scale_axon_diam=1;
    mout.scale_apic_diam=1;
    mout.scale_basal_diam=1;
    mout.scale_soma_area=1;
    mout.scale_basal_L=1;
    mout.min_myelinD=0.2; 
    mout.min_myelinL=20; 
    mout.max_myelin_order=0;    
elseif strcmp(nrn_model_ver,'max4rP60c')
    % adult, rat with myelinated axon, using Zhu 2000 scaling
    mout.nrn_model_ver = nrn_model_ver;   
    mout.myelinate_axon=1;
    mout.prune_axon=0;
    mout.replace_axon=0;
    mout.scale_axon_diam=1.322; % using Zhu 2000 L5 pc soma scaling
    mout.scale_apic_diam=1.248; % using Romand 2011 L5 PC growth
    mout.scale_basal_diam=1.133; % using Romand 2011 L5 PC growth
    mout.scale_soma_area=1.322; % using Zhu 2000 L5 pc soma scaling
    mout.scale_basal_L=1.17; % using Romand 2011 L5 PC growth
    mout.min_myelinD=0.2; 
    mout.min_myelinL=20; 
    mout.max_myelin_order=0;
elseif strcmp(nrn_model_ver,'maxH')
    % adult, human with myelinated axon using L2/3 PC rat BB: human Eyal 2018
    mout.nrn_model_ver = nrn_model_ver;   
    mout.myelinate_axon=1;
    mout.prune_axon=0;
    mout.replace_axon=0;
    mout.scale_axon_diam=2.453; % use soma scaling
    mout.scale_apic_diam=1.876;
    mout.scale_basal_diam=1.946;
    mout.scale_soma_area=2.453;
    mout.scale_basal_L=1.17; % using Romand 2011 L5 PC growth
    mout.min_myelinD=0.2;
    mout.min_myelinL=20; 
    mout.max_myelin_order=0;
elseif strcmp(nrn_model_ver,'umaxH')
    % adult, human with unmyelinated axon using L2/3 PC rat BB: human Eyal 2018
    mout.nrn_model_ver = nrn_model_ver;   
    mout.myelinate_axon=0;
    mout.prune_axon=0;
    mout.replace_axon=0;
    mout.scale_axon_diam=2.453; 
    mout.scale_apic_diam=1.876;
    mout.scale_basal_diam=1.946;
    mout.scale_soma_area=2.453;
    mout.scale_basal_L=1.17; 
    mout.min_myelinD=0.2; 
    mout.min_myelinL=20; 
    mout.max_myelin_order=0;
elseif strcmp(nrn_model_ver,'maxHlin') 
    % adult, human with myelinated axon using L2/3 PC rat BB: human Eyal 2018
    % replace axon reconstruction with linear, artificial myelinated axon
    mout.nrn_model_ver = nrn_model_ver;   
    mout.myelinate_axon=1;
    mout.prune_axon=0;
    mout.replace_axon=5; % calls replaceAxon() and addStraightMyelinatedAxon()
    mout.scale_axon_diam=2.453; % use soma scaling
    mout.scale_apic_diam=1.876;
    mout.scale_basal_diam=1.946;
    mout.scale_soma_area=2.453;
    mout.scale_basal_L=1.17; % using Romand 2011 L5 PC growth    
    mout.total_axonL = 1000; % um - total axon length
else
    error('Input %s not recognized',nrn_model_ver); 
end 
fprintf('Using params for nrn_model_ver: %s\n',nrn_model_ver)
    