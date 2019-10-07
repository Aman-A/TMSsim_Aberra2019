function NeuronPop = loadNeuronPop(nrn_pop_name,nrn_model_ver)
if nargin == 0
   nrn_pop_name = 'nrn_pop1';
   nrn_model_ver = 'maxH';
end
mat_dir = addPaths; 
layer_folder = fullfile(mat_dir,'output_data','layer_data'); 
nrn_model_ver_folder = fullfile(layer_folder,nrn_model_ver); 
NeuronPop_data = load(fullfile(nrn_model_ver_folder,...
                        sprintf('%s_%s.mat',nrn_pop_name,nrn_model_ver)));
NeuronPop = NeuronPop_data.NeuronPop; 