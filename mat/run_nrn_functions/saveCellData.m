function saveCellData(cell_id,nrn_model_ver,cell_model)
if nargin == 0
   cell_id = 1; 
   nrn_model_ver = 'maxH';
end
if nargin < 3 % no cell_model input    
    m.cell_id = cell_id; % 1 = L3, 2 = L5 pyr Mainen,3 = L3 asp ste, 4 = L4 spi ste       
    m.cell_model_name = cellModelNames(cell_id);   
    m = outputMorphParams(nrn_model_ver,m); 
    model_prefix = 'cell_dat';
    model_name = sprintf('%s_%s',model_prefix,nrn_model_ver);
    [simulation,cell_model,waveform,Efield] = loadParams(model_prefix,model_name,m);
    simulation.run_name = [m.cell_model_name '_' model_name];
end
mat_dir = addPaths; 
nrn_dir = fullfile(mat_dir,'../nrn');     
run_params_fold = fullfile(nrn_dir,'params',simulation.run_name); 
if exist(run_params_fold,'dir') == 0
    mkdir(run_params_fold); 
end
run_tmp_fold = fullfile(nrn_dir,'tmp',simulation.run_name); 
if exist(run_tmp_fold,'dir') == 0
    mkdir(run_tmp_fold); 
end
% Write params file
writeParamsHoc(run_params_fold,simulation,cell_model,Efield,waveform)                        
cd(nrn_dir); % launch init_run_stim.hoc from here
% launch neuron
nrn_file = fullfile(nrn_dir,'init_field_stim.hoc'); 
commands = sprintf(' -c "run_name=\\"%s"\\" -c "loadFiles()" -c "cell_chooser(cell_id,cell_model_name)" -c "save_cell_data()" -c "quit()"',simulation.run_name);                    
options = sprintf('-nopython -NSTACK %g -NFRAME %g ',100000,20000);

if ispc
    % Get the NEURON installation path to execute the code
    nrn_install_path = getenv('NEURONHOME');
    system([fullfile(nrn_install_path, 'bin', 'nrniv.exe') ' -nobanner ' options nrn_file commands]);
    %     system(['C:\nrn\bin\nrniv.exe -nobanner ' options nrn_file commands]);
elseif isunix
    system(['./special ' options nrn_file commands]);
end
    
cd(mat_dir); % launch init_run_stim.hoc from here
%% neuron data files
coords_file = fullfile(run_tmp_fold,sprintf('coordinates%g.txt',cell_id)); 
areas_file = fullfile(run_tmp_fold,sprintf('areas%g.txt',cell_id)); 
diams_file = fullfile(run_tmp_fold,sprintf('diams%g.txt',cell_id)); 
secnames_file = fullfile(run_tmp_fold,sprintf('secnames%g.txt',cell_id)); 
sectypes_file = fullfile(run_tmp_fold,sprintf('sectypes%g.txt',cell_id)); 
parentnames_file = fullfile(run_tmp_fold,sprintf('parentnames%g.txt',cell_id)); 
branchorders_file = fullfile(run_tmp_fold,sprintf('branchorders%g.txt',cell_id)); 
% Create matlab data folders
nrn_model_ver_fold = fullfile(mat_dir,'cell_data',nrn_model_ver);
coordinates_fold = [nrn_model_ver_fold filesep 'coordinates']; 
areas_fold = [nrn_model_ver_fold filesep 'areas']; 
diams_fold = [nrn_model_ver_fold filesep 'diams'];
secnames_fold = [nrn_model_ver_fold filesep 'secnames'];
comp_types_fold = [nrn_model_ver_fold filesep 'comp_types'];
sectypes_fold = [nrn_model_ver_fold filesep 'sectypes'];
parent_inds_fold = [nrn_model_ver_fold filesep 'parent_inds'];
branchorders_fold = [nrn_model_ver_fold filesep 'branchorders'];
graphs_fold = [nrn_model_ver_fold filesep 'graphs'];
if exist(nrn_model_ver_fold,'dir') == 0
    mkdir(nrn_model_ver_fold); % nrn_model_ver parent folder
    mkdir(coordinates_fold); % subdirectories
    mkdir(areas_fold);
    mkdir(diams_fold);
    mkdir(secnames_fold);
    mkdir(comp_types_fold);        
    mkdir(sectypes_fold);    
    mkdir(parent_inds_fold);  
    mkdir(branchorders_fold);
    mkdir(graphs_fold);
    fprintf('Created nrn_model_ver directories in cell_data/%s\n',nrn_model_ver);
end
% Save data to .mat files
%% Save coordinates (um)
C = dlmread(coords_file); 
save(fullfile(coordinates_fold,sprintf('coordinates%g.mat',cell_id)),'C'); 
%% Save areas (um2)
areas = dlmread(areas_file); 
save(fullfile(areas_fold,sprintf('areas%g.mat',cell_id)),'areas'); 
%% Save diams (um)
diams = dlmread(diams_file); 
save(fullfile(diams_fold,sprintf('diams%g.mat',cell_id)),'diams'); 
%% Save secnames
fid = fopen(secnames_file,'r'); 
secnames = textscan(fid,'%s'); 
fclose(fid);
secnames = secnames{1}; 
secnames{1} = [secnames{1} '.v(0.50)']; % add to conform to format        
save(fullfile(secnames_fold,sprintf('secnames%g.mat',cell_id)),'secnames'); 
%% Save comp_types (using secnames)
axon_comps = cellfun(@(x) ~isempty(x),strfind(secnames,'axon')); 
node_comps = 2*cellfun(@(x) ~isempty(x),strfind(secnames,'Node')); 
myelin_comps = 3*cellfun(@(x) ~isempty(x),strfind(secnames,'Myelin')); 
unmyelin_comps = 4*cellfun(@(x) ~isempty(x),strfind(secnames,'Unmyelin')); 
basal_comps = 5*cellfun(@(x) ~isempty(x),strfind(secnames,'dend')); 
apical_comps = 6*cellfun(@(x) ~isempty(x),strfind(secnames,'apic')); 
comp_types = axon_comps+node_comps+myelin_comps+unmyelin_comps+basal_comps+apical_comps;
save(fullfile(comp_types_fold,sprintf('comp_types%g.mat',cell_id)),'comp_types'); 
%% Save sectypes 
sectypes = dlmread(sectypes_file); 
save(fullfile(sectypes_fold,sprintf('sectypes%g.mat',cell_id)),'sectypes'); 
%% Save parent_inds (using parentnames & secnames)
numComp = size(C,1); 
parent_inds = zeros(numComp-1,1); 
fid = fopen(parentnames_file,'r'); 
parentnames = textscan(fid,'%s'); 
fclose(fid); 
parentnames = parentnames{1}; 
for i = 1:numComp-1
    parent_inds(i) = find(strcmp(secnames,parentnames{i}));
end
save(fullfile(parent_inds_fold,sprintf('parent_inds%g.mat',cell_id)),'parent_inds'); 
%% Save branchorders
branchorders = dlmread(branchorders_file);
save(fullfile(branchorders_fold,sprintf('branchorders%g.mat',cell_id)),'branchorders'); 
%% Delete .txt files
rmdir(run_params_fold,'s'); 
rmdir(run_tmp_fold,'s'); 
fprintf('Saved cell data to %s\n',nrn_model_ver_fold); 
end