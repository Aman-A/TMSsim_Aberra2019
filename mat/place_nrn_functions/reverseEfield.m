function reverseEfield(layer_set_num,nrn_pop_name,nrn_model_ver,Efield_solution)
if nargin == 0
   layer_set_num = 1; 
   nrn_pop_name = 'nrn_pop1';
   nrn_model_ver = 'maxH';
   Efield_solution = 'M1_PA_MCB70';   
end
mat_dir = addPaths;
layer_set_name = sprintf('layer_set_%g',layer_set_num); 
output_dir = fullfile(mat_dir,'output_data','nrn_efields',layer_set_name);
Efield_solution_r = [Efield_solution '_r']; % denote reversed E-field direction
newEfield_folder = fullfile(output_dir,Efield_solution_r); 
Efield_folder = fullfile(output_dir,Efield_solution); 
if exist(newEfield_folder,'dir') == 0
    mkdir(newEfield_folder);
    fprintf('Making new directory for reversed E-field %s\n',newEfield_folder);
else
    fprintf('Adding to existing reversed E-field folder %s\n',newEfield_folder);
end
nrn_pop_name_full = [nrn_pop_name '_' nrn_model_ver]; % specific population E-field
newEfield_pop_fold = fullfile(newEfield_folder,nrn_pop_name_full); 
Efield_pop_fold = fullfile(Efield_folder,nrn_pop_name_full); 
if exist(newEfield_pop_fold,'dir') == 0
    fprintf('Creating output directory %s\n',newEfield_pop_fold);
    mkdir(newEfield_pop_fold)
else    
    fprintf('Saving to existing output directory %s\n',newEfield_pop_fold);
end
d = dir(Efield_pop_fold);
E_files = d(arrayfun(@(x) ~strcmp(x.name(1),'.'),d)); % remove hidden files
% extract just .mat files for interpolation in MATLAB
E_files = E_files(arrayfun(@(x) strcmp(x.name(end-2:end),'mat'),E_files));
num_cells = length(E_files);
fprintf('Reversing orientation of field for %g cells\n LS: %g, E-field: %s\n',num_cells,layer_set_num, Efield_solution);
for i = 1:num_cells
    cell_file = E_files(i).name;
    % get cell_id from file name
    [~,endI] = regexp(cell_file,'cell');
    cell_id = str2double(cell_file(endI+1:end-4)); % id is after 'cell' and up to '.mat'
    % get original E-field data for cell
    Ecell_data = load(fullfile(Efield_pop_fold,cell_file));
    Ecell_old = Ecell_data.Ecell;
    num_positions = size(Ecell_old,1);
    Ecell = cell(num_positions,1); % initialize new reversed data
    for j = 1:num_positions % loop through each cell position, flip E-field
        Ecell{j} = [Ecell_old{j}(1,:); -1*Ecell_old{j}(2:end,:)];
    end       
    save(fullfile(newEfield_pop_fold,sprintf('cell%g.mat',cell_id)),'Ecell');    
    fprintf('Saved cell%g - %g of %g\n',cell_id,i,num_cells);
end
fprintf('Finished\n');