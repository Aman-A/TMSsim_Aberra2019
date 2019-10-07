function interpEfield(layer_set_num,nrn_pop_name,nrn_model_ver,Efield_solution,over_write)
%interpEfield Interpolates E-field at all neuron compartments of population
%
%   interpEfield(layer_set_num,nrn_pop_name,nrn_model_ver,Efield_solution)
%
%   Inputs
%   ------
%   layer_set_num : integer 
%                   specify layer set to populate with neurons                       
%   nrn_pop_name : string
%                  name of neuron population (set of neuron models in
%                  layer and their azimuthal orientations)
%   nrn_model_ver : string
%                   name corresponding to set of neuron model parameters,
%                   specified in outputMorphParams
%   Efield_solution : string
%                     name of .mat file with E-field vectors from FEM
%   over_write : 1 or 0
%               set to 1 to overwrite existing data, otherwise skips
%   Examples
%   ---------- 
%   1) 
%   layer_set_num = 1; 
%   nrn_pop_name = 'nrn_pop1';
%   nrn_model_ver = 'maxH';
%   Efield_solution = 'M1_PA_MCB70';
%   interpEfield(layer_set_num,nrn_pop_name,nrn_model_ver,Efield_solution)
%   2)
%   interpEfield(); % uses defaults (same as above)
if nargin == 0
   layer_set_num = 1;
   nrn_pop_name = 'nrn_pop1';
   nrn_model_ver = 'maxH';
   Efield_solution = 'M1_PA_MCB70';   
   over_write = 0; 
end
% Get cell coordinate file names
mat_dir = addPaths;  
output_folder = fullfile(mat_dir,'output_data');
input_folder = fullfile(mat_dir,'input_data');
layer_data_folder = fullfile('output_data','layer_data');
layer_set_name = sprintf('layer_set_%g',layer_set_num);
layers_file = load(fullfile(mat_dir,layer_data_folder,sprintf('%s.mat',layer_set_name)));
layers = layers_file.layers;
% cell coordinates placed in ROI
nrn_pop_name_full = [nrn_pop_name '_' nrn_model_ver];
pop_folder = fullfile(mat_dir,layer_data_folder,nrn_model_ver); 
NeuronPop_data = load(fullfile(pop_folder,[nrn_pop_name_full '.mat'])); 
NeuronPop = NeuronPop_data.NeuronPop; 
phis = NeuronPop.phis;
cell_ids = NeuronPop.cell_ids; 
num_layers = length(layers); 
num_cells = sum(cellfun(@length,cell_ids));
%% Load Efield
Efield_folder = fullfile(input_folder,'fem_efield_data'); 
Efield_ROI_folder = fullfile(Efield_folder,'ROI');
Efield_ROI_file = fullfile(Efield_ROI_folder,[Efield_solution '.mat']);
if exist(Efield_ROI_file,'file') == 0
    Efield_file = fullfile(Efield_folder,[Efield_solution '.mat']); 
    fprintf('Loading %s full Efield solution\n',Efield_file);
    Edata = load(Efield_file); 
    E = Edata.E;  
    if isfield(Edata,'E_wm') 
        E = [E;Edata.E_wm]; 
        fprintf('Incorporating E-field in WM\n'); 
    end
    MeshROI_data = load(fullfile(mat_dir,layer_data_folder,'MeshROI.mat'));
    MeshROI = MeshROI_data.MeshROI;    
    fprintf('Loaded MeshROI\n');     
    ROI = MeshROI.ROI;    
    % apply shift if necessary    
    % Use ROI and Noise to filter E-solution to just points of interest
    E = getEROI(E,ROI);
    if exist(Efield_ROI_folder,'dir') == 0
        mkdir(Efield_ROI_folder);
        fprintf('Made directory %s\n',Efield_ROI_folder);
    else
       fprintf('Saving to %s\n',Efield_ROI_folder);
    end
    save(Efield_ROI_file,'E');
    fprintf('Saved %s Efield in ROI\n',Efield_ROI_file)
else
    fprintf('Loading already created E-ROI file %s\n',Efield_ROI_file)
    Edata = load(Efield_ROI_file);
    E = Edata.E;
end  
Epts = E(:,1:3);
Evec = E(:,4:6);
Emag = sqrt(Evec(:,1).^2+Evec(:,2).^2+Evec(:,3).^2);
% remove outliers
Epts = Epts(Emag < 450,:); % extract points below 450 V/m
Evec = Evec(Emag < 450,:); 
% remove duplicate points
[Epts,uinds,~] = unique(Epts,'rows');
Evec = Evec(uinds,:);
%% Set up folders for saving E data
output_dir = fullfile(output_folder,'nrn_efields',layer_set_name);
if exist(output_dir,'dir') == 0
    mkdir(output_dir);
end
Efield_fold = fullfile(output_dir,Efield_solution);
if exist(Efield_fold,'dir') == 0
    mkdir(Efield_fold);
end
Efield_pop_fold = fullfile(Efield_fold,nrn_pop_name_full);
if exist(Efield_pop_fold,'dir') == 0
    fprintf('Creating output directory %s\n',Efield_pop_fold);
    mkdir(Efield_pop_fold)
else
    fprintf('Saving to existing output directory %s\n',Efield_pop_fold);    
end
%% Set up parpool
numCPUs = feature('numCores');    
fprintf('Number of CPUs available = %g\n',numCPUs);  
if numCPUs > 1 
    if numCPUs > 4
       pc_storage_dir = fullfile(mat_dir,'pc_storage',getenv('SLURM_JOB_ID'));        
       mkdir(pc_storage_dir);
       pc = parcluster('local');
       pc.JobStorageLocation =  pc_storage_dir;
    else % use default
        pc = parcluster('local');
    end
   poolobj = parpool(pc,numCPUs);
end
%% Interpolate E-field at neuron model compartment coordinates
cell_coord_fold = fullfile(mat_dir,'cell_data',nrn_model_ver,'coordinates'); 
for i = 1:num_layers
    num_cells_in_layer = length(cell_ids{i}); 
    numPositions = layers(i).num_elem; 
    % Get cell origins and normals
    cell_originsi = layers(i).cell_origins;
    cell_normalsi = layers(i).cell_normals; 
    for j = 1:num_cells_in_layer
        cell_id = cell_ids{i}(j); 
        Ecell_fileij = fullfile(Efield_pop_fold,sprintf('cell%g.mat',cell_id));
        if ~exist(Ecell_fileij,'file') || over_write 
            %  get azimuthal rotations for jth cell in layer 
            phisij = phis{i}{j}; % numPositions x 1 vector of azimuthal rotations         
            cell_coord_file = fullfile(cell_coord_fold,sprintf('coordinates%g.mat',cell_id)); 
            Cdata = load(cell_coord_file); 
            C = Cdata.C*1e-3; % coordinates of cell (local coordinates) in mm                         
            Ecell = cell(numPositions,1);
            fprintf('Looping through all positions...\n');
            tic;
            parfor k = 1:numPositions
                %             fprintf('Position %g\n',j)
                % Extract coordinates of jth position
                %             Cij = Call(1+(j-1)*numComp:j*numComp,:); % ith cell, jth position coordinates
                Cij = placeCell(cell_originsi(k,:),cell_normalsi(k,:),C,phisij(k)); % get coordinates of all compartments for kth cell (in mm)
                inds = knnsearch(Epts,Cij,'k',10); % find 10 nearest points in E
                unique_inds = unique(inds); % extract unique points
                pts_near_Cij = Epts(unique_inds,:); % extract coordinates
                E_near_i = Evec(unique_inds,:);
                % make scattered interpolants for each component of E
                Ex_int = scatteredInterpolant(pts_near_Cij(:,1),pts_near_Cij(:,2),pts_near_Cij(:,3),E_near_i(:,1),'linear');
                Ey_int = scatteredInterpolant(pts_near_Cij(:,1),pts_near_Cij(:,2),pts_near_Cij(:,3),E_near_i(:,2),'linear');
                Ez_int = scatteredInterpolant(pts_near_Cij(:,1),pts_near_Cij(:,2),pts_near_Cij(:,3),E_near_i(:,3),'linear');

                Eint = [Ex_int(Cij(:,1),Cij(:,2),Cij(:,3)),...
                    Ey_int(Cij(:,1),Cij(:,2),Cij(:,3)),...
                    Ez_int(Cij(:,1),Cij(:,2),Cij(:,3))]; % get interpolated field components
                if size(Eint) ~= size(Cij)
                    error('E-field has different number of elements from cell-coordinates');
                end
                % Add to full array of Efield vectors
                Ecell{k} = [cell_normalsi(k,:);Eint]; % place cell normal at top of array
                %             hold on;
                %             scale_factor = 0.1;
                %             plot3(pts(j,1),pts(j,2),pts(j,3),'Marker','o','MarkerSize',40,'Color','g')
                %             quiver3(pts_near_Cij(:,1),pts_near_Cij(:,2),pts_near_Cij(:,3),E_near_i(:,1)*scale_factor,E_near_i(:,2)*scale_factor,E_near_i(:,3)*scale_factor,'k','Autoscale','off')
                %             quiver3(pts(j,1),pts(j,2),pts(j,3),Eint(j,1)*scale_factor,Eint(j,2)*scale_factor,Eint(j,3)*scale_factor,'r','Autoscale','off');
                %             drawnow;
            end
            fprintf('Done\n');
            toc                              
            % Save cell array in .mat file
            save(Ecell_fileij,'Ecell');        
            fprintf('L%g of %g: Saved cell %g of %g (%g total cells)\n',i,num_layers,j,num_cells_in_layer,num_cells);
        else
           fprintf('cell%g.mat exists, over_write = %g\n',cell_id,over_write);  
        end
    end       
end
if exist('poolobj','var')
    delete(poolobj)
end
if exist(pc_storage_dir,'dir')
   rmdir(pc_storage_dir,'s'); % delete parfor temp files if necessary 
end
end