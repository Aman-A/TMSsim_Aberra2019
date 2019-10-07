% Generates table of cell model names and IDs and saves to .mat file
function writeCellModelTable()  
table_file_name = 'cell_models.mat';
mat_dir = addPaths; 
cell_data_dir = fullfile(mat_dir,'cell_data'); 
cell_dir = fullfile(mat_dir,'../nrn/cells');
d = dir(cell_dir); 
d = d(~strcmp({d.name},'.') & ~strcmp({d.name},'..') & ~strcmp({d.name},'.DS_Store')); 
cell_model_names = {d.name};
numCells = length(cell_model_names); 
layer = cell(numCells,1); 
mtype = cell(numCells,1); 
layer_mtype = cell(numCells,1); 
etype = cell(numCells,1); 
clone_num = zeros(numCells,1);
for i = 1:numCells
    [layeri,mtypei,layer_mtypei,etypei,clone_numi] = cellNameParser(cell_model_names{i}); 
    layer{i} = layeri; 
    mtype{i} = mtypei;
    layer_mtype{i} = layer_mtypei; 
    etype{i} = etypei;
    clone_num(i) = clone_numi; 
end
cell_id = (1:numCells)'; % index each cell from 0 to numCells-1
T = table(cell_id,layer,mtype,layer_mtype,etype,clone_num,'RowNames',cell_model_names); 
save(fullfile(cell_data_dir,table_file_name),'T'); 
fprintf('Saved table with %g cell models to %s\n',numCells,table_file_name); 