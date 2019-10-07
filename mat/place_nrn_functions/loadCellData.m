function cell_data = loadCellData(cell_id,nrn_model_ver,varargin)
% cellDataLoader loads cell morphology data for corresponding cell_id and
% nrn_model_ver
%
% cellDataLoader(cell_id,nrn_model_ver,varargin)
% 
% Optional inputs
% --------------
% 'coordinates' or 'C'
% 'comp_types'
% 'parent_inds'
% 'areas'
% 'branch_orders'
% 'diams'
% 'secnames'
% 'sectypes'
% 'graph'
if nargin < 2
   nrn_model_ver = 'maxH'; % default adult, myelinated human model version
end
if isempty(varargin)  % no data files specified
   % load files necessary for plotting
   data_files = {'coordinates','comp_types','parent_inds'}; 
else
    data_files = varargin; 
end
mat_dir = addPaths; 
cell_data_fold = fullfile(mat_dir,'cell_data',nrn_model_ver); 
cell_data = struct(); 
for i = 1:length(data_files)
   data_type = data_files{i};
   if strcmp(data_type,'C') 
      data_type = 'coordinates'; % replace with coordinates
   end
   data_file_path = fullfile(cell_data_fold,data_type,sprintf('%s%g.mat',data_type,cell_id));
   if exist(data_file_path,'file')
       data = load(data_file_path);
       var_name = fieldnames(data);  
       cell_data.(var_name{1}) = data.(var_name{1});
   end
end
if isempty(fieldnames(cell_data))
   error('No data was loaded, check input arguments');  
end
end
