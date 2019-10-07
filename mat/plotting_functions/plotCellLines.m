function C = plotCellLines(cell_id,nrn_model_ver,varargin)
%plotCellLines plot cell morphology using lines (edges of surface)
%
%   C = plotCellLines(cell_id,nrn_model_ver,varargin);
%
%   Usage Notes
%   ------
%   Input cell_id and nrn_model_ver to plot neuron with default coordinates
%   at origin (in µm) or to plot population of same cell, load cell data
%   externally and input with input origin, normal vector, and azimuthal
%   rotation
%
%   Inputs
%   ------
%   cell_id : integer 
%             cell_id corresponding to specific cell_model_name      
%   nrn_model_ver : string
%                   name corresponding to set of neuron model parameters,
%                   specified in outputMorphParams
%
%   Optional Inputs
%   ---------------
%   cell_data : struct with following fields (also output by
%   loadCellData.m)
%       C : num_comp x 3 array
%           coordinates of compartments for cell (default is to load from
%           cell_data/<nrn_model_ver>)
%       comp_types : num_comp x 1 vector
%                   compartment type of each compartment (default is to load):
%                   0 - soma, 1 - axon, 2 - node, 3 - myelin, 4 - unmyelin, 
%                   5 - basal, 6 - apic
%                   used for coloring if vals is empty
%       parent_inds : num_comp-1 x 1 vector
%                   index of parent compartment for each compartment,
%                   starting from 2nd compartment (1st coordinate is always
%                   root soma)
%   vals : num_comp x 1 vector
%          values for coloring neurons. (default is to use comp_types).
%          Otherwise can be used to plot potential/E-field distribution on
%          compartments, or any other set of values. 
%   cell_origin : 3 x 1 vector
%                 position to shift cell to in mm (default is no shift)
%   cell_normal : 3 x 1 vector
%                 unit vector corresponding to somatodendritic axis/element
%                 normal in layer mesh for rotating cell (default is no
%                 rotation)
%   phi : scalar
%         azimuthal rotation of cell model (default is 0)
%   lw : scalar
%        line width of lines used to plot morphology (default is 0.5)
%   Examples
%   ---------- 
%   1) **return coordinates and plot default cell position/orientation**
%   cell_id = 6; % 'L23_PC_cADpyr229_1'
%   nrn_model_ver = 'maxH'; 
%   C = plotCellLines(cell_id,nrn_model_ver); 
%   
%   2) ** pre-load cell data and plot at default cell position/orientation**
%   cell_id = 6; % 'L23_PC_cADpyr229_1'
%   nrn_model_ver = 'maxH'; 
%   cell_data = loadCellData(cell_id,nrn_model_ver); % loads C, comp_types, and parent_inds
%   figure; plotCellLines(cell_id,nrn_model_ver,'cell_data',cell_data); 
%   
%   3) ** pre-load cell_data and input values for coloring compartments
%   cell_data = loadCellData(cell_id,nrn_model_ver); % loads C, comp_types, and parent_inds
%   vals = cell_data.C(:,3); % use z coordinate
%   figure;
%   plotCellLines(cell_id,nrn_model_ver,'cell_data',cell_data,'vals',vals);
%   
%   4) ** plot specific cell within NeuronPop ** 
%   cell_id = 6; 
%   nrn_model_ver = 'maxH';
%   layer_set_num = 1; 
%   l_ind = 2; % layer index
%   c_ind = 1; % cell index within layer
%   pos_ind = 1500; % position index
%   layers = loadLayers(layer_set_num);
%   cell_data = loadCellData(cell_id,nrn_model_ver); % loads C, comp_types, and parent_inds
%   cell_data.C = cell_data.C*1e-3; % convert to mm
%   cell_origin = layers(l_ind).cell_origins(pos_ind,:);
%   cell_normal = layers(l_ind).cell_normals(pos_ind,:); 
%   nrn_pop_name = 'nrn_pop1';
%   NeuronPop = loadNeuronPop(nrn_pop_name,nrn_model_ver); 
%   phi = NeuronPop.phis{l_ind}{c_ind}(pos_ind);
%   figure;
%   plotCellLines(cell_id,nrn_model_ver,'cell_data',cell_data,'cell_origin',...
%                   cell_origin,'cell_normal',cell_normal,'phi',phi);
%   axis equal; view([0 0]); axis tight; 
if nargin == 1
   nrn_model_ver = 'maxH'; 
end
in.cell_data.C = [];  % default load
in.cell_data.comp_types = [];
in.cell_data.parent_inds = [];
in.comp_types = [];
in.parent_inds = [];
in.vals = []; % place holder, replace with comp_types later if not changed by user 
in.cell_origin = []; % no shift by default
in.cell_normal = []; % no rotation by default
in.phi = []; % no rotation by default
in.lw = 0.5;
in = sl.in.processVarargin(in,varargin);
cell_data_files = fieldnames(in.cell_data);
% Load non-input morphology data
load_data_files = cell_data_files(cellfun(@(x) isempty(in.cell_data.(x)),cell_data_files)); 
if ~isempty(load_data_files)
    in.cell_data = loadCellData(cell_id,nrn_model_ver,load_data_files{:});
end
if isempty(in.vals) % Use compartment types for coloring
   in.vals = in.cell_data.comp_types; % use compartment types 
   in.vals(in.vals==2|in.vals==3|in.vals==4) = 1; % set all axonal to 1
   in.vals(in.vals==5) = 2; % apical dend blue
   in.vals(in.vals==6) = 3; % basal dendrite green
   use_rgb_cmap = 1;
   use_sing_color = 0;
elseif length(in.vals) == 1 && isnumeric(in.vals) % single value to color morphology
    in.vals = in.vals*ones(length(in.cell_data.C),1);
    use_rgb_cmap = 0;
    use_sing_color = 0;
elseif length(in.vals) == 3 || ischar(in.vals) % color as [r,g,b] or 'color name'
    use_sing_color = 1;
    use_rgb_cmap = 0; 
else % use input color values
   use_rgb_cmap = 0; 
   use_sing_color = 0;
end
if ~isempty([in.cell_origin,in.cell_normal,in.phi]) % shift/rotate cell accordingly    
    C = placeCell(in.cell_origin,in.cell_normal,in.cell_data.C,in.phi); % converts to mm
else
    C = in.cell_data.C; % use loaded/input C
end
% colors = ['r','r','r','r','g','b'];
hold on;
% compartment sec types 
% 0 - soma, 1 - axon, 2 - node, 3 - myelin, 4 - unmyelin, 5 - basal, 6 - apic
num_comp = size(C,1);
Cp = zeros(3,(num_comp-1)*3); % 2 coordinates per compartment (not incl soma)
% make coordinate array for plotting
Cp(:,1:3:end) = C(in.cell_data.parent_inds,:)';
Cp(:,2:3:end) = C(2:end,:)';
Cp(:,3:3:end) = nan(3,num_comp-1);
% for i = 2:numComp
%    c0 = C(parent_inds(i-1),:)';
%    ci = C(i,:)';
%    Cp(:,(i-1)*3:(i-1)*3+2) = [c0,ci,[nan; nan; nan]]; 
%    valsp((i-1)*3:(i-1)*3+2) = [vals(parent_inds(i-1)) vals(i) nan]; 
% end
if use_sing_color
    surface(gca,[Cp(1,:);Cp(1,:)],[Cp(2,:);Cp(2,:)],[Cp(3,:);Cp(3,:)],...
      'FaceColor','flat','EdgeColor',in.vals,'LineWidth',in.lw);
else
    valsp =  zeros(1,(num_comp-1)*3);
    valsp(1:3:end) = in.vals(in.cell_data.parent_inds);
    valsp(2:3:end) = in.vals(2:end);
    valsp(3:3:end) = nan(1,num_comp-1);
    surface(gca,[Cp(1,:);Cp(1,:)],[Cp(2,:);Cp(2,:)],[Cp(3,:);Cp(3,:)],[valsp;valsp],...
      'FaceColor','flat','EdgeColor','flat','LineWidth',in.lw); 
end
if use_rgb_cmap % use r,g,b for coloring axon, basal, apical dendrs
   colormap(gca,[1 0 0;0 1 0;0 0 1]); 
   caxis(gca,[0 4]);  
end
%  axis off; axis equal; view([0 0]); axis tight; % preferred viewing options
end
