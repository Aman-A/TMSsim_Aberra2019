% rotate azimuthal angles (phis) of all cells by angle (°) (CCW)
% phis should be numLayers x 1 cell array of numClones x 1 cell arrays of
% numPositions x 1 vectors of cell angles
function new_phis = rotateNeuronPopPhis(phis,angle)
if nargin < 2
   mode = 'rand'; % rotate each cell randomly
else
   mode = 'uniform'; % rotate each cell uniformly with angle (°) 
end
num_layers = length(phis);
new_phis = cell(num_layers,1); 
for i= 1:num_layers    
   if strcmp(mode,'uniform')
       new_phis{i} = cellfun(@(x) x+angle,phis{i},'UniformOutput',0);
   else
       angles = rand(phis{i}{1},1)*360;
       new_phis{i} = cellfun(@(x) x+angles,phis{i},'UniformOutput',0);
   end
end

end