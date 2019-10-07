% input cell_name string outputs properties based on name
% Ex: cell_name = 'L1_DAC_cNAC187_1';
% layer = 'L1'
% mtype = 'DAC'
% etype = 'cNAC187"
% clone_num = 1
function [layer,mtype,layer_mtype,etype,clone_num] = cellNameParser(cell_name)
delims = regexp(cell_name,'_'); 
if length(delims) == 3
    layer = cell_name(1:delims(1)-1);
    mtype = cell_name(delims(1)+1:delims(2)-1);
    layer_mtype = [layer '_' mtype]; 
    etype = cell_name(delims(2)+1:delims(3)-1);
    clone_num = str2double(cell_name(delims(3)+1:end));
elseif length(delims) == 4
    layer = cell_name(1:delims(1)-1);
    mtype = cell_name(delims(1)+1:delims(3)-1);
    layer_mtype = [layer '_' mtype]; 
    etype = cell_name(delims(3)+1:delims(4)-1);
    clone_num = str2double(cell_name(delims(4)+1:end));
end
