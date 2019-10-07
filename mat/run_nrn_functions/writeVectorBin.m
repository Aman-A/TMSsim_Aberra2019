% input target directory, and column vectors 
% writes vector to text file in binary format with name, varname.txt
function writeVectorBin(direc,varargin)
    file_ext = '.txt';
    for i = 1:length(varargin)       
       file_name = [direc filesep inputname(i+1) file_ext];
       vector = varargin{i};
       fid = fopen(file_name,'w'); % index 2, 4, 6, etc. of arguments (file name)       
       header = [];
       header(2) = 4; % for double format
       header(1) = length(vector);
       fwrite(fid,header,'int32');
       fwrite(fid,vector,'double'); % index 1, 3, 5, etc. of arguments (vector)
       fclose(fid); 
    end
end