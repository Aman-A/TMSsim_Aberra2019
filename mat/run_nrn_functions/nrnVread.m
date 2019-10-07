% for reading vector vwrite data from NEURON
function [data] = nrnVread(file_name)

    data= [];
    fid = fopen(file_name,'r');
    if fid ==-1
        fprintf('file unable to be opened\n')
        return;
    else
    precision = 'double';
    [data,cnt] = fread(fid,precision);
    if length(data)>= 2
        fprintf('Successfully read %d samples from %s\n',cnt-1,file_name);
        data = data(2:end); % remove extraneous number at beginning
    else
        fprintf('Empty vector in %s\n',file_name);
    end
    
    fclose(fid);
    end
end