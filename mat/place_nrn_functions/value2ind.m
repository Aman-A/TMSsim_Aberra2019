function val_ind = value2ind(t0,time)
%VALUE2IND Takes an input t0 value and time vector, finds closest value in time vector to
%input t0, returns the index of that value in the time vector
%picks out first column if time input is actually a matrix, accomodates
%for input of 'raw' array from event structure
%
% AUTHOR : Aman Aberra
    if isrow(time)
       time = time'; 
    end
    time = time(:,1);
    diff = abs(time-t0);
    val_ind = find(diff == min(diff));
    if length(val_ind)> 1
        val_ind = val_ind(2); %ensures only one value is chosen, if in between, round up            
    end
end