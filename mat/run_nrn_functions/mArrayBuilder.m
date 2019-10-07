% generates num_run x 1 struct array for inputting into parfor loop
% generates struct array with all static fields initialized to the same
% value in each element of array, while varying fields are set to their
% correct value
% input single 1x1 struct with static fields already initialized, and
% num_run x 1 vectors with each element corresponding to the ith run's
% parameter setting, and the name of the vector being identical to the name
% of the field
% all vectors must be of size num_runs x 1, with each index i corresponding
% each field must contain variables of type double
% to parameters in ith run
function m_array = mArrayBuilder(m,varargin)    
    num_param_vects = length(varargin);
    if any(cellfun(@isvector,varargin)==0)
        error('Variable parameters must be entered as vectors, i.e. Nx1 or 1xN arrays')
    elseif num_param_vects > 1 && any(diff(cellfun('length',varargin)))
        error('Unequal number of elements in parameter arrays')    
    end 
    static_param_fields = fieldnames(m);
    new_param_fields = cell(num_param_vects,1);
    % grab field names from inputs, should correspond to fields in load_params/load_params_thresh
    for i = 1:num_param_vects
        new_param_fields{i} = inputname(i+1);
    end    
    % initialize eval_string with first field entry (also allows for
    % same string to be added with correct formatting)
    if isa(m.(static_param_fields{1}),'double')
        eval_string = sprintf('m_array = struct(''%s'',%d',static_param_fields{1},m.(static_param_fields{1}));
    elseif ischar(m.(static_param_fields{1}))
        eval_string = sprintf('m_array = struct(''%s'',''%s''',static_param_fields{1},m.(static_param_fields{1}));
    end
    for i = 2:length(static_param_fields)
        if isa(m.(static_param_fields{i}),'double')
            eval_string = [eval_string sprintf(',''%s'',%d',static_param_fields{i},m.(static_param_fields{i}))]; 
        elseif ischar(m.(static_param_fields{i}))
            eval_string = [eval_string sprintf(',''%s'',''%s''',static_param_fields{i},m.(static_param_fields{i}))];
        else
           error('field is not of type double'); % must be double to work with chosen formatspec identifier
        end
    end    
    for n = 1:length(new_param_fields)        
        if isrow(varargin{n})
           varargin{n} = varargin{n}'; % make sure it's a column vector for mat2cell to work
        end
        eval(sprintf('%s=varargin{%g};',new_param_fields{n},n)); % create variable with correct field name in workspace
        eval_string = [eval_string sprintf(',''%s'',mat2cell(%s,ones(1,length(%s)))',new_param_fields{n},new_param_fields{n},new_param_fields{n})];
    end
    eval_string = [eval_string ');']; % close struct initialization string
    eval(eval_string) % generates m_array for output
end