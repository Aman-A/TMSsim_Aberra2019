% input an axis handle, the model_prefix, cell_model_name, vector of
% thresholdEs, vector of theta and phi angles, cell array of section at
% which spike initiated, and corresponding vector of section types
% plots the cell morphology with a vector indicating direction of E-field
% for minimum threshold and a marker on the point of initiation, with the
% marker type corresponding to the section type
%
% Aman Aberra

function plotVeFunc(ax1,cellid,threshEs,thetas,phis,init_inds,nrn_model_ver)    
    expand_ax = 0;
%     x_arrow_shift = -200;
%     z_arrow_shift = 170;    
    [~,~,layer_mtype,~,clone_num] = cellNameParser(cellModelNames(cellid));
    if strcmp(layer_mtype,'L1_NGC-DA')
        xlimits = [-400 400]; zlimits = [-200 200];
        arrow_size = 316;
    elseif strcmp(layer_mtype,'L23_PC')            
        if clone_num == 2
            xlimits = [-1500 1100]; zlimits = [-800 750]; 
        else
            xlimits = [-1500 1100]; zlimits = [-1200 350]; 
        end
        arrow_size = 1000;
    elseif strcmp(layer_mtype,'L4_LBC')            
        xlimits = [-700 500]; zlimits = [-800 800]; 
        arrow_size = 620; 
    elseif strcmp(layer_mtype,'L5_TTPC2')
        xlimits = [-1000 1000]; zlimits = [-1000 1500]; 
        arrow_size = 850;
    elseif strcmp(layer_mtype,'L6_TPC')
        xlimits = [-1000 1000]; 
        if clone_num == 1
            zlimits = [-720 1500];
        else
            zlimits = [-750 1500];
        end
        arrow_size=600;
        expand_ax = 1;
    else
        xlimits = [-1000 1000]; zlimits = [-1000 1500]; 
        arrow_size = 850;
    end           
    mark_size = 10;    
    plotInitPoint(ax1,cellid,threshEs,thetas,phis,init_inds,arrow_size,...
                expand_ax,xlimits,zlimits,mark_size,nrn_model_ver)
end