% Plots desired layer to existing figure (gcf)
function p_all = plotLayers(layers,layer_nums,shift_dir,color,normals_on)    
    if nargin < 2
       layer_nums = 1:length(layers);  
    end
    if nargin < 3
        shift_dir = [0,-35,0]; % vector of direction of shift
    end
    if nargin < 4
        if length(layer_nums) == 1
            color = [0.9 0.9 0.9];
        else
           color = jet(length(layer_nums));            
        end
    end
    if nargin < 5
        normals_on = 0;
    end
    num_plot_layers = length(layer_nums); 
    p_all = cell(num_plot_layers,1); 
    for i = 1:num_plot_layers
        layer_num = layer_nums(i); 
        cell_origins = layers(layer_num).cell_origins;
        cell_normals = layers(layer_num).cell_normals;
        layer_surf = layers(layer_num).surface;
        layer_surf.vertices = layer_surf.vertices + (i-1)*repmat(shift_dir,size(layer_surf.vertices,1),1);
        hold on;
        p = patch(layer_surf);
        p.FaceColor = color(i,:); 
        p.EdgeColor = 'none'; 
        p.FaceLighting = 'gouraud';
        if normals_on
           quiver3(cell_origins(:,1),cell_origins(:,2),cell_origins(:,3),cell_normals(:,1),cell_normals(:,2),cell_normals(:,3),'k')
        end
        p_all{i} = p;
    end
    if length(findobj(gca,'Type','light')) < 1 % check if light exists
        camlight; 
    end
    axis equal; axis tight;
    xlabel('x'); ylabel('y'); zlabel('z');
%     view([ -137.5 28]);   
    view([-77.9 48.4])
    if length(p_all) == 1
       p_all = p_all{1}; % return patch object instead of cell array 
    end
end