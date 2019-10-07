% Plot E-field magnitude on layer surface
function E_plot = plotLayersE(layersE,layer_num,mode,normE,color_map,cam_view,shift_dir,vecs_on)    
    if nargin < 8
       vecs_on = 0;  
    end
    if nargin < 7 || isempty(shift_dir)
       shift_dir = [0,-35,0]; % vector of direction of shift -14 for S116 
    end
    if nargin < 6 || isempty(cam_view)
       cam_view = [ -137.5 28]; 
    end
    if nargin < 5 
       color_map = 'jet'; 
    elseif isempty(color_map)
       color_map = 'jet';          
    end
    if nargin < 4
       normE = 0; % normalize E to max 
    end
    if nargin < 3
       mode = 'mag'; % plot magnitude - 'mag', surface normal - 'norm', or tangent, 'tan'
    end
    hold on;   
    layersE = layersE(layer_num); % extract layers of interest
    num_layers = length(layersE);  
    if strcmp(mode,'mag')
        EfieldLoc = cell(num_layers,1); 
        [EfieldLoc{:}] = layersE.EfieldLoc; % put in format of threshEs
        Emag = cellfun(@(x) x(:,4), EfieldLoc,'UniformOutput',0); % E-field magnitude on layer            
        E_plot = Emag; % E-value to plot
        title_str = 'E magnitude';         
    elseif strcmp(mode,'norm')
        Efield = cell(num_layers,1); cell_normals = cell(num_layers,1);     
        [Efield{:}] = layersE.Efield; % put in format of threshEs
        [cell_normals{:}] = layersE.cell_normals; 
        Enorm = cellfun(@(x,y) dot(x(:,4:6),y,2), Efield,cell_normals,'UniformOutput',0); 
        E_plot = Enorm; 
        title_str = 'E normal'; 
    end    
    cell_ids = mat2cell((1:num_layers)',ones(1,num_layers));
    cell_model_names = cellModelNames(1:num_layers); % place-holder cell_model_names 
    E_plot = getAllVertexCData(E_plot,layersE,cell_ids,cell_model_names);                    
    if normE
       globalmaxE = max(cellfun(@(x) max(abs(x)), E_plot)); 
       E_plot = cellfun(@(x) x/globalmaxE, E_plot,'UniformOutput',0); 
       fprintf('Plotting %s normalized to max\n',title_str);
       title_str = [title_str ' (normalized to max)'];
    else
        title_str = [title_str ' (V/m)']; 
    end
    num_plot_layers = length(layer_num); 
    for i = 1:num_plot_layers
        layer_surf = layersE(i).surface;
        layer_surf.vertices = layer_surf.vertices + (i-1)*repmat(shift_dir,size(layer_surf.vertices,1),1);
        hold on;
        p = patch(layer_surf);
        p.FaceVertexCData = E_plot{i}; % extract vector 
        p.FaceColor = 'interp'; 
        p.EdgeColor = 'none';   
        if vecs_on
            cell_origins = layersE(i).cell_origins + (i-1)*repmat(shift_dir,size(layersE(i).cell_origins,1),1); 
            Evecs = layersE(i).Efield(:,4:6);
            Evecs = Evecs./vmag(Evecs); 
           quiver3(cell_origins(:,1),cell_origins(:,2),cell_origins(:,3),Evecs(:,1),Evecs(:,2),Evecs(:,3),'k'); 
        end
    end
    title(title_str); 
    axis equal; axis tight;
    view(cam_view); 
    colormap(color_map); 
%     colorbar;
    xlabel('x'); ylabel('y'); zlabel('z');
    camlight; lighting gouraud;
end