function plotFig3d(save_fig)
%PLOTFIG3D Generates threshold vs. Efield correlation plots for Fig. 3d 
%  
% AUTHOR    : Aman Aberra 
%%
if nargin==0
   save_fig = 0;
end
mat_dir = addPaths; 
% Simulation settings
nrn_model_ver = 'maxH';
mode = 1; % monophasic MagProX100 pulse
layer_set_num = 1;
Efield_name = 'M1_PA_MCB70'; 
cell_ids = {1:5;6:10;11:15;16:20;21:25};  
nrn_pop = 'nrn_pop1-nrn_pop6_all'; 
model_prefix = sprintf('tms_%s_w%g_ls_%g_E_%s_P_%s',nrn_model_ver,mode,...
                            layer_set_num,Efield_name,nrn_pop); 
%% Plot settings
fig_size = [17 8.5]; % cm 
font_size = 8; 
text_fontsize = 9; 
font_name = 'Times';
msize = 4; % marker size
lw = 1.5; % line width
%% Load data
layersE = loadLayers(layer_set_num,2,Efield_name); 
data_fold = fullfile(mat_dir,'nrn_sim_data');
data_struct = load(fullfile(data_fold,model_prefix)); 
threshEs = data_struct.threshEs; 
cell_model_names = data_struct.cell_model_names;
num_layers = length(layersE); 
EfieldLoc = cell(num_layers,1); 
[EfieldLoc{:}] = layersE.EfieldLoc; % put in format of threshEs
Emag = cellfun(@(x) x(:,4), EfieldLoc,'UniformOutput',0); % E-field magnitude on layer        
% normal E 
Efield = cell(num_layers,1); cell_normals = cell(num_layers,1);     
[Efield{:}] = layersE.Efield; % put in format of threshEs
[cell_normals{:}] = layersE.cell_normals; 
Enorm = cellfun(@(x,y) dot(x(:,4:6),y,2), Efield,cell_normals,'UniformOutput',0); 
%% Analysis
threshEs_med = cell(num_layers,1); 
inEinds = cell(num_layers,1); 
outEinds = cell(num_layers,1); 
for i = 1:num_layers
    [~,cell_inds] = intersect(cell_model_names,cellModelNames(cell_ids{i}));
    threshEs_med{i} =  median(cell2mat(threshEs(cell_inds)),2);            
    % Normalize E-mag and norm to max
    Emag{i} = Emag{i}/max(Emag{i}); 
    Enorm{i} = Enorm{i}/max(abs(Enorm{i}));    
    % get elements where normal component is inward/outward
    inEinds{i} = Enorm{i} < 0; 
    outEinds{i} = Enorm{i} > 0; 
    % use inverse threshold
    threshEs_med{i} = 1./threshEs_med{i};  
end
%% Regression
Emag_fits = zeros(num_layers,2); 
Emag_rsq = zeros(num_layers,1); 
Enorm_fits = zeros(num_layers,2); % inE (if sep_inoutE = 1)
Enorm_rsq = zeros(num_layers,1); 
Enorm_fits_out = zeros(num_layers,2); % outE
Enorm_rsq_out = zeros(num_layers,1); 
for i = 1:num_layers
    % magnitude
    Emag_fits(i,:) = polyfit(Emag{i},threshEs_med{i},1);
    yfit = polyval(Emag_fits(i,:),Emag{i});
    Emag_rsq(i) = 1 - sum((yfit - threshEs_med{i}).^2)/((length(yfit)-1)*var(threshEs_med{i}));
    % normal
    Enorm_fits(i,:) = polyfit(Enorm{i}(inEinds{i}),threshEs_med{i}(inEinds{i}),1);
    yfit2 = polyval(Enorm_fits(i,:),Enorm{i}(inEinds{i}));
    Enorm_rsq(i) = 1 - sum((yfit2 - threshEs_med{i}(inEinds{i})).^2)/((length(yfit2)-1)*var(threshEs_med{i}(inEinds{i})));
    % outE
    Enorm_fits_out(i,:) = polyfit(Enorm{i}(outEinds{i}),threshEs_med{i}(outEinds{i}),1);
    yfit3 = polyval(Enorm_fits_out(i,:),Enorm{i}(outEinds{i}));
    Enorm_rsq_out(i) = 1 - sum((yfit3 - threshEs_med{i}(outEinds{i})).^2)/((length(yfit3)-1)*var(threshEs_med{i}(outEinds{i})));
end
%% Plot
fig = figure('Color','w'); 
fig.Units = 'centimeters';
fig.Position(1:2) = [2 13.41]; 
fig.Position(3:4) = fig_size; 
xstart = 0.06; 
axwidth = 0.7/num_layers; 
axxgap = (0.95-axwidth*num_layers)/5;
ystart = 0.15; 
axheight = .15; 
axygap = 0.15;
textx = 0.55; 
texty = 0.125;

for i = 1:num_layers  
   % *** MAGNITUDE ***
   axi = axes('Position',[xstart + (axwidth+axxgap)*(i-1),ystart+axheight*2+axygap*2,axwidth,axheight]);    
   plot(axi,Emag{i},threshEs_med{i}/max(threshEs_med{i}),'k.','MarkerSize',msize); % normalize to max
   hold on;
   plot(axi,[min(Emag{i}),max(Emag{i})],polyval(Emag_fits(i,:),[min(Emag{i}),max(Emag{i})])/max(threshEs_med{i}),...
       'g','LineWidth',lw);
   ylim(axi,[0 1]);
   box(axi,'off'); 
   grid(axi,'off');
   axi.FontName = font_name;
   axi.FontSize = font_size;
   axlims = axis(axi);
   text(axlims(1)+diff(axlims(1:2))*(textx),axlims(3)+diff(axlims(3:4))*texty,num2str(Emag_rsq(i),'$$R^{2} = %.3f$$'),...
       'FontSize',text_fontsize,'FontName',font_name,'FontWeight','Bold','Interpreter','latex');
   xlabel(axi,'$$|\vec{E}|$$ (norm.)','Interpreter','latex');
   xlim(axi,[0 1]);
   axi.YColor = 'k'; axi.XColor = 'k';
   % *** NORMAL ***
   % in
   axi2 = axes('Position',[xstart + (axwidth+axxgap)*(i-1),ystart+axheight+axygap,axwidth,axheight]);       
   % points
   plot(axi2,Enorm{i}(inEinds{i}),threshEs_med{i}(inEinds{i})/max(threshEs_med{i}(inEinds{i})),'k.','MarkerSize',msize);
   hold on;
   % line
   plot(axi2,[min(Enorm{i}(inEinds{i})) max(Enorm{i}(inEinds{i}))],...
       polyval(Enorm_fits(i,:),[min(Enorm{i}(inEinds{i})) max(Enorm{i}(inEinds{i}))])/max(threshEs_med{i}(inEinds{i})),...
       'b','LineWidth',lw); 
   ylim(axi2,[0 1]);     
   box(axi2,'off'); 
   grid(axi2,'off'); 
   axi2.FontName = font_name;
   axi2.FontSize = font_size;  
   axi2.YAxisLocation = 'right';
   axlims2 = axis(axi2);   
   % in text
   text(axlims2(1)+diff(axlims2(1:2))*(-0.05),axlims2(3)+diff(axlims2(3:4))*(texty),num2str(Enorm_rsq(i),'$$R^{2} = %.3f$$'),...
   'FontSize',text_fontsize,'FontName',font_name,'FontWeight','Bold','Interpreter','latex','Color','k');    
   axi2.YColor = 'k'; axi2.XColor = 'k';
   xlabel(axi2,'$$\vec{E} \cdot \hat{n}$$ (norm.)','Interpreter','latex');
   xlim(axi2,[-1 0]);      
   % out
   axi3 = axes('Position',[xstart + (axwidth+axxgap)*(i-1),ystart,axwidth,axheight]);    
   % points
   plot(axi3,Enorm{i}(outEinds{i}),threshEs_med{i}(outEinds{i})/max(threshEs_med{i}(outEinds{i})),...
       'k.','MarkerSize',msize);
   hold on;
   % line
   plot(axi3,[min(Enorm{i}(outEinds{i})),max(Enorm{i}(outEinds{i}))],...
       polyval(Enorm_fits_out(i,:),[min(Enorm{i}(outEinds{i})),max(Enorm{i}(outEinds{i}))])/max(threshEs_med{i}(outEinds{i})),'r','LineWidth',lw);
   ylim(axi3,[0 1]);    
   box(axi3,'off'); 
   grid(axi3,'off'); 
   axi3.FontName = font_name;
   axi3.FontSize = font_size;      
   axlims2 = axis(axi3);
   % out text
   text(axlims2(1)+diff(axlims2(1:2))*textx,axlims2(3)+diff(axlims2(3:4))*(texty),num2str(Enorm_rsq_out(i),'$$R^{2} = %.3f$$'),...
           'FontSize',text_fontsize,'FontName',font_name,'FontWeight','Bold','Interpreter','latex','Color','k');   
   axi3.YColor = 'k'; axi3.XColor = 'k';
   xlabel(axi3,'$$\vec{E} \cdot \hat{n}$$ (norm.)','Interpreter','latex');
   xlim(axi3,[0 1]);      
end
%%
if save_fig   
   fig_name = 'Fig3d';
   fig_fold = fullfile(mat_dir,'figures');
   savefig(fig,fullfile(fig_fold,[fig_name '.fig']));   
   print(fig,fullfile(fig_fold,[fig_name '.eps']),'-depsc2','-cmyk');   
end
