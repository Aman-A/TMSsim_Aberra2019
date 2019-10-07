function fig = plotMWmorphGrid(threshEs,thetas,phis,init_inds,cell_ids,nrn_model_ver,mw_settings)
% Makes grid with morphologies next to threshold maps, used by plotFig1d
%% Generate grid and axes handles
fig = figure('units','normalized','outerposition',[0 0 1 1],'Color',[1 1 1]); 
num_rows = 5; num_cols = 10; ax_handles = cell(num_rows,num_cols);
fprintf('Plotting with %g rows and %g cols, %g entries\n',num_rows,num_cols,num_rows*num_cols);
for r = 1:num_rows
    %for c = 1:num_cols
    for c = 1:num_cols/2
        width = 7/(10*num_cols/2); height = 1/num_rows;
        width2 = 3/(10*num_cols/2); height2 = height;
        c1 = 2*c-1;
        ax_handles{r,c1} = axes('Position',[(width+width2)*(c-1),height*(num_rows-r),width,height]);
        ax_handles{r,c1+1} = axes('Position',[(width+width2)*(c-1)+width,height2*(num_rows-r),width2,height2 ]);
        % Map 1 axis settings
%         box(ax_handles{r,c1},'on');
        axis(ax_handles{r,c1},'off');
        hold(ax_handles{r,c1},'all');
        % Cell axis settings
%         box(ax_handles{r,c1+1},'on');
        axis(ax_handles{r,c1+1},'off');
        hold(ax_handles{r,c1+1},'all');
    end
end
%% Plot to grid
num_cells = length(cell_ids);
for i = 1:num_cells
   cell_id = cell_ids(i);
   [c,r] = ind2sub([num_cols,num_rows],2*i-1); % index runs across rows, then down columns
   [c2,r2] = ind2sub([num_cols,num_rows],2*i); 
   axi = ax_handles{r,c};
   mw_settings.ax = axi; 
   plotThreshMapMW(threshEs{i},thetas,phis,mw_settings);
   caxis(axi,'manual'); 
   axi2 = ax_handles{r2,c2};
   plotVeFunc(axi2,cell_id,threshEs{i},thetas,phis,init_inds{i},nrn_model_ver); 
end
end