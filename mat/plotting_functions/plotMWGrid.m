function plotMWGrid(ax_handles,threshEs,thetas,phis,num_rows,num_cols,mw_settings)
% PLOTMWGRID Makes grid of threshold maps, used by SuppFig5_6
%% Plot to grid
num_cells = length(threshEs);
for i = 1:num_cells   
   [c,r] = ind2sub([num_cols,num_rows],i); % go down coloumn, then
   axi = ax_handles{r,c};
   mw_settings.ax = axi; 
   plotThreshMapMW(threshEs{i},thetas,phis,mw_settings);
   caxis(axi,'manual');       
end
end