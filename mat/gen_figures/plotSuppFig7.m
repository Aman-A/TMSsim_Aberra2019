function plotSuppFig7(save_figs)
%PLOTSUPPFIG7 Plots E-field magnitude, normal component, and thresholds
%from Fig. 3 from alternate view for Supp. Fig. 6
%  
% AUTHOR    : Aman Aberra 
if nargin==0
   save_figs = 0;
end
plot_figs = [1 2]; 
fig_names = {'SuppFig7a','SuppFi76b'};
ax_view = [90 70]; % view from other side (facing posterior)
plotFig3ab(save_figs,plot_figs,fig_names,ax_view); % E mag (a) and E norm (b)
fig_name2 = 'SuppFig7c'; 
plotFig3c(save_figs,fig_name2,ax_view); % median thresholds (c)
