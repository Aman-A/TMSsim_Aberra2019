% open figure
% adds fill for all boxes for plotFig5_Sommerdata
function addBoxPlotFillAll(colors,num_subdiv,num_total)
% ADDBOXPLOTFILLALL adds filled in patch elements for box plots of gcf
if nargin < 2
   num_subdiv = 5;
   num_total = 30; 
end
ax = gca;
b = ax.Children(end); 
bc = b.Children; 
% box objects are 61:90
% colors = flipud(jet(5)); colors(2,:) = [0.5961 0.3059 0.6392]; 
% colors = flipud(colors);
boxes = bc((num_total+1):num_total*2); 
xData= reshape([boxes(:).XData],5,num_total)';
yData= reshape([boxes(:).YData],5,num_total)';
colors = flipud(colors); % flip to preserve original order
colors_all = repmat(colors,num_subdiv+1,1); 
% median lines
medlines = bc(1:num_total); 
% medlines = bc(31:60); 
[medlines.Color] = deal('w'); 
ps = cell(num_subdiv+1,1); 
for i = 1:num_subdiv
    ps{i} = patch(xData(i:num_subdiv:end,:)',yData(i:num_subdiv:end,:)','k'); 
    ps{i}.FaceVertexCData = repmat(colors_all(i,:),5*(num_total/num_subdiv),1); 
    ps{i}.FaceColor = colors_all(i,:);
    ps{i}.EdgeColor = 'flat';
end
ax.Children = [ax.Children(end);ax.Children(1:end-1)]; 
