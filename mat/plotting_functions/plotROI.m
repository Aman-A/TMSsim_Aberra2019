% Plot rectangle with patches of ROI = [xmin xmax ymin ymax zmin zmax]
% plots on existing axis (gca)
function plotROI(ROI)
    size = abs([ROI(2)-ROI(1),ROI(4)-ROI(3),ROI(6)-ROI(5)]);
    origin = [ROI(1), ROI(3), ROI(5)];
    x=([0 1 1 0 0 0;1 1 0 0 1 1;1 1 0 0 1 1;0 1 1 0 0 0]-0.5)*size(1)+origin(1)+size(1)/2;
    y=([0 0 1 1 0 0;0 1 1 0 0 0;0 1 1 0 1 1;0 0 1 1 1 1]-0.5)*size(2)+origin(2)+size(2)/2;
    z=([0 0 0 0 0 1;0 0 0 0 0 1;1 1 1 1 0 1;1 1 1 1 0 1]-0.5)*size(3)+origin(3)+size(3)/2;
    hold on;
    for i=1:6
        h=patch(x(:,i),y(:,i),z(:,i),'w');
        set(h,'edgecolor','k','LineWidth',4,'FaceAlpha',0)
    end
%     plot3(repmat(ROI(1:2),4,1),repmat(ROI(2:3)',2,2),[repmat(ROI(5),2,2);repmat(ROI(6),2,2)],'k','LineWidth',4);
%     plot3(repmat(ROI(1:2)',2,2),repmat(ROI(2:3),4,1),[repmat(ROI(5),2,2);repmat(ROI(6),2,2)],'k','LineWidth',4);
%     plot3([repmat(ROI(1),2,2);repmat(ROI(2),2,2)],repmat(ROI(2:3)',2,2),repmat(ROI(5:6),4,1),'k','LineWidth',4);
end