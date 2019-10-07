% Remove points in E outside of ROI 
% E is Npoints x 6 [x y z Ex Ey Ez;...]
% ROI is [xmin xmax ymin ymax zmin zmax;...], each row is
% different box 
function Ef = getEROI(E,ROI)   
    numROIs = size(ROI,1); 
    for i = 1:numROIs
        Ef = clipPoints3d(E,ROI(i,:));
    end          
end