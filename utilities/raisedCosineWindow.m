function img = raisedCosineWindow(xgrid,ygrid,minrad,maxrad,xc,yc)
% img = raisedCosineWindow(xgrid,ygrid,minrad,maxrad)

if ~exist('xc','var')
    xc = 0;
end
if ~exist('yc','var')
    yc = 0;
end

rad = sqrt((xgrid-xc).^2+(ygrid-yc).^2);
img = zeros(size(xgrid));
img(rad<=minrad) = 1;
id = rad>minrad & rad<=maxrad;
img(id) = (cos(pi*(rad(id)-minrad)/(maxrad-minrad))+1)/2;
