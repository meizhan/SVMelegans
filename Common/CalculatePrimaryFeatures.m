function [features, featurenames]=CalculatePrimaryFeatures(BWimage, calibration);
% [features, featurenames]=CalculatePrimaryFeatures(BWimage, calibration)
% calculates first layer features of 8-connected binary particles in the 
% input binarized image BWimage.
% The input calibration specifies the pixel to micron conversion factor 
% for the image (specified in microns/pixel).

%The output features is a vector with each particle in the image BWimage as
%a row and each feature demonted by the output featurenames as a column.

featurenames={'Area', 'Perimeter', 'Bounding Rectangle Width',...
    'Bounding Rectangle Height', 'ConvexArea', 'Maximum Feret Diameter',...
    'Maximum Feret Diameter Orientation', 'Hu Moment 1', 'Hu Moment 2',...
    'Hu Moment 3', 'Hu Moment 4', 'Hu Moment 5', 'Hu Moment 6', 'Hu Moment 7'};

BWimage=BWimage==1;

stats=regionprops(BWimage, 'Area', 'Perimeter', 'BoundingBox', 'ConvexArea', 'Extrema', 'Image');
if length(stats)==0  %Return no features if there are no particles
    features=[];
    return
end

BoundingBox=vertcat(stats.BoundingBox);
RectWidth=BoundingBox(:,4)*calibration;  %width of rectangular bounding box
RectHeight=BoundingBox(:,3)*calibration; %height of rectangular bounding box
Area=vertcat(stats.Area)*calibration^2;  %area of particle
Perim=vertcat(stats.Perimeter)*calibration; %perimeter of particle
ConvexArea=vertcat(stats.ConvexArea)*calibration^2; %convex area of particle


HuMo=[]; % initialize moments
FeretDiam=[];
FeretOrient=[];
for p=1:length(stats)   %calculate 7 Hu moments
    Image=stats(p).Image;
    HuMo(p,:)=HuInvMoments(Image); %Hu moments
    
    extrema=stats(p).Extrema;
    ex_x=extrema(:,1);
    ex_y=extrema(:,2);
    [X,Y]=meshgrid(ex_x, ex_y);
    dist=sqrt((X-X').^2+(Y-Y').^2);
    FeretDiam(p,:)=max(dist(:))*calibration;  %calculate maximum feret diameter from maximum extrema distance
    [a, b]=find(dist==max(dist(:)));
    feret_x1=ex_x(a);
    feret_y1=ex_y(a);
    feret_x2=ex_x(b);
    feret_y2=ex_y(b);
    orient=-atand((feret_x1-feret_x2)./(feret_y1-feret_y2)); %calculate maximum feret diameter orientation
    FeretOrient(p,:)=mean(orient);
end

features=[Area Perim RectHeight RectWidth ConvexArea FeretDiam FeretOrient HuMo];
