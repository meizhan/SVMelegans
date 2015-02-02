function [features2 featurenames2 setlabels]=CalculateSecondaryFeaturesCellID2Pair(BWimageL1, MPimage, calibration);
% [features2 featurenames2 setlabels]=
% CalculateSecondaryFeaturesCellID2Pair(BWimageL1, MPimage, calibration);
% calculates the second-layer relational features between potential cells
% identified in layer one of classification in the binary input image
% BWimageL1. CalculateSecondaryFeaturesCellID2Pair calculates relational
% features for two pairs of bilaterally symmetric cells (the ASI and ASJ
% pairs of neurons in C. elegans).
% The MPimage input is the maximally projected fluorescent image used to
% calculate the intensities of the cell pair candidates.
% The calibration input specifies the pixel to micron conversion factor 
% for the image (specified in microns/pixel).

% The features2 output is a vector with each potential cell tetrad in BWimageL1
% as a row and each feature described by the featurenames2 output as a column.
% The setlabels output identifies the cell particles in the sets corresponding 
% with features2. 8-connected particles are identified in setlabels numbered 
% from left to right based on the input image BWimageL1. The first two
% columns in setlabls correspond to the ASI pair and the second two columns in
% setID correspond to the ASJ pair.

BwimageL1=BWimageL1==1;
BWimageL1_lab=bwlabel(BWimageL1);
featurenames2={'Pair 1 X 1', 'Pair 1 X 2', 'Pair 1 Y1', 'Pair 1 Y2'...
    'Set 1 Norm Intensity 1', 'Set 1 Norm Intensity 2',...
    'Pair 2 X1', 'Pair 2 X2', 'Pair 2 Y1', 'Pair 2 Y2',...
    'Set 2 Norm Intensity 1', 'Set 2 Norm Intensity 2'};
num_particles=max(BWimageL1_lab(:));

if num_particles<4
    features2=[];
    setlabels=[];
    return;
end

stats=regionprops(BWimageL1, 'Centroid');
centroid=vertcat(stats.Centroid);

%identify potential tetrads of cells:
pairs=nchoosek(1:num_particles, 2);
sets=nchoosek(1:size(pairs,1),2);
pair1=pairs(sets(:,1),:); %all potential candidates for ASI pair
pair2=pairs(sets(:,2),:); %all potential candidates for ASJ pair
%identify pairings where individual cell candidates are doubly labelled:
compn=[1 1; 1 2; 2 1; 2 2];
doub=sum(pair1(:,compn(:,1))==pair2(:,compn(:,2)),2);
pair1=pair1(doub==0,:);
pair2=pair2(doub==0,:);
%Compose final pair sets:
pair1=[pair1; pair2];
pair2=[pair2; pair1(1:size(pair2,1),:)];
setlabels=[pair1 pair2];

%Calculate cell locations relative to set centroid
x=centroid(:,1);
y=centroid(:,2);
setx=mean(x(setlabels),2);
x1=(x(pair1)-[setx setx])*calibration;
x2=(x(pair2)-[setx setx])*calibration;
sety=mean(y(setlabels),2);
y1=(y(pair1)-[sety sety])*calibration;
y2=(y(pair2)-[sety sety])*calibration;

%calculate mean cell candidate intensities:
for i=1:num_particles
    I(i)=mean(MPimage(BWimageL1_lab==i)); 
end

%Calculate cell intensities normalized to the maximum intensity:
max_I=max(I);
Inten1=I(pair1)/max_I;
Inten2=I(pair2)/max_I;

%compose feature vector:
features2=[x1, y1, Inten1, x2, y2, Inten2];



