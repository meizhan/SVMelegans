function [features2 featurenames2 pairID]=CalculateSecondaryFeaturesCellID(BWimageL1, MPimage, calibration);
% [features2 featurenames2 pairID]=
% CalculateSecondaryFeaturesCellID(BWimageL1, MPimage, calibration)
% calculates the second-layer relational features between potential cells
% identified in layer one of classification in the binary input image
% BWimageL1. CalculateSecondaryFeaturesCellID calculates relational
% features for one pair of bilaterally symmetric cells.
% The MPimage input is the maximally projected fluorescent image used to
% calculate the intensities of the cell pair candidates.
% The calibration input specifies the pixel to micron conversion factor 
% for the image (specified in microns/pixel).

% The features2 output is a vector with each potential cell pair in BWimageL1
% as a row and each feature described by the featurenames2 output as a column.
% The pairID output identifies the cell particles that form the pairs 
% corresponding with features2. 8-connected particles are identified in 
% pairID numbered from left to right based on the input image BWimageL1.

BwimageL1=BWimageL1==1;
BWimageL1_lab=bwlabel(BWimageL1);
featurenames2={'Delta X', 'Delta Y', 'Norm Intensity 1', 'Norm Inten 2'};
num_particles=max(BWimageL1_lab(:));

% return no features if there are not enough particles to form a pair:
if num_particles<2
    features2=[];
    pairID=[];
    return;
end

stats=regionprops(BWimageL1, 'Centroid');
centroid=vertcat(stats.Centroid);
pairID=nchoosek(1:num_particles, 2); % identify the potential cell pairs
%Calculate the relative distances between each potential pair:
delta_x=abs(centroid(pairID(:,1),1)-centroid(pairID(:,2),1))*calibration;
delta_y=abs(centroid(pairID(:,1),2)-centroid(pairID(:,2),2))*calibration;


%Calculate the mean intensities of each potential cell:
for i=1:num_particles
    I(i)=mean(MPimage(BWimageL1_lab==i));
end

%Calculate the intensities of each potential pair of cells normalized to
%the maximum intensity in the image:
max_I=max(I);
Inten=I(pairID)/max_I;

%form feature vector:
features2=[delta_x, delta_y, Inten];


