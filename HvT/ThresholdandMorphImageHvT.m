function [BWimage BWimagefull]=ThresholdandMorphImageHvT(MP, r, k, calibration);
% [BWimage BWimagefull]=ThresholdandMorphImageHvT(MP, r, k, calibration)
% binarizes brightfield images or minimum intensity projections as the
% first step in pharygngeal grinder identification.
% The BWimage output is a binary image containing candidates for grinder
% particles.
% The BWimagefull output is the binary image resulting from direct Niblack
% segmentation prior to cleaning up the binary image using hole-filling,
% morphological and filtering operations.
% MP- Input brightfield or minimum projection image
% r- filter width (in um) for Niblack local thresholding
% k- Niblack parameter specifying the stringency of the threshold
% calibration - length scale calibration metric for input image (um/pixel)

r=floor(r/calibration); %calculate search radius in pixels
se1=strel('disk', ceil(0.25/calibration));  %structuring element for morphological opening
se2=strel('disk', ceil(0.5/calibration));

BWimagefull=niblack2d(MP, r, k, 1); %generate binary image

MinArea=round(12.5/calibration^2);    %Calculate rough particle filtration thresholds
MaxArea=round(150/calibration^2);

BWimage=imfill(BWimagefull, 'holes'); %fill holes in image

BWimage=imerode(BWimagefull,se1);     %erode image to disconnect structures

%clear small particles and dilate
BWimage=bwareaopen(BWimage,MinArea); 
BWimage=imdilate(BWimage, se2); %dilate to complete morphological open operation
BWimage=imfill(BWimage, 'holes');

%final particle filter:
BWimage=imclearborder(bwareaopen(BWimage,3*MinArea) & ~bwareaopen(BWimage, MaxArea)); 

