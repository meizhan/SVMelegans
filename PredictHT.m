function [HTresult BW0 BW1 BW2 BW3]=PredictHT(MPimage, calibration)
% [HTresult BW0 BW1 BW2 BW3]=PredictHT(MPimage, calibration)
% determines whether an input brightfield image, MPimage, shows the head of
% a C. elegans specimen via detection of the pharyngeal grinder.
% The calibration input species the pixel to micron conversion factor for
% the input image (in um/pixel)
%
% The input image (MPimage) is assumed to be oriented with the
% anterior-posterior axis aligned vertically.
%
% The HTresult output is 1 if the pharyngeal grinder is detected and 0 if
% it is not. The other outputs show intermediate outputs of the detection
% procedure:
% BW0 shows the binary image calculated from MPimage using Niblack
% semgentation
% BW1 shows the binary image, BW0, after preliminary filtration of spurious
% binary regions too small or large to be the grinder particle
% BW2 shows the binary image containing candidates for grinder after the
% first layer of SVM classification
% BW3 shows the binary image showing grinder particle(s) after the second
% layer of SVM classification.
% 
% Dependencies: 
% Use addpath(genpath(CODEFOLDER)), to add all dependencies to the path.
% ThresholdandMorphImageHvT (in HvT subfolder)
% CalculatePrimaryFeatures (in Common subfolder)
% svmpredict (in Common subfolder, from the LIBSVM 3.17 package)
% HvTSVMLayer1.mat (in HvT subfolder)
% HvTSVMLayer2.mat (in HvT subfolder)

load HvTSVMLayer1.mat;
load HvTSVMLayer2.mat;

%Threshold input image to generate candidates for the grinder particle:
[BW1 BW0]=ThresholdandMorphImageHvT(double(MPimage), 15, 0.75, calibration);
BW1_lab=bwlabel(BW1);
%Calculate features for the first layer of classification:
[featuresL1 featurenamesL1]=CalculatePrimaryFeatures(BW1, calibration);
for i=1:(size(featuresL1,2)) %scale features
    featuresL1(:,i)=(featuresL1(:,i)+shiftL1(i))*scaleL1(i);
end
%Perform the first layer of classification:
resultL1=svmpredict(ones([size(featuresL1,1),1]), featuresL1, svmStructL1);
BW2=ismember(BW1_lab,find(resultL1==1));
BW2_lab=bwlabel(BW2);
%perform the second layer of classification:
[featuresL2 featurenamesL2]=CalculateSecondaryFeaturesHvT(BW2, BW0, calibration);
for i=1:(size(featuresL2,2)) %scale features
    featuresL2(:,i)=(featuresL2(:,i)+shiftL2(i))*scaleL2(i);
end
resultL2=svmpredict(ones([size(featuresL2,1),1]), featuresL2, svmStructL2);
BW3=ismember(BW2_lab, find(resultL2==1));
%generate binary head versus tail determination based on the presence of a
%pharyngeal particle:
HTresult=sum(resultL2)~=0;

