function [BW3 BW0 BW1 BW2 CellID]=CellIDOnePair(MPimage, calibration)
% [BW3 BW0 BW1 BW2 CellID]=CellIDOnePair(MPimage, calibration) identifies
% bilaterally symmetric pairs of cells in the fluorescent maximum
% projection input image MPimage with the pixel to micron conversion factor
% specified by calibration (specified in microns/pixel).
%
% The input image (MPimage) is assumed to be oriented with the
% anterior-posterior axis aligned vertically.
%
% The first output, BW3, outputs a binary image with the identified pair
% (if any) shown as two binary regions. The other outputs show intermediate
% outputs:
% BW0 shows the binary image calculated from MPimage using Niblack
% semgentation
% BW1 shows the binary image, BW0, after preliminary filtration of spurious
% binary regions too small or large to be cells
% BW2 shows the binary image containing candidates for cell pairs after the
% first layer of SVM classification
% CellID is a vector identifying the regions identified as cells in a cell
% pair as referenced to BW1 (numbered from left to right).
% 
% Dependencies: 
% Use addpath(genpath(CODEFOLDER)), to add all dependencies to the path.
% ThresholdandMorphImageCellID (in CellID subfolder)
% CalculatePrimaryFeatures (in Common subfolder)
% svmpredict (in Common subfolder, from the LIBSVM 3.17 package)
% CellIDSVMLayer1.mat (in CellID subfolder)
% CellIDSVMLayer2.mat (in CellID subfolder)

load CellIDSVMLayer1.mat
load CellIDSVMLayer2.mat

%Threshold input image:
[BW1 BW0]=ThresholdandMorphImageCellID(double(MPimage), 5, 0.85, calibration);


%Calculate and scale features for the first layer of classification:
BW1_lab=bwlabel(BW1);
[featuresL1 featurenamesL1]=CalculatePrimaryFeatures(BW1, calibration);
for i=1:(size(featuresL1,2)) %scale features
    featuresL1(:,i)=(featuresL1(:,i)+shiftL1(i))*scaleL1(i);
end
%Perform the first layer of classification:
resultL1=svmpredict(ones([size(featuresL1,1),1]), featuresL1, svmStructL1);
resultL1ind=find(resultL1==1);
BW2=ismember(BW1_lab,resultL1ind);

%Calculate and scale features for the second layer of classification:
[featuresL2 featurenamesL2 pairs]=CalculateSecondaryFeaturesCellID(BW2, MPimage, calibration);
pairs=resultL1ind(pairs);
if size(pairs,2)~=2
    pairs=reshape(pairs, [size(pairs,2) 2]);
end
for i=1:(size(featuresL2,2)) %scale features
    featuresL2(:,i)=(featuresL2(:,i)+shiftL2(i))*scaleL2(i);
end
%Perform the second layer of classification:
[resultL2, junk, probL2]=svmpredict(ones([size(featuresL2,1),1]), featuresL2, svmStructL2, ['-b 1']);
%Restrict results to the most likely cell pair:
probL2=probL2(:,1);
if sum(resultL2==1)>1
    resultL2prob=-1*ones(size(resultL2));
    T=find(probL2==max(probL2));
    resultL2prob(T)=1;
else
    resultL2prob=resultL2;
end

%Identify cell pair and generate binary cell pair image:
CellID=pairs(resultL2prob==1, :);
if length(CellID)>1
    BW3=ismember(BW1_lab, CellID);
else
    BW3=zeros(size(BW1));
end



