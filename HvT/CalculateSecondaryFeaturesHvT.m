function [features2 featurenames2]=CalculateSecondaryFeaturesHvT(BWimageL1, BWimagefull, calibration);
% [features2 featurenames2]=
% CalculateSecondaryFeaturesHvT(BWimageL1, BWimagefull, calibration);
% calculates the second-layer regional features of the putative grinder
% particles identified in layer one of classification in the binary input image
% BWimageL1. 
% The BWimagefull input is the binarized image resulting directly from
% initial Niblack segmentation and is used to calculate the regional
% characteristics around the putatative grinder particle.
% The calibration input specifies the pixel to micron conversion factor 
% for the image (specified in microns/pixel).

% The features2 output is a vector with each putative grinder particle in BWimageL1
% as a row and each feature described by the featurenames2 output as a column.

BWimageL1=BWimageL1~=0; %Convert into binary image
stats=regionprops(BWimageL1, 'centroid');
centroid=vertcat(stats.Centroid);
BWimageL1_label=bwlabel(BWimageL1); 

[M,N]=size(BWimagefull);
[x,y]=meshgrid(1:N, 1:M); %index for calculating regional masks

features2=[];

%construct bins for radial and angular regional properties:
radius=25;
radiuspix=round(radius/calibration);
Rbin=[0 7.5 12.5 17.5 20 radius]; %bins for radius
numRbin=length(Rbin)-1;
Angbin=linspace(-180, 180, 5); %bins for angles
numAbin=length(Angbin)-1;

%calcualte regional properties for each candidate grinder particle:
for i=1:length(stats)
    BWpart=BWimageL1_label==i;
    %regional mask:
    mask=((x-centroid(i,1)).^2+(y-centroid(i,2)).^2)<radiuspix^2; %mask for inner circle
    %masked regional binary image:
    BWlocal=BWimagefull&mask&~BWpart;
    
    %calculate properties of the particles in the regional image:
    [pixlisty pixlistx]=find(BWlocal);
    statind=regionprops(BWlocal, 'Area', 'Centroid', 'Perimeter');
    centroidind=vertcat(statind.Centroid);
    areasind=vertcat(statind.Area);
    perimind=vertcat(statind.Perimeter);
    circind=perimind./(2*sqrt(pi*areasind));    
    
    if length(pixlistx)>0
        %locate regional pixels relative to the region center:
        pixdiff=[pixlistx-centroid(i,1) pixlisty-centroid(i,2)];
        angles=atan2d(pixdiff(:,2),pixdiff(:,1));
        r=sqrt(sum(pixdiff.^2,2));
        %locate regional particle centroids relative to the region center:
        pixinddiff=[centroidind(:,1)-centroid(i,1) centroidind(:,2)-centroid(i,2)];
        rc=sqrt(sum(pixinddiff.^2,2));
        anglesc=atan2d(pixinddiff(:,2),pixinddiff(:,1));
        
        %Calculate radial distributions of pixels in the region:
        [AreaDistRad]=histc(r*calibration, Rbin);
        AreaDistRad=reshape(AreaDistRad(1:end-1), [1 numRbin]);
        AreaDistRad=AreaDistRad*calibration^2./(pi*diff(Rbin).^2);
        %Calculate radial distributions of numbers of particles in the
        %region:
        [NumDistRad id]=histc(rc*calibration, Rbin);
        NumDistRad=reshape(NumDistRad(1:end-1), [1 numRbin])./(pi*diff(Rbin).^2);
        %Calculate radial distributions of mean regional particle size,
        %and circularity:
        for j =1:numRbin
            a=areasind(id==j);
            c=circind(id==j);
            if ~isnan(a)
                MeanSizeDistRad(j)=mean(a)*calibration^2;
                MeanCircDistRad(j)=mean(c);
            else
                MeanSizeDistRad(j)=0;
                MeanCircDistRad(j)=0;
            end
        end

        %Calculate angular distributions of pixels in the region:
        AreaDistAng=histc(angles,Angbin);
        AreaDistAng=reshape(AreaDistAng(1:end-1)/(pi*radiuspix.^2/4), [1 numAbin]);
        %Calculate angular distributions of numbers of particles in the
        %region:
        [NumDistAng id]=histc(anglesc,Angbin);
        NumDistAng=reshape(NumDistAng(1:end-1)/(pi*radiuspix.^2/4), [1 numAbin]);
        %Calculate angular distributions of mean regional particle size:
        for j=1:numAbin
            a=areasind(id==j);
            if ~isnan(a)
                MeanSizeDistAng(j)=mean(a)*calibration^2;
            else
                MeanSizeDistAng(j)=0;
            end
        end
        %calculate normalized moment of inertia for the region:
        NormMI=sum(r.^2/(length(r)*radiuspix^2));
    else
        %return zeros if no particles are found in the local region:
        AreaDistRad=zeros([1 numRbin]);
        NumDistRad=AreaDistRad;
        MeanSizeDistRad=AreaDistRad;
        MeanPerimDistRad=AreaDistRad;
        MeanCircDistRad=AreaDistRad;
        AreaDistAng=zeros([1 length(Angcent)]);
        NumDistAng=AreaDistAng;
        MeanSizeDistAng=AreaDistAng;
        NormMI=0;
    end
    %construct feature vector:
    features2(i,:)=[AreaDistRad NumDistRad MeanSizeDistRad MeanCircDistRad AreaDistAng NumDistAng MeanSizeDistAng NormMI];
end

%name features:
for i=1:numRbin
    featurenames2{i}=['Area Radial Distro ' num2str(i)];
    featurenames2{numRbin+i}=['Number Radial Distro ' num2str(i)];
    featurenames2{2*numRbin+i}=['Mean Size Radial Distro ' num2str(i)];
    featurenames2{3*numRbin+i}=['Mean Circularity Radial Distro ', num2str(i)];
end

for i=1:numAbin
    anglefeatnames{i}=['Area Angular Distro ', num2str(i)];
    anglefeatnames{numAbin+i}=['Number Angular Distro ' num2str(i)];
    anglefeatnames{2*numAbin+i}=['Mean Size Angular Distro ' num2str(i)];
end
featurenames2=[featurenames2 anglefeatnames {'Normalized Moment of Inertia'}];




