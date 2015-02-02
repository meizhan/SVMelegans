function X_bin= niblack2d(X, filt_radius, k, darkflag);
% X_bin= niblack2d(X, filt_radius, k, darkflag) provides the binary image
% X_bin generated via local thresholding with a Niblack filter with a 
% rectangular window defined by filt_radius
% The input parameter k, specifies the stringency of the threshold
% The input parameter darkflag specifies whether the thresholding procedure
% identifies dark (darkflag==1) or bright (darkflag==0) regions of the image

filt=ones(2*filt_radius+1);

local_mean=imfilter(X,double(filt), 'symmetric')/sum(filt(:));

local_std=stdfilt(X,filt);

if darkflag==0
    X_bin=X>=(local_mean+k*local_std);
else
    X_bin=X<=(local_mean-k*local_std);
end