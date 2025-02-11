function out = wrapper_contour(y,sigma)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% wrapper file for contourlet
% Please download contourlet package at
% http://minhdo.ece.illinois.edu/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ref = contour_denoise(y,sigma);
out = contour_denoise(y,sigma,ref);