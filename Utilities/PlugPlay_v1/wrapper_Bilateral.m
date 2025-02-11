function out = wrapper_Bilateral(in,sigma)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% out = wrapper_Bilateral(in,sigma)
% performs bilateral filter
% 
% Require RF package
% 
% Available in ./RF/
%
% Xiran Wang and Stanley Chan
% Copyright 2016
% Purdue University, West Lafayette, In, USA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sigmaA(1) = 3;
sigmaA(2) = sigma;
ref = bfilter2(in,3,sigmaA,in,sigma);
out = bfilter2(in,3,sigmaA,ref,sigma);

end