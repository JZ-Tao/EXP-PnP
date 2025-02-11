%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: 
%           Filter-based Quality with No Reference (FQNR): Spectral distortion index.
% 
% Interface:
%           [Dl,Dl_map] = D_lambda_F(I_F, I_MS,ratio,sensor,S)
%
% Inputs:
%           I_F:       Pansharpened image;
%           I_MS:      MS image resampled to panchromatic scale;
%           ratio:     Resolution ratio;
%           sensor:    Type of sensor;
%           S:         Block size (optional); Default value: 32.
% 
% Outputs:
%           Dl:        D_lambda index FQNR.
%           Dl_map:    D_lambda map FQNR.
% 
% Reference:
%           [Khan09]   M. M. Khan, L. Alparone, and J. Chanussot, "Pansharpening quality assessment using the modulation transfer functions of instruments,"
%                      IEEE Transactions on Geoscience and Remote Sensing, vol. 47, no. 11, pp. 3880-3891, 2009.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Dl, Dl_map] = D_lambda_F(I_F, I_MS, ratio, sensor, S)

% Input check
if (size(I_F,1) ~= size(I_MS,1) || size(I_F,2) ~= size(I_MS,2))
    error('The two images must have the same dimensions!')
end

[N,M,~] = size(I_F);
if (rem(N,S) ~= 0)
    error('Number of rows must be multiple of the block size!')
end
if (rem(M,S) ~= 0)
    error('Number of columns must be multiple of the block size!')
end

% MTF-based spatial degradation of the fused imagery
fused_degraded = MTF(I_F, sensor, ratio);

% Q2n between the expanded MS and degraded fused MS
[Q2n_index, Q2n_map] = q2n(I_MS, fused_degraded, S, S);

% Spectral distortion index
Dl = 1-Q2n_index;

% Spectral distortion map
Dl_map = 1-Q2n_map;

end