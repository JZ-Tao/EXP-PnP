%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: 
%           Filter-based Quality with No Reference (FQNR) index. 
% 
% Interface:
%           [FQNR_index,Dl_F,Ds_F,FQNR_map,Dl_F_map,Ds_F_map] = FQNR(I_F,I_MS,I_MS_LR,I_PAN,ratio,S,sensor,alpha,beta)
%
% Inputs:
%           I_F:             Pansharpened image;
%           I_MS:            MS image resampled to panchromatic scale;
%           I_MS_LR:         Original MS image;
%           I_PAN:           Panchromatic image;
%           ratio:           Scale ratio between MS and PAN. Pre-condition: Integer value;
%           S:               Block size;
%           sensor:          String for type of sensor (e.g.,'QB','IKONOS,'GeoEye1','WV2','WV3' or 'none');
%           alpha            Spectral quality weight 
%           beta             Spatial quality weight 
% 
% Outputs:
%           FQNR_index:      FQNR index;
%           Dl_F:            D_lambda FQNR index;
%           Ds_F:            D_s FQNR index;
%           FQNR_map:        FQNR map;
%           Dl_F_map:        D_lambda FQNR map;
%           Ds_F_map:        D_s FQNR map.
% 
% References:
%           [Khan09]         M. M. Khan, L. Alparone, and J. Chanussot, "Pansharpening quality assessment using the modulation transfer functions of instruments,"
%                            IEEE Transactions on Geoscience and Remote Sensing, vol. 47, no. 11, pp. 3880-3891, 2009.
%           [Vivone15]       G. Vivone, L. Alparone, J. Chanussot, M. Dalla Mura, A. Garzelli, G. Licciardi, R. Restaino, and L. Wald, �A Critical Comparison Among Pansharpening Algorithms�, 
%                            IEEE Transactions on Geoscience and Remote Sensing, vol. 53, no. 5, pp. 2565�2586, May 2015.
%           [Arienzo22]      A. Arienzo, G. Vivone, A. Garzelli, L. Alparone and J. Chanussot, "Full Resolution Quality Assessment of Pansharpening: Theoretical and hands-on Approaches", 
%                            IEEE Geoscience and Remote Sensing Magazine, 10(2):2-35, 2022.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [FQNR_index, Dl_F_index, Ds_F_index, FQNR_map, Dl_F_map, Ds_F_map] = FQNR(I_F, I_MS, I_MS_LR, I_PAN, ratio, S, sensor, alpha, beta)

% Spectral distortion
[Dl_F_index, Dl_F_map] = D_lambda_F2(I_F, I_MS, ratio, sensor, S);

% Spatial distortion
[Ds_F_index, Ds_F_map] = D_s_F2(I_F, I_PAN, I_MS_LR, ratio, sensor, S); 

% FQNR Quality index
FQNR_index = ((1-Dl_F_index)^alpha)*((1-Ds_F_index)^beta);

% Resize
Dl_F_map=imresize(Dl_F_map, [size(Ds_F_map,1),size(Ds_F_map,2)], 'cubic');

% FQNR Quality map
FQNR_map = ((1-Dl_F_map).^alpha).*((1-Ds_F_map).^beta);

end