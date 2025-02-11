%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: 
%           Regression-based Quality with No Reference (RQNR) index. 
% 
% Interface:
%           [RQNR_index,Dl_F_index,Ds_R_index,RQNR_map,Dl_F_map,Ds_R_map] = RQNR(I_F,I_MS,I_PAN,ratio,S,sensor,alpha,beta)
%
% Inputs:
%           I_F:             Pansharpened image;
%           I_MS:            MS image resampled to panchromatic scale;
%           I_PAN:           Panchromatic image;
%           ratio:           Scale ratio between MS and PAN. Pre-condition: Integer value;
%           S:               Block size;
%           sensor:          String for type of sensor (e.g.,'QB','IKONOS,'GeoEye1','WV2','WV3' or 'none');
%           alpha            Spectral quality weight;
%           beta             Spatial quality weight.
%
% Outputs:
%           RQNR_index:      RQNR index
%           Dl_F:            D_lambda FQNR index;
%           Ds_R:            D_s RQNR index;
%           RQNR_map:        RQNR map;
%           Dl_F_map:        D_lambda FQNR map;
%           Ds_R_map:        D_s RQNR map.
% 
% References:
%           [Khan09]         M. M. Khan, L. Alparone, and J. Chanussot, "Pansharpening quality assessment using the modulation transfer functions of instruments", 
%                            IEEE Trans. Geosci. Remote Sens., vol. 11, no. 47, pp. 3880�3891, Nov. 2009.
%           [Vivone15]       G. Vivone, L. Alparone, J. Chanussot, M. Dalla Mura, A. Garzelli, G. Licciardi, R. Restaino, and L. Wald, �A Critical Comparison Among Pansharpening Algorithms�, 
%                            IEEE Transactions on Geoscience and Remote Sensing, vol. 53, no. 5, pp. 2565�2586, May 2015.
%           [Alparone18]     L. Alparone, A. Garzelli, and G. Vivone, "Spatial consistency for fullscale assessment of pansharpening",
%                            Proc. IGARSS, 2018, pp. 5132�5134.
%           [Arienzo22]      A. Arienzo, G. Vivone, A. Garzelli, L. Alparone and J. Chanussot, "Full Resolution Quality Assessment of Pansharpening: Theoretical and hands-on Approaches", 
%                            IEEE Geoscience and Remote Sensing Magazine, 10(2):2-35, 2022.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [RQNR_index, Dl_F, Ds_R, RQNR_map, Dl_F_map, Ds_R_map] = RQNR(I_F, I_MS, I_PAN, ratio, S, sensor, alpha, beta)

% Spectral distortion (the same of FQNR)
[Dl_F, Dl_F_map] = D_lambda_F(I_F, I_MS, ratio, sensor, S);

% Spatial distortion 
[Ds_R, Ds_R_map] = D_s_R(I_F, I_PAN,S);

% RQNR Quality index
RQNR_index = ((1-Dl_F)^alpha)*((1-Ds_R)^beta);

% Resize
Dl_F_map=imresize(Dl_F_map, [size(Ds_R_map, 1),size(Ds_R_map, 2)], 'cubic');

% RQNR Quality map
RQNR_map = ((1-Dl_F_map).^alpha).*((1-Ds_R_map).^beta);

end