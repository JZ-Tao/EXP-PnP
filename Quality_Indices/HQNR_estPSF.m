%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: 
%           Hybrid Quality with No Reference (HQNR) index. 
% 
% Interface:
%           [HQNR_index,Dl_F,Ds_Q,HQNR_map,Dl_F_map,Ds_Q_map] = HQNR(I_F,I_MS,I_MS_LR,I_PAN,ratio,S,sensor,alpha,beta)
%
% Inputs:
%           I_F:             Pansharpened image;
%           I_MS:            MS image resampled to panchromatic scale;
%           I_MS_LR:         Original MS image;
%           I_PAN:           Panchromatic image;
%           ratio:           Scale ratio between MS and PAN. Pre-condition: Integer index.
%           S:               Block size;
%           sensor:          String for type of sensor (e.g.,'QB','IKONOS,'GeoEye1','WV2','WV3' or 'none');
%           alpha            Spectral quality weight;
%           beta             Spatial quality weight.
% 
% Outputs:
%           HQNR_index:      HQNR index;
%           Dl_F:            D_lambda FQNR index;
%           Ds_Q:            D_s QNR index;
%           HQNR_map:        HQNR map;
%           Dl_F_map:        D_lambda FQNR map;
%           Ds_Q_map:        D_s QNR map.
% 
% References:
%           [Alparone08]     L. Alparone, B. Aiazzi, S. Baronti, A. Garzelli, F. Nencini, and M. Selva, "Multispectral and panchromatic data fusion assessment without reference,"
%                            Photogrammetric Engineering and Remote Sensing, vol. 74, no. 2, pp. 193�200, February 2008. 
%           [Khan09]         M. M. Khan, L. Alparone, and J. Chanussot, "Pansharpening quality assessment using the modulation transfer functions of instruments", 
%                            IEEE Trans. Geosci. Remote Sens., vol. 11, no. 47, pp. 3880�3891, Nov. 2009.
%           [Aiazzi14]       B. Aiazzi, L. Alparone, S. Baronti, R. Carl�, A. Garzelli, and L. Santurri, 
%                            "Full scale assessment of pansharpening methods and data products", in SPIE Remote Sensing, pp. 924 402 � 924 402, 2014.
%           [Arienzo22]      A. Arienzo, G. Vivone, A. Garzelli, L. Alparone and J. Chanussot, "Full Resolution Quality Assessment of Pansharpening: Theoretical and hands-on Approaches", 
%                            IEEE Geoscience and Remote Sensing Magazine, 10(2):2-35, 2022.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [HQNR_index, Dl_F, Ds_Q, HQNR_map, Dl_F_map, Ds_Q_map] = HQNR_estPSF(I_F, I_MS, I_MS_LR, I_PAN, ratio, S, sensorInf, alpha, beta)

% Spectral distortion (the same of FQNR)
[Dl_F, Dl_F_map] = D_lambda_F_estPSF(I_F, I_MS, ratio, sensorInf, S);

% Spatial distortion (the same of QNR)

[Ds_Q, Ds_Q_map] = D_s_Q_est(I_F, I_MS, I_MS_LR, I_PAN, ratio, sensorInf, S);

% HQNR Quality index
HQNR_index = ((1-Dl_F)^alpha)*((1-Ds_Q)^beta);

% Resizing at the Pan scale
Dl_F_map = imresize(Dl_F_map, [size(I_PAN,1) size(I_PAN,2)], 'cubic');
Ds_Q_map = imresize(Ds_Q_map, [size(I_PAN,1) size(I_PAN,2)], 'cubic');

% HQNR Quality map
HQNR_map = ((1-Dl_F_map).^alpha).*((1-Ds_Q_map).^beta);

end