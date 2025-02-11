%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: 
%           Filter-based Quality with No Reference (FQNR): Spatial distortion index.
% 
% Interface:
%           [Ds_F_index,Ds_F_map] = D_s_F(I_F,I_PAN,I_MS_LR,ratio,sensor,S)
%
% Inputs:
%           I_F:         Pansharpened image;
%           I_PAN:       Panchromatic image;
%           I_MS_LR:     Original MS image;
%           ratio:       Scale ratio between MS and PAN. Pre-condition: Integer value;
%           sensor:      String for type of sensor (e.g.,'QB','IKONOS,'GeoEye1','WV2','WV3' or 'none');
%           S:           Block size (optional); Default value: 32;
% 
% Outputs:
%           Ds_F_index:  D_s FQNR index;
%           Ds_F_map:    D_s FQNR map;
% 
% References:
%           [Khan09]     M. M. Khan, L. Alparone, and J. Chanussot, "Pansharpening quality assessment using the modulation transfer functions of instruments,"
%                        IEEE Transactions on Geoscience and Remote Sensing, vol. 47, no. 11, pp. 3880-3891, 2009.
%           [Vivone15]   G. Vivone, L. Alparone, J. Chanussot, M. Dalla Mura, A. Garzelli, G. Licciardi, R. Restaino, and L. Wald, 揂 Critical Comparison Among Pansharpening Algorithms�, 
%                        IEEE Transactions on Geoscience and Remote Sensing, vol. 53, no. 5, pp. 2565�2586, May 2015.
%           [Arienzo22]  A. Arienzo, G. Vivone, A. Garzelli, L. Alparone and J. Chanussot, "Full Resolution Quality Assessment of Pansharpening: Theoretical and hands-on Approaches", 
%                        IEEE Geoscience and Remote Sensing Magazine, 10(2):2-35, 2022.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Ds_F_index, Ds_F_map] = D_s_F_estPSF(I_F, I_PAN, I_MS_LR, ratio, sensorInf, S)

if (size(I_F) ~= size(I_MS_LR))
    error('The two images must have the same dimensions')
end

[N, M, Nb] = size(I_F);

if (rem(N,S) ~= 0)
    error('number of rows must be multiple of the block size')
end

if (rem(M,S) ~= 0)
    error('number of columns must be multiple of the block size')
end
sensorInf_est = sensorInf;
sensorInf_est.PSF_G = sensorInf_est.estPSF_G;
% Preallocation
Q_map_HR = zeros(size(I_F));
Q_map_LR = zeros(size(I_MS_LR));
    
    %%%%%%%%%% HR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % MTF filtering of fused MS
    %[I_F_degraded] = MTF(I_F, sensor, ratio);
    [I_F_degraded] = MTF_conv_sample(I_F,sensorInf_est,ratio,0);
    % Details extraction
    details_fused = I_F - I_F_degraded;
    
    % Ideal filtering of PAN
    I_PAN_degraded = idfilter23tap(I_PAN, ratio);
    % Details extraction
    details_pan = I_PAN - I_PAN_degraded;
    
    Q_det_HR = zeros(1,size(details_fused, 3));
    for ii = 1 : size(details_fused, 3)
        [Q_det_HR(ii), Q_map_HR(:,:,ii)] = img_qi_mod(details_fused(:,:,ii), details_pan, S);
        Q_back_HR = Q_map_HR(:,:,ii);
        % Clipping negative values
        Q_back_HR(Q_back_HR<0) = 0;
        Q_map_HR(:,:,ii) = Q_back_HR;
        Q_det_HR(ii) = mean2(Q_map_HR(:,:,ii));
    end

    %%%%%%%%%% LR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % MTF-filtering of PAN
    %panLR = MTF_PAN(I_PAN, sensor, ratio);
    %panLR = imresize(panLR, 1./ratio, 'nearest');
    panLR = imresize(I_PAN, 1./ratio);
    % MTF-filtering of original MS
    %[I_MS_degraded] = MTF(I_MS_LR, sensor, ratio);
    [I_MS_degraded] = MTF_conv_sample(I_MS_LR,sensorInf_est,ratio,0);
    % MS details extraction
    details_ms = I_MS_LR - I_MS_degraded;

    % Ideal filtering of degraded PAN
    pan_LR_degraded = idfilter23tap(panLR, ratio);
    % PAN details extraction
    details_pan_LR = panLR - pan_LR_degraded;

    Q_det_LR = zeros(1, size(details_ms,3));
    for ii = 1 : size(details_ms, 3)
        [Q_det_LR(ii), Q_map_LR(:,:,ii)] = img_qi_mod(details_ms(:,:,ii), details_pan_LR, round(S./ratio));
        Q_back_LR = Q_map_LR(:,:,ii);
        % Clipping of negative values
        Q_back_LR(Q_back_LR<0) = 0;
        Q_map_LR(:,:,ii) = Q_back_LR;
        Q_det_LR(ii) = mean2(Q_map_LR(:,:,ii));
    end
    
    %%%%%%%%%% Computation of FQNR spatial distortion %%%%%%%%%%%%%%%%%%%%%

    Ds_F_index = 0;
    Q_map_LR = imresize(Q_map_LR, [size(Q_map_HR,1) size(Q_map_HR,2)], 'cubic');
    Ds_F_map = zeros(size(Q_map_HR,1), size(Q_map_HR,2));
    
    for ii=1:size(details_ms, 3)     
        Ds_F_map = Ds_F_map + abs(Q_map_HR(:,:,ii)-Q_map_LR(:,:,ii));
        Ds_F_index = Ds_F_index + abs(Q_det_HR(ii)-Q_det_LR(ii));        
    end
    
    Ds_F_index = Ds_F_index/size(I_MS_LR,3);
    Ds_F_map   = Ds_F_map./size(I_MS_LR,3);
    
end