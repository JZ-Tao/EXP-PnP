%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: 
%           Quality with No Reference (QNR): Spatial distortion index.
% 
% Interface:
%           [Ds_Q_index,Ds_Q_map] = D_s_Q(I_F,I_MS,I_MS_LR,I_PAN,ratio,S)
%
% Inputs:
%           I_F:          Pansharpened image;
%           I_MS:         MS image resampled to panchromatic scale;
%           I_MS_LR:      Original MS image;
%           I_PAN:        Panchromatic image;
%           ratio:        Scale ratio between MS and PAN. Pre-condition: Integer value;
%           S:            Block size (optional); Default value: 32;
% 
% Outputs:
%           D_s_Q_index:  D_s QNR index;
%           D_s_Q_map:    D_s QNR map;
% 
% References:
%           [Alparone08]  L. Alparone, B. Aiazzi, S. Baronti, A. Garzelli, F. Nencini, and M. Selva, "Multispectral and panchromatic data fusion assessment without reference,"
%                         Photogrammetric Engineering and Remote Sensing, vol. 74, no. 2, pp. 193–200, February 2008. 
%           [Vivone15]       G. Vivone, L. Alparone, J. Chanussot, M. Dalla Mura, A. Garzelli, G. Licciardi, R. Restaino, and L. Wald, “A Critical Comparison Among Pansharpening Algorithms”, 
%                            IEEE Transactions on Geoscience and Remote Sensing, vol. 53, no. 5, pp. 2565–2586, May 2015.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Ds_index, Ds_map] = D_s_Q(I_F, I_MS, I_MS_LR, I_PAN, ratio, S)

% Flag for choosing the computation domain
%  (if 0: domain of expanded MS; 
%   if 1: domain of original MS as proposed in the original QNR paper.)
flag_orig_paper = 0;

if (size(I_F) ~= size(I_MS))
    error('The two images must have the same dimensions')
end

[N, M, Nb] = size(I_F);

if (rem(N,S) ~= 0)
    error('number of rows must be multiple of the block size')
end

if (rem(M,S) ~= 0)
    error('number of columns must be multiple of the block size')
end

if flag_orig_paper == 0
    %%% Opt.1  Domain of expanded MS; 
    pan_filt = interp23tap(imresize(I_PAN,1./ratio),ratio);
else
    %%% Opt.2  Domain of original MS as proposed in the original QNR paper.
    pan_filt = imresize(I_PAN,1./ratio);
end

Ds_index = 0;

for i = 1:Nb
        band1 = I_F(:,:,i);
        band2 = I_PAN;
        fun_uqi = @(bs) uqi(bs.data,...
            band2(bs.location(1):bs.location(1)+S-1,...
            bs.location(2):bs.location(2)+S-1));
        Qmap_high = blockproc(band1,[S S],fun_uqi);        
        Q_high = mean2(Qmap_high);
        if i==1
               Ds_map = zeros(size(Q_high));
            else
        end
        if flag_orig_paper == 0
            %%% Opt.1 (Domain of expanded MS)
            band1 = I_MS(:,:,i);
            band2 = pan_filt;
            fun_uqi = @(bs) uqi(bs.data,...
            band2(bs.location(1):bs.location(1)+S-1,...
            bs.location(2):bs.location(2)+S-1));
            Qmap_low = blockproc(band1,[S S],fun_uqi);
            
        else
            %%%% Opt.2 (Domain of original MS)
            band1 = I_MS_LR(:,:,i);
            band2 = pan_filt;
            fun_uqi = @(bs) uqi(bs.data,...
            band2(bs.location(1):bs.location(1)+S/ratio-1,...
            bs.location(2):bs.location(2)+S/ratio-1));
            Qmap_low = blockproc(band1,[S/ratio S/ratio],fun_uqi);
        end        
        Q_low = mean2(Qmap_low);
        Ds_index= Ds_index + abs(Q_high-Q_low);
        Ds_map= Ds_map + abs(Qmap_high-Qmap_low);
end
Ds_index = (Ds_index/Nb);
Ds_map= (Ds_map/Nb);

end

% Add. function: UIQI on x and y images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Q = uqi(x,y)

x = double(x(:));
y = double(y(:));
mx = mean(x);
my = mean(y);
C = cov(x,y);

Q = 4 * C(1,2) * mx * my / (C(1,1)+C(2,2)) / (mx^2 + my^2);  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end