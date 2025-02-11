%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: 
%           Regression-based Quality with No Reference (RQNR): Spatial distortion index.
% 
% Interface:
%           [Ds_R_index, Ds_R_map] = D_s_R(I_F, I_PAN, S)
%
% Inputs:
%           I_F:         Pansharpened image;
%           I_PAN:       Panchromatic image;
%           S:           Block size (suggested 32).
% 
% Outputs:
%           Ds_R_index:  D_s RQNR index;
%           Ds_R_map:    D_s RQNR map;
% 
% References:
%           [Khan09]     M. M. Khan, L. Alparone, and J. Chanussot, "Pansharpening quality assessment using the modulation transfer functions of instruments", 
%                        IEEE Trans. Geosci. Remote Sens., vol. 11, no. 47, pp. 3880–3891, Nov. 2009.
%           [Vivone15]   G. Vivone, L. Alparone, J. Chanussot, M. Dalla Mura, A. Garzelli, G. Licciardi, R. Restaino, and L. Wald, “A Critical Comparison Among Pansharpening Algorithms”, 
%                        IEEE Transactions on Geoscience and Remote Sensing, vol. 53, no. 5, pp. 2565–2586, May 2015.
%           [Alparone18] L. Alparone, A. Garzelli, and G. Vivone, "Spatial consistency for fullscale assessment of pansharpening",
%                        Proc. IGARSS, 2018, pp. 5132–5134.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Ds_R_index, Ds_R_map] = D_s_R(I_F, I_PAN, S)
    
% Multivariate linear regresssion between fused MS and original PAN
[~,~,CD,CD_map,~,~] = LSR(I_F, I_PAN, S);
    
% Spatial distortion index
Ds_R_index = 1 - CD;

% Spatial distortion map
Ds_R_map = 1 - CD_map;
    
end