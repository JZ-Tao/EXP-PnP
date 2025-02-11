%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: 
%           Least Squares Regression.
% 
% Interface:
%           [w,res,cd,cd_map,varPAN,err_reg] = LSR(I_MS, I_PAN, S)
%
% Inputs:
%           I_MS:        Pansharpened image;
%           I_PAN:       Panchromatic image;
%           S:           Block size (optional); Default value: 32;
% 
% Outputs:
%           w:           Estimated regression coefficients
%           res:         Regression residual value
%           cd:          Coefficient of determination
%           cd_map:      Local map of the coefficient of determination
%           var_pan:     Variance of PAN
%           err_reg:     Spatially-variant regression residual
% 
% References:
%           [Khan09]     M. M. Khan, L. Alparone, and J. Chanussot, "Pansharpening quality assessment using the modulation transfer functions of instruments", 
%                        IEEE Trans. Geosci. Remote Sens., vol. 11, no. 47, pp. 3880–3891, Nov. 2009.
%           [Vivone15]   G. Vivone, L. Alparone, J. Chanussot, M. Dalla Mura, A. Garzelli, G. Licciardi, R. Restaino, and L. Wald, “A Critical Comparison Among Pansharpening Algorithms”, 
%                        IEEE Transactions on Geoscience and Remote Sensing, vol. 53, no. 5, pp. 2565–2586, May 2015.
%           [Alparone18] L. Alparone, A. Garzelli, and G. Vivone, "Spatial consistency for fullscale assessment of pansharpening",
%                        Proc. IGARSS, 2018, pp. 5132–5134.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [w,res,cd,cd_map,var_pan,err_reg] = LSR(I_F,I_PAN,S)

if nargin < 3, S=32; end

% Vectorization of PAN and fused MS
IHc = reshape(I_PAN,[numel(I_PAN) 1]);
ILRc = reshape(I_F,[size(I_F,1)*size(I_F,2) size(I_F,3)]);

% Multivariate linear regression
w = ILRc\IHc;
alpha(1,1,:) = w;

% Fitted Least squares intensity
I_R = sum(I_F .* repmat(alpha,[size(I_F,1) size(I_F,2) 1]),3);

% Space-varying least squares error
err_reg = I_PAN(:) - I_R(:);

% Regression residual value
res = sqrt(mean(err_reg.^2));

% Coefficient of determination
cd = 1 - (var(err_reg)./var(I_PAN(:)));

% Variance of PAN
var_pan= var(I_PAN(:));
err_reg=reshape(err_reg,[size(I_PAN,1),size(I_PAN,2)]);

% Local map of the coefficient of determination
if rem(S,2) == 0
    s = S -1;
else
    s = S;
end

err_var_map = stdfilt(reshape(err_reg, [size(I_PAN,1),size(I_PAN,2)]), ones(s,s)).^2;
pan_var_map = stdfilt(I_PAN, ones(s,s)).^2;

% Clipping
err_var_map(err_var_map<0)=0;
pan_var_map(pan_var_map<0)=0;
        
cd_map = 1-(err_var_map./pan_var_map);

% Clipping 
cd_map(cd_map>1) = 1;
cd_map(cd_map<0) = 0;

end