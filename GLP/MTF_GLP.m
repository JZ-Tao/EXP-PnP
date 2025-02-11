%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: 
%           MTF_GLP fuses the upsampled MultiSpectral (MS) and PANchromatic (PAN) images by 
%           exploiting the Modulation Transfer Function - Generalized Laplacian Pyramid (MTF-GLP) algorithm. 
% 
% Interface:
%           I_Fus_MTF_GLP = MTF_GLP(I_PAN,I_MS,sensor,tag,ratio)
%
% Inputs:
%           I_PAN:              PAN image;
%           I_MS:               MS image upsampled at PAN scale;
%           sensor:             String for type of sensor (e.g. 'WV2','IKONOS');
%           tag:                Image tag. Often equal to the field sensor. It makes sense when sensor is 'none'. It indicates the band number
%                               in the latter case;
%           ratio:              Scale ratio between MS and PAN. Pre-condition: Integer value.
%
% Outputs:
%           I_Fus_MTF_GLP:      MTF_GLP pansharpened image.
% 
% References:
%           [Aiazzi02]          B. Aiazzi, L. Alparone, S. Baronti, and A. Garzelli, “Context-driven fusion of high spatial and spectral resolution images based on
%                               oversampled multiresolution analysis,?IEEE Transactions on Geoscience and Remote Sensing, vol. 40, no. 10, pp. 2300?312, October
%                               2002.
%           [Aiazzi06]          B. Aiazzi, L. Alparone, S. Baronti, A. Garzelli, and M. Selva, “MTF-tailored multiscale fusion of high-resolution MS and Pan imagery,?
%                               Photogrammetric Engineering and Remote Sensing, vol. 72, no. 5, pp. 591?96, May 2006.
%           [Vivone14a]         G. Vivone, R. Restaino, M. Dalla Mura, G. Licciardi, and J. Chanussot, “Contrast and error-based fusion schemes for multispectral
%                               image pansharpening,?IEEE Geoscience and Remote Sensing Letters, vol. 11, no. 5, pp. 930?34, May 2014.
%           [Vivone14b]         G. Vivone, L. Alparone, J. Chanussot, M. Dalla Mura, A. Garzelli, G. Licciardi, R. Restaino, and L. Wald, “A Critical Comparison Among Pansharpening Algorithms? 
%                               IEEE Transaction on Geoscience and Remote Sensing, 2014. (Accepted)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function I_Fus_MTF_GLP = MTF_GLP(I_PAN,I_MS,sensorInf,ratio)

imageHR = double(I_PAN);
I_MS = double(I_MS);

%%% Equalization
nBands = size(I_MS,3);
if size(imageHR,3) ~= nBands
    imageHR = repmat(imageHR,[1 1 nBands]);
end

for ii = 1 : nBands
   imageHR(:,:,ii) = (imageHR(:,:,ii) - mean2(imageHR(:,:,ii))).*(std2(I_MS(:,:,ii))./std2(imageHR(:,:,ii))) + mean2(I_MS(:,:,ii));
end

PAN_LP = MTF_conv_sample(imageHR,sensorInf,ratio,2);

I_Fus_MTF_GLP = I_MS + imageHR - PAN_LP;

end