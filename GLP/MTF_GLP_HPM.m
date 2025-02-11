%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: 
%           MTF_GLP_HPM fuses the upsampled MultiSpectral (MS) and PANchromatic (PAN) images by 
%           exploiting the Modulation Transfer Function - Generalized Laplacian Pyramid (MTF-GLP) with High Pass Modulation (HPM) injection model algorithm. 
% 
% Interface:
%           I_Fus_MTF_GLP_HPM = MTF_GLP_HPM(I_PAN,I_MS,sensor,tag,ratio)
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
%           I_Fus_MTF_GLP_HPM:  MTF_GLP_HPM pansharpened image.
% 
% References:
%           [Aiazzi03]          B. Aiazzi, L. Alparone, S. Baronti, A. Garzelli, and M. Selva, “An MTF-based spectral distortion minimizing model for Pan-sharpening
%                               of very high resolution multispectral images of urban areas,?in Proceedings of URBAN 2003: 2nd GRSS/ISPRS Joint Workshop on
%                               Remote Sensing and Data Fusion over Urban Areas, 2003, pp. 90?4.
%           [Aiazzi06]          B. Aiazzi, L. Alparone, S. Baronti, A. Garzelli, and M. Selva, “MTF-tailored multiscale fusion of high-resolution MS and Pan imagery,?
%                               Photogrammetric Engineering and Remote Sensing, vol. 72, no. 5, pp. 591?96, May 2006.
%           [Vivone14a]         G. Vivone, R. Restaino, M. Dalla Mura, G. Licciardi, and J. Chanussot, “Contrast and error-based fusion schemes for multispectral
%                               image pansharpening,?IEEE Geoscience and Remote Sensing Letters, vol. 11, no. 5, pp. 930?34, May 2014.
%           [Vivone14b]         G. Vivone, L. Alparone, J. Chanussot, M. Dalla Mura, A. Garzelli, G. Licciardi, R. Restaino, and L. Wald, “A Critical Comparison Among Pansharpening Algorithms? 
%                               IEEE Transaction on Geoscience and Remote Sensing, 2014. (Accepted)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function I_Fus_MTF_GLP_HPM = MTF_GLP_HPM(I_PAN,I_MS,sensorInf,ratio)

imageHR = double(I_PAN);
I_MS = double(I_MS);

%%% Equalization
imageHR = repmat(imageHR,[1 1 size(I_MS,3)]);

for ii = 1 : size(I_MS,3)    
  imageHR(:,:,ii) = (imageHR(:,:,ii) - mean2(imageHR(:,:,ii))).*(std2(I_MS(:,:,ii))./std2(imageHR(:,:,ii))) + mean2(I_MS(:,:,ii));  
end
PAN_LP = MTF_conv_sample(imageHR,sensorInf,ratio);


I_Fus_MTF_GLP_HPM = I_MS .* (imageHR ./ (PAN_LP + eps));

end