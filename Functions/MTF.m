%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: 
%           MTF filters the image I_MS using a Gaussin filter matched with the Modulation Transfer Function (MTF) of the MultiSpectral (MS) sensor. 
% 
% Interface:
%           I_Filtered = MTF(I_MS,sensor,tag,ratio)
%
% Inputs:
%           I_MS:           MS image;
%           sensor:         String for type of sensor (e.g. 'WV2', 'IKONOS');
%           tag:            Image tag. Often equal to the field sensor. It makes sense when sensor is 'none'. It indicates the band number
%                           in the latter case;
%           ratio:          Scale ratio between MS and PAN.
%
% Outputs:
%           I_Filtered:     Output filtered MS image.
% 
% References:
%           [Aiazzi06]      B. Aiazzi, L. Alparone, S. Baronti, A. Garzelli, and M. Selva, “MTF-tailored multiscale fusion of high-resolution MS and Pan imagery,?
%                           Photogrammetric Engineering and Remote Sensing, vol. 72, no. 5, pp. 591?96, May 2006.
%           [Lee10]         J. Lee and C. Lee, “Fast and efficient panchromatic sharpening,?IEEE Transactions on Geoscience and Remote Sensing, vol. 48, no. 1,
%                           pp. 155?63, January 2010.
%           [Vivone14]      G. Vivone, L. Alparone, J. Chanussot, M. Dalla Mura, A. Garzelli, G. Licciardi, R. Restaino, and L. Wald, “A Critical Comparison Among Pansharpening Algorithms? 
%                           IEEE Transaction on Geoscience and Remote Sensing, 2014. (Accepted)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function I_Filtered = MTF(I_MS,sensor,ratio)

I_MS_LP = zeros(size(I_MS));
nBands = size(I_MS,3);

PSF_G = Blur_Kernel(nBands, sensor, ratio);   
for ii = 1 : nBands
    h = PSF_G(:,:,ii);
    I_MS_LP(:,:,ii) = imfilter(I_MS(:,:,ii),real(h),'circular');
end

I_Filtered= double(I_MS_LP);

end