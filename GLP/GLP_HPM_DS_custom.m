function I_GLP_HRI = GLP_HPM_DS_custom(I_MS, I_PAN, sensorInf, ratio)

if size(I_PAN,3) == 1
    imageHR = repmat(I_PAN,[1 1 size(I_MS,3)]);
else
    imageHR = I_PAN;
    I_PAN = imageHR(:,:,1);
end
I_PAN_Filt = MTF_conv_sample(imageHR,sensorInf,ratio);

gHRI = zeros(1,size(I_MS,3));
gLRI = zeros(1,size(I_MS,3));
g = zeros(1,size(I_MS,3));
crop = 0;
mu = 0.95;
for ii = 1 : size(I_MS,3)
    imPAN = I_PAN(1+crop:end-crop, 1+crop:end-crop,:); 
    imMSB = I_MS(1+crop:end-crop, 1+crop:end-crop,ii);
    imPANLRB = I_PAN_Filt(1+crop:end-crop, 1+crop:end-crop,ii);
        
    CMSPAN = cov(imMSB(:), imPAN(:));
    CMSPANLR = cov(imMSB(:), imPANLRB(:));
    CPANPANLR = cov(imPANLRB(:), imPAN(:));
    gHRI(ii) = CMSPAN(1,2)./CPANPANLR(1,2);
    gLRI(ii) = CMSPANLR(1,2)./CPANPANLR(1,2);
    g(ii) = mu.*gHRI(ii) + (1-mu).*gLRI(ii);
end

%I_GLP_HRI = Fusion_General(imageHR,I_MS,gHRI,I_PAN_Filt);
imageHR = double(I_PAN);
I_MS = double(I_MS);

%%% Equalization
imageHR = repmat(imageHR,[1 1 size(I_MS,3)]);
I_GLP_HRI = zeros(size(I_MS));
% for ii = 1 : size(I_MS,3)    
%   imageHR(:,:,ii) = (imageHR(:,:,ii) - mean2(imageHR(:,:,ii))).*(std2(I_MS(:,:,ii))./std2(imageHR(:,:,ii))) + mean2(I_MS(:,:,ii));  
% end
PAN_LP = MTF_conv_sample(imageHR,sensorInf,ratio);

for ii = 1 : size(I_MS,3)  
    I_GLP_HRI(:,:,ii) = I_MS(:,:,ii) .* ((imageHR(:,:,ii) - mean2(imageHR(:,:,ii)) + mean2(I_MS(:,:,ii))./g(ii)) ./ ((PAN_LP(:,:,ii) - mean2(PAN_LP(:,:,ii))+ mean2(I_MS(:,:,ii))./g(ii))+ eps));
end

end