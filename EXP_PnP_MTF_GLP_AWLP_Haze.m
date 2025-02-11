function I_Fus_MTF_GLP_AWLP = EXP_PnP_MTF_GLP_AWLP_Haze(I_MS_LR,I_PAN,I_MS,ratio,decimation, sensorInf,EXP_Opts)

vCorr = zeros(1, size(I_MS,3));
for ii = 1 : size(I_MS,3)
    b = I_MS(:,:,ii);
    vCorr(ii) = min(b(:));
end

v3Corr = zeros(1,1,size(I_MS,3));
v3Corr(1,1,:) = vCorr;
LCorr = repmat(v3Corr, [size(I_MS,1) size(I_MS,2)]);

imageHR = double(I_PAN);
I_MS = double(I_MS);

%%% Intensity
imageHR_LP = LPfilterGauss(imageHR,ratio);

h = estimation_alpha(I_MS,imageHR_LP,'global');
alpha(1,1,:) = h;
I = sum((I_MS - LCorr) .* repmat(alpha,[size(I_MS,1) size(I_MS,2) 1]),3); 

imageHR = (imageHR - mean2(LPfilterGauss(imageHR,ratio))).*(std2(I)./std2(LPfilterGauss(imageHR,ratio))) + mean2(I);  

IntensityRatio = (I_MS - LCorr) ./(repmat(I,[1 1 size(I_MS,3)]) + eps);

PAN_LP = LPfilterGauss(imageHR,ratio);
imageHR = IntensityRatio.*repmat(imageHR,[1 1 size(I_MS,3)]);
PAN_LP = IntensityRatio.*repmat(PAN_LP,[1 1 size(I_MS,3)]);

if decimation
    I_PAN_in = downsampleWrapper(PAN_LP, ratio, sensorInf.downsampling);
    I_MS_in = I_MS_LR;
else
    EXP_Opts.need_interp = 0;
    I_PAN_in = PAN_LP;
    I_MS_in = I_MS;
end
% 
% if decimation
%     PAN_LR = downsampleWrapper(PAN_LP, ratio, sensorInf.downsampling);
%     Res = I_MS_LR - PAN_LR;
% else
%     EXP_Opts.need_interp = 0;
%     Res = I_MS - PAN_LP;
% end

if EXP_Opts.matching_blur
    tap = 41;
    if decimation
        r = 1;
    else 
        r = ratio;
    end
    [I_MS_in, I_PAN_in] = matching_blur_level(I_MS_in, I_PAN_in, sensorInf.mGNyq_e, sensorInf.mGNyq, sensorInf, r, tap);
end
Res = I_MS_in - I_PAN_in;
Res_HR = Res_EXP(Res, ratio, sensorInf, EXP_Opts);

I_Fus_MTF_GLP_AWLP = Res_HR + imageHR;

end