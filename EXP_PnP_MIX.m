function I_Fus = EXP_PnP_MIX(I_MS_LR, I_PAN, sensorInf, ratio, base_method, EXP_Opts)
% Basic forms of mixed methods: I_Fus = I_MS_deb + IntensityRatio_deb .* repmat(imageHR - PAN_LP_deb,[1 1 size(I_MS,3)]);

I_MS = interpWrapper(I_MS_LR,ratio,sensorInf.upsampling);

% Histogram matching for PAN and PAN_LP generation
switch(base_method)
    case 'AWLP'
        decimation = 0;
        imageHR = repmat(I_PAN,[1 1 size(I_MS,3)]);
        for ii = 1 : size(I_MS,3)    
            imageHR(:,:,ii) = (imageHR(:,:,ii) - mean2(imageHR(:,:,ii))).*(std2(I_MS(:,:,ii))./std2(imageHR(:,:,ii))) + mean2(I_MS(:,:,ii));  
        end
        I_PAN_LP = genPAN_LP(imageHR, ratio, 'AWT', sensorInf);
    case 'AWLP-Dec'
         I_Fus = EXP_PnP_AWLP(I_MS_LR,I_MS,I_PAN,ratio,sensorInf,EXP_Opts);
         return;
        decimation = 1;
%         if strcmp(EXP_Opts.type, 'SR')
%             decimation = 1;
%         else
%             decimation = 0;
%         end
        imageHR = repmat(I_PAN,[1 1 size(I_MS,3)]);
        for ii = 1 : size(I_MS,3)    
            imageHR(:,:,ii) = (imageHR(:,:,ii) - mean2(imageHR(:,:,ii))).*(std2(I_MS(:,:,ii))./std2(imageHR(:,:,ii))) + mean2(I_MS(:,:,ii));  
        end
        I_PAN_LP = genPAN_LP(imageHR, ratio, 'AWT', sensorInf);
    case 'AWLP-H'
        decimation = 1;
         I_Fus = EXP_PnP_MTF_GLP_AWLP_Haze(I_MS_LR,I_PAN,I_MS,ratio,decimation, sensorInf,EXP_Opts);
         return;
%         decimation = 1;
%         Lp_MS = repmat(reshape(min(im2mat(I_MS),[],2), [1,1,size(I_MS,3)]), [size(I_MS,1) size(I_MS,2)]);
%         imageHR = double(I_PAN);
%         imageHR_LP = genPAN_LP(imageHR, ratio, 'Gauss', sensorInf);
%         alpha(1,1,:) = estimation_alpha(I_MS,imageHR_LP,'global');
%         I = sum((I_MS - Lp_MS) .* repmat(alpha,[size(I_MS,1) size(I_MS,2) 1]),3); 
%         imageHR = (imageHR - mean2(LPfilterGauss(imageHR,ratio))).*(std2(I)./std2(LPfilterGauss(imageHR,ratio))) + mean2(I);  
%         I_PAN_LP = genPAN_LP(imageHR, ratio, 'Gauss', sensorInf);
    case 'AWLP-H-noDec'
%         I_Fus = EXP_PnP_MTF_GLP_AWLP_Haze(I_MS_LR,I_PAN,I_MS,ratio,decimation, sensorInf,EXP_Opts);
%         return;
        decimation = 0;
        Lp_MS = repmat(reshape(min(im2mat(I_MS),[],2), [1,1,size(I_MS,3)]), [size(I_MS,1) size(I_MS,2)]);
        imageHR = double(I_PAN);
        imageHR_LP = genPAN_LP(imageHR, ratio, 'Gauss', sensorInf);
        alpha(1,1,:) = estimation_alpha(I_MS,imageHR_LP,'global');
        I = sum((I_MS - Lp_MS) .* repmat(alpha,[size(I_MS,1) size(I_MS,2) 1]),3); 
        imageHR = (imageHR - mean2(LPfilterGauss(imageHR,ratio))).*(std2(I)./std2(LPfilterGauss(imageHR,ratio))) + mean2(I);  
        I_PAN_LP = genPAN_LP(imageHR, ratio, 'Gauss', sensorInf);
end

% If the model does not do decimation by default, PAN does not have
% LR-scale images and cannot directly use need_interp == 0.
% Either input HR variables directly from outside and mark them as not
% interpolated, or force decimation.
% {need_interp, decimation}.
% {0, 0}: Pass in LR image. need_interp stays at 0, decimation replaces
% with 1.
% {0, 1}: Pass in LR image. need_interp stays at 0.
% {1, 0}: Deb method or less SR method. Pass in HR image, need_interp set
% to 0 (or decimation set to 1, pass in LR image).
% {1, 1}: Deb method or less SR method. Pass in LR image, need_interp held
% to 1.

if ~(EXP_Opts.need_interp || decimation)
    warning(['The selected PnP optimizer (' EXP_Opts.exp_method ') does not work for the current base method (' base_method '). Change to decimation']);
    decimation = 1;
else
    EXP_Opts.need_interp = EXP_Opts.need_interp && decimation;
end
I_PAN_LR = downsampleWrapper(I_PAN_LP, ratio, sensorInf.downsampling);
if decimation
    I_PAN_in = I_PAN_LR;
    I_MS_in = I_MS_LR;
else
    I_PAN_in = I_PAN_LP;
    I_MS_in = I_MS;
end

if EXP_Opts.matching_blur
    tap = 41;
    if decimation
        r = 1;
    else 
        r = ratio;
    end
    [I_MS_in, I_PAN_in] = matching_blur_level(I_MS_in, I_PAN_in, sensorInf.mGNyq_e, sensorInf.mGNyq, sensorInf, r, tap);
end
sensorInf_P = sensorInf;

I_MS_deb = Res_EXP(I_MS_in, ratio, sensorInf, EXP_Opts);


PAN_LP_deb = Res_EXP(I_PAN_in, ratio, sensorInf_P, EXP_Opts);



% I-component generation
switch(base_method)
    case {'AWLP', 'AWLP-Dec'} 
        SumImage=sum(I_MS_deb,3)/size(I_MS,3);
        IntensityRatio = zeros(size(I_MS_deb));
        for i=1:size(I_MS,3)
            IntensityRatio(:,:,i)=I_MS_deb(:,:,i)./(SumImage+eps);
        end
    case {'AWLP-H','AWLP-H-noDec'}
        Lp_MS_deb = repmat(reshape(min(im2mat(I_MS_deb),[],2), [1,1,size(I_MS_deb,3)]), [size(I_MS_deb,1) size(I_MS_deb,2)]);
        %%% Intensity
        I_deb = sum((I_MS_deb - Lp_MS_deb) .* repmat(alpha,[size(I_MS_deb,1) size(I_MS_deb,2) 1]),3); 
        IntensityRatio = (I_MS_deb - Lp_MS_deb) ./(repmat(I_deb,[1 1 size(I_MS_deb,3)]) + eps);
end

I_Fus = I_MS_deb + IntensityRatio .* (imageHR - PAN_LP_deb);

end