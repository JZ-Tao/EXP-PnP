function I_Fus = EXP_PnP_ADD_General(I_MS_LR, I_PAN, sensorInf, ratio, base_method, EXP_Opts)

n_band = size(I_MS_LR,3);
imageHR = repmat(I_PAN,[1 1 n_band]);
I_MS = interpWrapper(I_MS_LR,ratio,sensorInf.upsampling);

imageHR_LP = MTF_conv_sample(imageHR,sensorInf,ratio);
% Histogram matching for PAN and PAN_LP generation
switch(base_method)
    case 'GLP' 
        I_PAN_mod = PanHistEqualization(I_MS, imageHR, imageHR_LP, -1);
        LP_method = 'MTF';
        decimation = 1;
        I_PAN_LP = genPAN_LP(I_PAN_mod, ratio, LP_method, sensorInf);
    case 'GLP-HS'
        % LR PAN
        mtf = sensorInf.mGNyq;%sensor_mtf(sensorInf.sensor);
        pL = pandeg(I_PAN,ratio,mtf);
       
        mth = 'HS-WLS';
        samp = 1;
        
        if samp == 1
            intp_mth = 'nearest';
            pL = imresize(pL,1/ratio,intp_mth);
            pL = interpWrapper(pL,ratio,sensorInf.upsampling);
        end
        % HM preprocessing
        [pL,a,b] = panhm2(pL,I_MS,'std');
        I_PAN_mod = a.*I_PAN+b;
        
        % HS estimation of injection gains 
        g = gainest_hs(I_MS,I_PAN_mod,pL,mth);
        I_PAN_mod = g.*I_PAN_mod;
        decimation = 1;
        I_PAN_LP = pandeg(I_PAN_mod,ratio,mtf);
    case 'GLP-REG'
        gHRI = zeros(1,size(I_MS,3));
        crop = 0;
        for ii = 1 : size(I_MS,3)
            imPAN = I_PAN(1+crop:end-crop, 1+crop:end-crop,:); 
            imMSB = I_MS(1+crop:end-crop, 1+crop:end-crop,ii);
            imPANLRB = imageHR_LP(1+crop:end-crop, 1+crop:end-crop,ii);

            CMSPAN = cov(imMSB(:), imPAN(:));    
            CPANPANLR = cov(imPANLRB(:), imPAN(:));
            gHRI(ii) = CMSPAN(1,2)./CPANPANLR(1,2);
        end
        coeff(1,1,:) = gHRI;
        G_HRI = repmat(coeff,[size(I_PAN,1) size(I_PAN,2) 1]);
        I_PAN_mod = G_HRI.*imageHR;
        LP_method = 'MTF';
        decimation = 1;
        I_PAN_LP = genPAN_LP(I_PAN_mod, ratio, LP_method, sensorInf);
    case 'ATWT'
        I_PAN_mod = PanHistEqualization(I_MS, imageHR, [], -1);
        LP_method = 'WT';
        decimation = 0; %%%%%%%%%%%%%%%%%%%%%%
        I_PAN_LP = genPAN_LP(I_PAN_mod, ratio, LP_method, sensorInf);

    case 'ATWT-Dec'
        I_PAN_mod = PanHistEqualization(I_MS, imageHR, [], -1);
        LP_method = 'WT';
        decimation = 1; %%%%%%%%%%%%%%%%%%%%%%
        I_PAN_LP = genPAN_LP(I_PAN_mod, ratio, LP_method, sensorInf);
%         decimation = 1;
%         I_Fus = EXP_PnP_ATWT(I_MS_LR,I_MS,I_PAN,ratio,sensorInf,EXP_Opts,decimation);
%         return;
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

if EXP_Opts.anti_aliasing

    if decimation
        I_PAN_in = interpWrapper(I_PAN_in,ratio,sensorInf.upsampling);
    else
        I_PAN_in = I_PAN_LP;
    end
    [I_MS_HR, I_PAN_HR] = matching_blur_level(I_MS, I_PAN_in, sensorInf.mGNyq_e, sensorInf.mGNyq, sensorInf, ratio);
    Res_HR = I_MS_HR - I_PAN_HR;
else
    Res_in = I_MS_in - I_PAN_in;
    Res_HR = Res_EXP(Res_in, ratio, sensorInf, EXP_Opts);
end
I_Fus = I_PAN_mod + Res_HR;

end