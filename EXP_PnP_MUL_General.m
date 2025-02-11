function I_Fus = EXP_PnP_MUL_General(I_MS_LR, I_PAN, sensorInf, ratio, base_method, EXP_Opts)


n_band = size(I_MS_LR,3);
imageHR = repmat(I_PAN,[1 1 n_band]);
I_MS = interpWrapper(double(I_MS_LR),ratio,sensorInf.upsampling);


% offset

switch(base_method)
    case 'GLP-HPM-H'
        minMS = getMinMS(I_MS);

        LP_method = 'MTF';
        I_PAN_LP = genPAN_LP(I_PAN, ratio, LP_method, sensorInf);

        w = estimation_alpha(cat(3,ones(size(I_PAN_LP)),I_MS),I_PAN_LP,'global');
        Lp_p = w' * [1;squeeze(minMS)]; 
        Lp_MS = minMS;

    case 'BT-H'
        minMS = getMinMS(I_MS);
        Lp_MS = minMS;
        Lp_p = zeros([1 1 size(I_MS,3)]);
    case 'GLP-HPM-R'
        %%% Low resolution PAN image
        I_PAN_LP = MTF_conv_sample(I_PAN,sensorInf,ratio,2);
        Lp_p = zeros([1 1 size(I_MS,3)]);
        for ii = 1:size(I_MS,3)
            %%%% Regression coefficients
            MSB = I_MS(:,:,ii);
            C = cov(MSB(:),I_PAN_LP(:));
            g = C(1,2)./C(2,2);
            Lp_p(:,:,ii) = mean(I_PAN(:)) - mean(MSB(:))./g;
        end
        Lp_MS = zeros([1 1 size(I_MS,3)]);
    otherwise
        Lp_MS = zeros([1 1 size(I_MS,3)]);
        Lp_p = zeros([1 1 size(I_MS,3)]);
end

% Histogram matching for PAN and PAN_LP generation
switch(base_method)
   case 'MF'
        %imageHR_LP = MTF_conv_sample(imageHR,sensorInf,ratio);
        I_PAN_mod = PanHistEqualization(I_MS, imageHR, [], -1);
        textse= [0 1 0; 1 1 1; 0 1 0];
        int_meth='bilinear';
        lev=ceil(log2(ratio))+1;
        P = Pyr_Dec(I_PAN_mod,textse,lev,int_meth);
        I_PAN_mod = P(:,:,:,1);
        I_PAN_LP = P(:,:,:,lev);
        decimation = 1;
    case 'GLP-HPM'
         imageHR_LP = MTF_conv_sample(imageHR,sensorInf,ratio);
         I_PAN_mod = PanHistEqualization(I_MS, imageHR, imageHR_LP, -1);
         LP_method = 'MTF';
         decimation = 1;
         I_PAN_LP = genPAN_LP(I_PAN_mod, ratio, LP_method, sensorInf);
    case 'FE-HPM'
        %%%%%%%%%%%%%%%%%%%%%%%%%% Parameters setting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tap = 25; % dimension of the support (at least 3*ratio) 
        num_iter = 5; % max number of iteration (at least 3; not sensitive)
        lambda = 10^5; % coefficient to weight the regularization term 
        mu = 10^5; % coefficient to weight the regularization term
        th = 10^(-3); % threshold on the kernel (it cuts to 0 values below threshold)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        imageHR_LP = MTF_conv_sample(imageHR,sensorInf,ratio);
        I_PAN_mod = PanHistEqualization(I_MS, imageHR, imageHR_LP, 1);
        PSF_l = FE(I_MS,I_PAN,ratio,tap,lambda,mu,th,num_iter,'Basic');

        PAN_LP = zeros(size(imageHR));
        for ii = 1 : n_band
            PAN_LP(:,:,ii) = imfilter(I_PAN_mod(:,:,ii),PSF_l,'replicate');
        end
        I_PAN_LP = double(PAN_LP);
        decimation = 1;

    case {'GLP-HPM-H', 'GLP-HPM-R'}
        I_PAN_mod = imageHR;
        LP_method = 'MTF';
        decimation = 1;
        I_PAN_LP = genPAN_LP(I_PAN_mod, ratio, LP_method, sensorInf);
    case 'SFIM'
        I_PAN_mod = double(PanHistEqualization(I_MS, imageHR, [], -1));
        LP_method = 'box';
        decimation = 0;
        I_PAN_LP = genPAN_LP(I_PAN_mod, ratio, LP_method, sensorInf);
    case 'GLP-HPM-DS'
        mu = 0.95;
        gHRI = zeros(1,size(I_MS,3));
        gLRI = zeros(1,size(I_MS,3));
        g = zeros(1,size(I_MS,3));
        crop = 0;
        I_PAN_Filt = MTF_conv_sample(imageHR,sensorInf,ratio);
        for ii = 1 : size(I_MS,3)
            imPAN = imageHR(1+crop:end-crop, 1+crop:end-crop,ii); 
            imMSB = I_MS(1+crop:end-crop, 1+crop:end-crop,ii);
            imPANLRB = I_PAN_Filt(1+crop:end-crop, 1+crop:end-crop,ii);

            CMSPAN = cov(imMSB(:), imPAN(:));
            CMSPANLR = cov(imMSB(:), imPANLRB(:));
            CPANPANLR = cov(imPANLRB(:), imPAN(:));
            gHRI(ii) = CMSPAN(1,2)./CPANPANLR(1,2);
            gLRI(ii) = CMSPANLR(1,2)./CPANPANLR(1,2);
            g(ii) = mu.*gHRI(ii) + (1-mu).*gLRI(ii);
        end
        %%% Equalization
        for ii = 1 : size(I_MS,3)   
            I_PAN_mod(:,:,ii) = imageHR(:,:,ii) - mean2(imageHR(:,:,ii)) + (mean2(I_MS(:,:,ii))./g(ii));
        end
        decimation = 1;
        LP_method = 'MTF';
        I_PAN_LP = genPAN_LP(I_PAN_mod, ratio, LP_method, sensorInf);

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
if strcmp(base_method, 'GLP-HPM-DS')
    for ii = 1 : size(I_MS,3)   
        I_PAN_LR(:,:,ii) = I_PAN_LR(:,:,ii) - mean2(I_PAN_LR(:,:,ii))+ mean2(I_MS_LR(:,:,ii))./g(ii);
    end
end

if decimation
    I_PAN_in = I_PAN_LR;
    I_MS_in = I_MS_LR;
else
    I_PAN_in = I_PAN_LP;
    I_MS_in = I_MS;
end

I_PAN_mod = bsxfun(@minus, I_PAN_mod, Lp_p);
I_PAN_in = bsxfun(@minus, I_PAN_in, Lp_p);
I_MS_in = bsxfun(@minus, I_MS_in, Lp_MS);

if EXP_Opts.matching_blur
    tap = 41;
    if decimation
        r = 1;
    else 
        r = ratio;
    end
    [I_MS_in, I_PAN_in] = matching_blur_level(I_MS_in, I_PAN_in, sensorInf.mGNyq_e, sensorInf.mGNyq, sensorInf, r, tap);
end
Res_in = I_MS_in./ (I_PAN_in + eps);
Res_HR = Res_EXP(Res_in, ratio, sensorInf, EXP_Opts);

I_Fus = bsxfun(@plus, (Res_HR .* I_PAN_mod), Lp_MS);


end