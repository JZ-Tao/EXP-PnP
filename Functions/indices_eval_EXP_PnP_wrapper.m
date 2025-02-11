function QI = indices_eval_EXP_PnP_wrapper(I_F,I_GT,I_MS,I_MS_LR, I_PAN, ratio,L,sensorInf,is_FS,Qblocks_size,flag_cut_bounds,dim_cut,th_values, using_estPSF)

if th_values
    I_F(I_F > 2^L-1) = 2^L-1;
    I_F(I_F < 0) = 0;
end
if is_FS
    if flag_cut_bounds
        I_PAN = I_PAN(1+dim_cut:end-dim_cut,1+dim_cut:end-dim_cut,:);
        I_F = I_F(1+dim_cut:end-dim_cut,1+dim_cut:end-dim_cut,:);
        I_MS = I_MS(1+dim_cut:end-dim_cut,1+dim_cut:end-dim_cut,:);
        I_MS_LR = I_MS_LR(1+floor(dim_cut/ratio):end-floor(dim_cut/ratio),1+floor(dim_cut/ratio):end-floor(dim_cut/ratio),:);
    end

    tag_case = 'unitary';
    [HQNR_alpha,HQNR_beta] = get_alpha_beta('HQNR',tag_case,sensorInf.sensor);
    [RQNR_alpha,RQNR_beta] = get_alpha_beta('RQNR',tag_case,sensorInf.sensor);
    if using_estPSF
        [QI.HQNR,QI.D_lambda_F,QI.D_S_Q, ~, ~, ~] = HQNR_estPSF(I_F, I_MS, I_MS_LR, I_PAN, ratio, Qblocks_size, sensorInf, HQNR_alpha, HQNR_beta);
        [QI.RQNR,~,QI.D_S_R, ~, ~, ~] = RQNR_estPSF(I_F, I_MS, I_PAN, ratio, Qblocks_size, sensorInf, RQNR_alpha, RQNR_beta);
    else
        [QI.HQNR,QI.D_lambda_F,QI.D_S_Q, ~, ~, ~] = HQNR(I_F, I_MS, I_MS_LR, I_PAN, ratio, Qblocks_size, sensorInf.sensor, HQNR_alpha, HQNR_beta);
        [QI.RQNR,~,QI.D_S_R, ~, ~, ~] = RQNR(I_F, I_MS, I_PAN, ratio, Qblocks_size, sensorInf.sensor, RQNR_alpha, RQNR_beta);
    end

    nb = size(I_GT, 3);
    tmp = 0;    
    for i=1:nb
         tmp = tmp + niqe(double(I_F(:,:,i)));
    end
    QI.NIQE = tmp/nb;
else
    if flag_cut_bounds
        I_GT = I_GT(1+dim_cut:end-dim_cut,1+dim_cut:end-dim_cut,:); %%%%%
        I_F = I_F(1+dim_cut:end-dim_cut,1+dim_cut:end-dim_cut,:);%%%%%%%
    end

    QI.Q = q2n(I_GT,I_F,Qblocks_size,Qblocks_size);
    QI.SAM = SAM(I_GT,I_F);
    QI.ERGAS = ERGAS(I_GT,I_F,ratio);
    QI.RMSE = RMSE(I_GT, I_F);
    psnr_out = PSNR(I_GT, I_F);
    QI.PSNR =psnr_out.ave;
    nb = size(I_GT, 3);
    tmp = 0;
    for i=1:nb
         tmp = tmp + ssim(double(I_F(:,:,i)), double(I_GT(:,:,i)),'DynamicRange',2^double(L)-1);
    end
    QI.SSIM = tmp/nb;
end
