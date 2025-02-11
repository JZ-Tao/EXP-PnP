function Res_HR = Res_EXP(Res, ratio, sensorInf, Opts)
exp_method = Opts.exp_method;

if ~isfield(Opts, 'sigma') 
    Opts.sigma = 0.0008; 
end
if ~isfield(Opts, 'global_lambda') 
    Opts.global_lambda = -1; 
end
if Opts.need_interp
    Res = interpWrapper(Res, ratio, sensorInf.upsampling);

end

PSF_G = sensorInf.PSF_G;

switch(Opts.normalize_type)
    case 0
        max_v = 1;
    case 1
        [M_Res, PS] = mapminmax(im2mat(Res),0,1);
        Res = mat2im(M_Res, size(Res,1));
    case 2
        L = Opts.L;

        max_v = max(2^L-1, 1);%max(max(I_MS(:)), max(I_PAN(:)));
        Res = Res./max_v;
        %sigma = sigma./max_v;        
    case 3

        tol = [0.0001 0.9999];
        for i=1:size(Res,3)
            b = double(Res(:,:,i));
            t(1) =  quantile(b(:),tol(1));
            t(2)  = quantile(b(:),tol(2));       
            b(b<t(1))=t(1);
            b(b>t(2))=t(2);
%             b = (b-t(1))/(t(2)-t(1));
            Res(:,:,i) = reshape(b, size(Res,1), size(Res,2));
        end
        [M_Res, PS] = mapminmax(im2mat(Res),0,1);
        Res = mat2im(M_Res, size(Res,1));
end


hsize_b = size(PSF_G, 3);
n_band = size(Res, 3);
done = 0;
if hsize_b == 1
    PSF_G = repmat(PSF_G,[1 1 n_band]);
end
if strcmp(exp_method, 'Donothing') || strcmp(exp_method, 'EXP')
    Res_HR = Res;
    done = 1;
end

sensorInf_org = sensorInf;
if Opts.deb_PSF_e
    PSF = sensorInf.PSF_e;
    sensorInf.PSF_G = PSF;
else
    PSF = PSF_G;
end

if strcmp(exp_method, 'LMG')
    kernel_name = ['LMG_k_dat' num2str(Opts.dataset_idx) '_FS_' num2str(Opts.is_FS) '_sz' num2str(size(Res,1)) '.mat'];
    LMG_opts.prescale = 1; %% downsampling
    LMG_opts.xk_iter = 5; %% the iterations
    LMG_opts.gamma_correct = 1.0;
    LMG_opts.k_thresh = 20;
    
    if exist(kernel_name, 'file')
        load(kernel_name, 'kernel');
        LMG_opts.kernel = kernel;
        LMG_opts.kernel_size = size(kernel,1);
        new_kernel = 0;
    else
        new_kernel = 1;
        LMG_opts.kernel_size = 27;
    end
    debRef = Opts.I_MS(:,:,Opts.color_band_idx);
    debRef = debRef./(max(debRef(:)));
    [Res_HR, kernel] = LMG_deblur_wrapper(debRef,Res,LMG_opts);
    if new_kernel
        save(kernel_name, 'kernel');
    end
    done = 1;
end
if strcmp(exp_method, 'GraphBID')
    debRef = Opts.I_MS(:,:,Opts.color_band_idx);
    debRef = debRef./(max(debRef(:)));
    Res_HR = GraphBID_wrapper(debRef,Res);
    done = 1;
end

if ~done && strcmp(Opts.type, 'Deblur') && ~strcmp(exp_method, 'PnP-ADMM_Deblur')
    Res_HR = MSDeblur(Res, PSF, Opts.L, exp_method);
    done = 1;
end
if strcmp(exp_method, 'PnP-ADMM-SR') || strcmp(exp_method, 'PnP-ADMM-Deblur')
    denoiser_name = Opts.ADMM.denoiser;%'FFDNet'; 
    if ~iscell(denoiser_name)
        if strcmp(version('-release'), '2014a')
            if strcmp(denoiser_name, 'FFDNet') || strcmp(denoiser_name, 'CNN')
                    denoiser_name = 'TV';
                    warning('FFDNet or CNN denoiser only supported in matlab r2015b or later. Changed to TV');
            end 
        end
        switch denoiser_name
           case 'MCWNNM'
                denoiser = @wrapper_MCWNNM;
                %lambda = 0.00001;
            case 'WNNM'
                denoiser = @wrapper_WNNM;
                %lambda = 0.00001;
            case 'RF'
                denoiser = @wrapper_RF;
                %lambda = 0.0003;
            case 'NLM'
                denoiser = @wrapper_NLM;
                %lambda = 0.005;
            case 'BM3D'
                denoiser = @wrapper_BM3D;
                %lambda = 0.000003;%0.000025;
            case 'TV'
                denoiser = @wrapper_TV;
                %lambda = 0.000002;%0.000025;%0.0001;
            case 'FFDNet'
                %lambda = 0.000001;%0.000005;%0.00001;
                denoiser = @wrapper_FFDNet;
            case 'IRCNN'
                %lambda = 0.00002;
                denoiser = @wrapper_CNN;
            case 'Bilateral'
                %lambda = 0.0001;
                denoiser = @wrapper_Bilateral;
            otherwise
                error('unknown denoiser \n');
        end
    end

    if Opts.global_lambda ~= -1
        lambda = Opts.global_lambda;
    else
        lambda = 5e-6;
    end
    %lambda = 5e-6;
    ADMM_opts.method = denoiser_name;
    %optional parameters
    ADMM_opts.rho     = 1;
    ADMM_opts.gamma   = 1;
    ADMM_opts.max_itr = 25;
    ADMM_opts.print   = true;
    ADMM_opts.GT_mode = Opts.GT_mode;
    ADMM_opts.step_gamma = Opts.step_gamma;
    
    Res_LP = Res;
  
    if strcmp(exp_method, 'PnP-ADMM-Deblur') % Deblur

        max(Res_LP(:))
        min(Res_LP(:))
        Res_HR = PnP_ADMM_deblur(Res_LP,sensorInf,lambda,denoiser,ADMM_opts);
        
        
    else

        Res_HR = PnP_ADMM_SR(Res,sensorInf,ratio,lambda,denoiser,ADMM_opts); % PnP_ADMM_SR PAMP_SR
    end

    done = 1;
end


if strcmp(exp_method, 'FFDNet')

    if size(Res,3) == 3
        Res_HR = wrapper_FFDNet(Res, 0);
    else
        for i = 1:size(Res,3)
            Res_HR(:,:,i) = wrapper_FFDNet(Res(:,:,i), 0);
        end
    end
    done = 1;
end

if strcmp(exp_method, 'IRCNN-deb')
    if Opts.global_lambda ~= -1
        IRCNN_lambda = Opts.global_lambda;
    else
        IRCNN_lambda = 0.0005;
    end

    
    if size(Res,3) == 3
        Res_HR = IRCNN_deb_wrapper(Res, PSF(:,:,1), IRCNN_lambda);
    else
        for i = 1:size(Res,3)
            Res_HR(:,:,i) = IRCNN_deb_wrapper(Res(:,:,i), PSF(:,:,i), IRCNN_lambda);
        end
    end
    done = 1;
end
if strcmp(exp_method, 'IRCNN-SR')
    if Opts.global_lambda ~= -1
        IRCNN_lambda = Opts.global_lambda;
    else
        error('No IRCNN lambda');
        IRCNN_lambda = 0.0005;
    end
    Res_HR = zeros(size(imresize(Res, ratio)));
    %IRCNN_lambda = 7.5e-5;
    
    if size(Res,3) == 3
        Res_HR = IRCNN_SR_wrapper(Res, sensorInf, PSF(:,:,1), ratio, IRCNN_lambda, Opts);
    else
        for i = 1:size(Res,3)
            Res_HR(:,:,i) = IRCNN_SR_wrapper(Res(:,:,i), sensorInf, PSF(:,:,i), ratio, IRCNN_lambda, Opts);
        end
    end
    done = 1;
end


if strcmp(Opts.type, 'SR') && ~done
    SROpts.L = Opts.L;
    SROpts.normalize_type = 0;
    SROpts.upsampling = sensorInf.upsampling;
    if Opts.extra_deb == 1
        Res = MSDeblur(Res, Opts.PSF_LR, Opts.L, 'Hyper-Laplacian');
    end

    is_LR = (size(Res,1) ~= Opts.HR_size(1)) && (size(Res,2) ~= Opts.HR_size(2));
    Res_HR = SROPwrapper(Res, is_LR, exp_method, ratio, SROpts);
end

switch(Opts.normalize_type)
    case {1,3}
        Res_HR = mat2im(mapminmax('reverse', im2mat(Res_HR), PS), size(Res_HR,1));
    case {0,2}
        Res_HR = double(Res_HR).*max_v;
end




