function [func_handle, xls_prefix] = genFuncHandle(func_name, sensorInf, ratio, L, dat_tag)
    switch(func_name)
        case 'EXP'
            func_handle = @(I_MS_LR, I_MS,I_PAN,ratio) I_MS;
        case 'Indusion'
            func_handle = @(I_MS_LR, I_MS,I_PAN,ratio) Indusion(I_PAN,I_MS_LR,ratio);
        case 'PWMBF'
                %%%%%%%%%%%%%%%%%%%%%%%%%% Parameters setting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                r = 3;
                degrade=0;
                wavelet=1;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                func_handle = @(I_MS_LR, I_MS,I_PAN,ratio) PWMBF(I_PAN,I_MS_LR,ratio,r,wavelet,degrade);
        case 'STEM-MA'
            STEM_opts = init_STEM_options('MA', dat_tag);
            func_handle = @(I_MS_LR, I_MS,I_PAN,ratio) STEM_wrapper(I_MS_LR, I_PAN, ratio, sensorInf, L, STEM_opts);
            xls_prefix = 'M';
        case 'STEM-MS'
            STEM_opts = init_STEM_options('MS', dat_tag);
            func_handle = @(I_MS_LR, I_MS,I_PAN,ratio) STEM_wrapper(I_MS_LR, I_PAN, ratio, sensorInf, L, STEM_opts);
            xls_prefix = 'M';
        case 'BDSD-PC'
            func_handle = @(I_MS_LR, I_MS,I_PAN,ratio) BDSD_PC(I_MS,I_PAN,ratio,sensorInf);
            xls_prefix = 'C';

        case 'BT-H'
            func_handle = @(I_MS_LR, I_MS,I_PAN,ratio) BroveyRegHazeMin(I_MS,I_PAN,ratio);
            xls_prefix = 'C';
        case 'SFIM'
            func_handle = @(I_MS_LR, I_MS,I_PAN,ratio) SFIM(I_MS,I_PAN,ratio);
            xls_prefix = 'M';
        case 'IHS'
            func_handle = @(I_MS_LR, I_MS,I_PAN,ratio) IHS(I_MS,I_PAN);
            xls_prefix = 'C';
        case 'GS'
            func_handle = @(I_MS_LR, I_MS,I_PAN,ratio) GS(I_MS,I_PAN);
            xls_prefix = 'C';
        case 'Brovey'
            func_handle = @(I_MS_LR, I_MS,I_PAN,ratio) Brovey(I_MS,I_PAN);
            xls_prefix = 'C';
        case 'MF'
            func_handle = @(I_MS_LR, I_MS,I_PAN,ratio) MF_HG_Pansharpen(I_PAN,I_MS,ratio);
            xls_prefix = 'M';
        case 'GSA'
            func_handle = @(I_MS_LR, I_MS,I_PAN,ratio) GSA(I_MS,I_PAN,I_MS_LR,ratio);
            xls_prefix = 'C';
        case 'GS-Segm'
            %%%%%%%%%%%%%%%%%%%%%%%%%% Parameters setting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            PS_algorithm = 'GSA'; % Pansharpening algorithm 
            n_segm = 5; % Number of segments
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            func_handle = @(I_MS_LR, I_MS,I_PAN,ratio) GS_Segm(I_MS,I_PAN,gen_LP_image(PS_algorithm,I_MS,I_PAN,I_MS_LR,ratio,sensorInf.sensor), k_means_clustering(I_MS,n_segm));
            xls_prefix = 'C';
        case 'AWLP'
            func_handle = @(I_MS_LR, I_MS,I_PAN,ratio) AWLP(I_MS,I_PAN,ratio);
            xls_prefix = 'M';
        case 'AWLP-H'
            func_handle = @(I_MS_LR, I_MS,I_PAN,ratio) MTF_GLP_AWLP_Haze(I_PAN,I_MS,ratio,1,sensorInf);
            xls_prefix = 'M';
        case 'ATWT'
            func_handle = @(I_MS_LR, I_MS,I_PAN,ratio) ATWT(I_MS,I_PAN,ratio);
            xls_prefix = 'M';
        case 'ATWT-M2'
            func_handle = @(I_MS_LR, I_MS,I_PAN,ratio) ATWT_M3(I_MS,I_PAN,ratio);
            xls_prefix = 'M';
        case 'ATWT-M3'
            func_handle = @(I_MS_LR, I_MS,I_PAN,ratio) ATWT_M3(I_MS,I_PAN,ratio);
            xls_prefix = 'M';
        case 'PCA'
            func_handle = @(I_MS_LR, I_MS,I_PAN,ratio) PCA(I_MS,I_PAN);
            xls_prefix = 'C';
        case 'GFPCA'
            func_handle = @(I_MS_LR, I_MS,I_PAN,ratio) GFPCA(I_MS_LR,I_PAN, 3, 5, 0.001^2);
            xls_prefix = 'M';
        case {'MTF-GLP-HPM','GLP-HPM'}
            func_handle = @(I_MS_LR, I_MS,I_PAN,ratio) MTF_GLP_HPM(I_PAN, I_MS, sensorInf, ratio);
            xls_prefix = 'M';
        case 'MTF-GLP-HPM_R'
            func_handle = @(I_MS_LR, I_MS,I_PAN,ratio) MTF_GLP_HPM_R(I_MS,I_PAN,sensorInf.sensor,ratio);
            xls_prefix = 'M';
        case {'MTF-GLP-HPM-H','GLP-HPM-H'}
            func_handle = @(I_MS_LR, I_MS,I_PAN,ratio) MTF_GLP_HPM_Haze_min(I_MS,I_PAN,sensorInf.sensor,ratio,1);
            xls_prefix = 'M';
        case {'GLP-HPM-DS','MTF-GLP-HPM_DS'}
            func_handle = @(I_MS_LR, I_MS,I_PAN,ratio) GLP_HPM_DS_custom(I_MS, I_PAN, sensorInf, ratio);
            xls_prefix = 'M';
        case 'FE-HPM'
            %%%%%%%%%%%%%%%%%%%%%%%%%% Parameters setting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            tap = 25; % dimension of the support (at least 3*ratio) 
            num_iter_max = 5; % max number of iteration (at least 3; not sensitive)
            lambda = 10^5; % coefficient to weight the regularization term 
            mu = 10^5; % coefficient to weight the regularization term
            threshold = 10^(-3); % threshold on the kernel (it cuts to 0 values below threshold)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            func_handle = @(I_MS_LR, I_MS,I_PAN,ratio) FE_HPM(I_PAN,I_MS,ratio,max(tap,3*ratio),lambda,mu,threshold,max(num_iter_max,3),'Basic'); 
            xls_prefix = 'M';
        case {'MTF-GLP','GLP'}
            func_handle = @(I_MS_LR, I_MS,I_PAN,ratio) MTF_GLP(I_PAN,I_MS,sensorInf,ratio);
            xls_prefix = 'M';
        case 'C-GLP-Robust'
            nclusters = 5;
            %%% Multiplicative coefficient (selected as in the original paper)
            zeta = 1;
            %%% Threshold NDVI (selected as in the original paper)
            thNDVI = 0.5;
            %%% Select filter for PAN image
            flagMTFPAN = 0;
            func_handle = @(I_MS_LR, I_MS,I_PAN,ratio) C_GLP_Robust(I_MS, I_PAN, ratio, sensorInf.sensor, nclusters, zeta, thNDVI, flagMTFPAN);
        case 'MTF-GLP-CBD'
            blockSize = [32 32];
            threshold = -Inf;
            func_handle = @(I_MS_LR, I_MS,I_PAN,ratio) MTF_GLP_CBD(I_MS,I_PAN,ratio,blockSize,threshold,sensorInf);
            %func_handle = @(I_MS_LR, I_MS,I_PAN,ratio) MTF_GLP(I_PAN,I_MS,sensorInf,ratio);
            xls_prefix = 'M';
        case 'GS2GLP'
            func_handle = @(I_MS_LR, I_MS,I_PAN,ratio) GS2_GLP(I_MS,I_PAN,ratio,sensorInf.sensor);
            xls_prefix = 'M';
        case 'REG-GLP'
            func_handle = @(I_MS_LR, I_MS,I_PAN,ratio) GLP_Reg_FS(I_MS,I_PAN, sensorInf, ratio, 0);
            xls_prefix = 'M';
        case 'GLP-HS'    
            mth = 'HS-WLS';
            samp = 1;
            func_handle = @(I_MS_LR, I_MS,I_PAN,ratio) pansharp_hs(I_MS,I_PAN,mth,ratio,sensorInf,samp);
            xls_prefix = 'M';
        case 'PRACS'
            func_handle = @(I_MS_LR, I_MS,I_PAN,ratio) PRACS(I_MS,I_PAN,ratio);
            xls_prefix = 'C';
        case 'BDSD'
            func_handle = @(I_MS_LR, I_MS,I_PAN,ratio) BDSD(I_MS,I_PAN,ratio,size(I_MS,1),sensorInf);
            xls_prefix = 'C';
        case 'NonlinearIHS'
            P       = 5;      % patch size for partioning input images
            q       = 3;      % overlap between adjacent patches

            eta     = 1;      % balanced parameter of global synthesis[see equ.(17)]
            maxIter = 10;     % the number of iterations in global synthesis procedure
            mu      = maxIter^-1; % step size of the gradient descent[see equ.(17)]
            gamma   = 10^-9;  % gamma controls the magnitude on the edges[see eq.(5)]
            eps     = 10^-10; % eps enforces a nonzero denominator [see eq.(5)]
            func_handle = @(I_MS_LR, I_MS,I_PAN,ratio) NonlinearIHS(I_PAN,I_MS,P, q, ratio,maxIter, mu, eta, gamma, eps);
            xls_prefix = 'C';
        case 'RR'
            func_handle = @(I_MS_LR, I_MS,I_PAN,ratio) RRWrapper(I_MS_LR,I_MS,I_PAN,sensorInf);
        case 'TV'
            func_handle = @(I_MS_LR, I_MS,I_PAN,ratio) TVWrapper(I_MS_LR,I_MS,I_PAN,sensorInf);
        case 'ST-naive'
            ST_naive_opts.prep.hm_mode = 1;
            ST_naive_opts.prep.I_type = 'H';
            ST_naive_opts.prep.G_type = 'ratio-H';
            ST_naive_opts.fusion.clevels = 3;
            func_handle = @(I_MS_LR, I_MS,I_PAN,ratio) ST_naive_wrapper(I_MS_LR, I_PAN, ratio, sensorInf, L, ST_naive_opts);
        case 'PNN'
            NDxI_flag = false;
            func_handle = @(I_MS_LR, I_MS,I_PAN,ratio) PNN(I_MS_LR, I_PAN, sensorInf, L, NDxI_flag, I_MS);
        case 'PNN-plus'
            func_handle = @(I_MS_LR, I_MS,I_PAN,ratio) PNNplus(I_MS_LR, I_PAN, sensorInf, 0, L, [], I_MS);
        case 'PNN-plus-FT'
            func_handle = @(I_MS_LR, I_MS,I_PAN,ratio) PNNplus(I_MS_LR, I_PAN, sensorInf, 50, L, [], I_MS);
        case 'SR-D'
            %%%%%%%%%%%%%%%%%%%%%%%%%% Parameters setting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            TS = 7; % Tiling (dimensions of the patches are TS x TS)
            ol = 4; % Overlap (in pixels) between contiguous tile
            n_atoms = 10; % Max number of representation atoms (default value = 10)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            func_handle = @(I_MS_LR, I_MS,I_PAN,ratio) CS(I_MS,I_PAN,I_MS_LR,ratio,sensorInf.sensor,TS,ol,n_atoms);
        case 'MMP'
            func_handle = @(I_MS_LR, I_MS,I_PAN,ratio) mmp_wrapper(I_MS_LR,I_PAN,sensorInf);
        otherwise
            warning(['Invalid BP method name: ' func_name '.Set to EXP']);
            func_handle = @(I_MS_LR, I_MS,I_PAN,ratio) I_MS;
    end
