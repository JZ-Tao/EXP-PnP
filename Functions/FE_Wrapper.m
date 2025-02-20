function [PSF,GNyq] = FE_Wrapper(I_MS_LR, I_PAN, ratio, sensorInf, type)
if ~exist('type','var')
   type = 'MBFE';
end
tap = 25;
lambda = 10^5;
mu = 10^5;
num_iter = 10;
th = 0;
filtername = 'other';

if size(I_MS_LR, 1) ~= size(I_PAN, 1)
    I_MS = interpWrapper(I_MS_LR,ratio,sensorInf.upsampling);

else
    I_MS = I_MS_LR;
end

switch type
    case 'None'
        PSF = sensorInf.PSF_G;
        GNyq = getGNyqBySensor(sensorInf.sensor, size(I_MS_LR,3));
    case 'FE'
        PSF = FE(I_MS,I_PAN,ratio,tap,lambda,mu,th,num_iter,filtername);
        GNyq = MTF_GNyq_Est(PSF, ratio);
    case 'FE_MS'
        if size(I_PAN,3) == 1
            I_PAN = repmat(I_PAN, [1,1,size(I_MS,3)]);
        end
        % Hist match?
        for ii = 1:size(I_MS,3)
            I_PAN(:,:,ii) = (I_PAN(:,:,ii) - mean2(LPfilterGauss(I_PAN(:,:,ii),ratio)))...
                .*(std2(I_MS(:,:,ii))./std2(LPfilterGauss(I_PAN,ratio))) + mean2(I_MS(:,:,ii));  
        end
        PSF = FE_MS(I_MS,I_PAN,tap,lambda,mu,th,filtername);
        GNyq = MTF_GNyq_Est(PSF, ratio);
    case 'MBFE_H'
        %I_MS_hat = BDSD_PC(I_MS,I_PAN,ratio,sensorInf);
        %I_MS_hat = GSA(I_MS,I_PAN,I_MS_LR,ratio);
        I_MS_hat1 = MTF_GLP_HPM_Haze_min(I_MS,I_PAN,sensorInf.sensor,ratio,1);
        I_MS_hat2 = BroveyRegHazeMin(I_MS,I_PAN,ratio);
        I_MS_hat3 = MTF_GLP_AWLP_Haze(I_PAN,I_MS,ratio,1,sensorInf);
        PSF1 = FE_MS(I_MS,I_MS_hat1,tap,lambda,mu,th,filtername);
        PSF2 = FE_MS(I_MS,I_MS_hat2,tap,lambda,mu,th,filtername);
        PSF3 = FE_MS(I_MS,I_MS_hat3,tap,lambda,mu,th,filtername);
        PSF = (PSF1 + PSF2 + PSF3)./3;
        GNyq1 = MTF_GNyq_Est(PSF1, ratio);
        GNyq2 = MTF_GNyq_Est(PSF2, ratio);
        GNyq3 = MTF_GNyq_Est(PSF3, ratio);
        GNyq = (GNyq1 + GNyq2 + GNyq3)./3;
    case 'MBFE'
        %I_MS_hat = BDSD_PC(I_MS,I_PAN,ratio,sensorInf);
        %I_MS_hat = GSA(I_MS,I_PAN,I_MS_LR,ratio);
        %%%%%%%%%%%%%%%%%%%%%%%%%% Parameters setting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        TS = 7; % Tiling (dimensions of the patches are TS x TS)
        ol = 4; % Overlap (in pixels) between contiguous tile
        n_atoms = 10; % Max number of representation atoms (default value = 10)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         FT = 50;%%%%%%% 50:PNN_plus_FT
%         tmp = PNNplus(I_MS_LR, I_PAN, sensorInf, FT, sensorInf.L, [], I_MS);
%         if ~isempty(tmp)
%             I_MS_hat = tmp;
%         else
%             error('Invalid sensor model for PNN+. Skipping');
%         end
        I_MS_hat = BroveyRegHazeMin(I_MS,I_PAN,ratio);
        
        %I_MS_hat = CS(I_MS,I_PAN,I_MS_LR,ratio,sensorInf.sensor,TS,ol,n_atoms);
        
        PSF = FE_MS(I_MS,I_MS_hat,tap,lambda,mu,th,filtername);
        GNyq = MTF_GNyq_Est(PSF, ratio);
    case 'MBFE_BDSD-PC'
        I_MS_hat = BDSD_PC(I_MS,I_PAN,ratio,sensorInf);
        PSF = FE_MS(I_MS,I_MS_hat,tap,lambda,mu,th,filtername);
        GNyq = MTF_GNyq_Est(PSF, ratio);
    case 'MBFE_GSA'
        I_MS_hat = GSA(I_MS,I_PAN,I_MS_LR,ratio);
        PSF = FE_MS(I_MS,I_MS_hat,tap,lambda,mu,th,filtername);
        GNyq = MTF_GNyq_Est(PSF, ratio);
    case 'MBFE_APNN'
        FT = 50;%%%%%%% 50:PNN_plus_FT
        tmp = PNNplus(I_MS_LR, I_PAN, sensorInf, FT, sensorInf.L, [], I_MS);
        if ~isempty(tmp)
            I_MS_hat = tmp;
        else
            error('Invalid sensor model for PNN+. Skipping');
        end
        PSF = FE_MS(I_MS,double(I_MS_hat),tap,lambda,mu,th,filtername);
        GNyq = MTF_GNyq_Est(PSF, ratio);
    case 'MBFE_SR-D'
        TS = 7; % Tiling (dimensions of the patches are TS x TS)
        ol = 4; % Overlap (in pixels) between contiguous tile
        n_atoms = 10; % Max number of representation atoms (default value = 10)
        I_MS_hat = CS(I_MS,I_PAN,I_MS_LR,ratio,sensorInf.sensor,TS,ol,n_atoms);
        PSF = FE_MS(I_MS,I_MS_hat,tap,lambda,mu,th,filtername);
        GNyq = MTF_GNyq_Est(PSF, ratio); 
    case 'MBFE_BT-H'
        I_MS_hat = BroveyRegHazeMin(I_MS,I_PAN,ratio);
        PSF = FE_MS(I_MS,I_MS_hat,tap,lambda,mu,th,filtername);
        GNyq = MTF_GNyq_Est(PSF, ratio);
    case 'bruse'
        off = 20;
        I_MS = I_MS(1+off:end-off, 1+off:end-off,:);
        I_PAN = I_PAN(1+off:end-off, 1+off:end-off,:);
        I_MS_LR = I_MS_LR(1+off/ratio:end-off/ratio, 1+off/ratio:end-off/ratio,:);
        
        n_band = size(I_MS,3);
        GNyq = zeros(1,n_band);
        PSF = zeros(tap, tap, n_band);
        min_err = ones(1,n_band)*Inf;
        sensorInf_P = sensorInf;
        Lap_MS_LR = zeros(size(I_MS_LR));
        for i = 1:size(I_MS_LR,3)
            Lap_MS_LR(:,:,i)= imfilter(I_MS_LR(:,:,i),fspecial('sobel'));
        end
        

        %imageHR = GSA(I_MS,I_PAN,I_MS_LR,ratio);% spectral refine
        imageHR = BroveyRegHazeMin(I_MS,I_PAN,ratio);

        for iGNyq = 0.1:0.01:0.9
            iPSF = MTF_GNyq2PSF(iGNyq*ones(1, size(I_MS_LR,3)), tap, ratio);
            sensorInf_P.PSF_G = iPSF;

            imageHR_LR = MTF_conv_sample(imageHR, sensorInf_P, ratio, 1);

            Lap_imageHR_LR = imfilter(imageHR_LR,fspecial('sobel'));
            Lap_ERR = Lap_MS_LR - Lap_imageHR_LR;

            
            for ii = 1:n_band    
                t = norm(Lap_ERR(:,:,ii), 'fro');
                if t < min_err(ii)
                    min_err(ii) = t;
                    GNyq(ii) = iGNyq;
                    PSF(:,:,ii) = iPSF(:,:,1);
                end
            end
        end
    case 'bruse_P2P'
        off = 10;
        
        I_PAN = I_PAN(1+off:end-off, 1+off:end-off,1);
        I_PAN_LR0 = I_MS_LR(1+off/ratio:end-off/ratio, 1+off/ratio:end-off/ratio,1);
        min_err = Inf;
        sensorInf_P = sensorInf;

        Lap_I_PAN_LR0 = imfilter(I_PAN_LR0,fspecial('sobel'));
        for iGNyq = 0.1:0.01:0.9
            iPSF = MTF_GNyq2PSF(iGNyq, tap, ratio);
            sensorInf_P.PSF_G = iPSF;
            I_PAN_LR = MTF_conv_sample(I_PAN, sensorInf_P, ratio, 1);
            
            Lap_I_PAN_LR = imfilter(I_PAN_LR,fspecial('sobel'));
            Lap_ERR = Lap_I_PAN_LR0 - Lap_I_PAN_LR;
    
            t = norm(Lap_ERR, 'fro');
            if t < min_err
                min_err = t;
                GNyq = iGNyq;
                PSF = iPSF;
            end
            
        end
    case 'FE_LR'
        PSF = FE_LR(I_MS, I_PAN,20,tap,lambda,mu,th,filtername);
        GNyq = MTF_GNyq_Est(PSF, 1);
    case 'bruse_LR'
        n_band = size(I_MS,3);
        GNyq = zeros(1,n_band);
        PSF = zeros(tap, tap, n_band);
        min_err = ones(1,n_band)*Inf;
        sensorInf_P = sensorInf;
        off = 10;
        I_MS_LR = I_MS_LR(1+off:end-off, 1+off:end-off,:);
        I_PAN = I_PAN(1+off:end-off, 1+off:end-off,:);
        Lap_MS_LR = zeros(size(I_MS_LR));
        
        for i = 1:size(I_MS_LR,3)
            Lap_MS_LR(:,:,i)= imfilter(I_MS_LR(:,:,i),fspecial('sobel'));
        end
        
        for iGNyq = 0.1:0.01:0.9
            iPSF = MTF_GNyq2PSF(iGNyq*ones(1, size(I_MS_LR,3)), tap, 1);
            sensorInf_P.PSF_G = iPSF;

            imageHR_LP = MTF_conv_sample(I_PAN, sensorInf_P, 1, 0);
            Lap_imageHR_LP = imfilter(imageHR_LP,fspecial('sobel'));
            Lap_ERR = Lap_MS_LR - Lap_imageHR_LP;
            for ii = 1:n_band    
                t = norm(Lap_ERR(:,:,ii), 'fro');
                if t < min_err(ii)
                    min_err(ii) = t;
                    GNyq(ii) = iGNyq;
                    PSF(:,:,ii) = iPSF(:,:,1);
                end
            end
        end        
    case 'Hysure'
        intersection = cell(1,1);
        intersection{1} = 1:size(I_MS,3);
        p = size(I_MS,3);
        contiguous = intersection;
        lambda_R = 1e1;
        lambda_B = 1e1;
        hsize_h = 2*ratio-1;
        hsize_w = 2*ratio-1;
        shift = sensorInf.downsampling.offset; % 'phase' parameter in MATLAB's 'upsample' function
        blur_center = mod(ratio+1,2); % to center the blur kernel according to the simluated data
        [~, ~, B_est] = sen_resp_est(I_MS_LR, I_PAN, ratio, intersection, contiguous, ...
            p, lambda_R, lambda_B, hsize_h, hsize_w, shift, blur_center);
        PSF = MatrixToKernel(B_est,hsize_h,hsize_w);
        
        GNyq = MTF_GNyq_Est(PSF, ratio);
        PSF = repmat(PSF, [1,1,size(I_MS_LR,3)]);
        GNyq = repmat(GNyq, [1,1,size(I_MS_LR,3)]);

end