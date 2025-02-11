% A demo for EXP-PnP.
% EXP-PnP: General Extended Model for Detail-injected Pan-sharpening with Plug-and-Play Residual Optimization. 
% https://doi.org/10.1109/TGRS.2025.3539776
% https://github.com/JZ-Tao/EXP-PnP
clc;
clear;
close all;
file_path = matlab.desktop.editor.getActive;
cd(fileparts(file_path.Filename));

n_run_times = 1;
update_output_qi = 1;

% Input data and global parameter preparation

data_name = 'F-WV2_19';
using_estPSF_QI = 0;
load(['Datasets/' data_name '.mat']);
FE_method = 'bruse';
n_band = size(I_MS_LR,3);
[PSF_org, GNyq_org] = Blur_Kernel(n_band, sensorInf.sensor, ratio, 1);
sensorInf.PSF_G = PSF_org;
sensorInf.GNyq = getGNyqBySensor(sensorInf.sensor, n_band);
sensorInf.mGNyq = mean(sensorInf.GNyq);
if is_FS
    [PSF_e, GNyq_e] = FE_Wrapper(I_MS_LR,I_PAN,ratio, sensorInf,FE_method);
    if using_estPSF_QI
        sensorInf.PSF_G = PSF_e;
        sensorInf.mGNyq = mean(GNyq_e);
        sensorInf.GNyq = GNyq_e;
    end
    sensorInf.GNyq_e = GNyq_e;
    sensorInf.mGNyq_e = mean(GNyq_e);
    sensorInf.estPSF_G = PSF_e;
end


% You can uncomment and use the code below to generate the corresponding
% input image for the PAirMax data. Note that the original I_MS image of
% this data is upsampled differently than tap-23, so it needs to be
% regenerated.

% if ~is_FS
% 	FS_str = 'DS';
%     PAirMax_FS_str = 'RR';
% else
% 	FS_str = 'FS';   
%     PAirMax_FS_str = 'FR';
% end
% dataset_idx = 1;
% FE_method = 'bruse';
% 
% if dataset_idx >= 1 && dataset_idx <= 1+8 %% PAirMax\
%     PM_dat_folders = {'GE_Lond_Urb', 'GE_Tren_Urb', 'W2_Miam_Mix', 'W2_Miam_Urb', ...
%         'W3_Muni_Mix', 'W3_Muni_Nat', 'W3_Muni_Urb', 'W4_Mexi_Nat', 'W4_Mexi_Urb'};
%     PM_sensors = {'GeoEye1', 'GeoEye1', 'WV2', 'WV2', 'WV3', 'WV3', 'WV3', 'WV4', 'WV4'};
%     PM_idx = dataset_idx;
%     folder_path = ['image/PAirMax/' PM_dat_folders{PM_idx} '/' PAirMax_FS_str];
%     L_org = 11;
%     im_prepare = 'donothing';
%     sensor = PM_sensors{PM_idx};
%     if ~is_FS
%         I_GT = double(imread([folder_path '/GT.tif']));
%         I_MS_LR = double(imread([folder_path '/MS_LR.tif']));
%     else
%         I_GT = double(imread([folder_path '/MS_LR.tif']));
%         I_MS = double(imread([folder_path '/MS.tif']));
%     end
% 
%     I_PAN = double(imread([folder_path '/PAN.tif']));
%     if size(I_GT,3) == 4
%         color_band_idx = [3 2 1];
%     else
%         color_band_idx = [5 3 1];
%     end
%     offy = 222; offx = 222; % 20 20
%     if ~is_FS 
%         offy = 0; offx = 0;
%         I_PAN = I_PAN(1+offy*ratio:GT_size+offy*ratio, 1+offx*ratio:GT_size+offx*ratio);
%         I_GT = I_GT(1+offy*ratio:GT_size+offy*ratio, 1+offx*ratio:GT_size+offx*ratio,:);
%         I_MS_LR = I_MS_LR(1+offy:GT_size/ratio+offy, 1+offx:GT_size/ratio+offx, :);
%         %I_MS = I_MS(1+offy*ratio:GT_size+offy*ratio, 1+offx*ratio:GT_size+offx*ratio,:);
%         I_MS = interpWrapper(double(I_MS_LR),ratio,sensorInf.upsampling);
%     else
%         sensorInf.sensor = sensor;
%         [I_MS, I_GT, I_PAN] = FSregistration(I_GT, I_PAN, [GT_size, GT_size], [ceil(GT_size/2)+offy, ceil(GT_size/2)+offx], ratio, sensorInf);
%     end
% end
% 
% sensorInf.L = L;
% sensorInf.sensor = sensor;
% [PSF_org, GNyq_org] = Blur_Kernel(n_band, sensor, ratio, 1);
% sensorInf.PSF_G = PSF_org;
% sensorInf.GNyq = getGNyqBySensor(sensor, n_band);
% sensorInf.mGNyq = mean(sensorInf.GNyq);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Resize data
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% if is_FS
%     I_MS_LR = I_GT;
%     I_PAN = I_PAN_GT;
%     [PSF_e, GNyq_e] = FE_Wrapper(I_MS_LR,I_PAN,ratio, sensorInf,FE_method);
%     if using_estPSF_QI
%         sensorInf.PSF_G = PSF_e;
%         sensorInf.mGNyq = mean(GNyq_e);
%         sensorInf.GNyq = GNyq_e;
%     end
%     sensorInf.GNyq_e = GNyq_e;
%     sensorInf.mGNyq_e = mean(GNyq_e);
%     sensorInf.estPSF_G = PSF_e;
% else
%     if strcmp(im_prepare,'resize')
% 
%         I_MS_LR = MTF_conv_sample(I_GT, sensorInf, ratio, 1);
%         
%         if size(I_GT(:,:,1)) ~= size(I_PAN_GT)
%             I_PAN = imresize(I_PAN_GT, 1/ratio);
%         else
%             I_PAN = I_PAN_GT;
%         end
%     end
% end
test_MTF_GLP = 1;
test_MTF_GLP_HPM_H = 1;


n_methods = 0;
% PnP Method name. Add '-BP' at the end of the name for post-processing
% optimization using BP.
EXP_PnP_cases = {'EXP','Wiener','SRCNN','Wiener-BP'};
% Baseline method name. Add '+MB' at the end of the name to enable blur matching.
base_methods = {'GLP','GLP-HPM-H','GLP-HPM-H+MB'};
ADMM_global_lambda = 1e-5;
IRCNN_SR_g_lambda = 0.00006;
IRCNN_deb_g_lambda = 0.06;

[PSF_e, GNyq_e] = FE_Wrapper(I_MS_LR,I_PAN,ratio,sensorInf,FE_method);
sensorInf_e = sensorInf;
sensorInf_e.PSF_G = PSF_e;
sensorInf.PSF_e = PSF_e;

if ~is_FS
    PnP_BP_Opts = init_BP_options("BP_I", is_FS);
else
    PnP_BP_Opts = init_BP_options("BP_T", is_FS);
end
PnP_BP_Opts.n_reg = 5;
FE_method = 'bruse';

n_tests = length(EXP_PnP_cases);

deblur_lambda = 1.0;

[PSF_e, GNyq_e] = FE_Wrapper(I_MS_LR,I_PAN,ratio, sensorInf, FE_method);
global_lambda = ADMM_global_lambda;
sensorInf_exp = sensorInf;
tap = 41;
sensorInf_exp.GNyq_e = deblur_lambda.*GNyq_e+(1-deblur_lambda).*sensorInf.GNyq;
sensorInf_exp.PSF_e = MTF_GNyq2PSF(sensorInf_exp.GNyq_e, tap, ratio);
sensorInf_exp.mGNyq_e = mean(sensorInf_exp.GNyq_e);
for iii = 1:length(base_methods)
    base_method = base_methods{iii};
    for ii = 1:n_tests
        EXP_Opts = init_EXP_PnP_options_by_test_case(EXP_PnP_cases{ii});
        EXP_Opts.dataset_idx = dataset_idx;
        EXP_Opts.I_MS = I_MS;
        EXP_Opts.anti_aliasing = 0;
        EXP_Opts.matching_blur = 0;
        EXP_Opts.color_band_idx = color_band_idx;
        EXP_Opts.deb_PSF_e = 0;
        switch(EXP_Opts.exp_method)
            case {'PnP-ADMM-SR', 'MS-LapSRN', 'DPSR','LapSRN'}
                EXP_Opts.need_interp = 0;
                global_lambda = ADMM_global_lambda;
            case 'IRCNN-deb'
                global_lambda = IRCNN_SR_g_lambda;
            case 'IRCNN-SR'
                global_lambda = IRCNN_deb_g_lambda;
        end
        n_methods = n_methods+1;
        Method_{n_methods} = [base_method '_' EXP_PnP_cases{ii} '_mGNyq_' num2str(sensorInf.mGNyq)];
        t2=tic;
        n_reg = 15;
        EXP_Opts.HR_size = size(I_MS);
        EXP_Opts.L = L;
        EXP_Opts.threshold_on = 1;
        EXP_Opts.global_lambda = global_lambda;

        EXP_Opts.is_FS = is_FS;

        if is_FS
            EXP_Opts.extra_deb = 0;
        else
            EXP_Opts.extra_deb = 1;
        end
        EXP_Opts.color_band_idx = color_band_idx;    
        if EXP_Opts.extra_deb
            I_PAN_LR_bic = imresize(I_PAN, 1/ratio);
            I_PAN_LR = MTF_conv_sample(I_PAN, sensorInf_exp, ratio, 1);

            [PSF_LR, GNyq_LR] = FE_Wrapper(I_PAN_LR, I_PAN_LR_bic, 1, sensorInf_exp, 'bruse_LR');
            EXP_Opts.PSF_LR = PSF_LR;

        end
        for i = 1:n_run_times
            I_t = EXP_PnP_Wrapper(I_MS_LR, I_PAN, sensorInf_exp, ratio, base_method, EXP_Opts);
            if EXP_Opts.BP
                I_t = BP_Wrapper(I_t, I_MS_LR, I_PAN, ratio, sensorInf_exp, PnP_BP_Opts);
            end
            I_F{n_methods} = I_t;
        end
        Time_{n_methods}=toc(t2)/n_run_times;
        fprintf('Elaboration time %s: %.4f [sec]\n',Method_{n_methods},Time_{n_methods});
    end
end


%% MTF-GLP
if test_MTF_GLP
    n_methods = n_methods+1;
    Method_{n_methods} = 'GLP';
    
    t2=tic;
    for i = 1:n_run_times
        I_F{n_methods} = MTF_GLP(I_PAN,I_MS,sensorInf,ratio);
    end
    Time_{n_methods} = toc(t2)/n_run_times;
    fprintf('Elaboration time %s: %.4f [sec]\n',Method_{n_methods},Time_{n_methods});
    
   
end

%% MTF_GLP_HPM_H
if test_MTF_GLP_HPM_H
    n_methods = n_methods+1;
    Method_{n_methods} = 'GLP-HPM-H';
    
    t2=tic;
    for i = 1:n_run_times
        I_F{n_methods} = MTF_GLP_HPM_Haze_min(I_MS,I_PAN,sensorInf.sensor,ratio,1);
    end
    Time_{n_methods} = toc(t2)/n_run_times;
    fprintf('Elaboration time %s: %.4f [sec]\n',Method_{n_methods},Time_{n_methods});
end



dat_time = datestr(now,30);

if is_FS
    I_GT = [];
end

QI_ = cell(1, n_methods);
for i = 1:n_methods
    QI_{i} = indices_eval_EXP_PnP_wrapper(I_F{i}, I_GT, I_MS, I_MS_LR, I_PAN, ratio,...
        L, sensorInf, is_FS, Qblocks_size, flag_cut_bounds, dim_cut,thvalues,using_estPSF_QI);    
end
QI_matrix = QI2matrix(QI_,Method_,Time_,is_FS);
if update_output_qi
    disp( ['Writting xls of Dataset: ' data_name]);
    xlswrite(['Results/qi/' datestr(now,30) '_' data_name], QI_matrix, ['Dataset ' data_name]);
end
