function EXP_Opts = init_EXP_PnP_options_by_test_case(case_name)

% Some default value for common parameters
EXP_Opts.gain_method = 'RS';
%EXP_Opts.exp_method = case_name{ii};%'Wiener'; %'PnP_ADMM';
EXP_Opts.ADMM.subspace = 0;
EXP_Opts.ADMM.denoiser = 'TV';
EXP_Opts.debug.seperate = 0; %%%%%%%%%
EXP_Opts.normalize_type = 1;
EXP_Opts.extra_deb = 0;

EXP_Opts.step_gamma = 1;
EXP_Opts.GT_mode = 1;
EXP_Opts.need_interp = 1;
EXP_Opts.anti_aliasing = 0;
EXP_Opts.spatial_adaption = 0;

if endsWith(case_name, '-BP')
    EXP_Opts.BP = 1;
    case_name = strrep(case_name, '-BP', '');
else
    EXP_Opts.BP = 0;
end

if strcmp(case_name, 'LMG')
    EXP_Opts.exp_method = 'LMG';
    EXP_Opts.type = 'Deblur';
    EXP_Opts.debug.seperate = 0;
    EXP_Opts.need_interp = 1;
    return;
end

if strcmp(case_name, 'GraphBID')
    EXP_Opts.exp_method = 'GraphBID';
    EXP_Opts.type = 'Deblur';
    EXP_Opts.debug.seperate = 0;
    EXP_Opts.need_interp = 1;
    return;
end

if strcmp(case_name, 'Seperate-MS-LapSRN')
    EXP_Opts.exp_method = 'MS-LapSRN';
    EXP_Opts.type = 'SR';
    EXP_Opts.debug.seperate = 1;
    EXP_Opts.need_interp = 0;
    return;
end

if strcmp(case_name, 'IRCNN-SR')
    EXP_Opts.exp_method = 'IRCNN-SR';
    EXP_Opts.type = 'SR';
    EXP_Opts.debug.seperate = 0;
    EXP_Opts.need_interp = 0;
    return;
end

if strcmp(case_name, 'IRCNN-deb')
    EXP_Opts.exp_method = 'IRCNN-deb';
    EXP_Opts.type = 'Deblur';
    EXP_Opts.debug.seperate = 0;
    EXP_Opts.need_interp = 1;
    return;
end

if strcmp(case_name, 'FFDNet')
    EXP_Opts.exp_method = 'FFDNet';
    EXP_Opts.type = 'SR';
    EXP_Opts.debug.seperate = 0;
    EXP_Opts.need_interp = 1;
    return;
end

if strcmp(case_name, 'RealESRGAN')
    EXP_Opts.exp_method = 'RealESRGAN';
    EXP_Opts.type = 'SR';
    EXP_Opts.debug.seperate = 0;
    EXP_Opts.need_interp = 0;
    return;
end

if strcmp(case_name, 'RealESRNet')
    EXP_Opts.exp_method = 'RealESRNet';
    EXP_Opts.type = 'SR';
    EXP_Opts.debug.seperate = 0;
    EXP_Opts.need_interp = 0;
    return;
end

if strcmp(case_name, 'MS-LapSRN')
    EXP_Opts.exp_method = 'MS-LapSRN';
    EXP_Opts.type = 'SR';
    EXP_Opts.debug.seperate = 0;
    EXP_Opts.need_interp = 0;
    return;
end
if strcmp(case_name, 'TTST')
    EXP_Opts.exp_method = 'TTST';
    EXP_Opts.type = 'SR';
    EXP_Opts.debug.seperate = 0;
    EXP_Opts.need_interp = 0;
    return;
end

if strcmp(case_name, 'LapSRN')
    EXP_Opts.exp_method = 'LapSRN';
    EXP_Opts.type = 'SR';
    EXP_Opts.debug.seperate = 0;
    EXP_Opts.need_interp = 0;
    return;
end

if strcmp(case_name, 'NEDI')
    EXP_Opts.exp_method = 'NEDI';
    EXP_Opts.type = 'SR';
    EXP_Opts.debug.seperate = 0;
    EXP_Opts.need_interp = 0;
    return;
end

if strcmp(case_name, 'MS-LapSRN-deb')
    EXP_Opts.exp_method = 'MS-LapSRN';
    EXP_Opts.type = 'SR';
    EXP_Opts.debug.seperate = 0;
    EXP_Opts.extra_deb = 1;
    EXP_Opts.need_interp = 0;
    return;
end

if strcmp(case_name, 'DPSR')
    EXP_Opts.exp_method = 'DPSR';
    EXP_Opts.type = 'SR';
    EXP_Opts.debug.seperate = 0;
    return;
end

if strcmp(case_name, 'A-PLUS')
    EXP_Opts.exp_method = 'A-PLUS';
    EXP_Opts.type = 'SR';
    EXP_Opts.debug.seperate = 0;
    return;
end

if strcmp(case_name, 'FSRCNN')
    EXP_Opts.exp_method = 'FSRCNN';
    EXP_Opts.type = 'SR';
    EXP_Opts.debug.seperate = 0;
    EXP_Opts.need_interp = 0;
    return;
end

if strcmp(case_name, 'SRCNN')
    EXP_Opts.exp_method = 'SRCNN';
    EXP_Opts.type = 'SR';
    EXP_Opts.debug.seperate = 0;
    return;
end

if strcmp(case_name, 'SRCNN-deb')
    EXP_Opts.exp_method = 'SRCNN';
    EXP_Opts.type = 'SR';
    EXP_Opts.debug.seperate = 0;
    EXP_Opts.extra_deb = 1;
    return;
end

if strcmp(case_name, 'SRCNN-deb2')
    EXP_Opts.exp_method = 'SRCNN';
    EXP_Opts.type = 'SR';
    EXP_Opts.debug.seperate = 0;
    EXP_Opts.extra_deb = 2;
    return;
end
if strcmp(case_name, 'PnP-ADMM-SR-TV-BM3D')
    EXP_Opts.exp_method = 'PnP-ADMM-SR';
    EXP_Opts.ADMM.denoiser = {'TV','BM3D'};
    EXP_Opts.type = 'SR';
    EXP_Opts.debug.seperate = 0;
    return;
end

if strcmp(case_name, 'PnP-ADMM-Deblur-IRCNN')
    EXP_Opts.exp_method = 'PnP-ADMM-Deblur';
    EXP_Opts.ADMM.denoiser = 'IRCNN';
    EXP_Opts.type = 'Deblur';
    EXP_Opts.debug.seperate = 0;
    return;
end

if strcmp(case_name, 'PnP-ADMM-SR-IRCNN')
    EXP_Opts.exp_method = 'PnP-ADMM-SR';
    EXP_Opts.ADMM.denoiser = 'IRCNN';
    EXP_Opts.type = 'SR';
    EXP_Opts.debug.seperate = 0;
    
    return;
end


if strcmp(case_name, 'PnP-ADMM-Deblur-FFDNet')
    EXP_Opts.exp_method = 'PnP-ADMM-Deblur';
    EXP_Opts.ADMM.denoiser = 'FFDNet';
    EXP_Opts.type = 'Deblur';
    EXP_Opts.debug.seperate = 0;
    return;
end

if strcmp(case_name, 'PnP-ADMM-SR-FFDNet')
    EXP_Opts.exp_method = 'PnP-ADMM-SR';
    EXP_Opts.ADMM.denoiser = 'FFDNet';
    EXP_Opts.type = 'SR';
    %EXP_Opts.debug.seperate = 0;
    
    return;
end

if strcmp(case_name, 'PnP-ADMM-SR-TV-MCWNNM')
    EXP_Opts.exp_method = 'PnP-ADMM-SR';
    EXP_Opts.ADMM.denoiser = {'TV','MCWNNM'};
    EXP_Opts.type = 'SR';
    EXP_Opts.debug.seperate = 0;
    return;
end
if strcmp(case_name, 'PnP-ADMM-Deblur-TV-BM3D')
    EXP_Opts.exp_method = 'PnP-ADMM-Deblur';
    EXP_Opts.ADMM.denoiser = {'TV','BM3D'};
    EXP_Opts.type = 'Deblur';
    EXP_Opts.debug.seperate = 0;
    return;
end
if strcmp(case_name, 'PnP-ADMM-Deblur-TV-MCWNNM')
    EXP_Opts.exp_method = 'PnP-ADMM-Deblur';
    EXP_Opts.ADMM.denoiser = {'TV','MCWNNM'};
    EXP_Opts.type = 'Deblur';
    EXP_Opts.debug.seperate = 0;
    return;
end
if strcmp(case_name, 'PnP-ADMM-SR-MCWNNM')
    EXP_Opts.exp_method = 'PnP-ADMM-SR';
    EXP_Opts.ADMM.denoiser = 'MCWNNM';
    EXP_Opts.type = 'SR';
    EXP_Opts.debug.seperate = 0;
    return;
end
if strcmp(case_name, 'PnP-ADMM-Deblur-MCWNNM')
    EXP_Opts.exp_method = 'PnP-ADMM-Deblur';
    EXP_Opts.ADMM.denoiser = 'MCWNNM';
    EXP_Opts.type = 'Deblur';
    EXP_Opts.debug.seperate = 0;
    return;
end


if strcmp(case_name, 'PnP-ADMM-SR-WNNM')
    EXP_Opts.exp_method = 'PnP-ADMM-SR';
    EXP_Opts.ADMM.denoiser = 'WNNM';
    EXP_Opts.type = 'SR';
    EXP_Opts.debug.seperate = 0;
    return;
end
if strcmp(case_name, 'PnP-ADMM-Deblur-WNNM')
    EXP_Opts.exp_method = 'PnP-ADMM-Deblur';
    EXP_Opts.ADMM.denoiser = 'WNNM';
    EXP_Opts.type = 'Deblur';
    EXP_Opts.debug.seperate = 0;
    return;
end


if strcmp(case_name, 'Seperate-PnP-ADMM-Deblur-TV')
    EXP_Opts.exp_method = 'PnP-ADMM-Deblur';
    EXP_Opts.ADMM.denoiser = 'TV';
    EXP_Opts.type = 'Deblur';
    EXP_Opts.debug.seperate = 1;
    return;
end
if strcmp(case_name, 'Seperate-PnP-ADMM-SR-TV')
    EXP_Opts.exp_method = 'PnP-ADMM-SR';
    EXP_Opts.ADMM.denoiser = 'TV';
    EXP_Opts.type = 'SR';
    EXP_Opts.debug.seperate = 1;
    return;
end
if strcmp(case_name, 'Seperate-PnP-ADMM-SR-BM3D')
    EXP_Opts.exp_method = 'PnP-ADMM-SR';
    EXP_Opts.ADMM.denoiser = 'BM3D';
    EXP_Opts.type = 'SR';
    EXP_Opts.debug.seperate = 1;
    return;
end
if strcmp(case_name, 'Seperate-Wiener')
    EXP_Opts.exp_method = 'Wiener';
    EXP_Opts.type = 'Deblur';
    EXP_Opts.debug.seperate = 1;
    return;
end
if strcmp(case_name, 'PnP-ADMM-SR-TV')
    EXP_Opts.exp_method = 'PnP-ADMM-SR';
    EXP_Opts.ADMM.denoiser = 'TV';
    EXP_Opts.type = 'SR';
    return;
end

if strcmp(case_name, 'PnP-ADMM-SR-BM3D-0.0015')
    EXP_Opts.exp_method = 'PnP-ADMM-SR';
    EXP_Opts.ADMM.denoiser = 'BM3D';
    EXP_Opts.type = 'SR';
    EXP_Opts.lambda = 0.0015;
    return;
end
if strcmp(case_name, 'PnP-ADMM-SR-BM3D-0.001')
    EXP_Opts.exp_method = 'PnP-ADMM-SR';
    EXP_Opts.ADMM.denoiser = 'BM3D';
    EXP_Opts.type = 'SR';
    EXP_Opts.lambda = 0.001;
    return;
end
if strcmp(case_name, 'PnP-ADMM-SR-BM3D-0.000025')
    EXP_Opts.exp_method = 'PnP-ADMM-SR';
    EXP_Opts.ADMM.denoiser = 'BM3D';
    EXP_Opts.type = 'SR';
    EXP_Opts.lambda = 0.000025;
    return;
end
if strcmp(case_name, 'PnP-ADMM-SR-BM3D-0.0001')
    EXP_Opts.exp_method = 'PnP-ADMM-SR';
    EXP_Opts.ADMM.denoiser = 'BM3D';
    EXP_Opts.type = 'SR';
    EXP_Opts.lambda = 0.0001;
    return;
end
if strcmp(case_name, 'PnP-ADMM-SR-BM3D')
    EXP_Opts.exp_method = 'PnP-ADMM-SR';
    EXP_Opts.ADMM.denoiser = 'BM3D';
    EXP_Opts.type = 'SR';
    return;
end

if strcmp(case_name, 'PnP-ADMM-SR-FFDNET')
    EXP_Opts.exp_method = 'PnP-ADMM-SR';
    EXP_Opts.ADMM.denoiser = 'FFDNET';
    EXP_Opts.type = 'SR';
    return;
end
if strcmp(case_name, 'PnP-ADMM-Deblur-TV')
    EXP_Opts.exp_method = 'PnP-ADMM-Deblur';
    EXP_Opts.ADMM.denoiser = 'TV';
    EXP_Opts.type = 'Deblur';
    return;
end

if strcmp(case_name, 'PnP-ADMM-Deblur-BM3D')
    EXP_Opts.exp_method = 'PnP-ADMM-Deblur';
    EXP_Opts.ADMM.denoiser = 'BM3D';
    EXP_Opts.type = 'Deblur';
    return;
end

if strcmp(case_name, 'PnP-ADMM-Deblur-FFDNET')
    EXP_Opts.exp_method = 'PnP-ADMM-Deblur';
    EXP_Opts.ADMM.denoiser = 'FFDNET';
    EXP_Opts.type = 'Deblur';
    return;
end
if strcmp(case_name, 'Wiener')
    EXP_Opts.exp_method = 'Wiener';
    EXP_Opts.type = 'Deblur';
    return;
end

if strcmp(case_name, 'BM3D-SOS')
    EXP_Opts.exp_method = 'BM3D-SOS';
    EXP_Opts.type = 'Deblur';
    EXP_Opts.sigma = 0.005;
    return;
end
if strcmp(case_name, 'BM3D-SOS-0.1')
    EXP_Opts.exp_method = 'BM3D-SOS';
    EXP_Opts.type = 'Deblur';
    EXP_Opts.sigma = 0.1;
    return;
end
if strcmp(case_name, 'BM3D-SOS-0.01')
    EXP_Opts.exp_method = 'BM3D-SOS';
    EXP_Opts.type = 'Deblur';
    EXP_Opts.sigma = 0.01;
    return;
end
if strcmp(case_name, 'BM3D-SOS-0.001')
    EXP_Opts.exp_method = 'BM3D-SOS';
    EXP_Opts.type = 'Deblur';
    EXP_Opts.sigma = 0.001;
    return;
end
if strcmp(case_name, 'BM3D-SOS-0.0001')
    EXP_Opts.exp_method = 'BM3D-SOS';
    EXP_Opts.type = 'Deblur';
    EXP_Opts.sigma = 0.0001;
    return;
end

if strcmp(case_name, 'BM3D')
    EXP_Opts.exp_method = 'BM3D';
    EXP_Opts.type = 'Deblur';
    EXP_Opts.sigma = 0.005;
    return;
end
if strcmp(case_name, 'BM3D-0.1')
    EXP_Opts.exp_method = 'BM3D';
    EXP_Opts.type = 'Deblur';
    EXP_Opts.sigma = 0.1;
    return;
end
if strcmp(case_name, 'BM3D-0.01')
    EXP_Opts.exp_method = 'BM3D';
    EXP_Opts.type = 'Deblur';
    EXP_Opts.sigma = 0.01;
    return;
end
if strcmp(case_name, 'BM3D-0.001')
    EXP_Opts.exp_method = 'BM3D';
    EXP_Opts.type = 'Deblur';
    EXP_Opts.sigma = 0.001;
    return;
end
if strcmp(case_name, 'BM3D-0.0001')
    EXP_Opts.exp_method = 'BM3D';
    EXP_Opts.type = 'Deblur';
    EXP_Opts.sigma = 0.0001;
    return;
end
if strcmp(case_name, 'BM3D-0.00005')
    EXP_Opts.exp_method = 'BM3D';
    EXP_Opts.type = 'Deblur';
    EXP_Opts.sigma = 0.00005;
    return;
end

if strcmp(case_name, 'Hyper-Laplacian')
    EXP_Opts.exp_method = 'Hyper-Laplacian';
    EXP_Opts.type = 'Deblur';
    return;
end

if strcmp(case_name, 'DWDN')
    EXP_Opts.exp_method = 'DWDN';
    EXP_Opts.type = 'Deblur';
    
    return;
end

if strcmp(case_name, 'GF1') || strcmp(case_name, 'GF2')
    EXP_Opts.exp_method = case_name;
    EXP_Opts.type = 'SR';
    EXP_Opts.need_interp = 1;
    return;
end
if strcmp(case_name, 'EXP')
    EXP_Opts.exp_method = 'EXP';
    EXP_Opts.type = 'SR';
    
    EXP_Opts.need_interp = 1;
    return;
end

error(['Invalid case name:' case_name]);