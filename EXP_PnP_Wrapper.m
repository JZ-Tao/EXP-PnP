% EXP-PnP: General Extended Model for Detail-injected Pan-sharpening with Plug-and-Play Residual Optimization. 
% https://doi.org/10.1109/TGRS.2025.3539776
% https://github.com/JZ-Tao/EXP-PnP
function I_t = EXP_PnP_Wrapper(I_MS_LR, I_PAN, sensorInf, ratio, base_method, EXP_Opts)
if endsWith(base_method, '+MB')
    EXP_Opts.matching_blur = 1;
    base_method = strrep(base_method, '+MB', '');
end
if endsWith(base_method, '+AA')
    EXP_Opts.anti_aliasing = 1;
    base_method = strrep(base_method, '+AA', '');
end
switch(base_method)
    % Additive model
    case {'ATWT','ATWT-Dec', 'GLP','GLP-REG','GLP-HS','GLP-CBD'} 
        I_t = EXP_PnP_ADD_General(I_MS_LR, I_PAN, sensorInf, ratio, base_method, EXP_Opts);
    % Mmultiplicative model
    case {'GLP-HPM','GLP-HPM2', 'GLP-HPM-DS','GLP-HPM-H', 'GLP-HPM-R', 'SFIM','MF','FE-HPM'}
        I_t = EXP_PnP_MUL_General(I_MS_LR, I_PAN, sensorInf, ratio, base_method, EXP_Opts);
    % Mixed model
    case {'AWLP', 'AWLP-Dec','AWLP-H','AWLP-H-noDec'}
        I_t = EXP_PnP_MIX(I_MS_LR, I_PAN, sensorInf, ratio, base_method, EXP_Opts);
    case 'BDSD-PC'
        I_t = CS_EXP_BDSD_PC(EXP_Opts.I_MS,I_MS_LR,I_PAN,ratio,sensorInf,EXP_Opts);
    otherwise 
        base_func = genFuncHandle(base_method, sensorInf, ratio, EXP_Opts.L);
        EXP_Opts.base_func = base_func;
        I_t = EXP_PnP_CS_general(I_MS_LR,I_PAN,ratio,sensorInf,EXP_Opts);
end
if EXP_Opts.threshold_on
    max_v = 2^(EXP_Opts.L)-1;
    I_t(I_t<0) = 0; 
    I_t(I_t>max_v) = max_v; 
end