
function I_Fus = EXP_PnP_CS_general(I_MS_LR,I_PAN,ratio,sensorInf, EXP_Opts)
I_MS_en = Res_EXP(I_MS_LR, ratio, sensorInf, EXP_Opts); 
I_Fus = EXP_Opts.base_func(I_MS_LR, I_MS_en, I_PAN, ratio);
end