function [I_MS_adjusted, I_PAN_adjusted] = matching_blur_level(I_MS, I_PAN, GNyq_MS, GNyq_PAN, sensorInf, ratio, tap)
% Default parameter settings
if nargin < 7
    tap = 25;    % PSF size control
end
% Adjustments are made to the PAN only, with as few secondary changes to
% the MS as possible.
I_MS_adjusted = I_MS;
lambda = 1.0;
if GNyq_MS > GNyq_PAN
    GNyq_relative = (GNyq_PAN/GNyq_MS)*lambda;
else
    GNyq_relative = (GNyq_MS/GNyq_PAN)*lambda;
end
    
%[PSF_relative, GNyq_relative] = compute_relative_PSF_by_GNyq(GNyq_MS, GNyq_PAN, ratio, tap);
if GNyq_relative > 0.9 % Small variance, no adjustment
    I_PAN_adjusted = I_PAN;
else
    PSF_relative = MTF_GNyq2PSF(GNyq_relative, tap, ratio);
    if GNyq_PAN < GNyq_MS % Most sensors, e.g. WV series
        I_PAN_adjusted = MSDeblur(I_PAN, PSF_relative, 0, 'Hyper-Laplacian');
    %     EXP_Opts.need_interp = 0;
    %     sensorInf_P = sensorInf;
    %     sensorInf_P.PSF_G = PSF_relative;
    %     sensorInf_P.PSF_e = PSF_relative;
    %     I_PAN_adjusted = Res_EXP(I_PAN, 1, sensorInf_P, EXP_Opts);
        %I_PAN_adjusted = ringing_artifacts_removal(I_PAN, PSF_relative, 2e-3, 1e-3, 0);
    else % GF-2
        sensorInf.PSF_G = PSF_relative;
        I_PAN_adjusted = MTF_conv_sample(I_PAN,sensorInf,1,0);
    end
end
