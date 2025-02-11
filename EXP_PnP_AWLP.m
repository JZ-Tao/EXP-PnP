%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: 
%           AWLP fuses the upsampled MultiSpectral (MS) and PANchromatic (PAN) images by 
%           exploiting the Additive Wavelet Luminance Proportional (AWLP) algorithm.
% 
% Interface:
%           I_Fus_AWLP = AWLP(I_MS,I_PAN,ratio)
%
% Inputs:
%           I_MS:           MS image upsampled at PAN scale;
%           I_PAN:          PAN image;
%           ratio:          Scale ratio between MS and PAN. Pre-condition: Integer value.
%
% Outputs:
%           I_Fus_AWLP:     AWLP pasharpened image.
% 
% References:
%           [Otazu05]       X. Otazu, M. Gonz碼lez-Aud?cana, O. Fors, and J. N磚榥ez, 揑ntroduction of sensor spectral response into image fusion methods.
%                           Application to wavelet-based methods,?IEEE Transactions on Geoscience and Remote Sensing, vol. 43, no. 10, pp. 2376?385,
%                           October 2005.
%           [Alparone07]    L. Alparone, L. Wald, J. Chanussot, C. Thomas, P. Gamba, and L. M. Bruce, 揅omparison of pansharpening algorithms: Outcome
%                           of the 2006 GRS-S Data Fusion Contest,?IEEE Transactions on Geoscience and Remote Sensing, vol. 45, no. 10, pp. 3012?021,
%                           October 2007.
%           [Vivone14]      G. Vivone, L. Alparone, J. Chanussot, M. Dalla Mura, A. Garzelli, G. Licciardi, R. Restaino, and L. Wald, 揂 Critical Comparison Among Pansharpening Algorithms? 
%                           IEEE Transaction on Geoscience and Remote Sensing, 2014. (Accepted)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function I_Fus_AWLP = EXP_PnP_AWLP(I_MS_LR,I_MS,I_PAN,ratio,sensorInf,EXP_Opts)

[Height,Width,Bands]=size(I_MS);
I_Fus_AWLP=zeros(Height,Width,Bands,'double');



SumImage=sum(I_MS,3)/Bands;

IntensityRatio = zeros(size(I_MS));
for i=1:Bands
    IntensityRatio(:,:,i)=I_MS(:,:,i)./(SumImage+eps);
end

I_PAN = repmat(I_PAN,[1 1 size(I_MS,3)]);

for ii = 1 : size(I_MS,3)    
  I_PAN(:,:,ii) = (I_PAN(:,:,ii) - mean2(I_PAN(:,:,ii))).*(std2(I_MS(:,:,ii))./std2(I_PAN(:,:,ii))) + mean2(I_MS(:,:,ii));  
end

h=[1 4 6 4 1 ]/16;
g=[0 0 1 0 0 ]-h;
htilde=[ 1 4 6 4 1]/16;
gtilde=[ 0 0 1 0 0 ]+htilde;
h=sqrt(2)*h;
g=sqrt(2)*g;
htilde=sqrt(2)*htilde;
gtilde=sqrt(2)*gtilde;
WF={h,g,htilde,gtilde};

Levels = ceil(log2(ratio));
PAN_LP = zeros(size(I_PAN));
for i=1:Bands
    
    WT = ndwt2_working(I_PAN(:,:,i),Levels,WF);
        
    for ii = 2 : numel(WT.dec), WT.dec{ii} = zeros(size(WT.dec{ii})); end
    PAN_LP(:,:,i) = indwt2_working(WT,'c');
end

% for i=1:Bands
%     
% end
I_PAN = IntensityRatio.*I_PAN;
PAN_LP = IntensityRatio.*PAN_LP;
if strcmp(EXP_Opts.exp_method, 'EXP')
    Res_HR = I_MS - PAN_LP;
elseif strcmp(EXP_Opts.type, 'SR')
    PAN_LR = downsampleWrapper(PAN_LP, ratio, sensorInf.downsampling);
    
    Res = I_MS_LR - PAN_LR;
    Res_HR = Res_EXP(Res, ratio, sensorInf, EXP_Opts);
%     if EXP_Opts.anti_aliasing
%         I_MS_deb = I_MS;
%     else
%         I_MS_deb = Res_EXP(I_MS_LR, ratio, sensorInf, EXP_Opts);
%     end
%     PAN_LP_deb = Res_EXP(PAN_LR, ratio, sensorInf, EXP_Opts);
else
    
    EXP_Opts.need_interp = 0;
    Res = I_MS - PAN_LP;
    Res_HR = Res_EXP(Res, ratio, sensorInf, EXP_Opts);

%     if EXP_Opts.anti_aliasing
%         I_MS_deb = I_MS;
%     else
%         I_MS_deb = Res_EXP(I_MS, ratio, sensorInf, EXP_Opts);
%     end
%     PAN_LP_deb = Res_EXP(PAN_LP, ratio, sensorInf, EXP_Opts);
end

I_Fus_AWLP = Res_HR + I_PAN;


end