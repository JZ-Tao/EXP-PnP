function [I_MS_out,I_MS_LR_out, I_PAN_out] = FSregistration(I_MS_LR, I_PAN, size_MS, center_MS, ratio, sensorInf)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modified from:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   FULL RESOLUTION REGISTRATION TOOLBOX 
% 
% Please, refer to the following paper:
% G. Vivone, M. Dalla Mura, A. Garzelli, and F. Pacifici, "A Benchmarking Protocol for Pansharpening: Dataset, Pre-processing, and Quality Assessment", 
% IEEE Journal of Selected Topics in Applied Earth Observations and Remote Sensing, 2021.
% 
% % % % % % % % % % % % % 
% 
% Version: 1
% 
% % % % % % % % % % % % % 
% 
% Copyright (C) 2021
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('ratio', 'var')
    ratio = 4;
end
if ~exist('size_MS', 'var')
    size_MS = [80, 80];
end
if ~exist('center_MS', 'var')
    size_MS = [40+20, 40+20];
end
center_r_MS = center_MS(1);
center_c_MS = center_MS(2);
Wr = size_MS(1);
Wc = size_MS(2);
%% Cut area

%%% Extract MS area
I_MS_LR_little = I_MS_LR(center_r_MS-round(Wr/2)+1:center_r_MS+round(Wr/2),center_c_MS-round(Wc/2)+1:center_c_MS+round(Wc/2),:);
size(I_MS_LR_little)

%%% Calculate equivalent PAN area
I_PAN_little = I_PAN(center_r_MS*ratio-round(Wr*ratio/2)+3:center_r_MS*ratio+round(Wr*ratio/2)+2,center_c_MS*ratio-round(Wc*ratio/2)+3:center_c_MS*ratio+round(Wc*ratio/2)+2,:);
size(I_PAN_little)

%% Interpolation MS

I_MS = interp23tap(I_MS_LR_little,ratio);

I_PAN_little_LP = MTF_conv_sample(I_PAN_little, sensorInf, ratio, 2);

usfac = 20;
%%% Check sub-pixel registration between MS and PAN 
%output = round(dftregistration(fft2(mean(I_MS,3)),fft2(I_PAN_little),100));
%output = dftregistration(fft2(mean(I_MS,3)),fft2(I_PAN_little),100);

%%% Cut the PAN image to align it to the MS image (alignment to the 2.5/3 pixel) 
%I_PAN_little = I_PAN(center_r_MS*ratio-round(Wr*ratio/2)+3 - output(3):center_r_MS*ratio+round(Wr*ratio/2)+2 - output(3),center_c_MS*ratio-round(Wc*ratio/2)+3- output(4):center_c_MS*ratio+round(Wc*ratio/2)+2- output(4),:);

% vect_tag_interp = {'e_e', 'o_o', 'e_o', 'o_e'};
% tap = 44; % tap filter
% misals = zeros(2,numel(vect_tag_interp));
% usfac = 50;
% for ii = 1 : numel(vect_tag_interp)
%     %%% Interpolator odd or even on rows/columns checking the misalignments
%     tag_interp = vect_tag_interp{ii};
% 
%     %%% Interpolation
%     I_MS = interpGeneral(I_MS_LR_little,ratio,tap,tag_interp,1,1);
%     %%% Check sub-pixel registration between MS and PAN 
%     output = dftregistration(fft2(mean(I_MS,3)),fft2(I_PAN_little),usfac);
%     misals(:,ii) = output(3:4);
% end
% % Select the best combination of interpolators (row/column)
% [~,indmin] = min(mean(abs(misals),1));
% tag_interp = vect_tag_interp{indmin};
% %%% Interpolation
% I_MS = interpGeneral(I_MS_LR_little,ratio,tap,tag_interp,1,1);
% I_MS_out = I_MS;
usfac = 50;
[output, I_PAN_temp] = dftregistration(fft2(mean(I_MS,3)),fft2(I_PAN_little),usfac);
output(3:4)
%% Final products
I_MS_LR_out = I_MS_LR_little;
I_PAN_out = abs(ifft2(I_PAN_temp));

%%% Interpolation

I_MS_out = interp23tap(I_MS_LR_out,ratio);

% [PSF_e, GNyq_e] = FE_Wrapper(I_MS_LR_out,I_PAN_out,ratio, sensorInf,'bruse');
% I_MS_de = MSDeblur(I_MS_out, PSF_e);
% usfac = 50;
% [output, I_PAN_temp] = dftregistration(fft2(mean(I_MS_de,3)),fft2(I_PAN_little),usfac);
% output(3:4)
% I_PAN_out = abs(ifft2(I_PAN_temp));
%% Measure the final misalignments
% output = dftregistration(fft2(mean(I_MS_out,3)),fft2(I_PAN_out),100);
% output(3:4)
