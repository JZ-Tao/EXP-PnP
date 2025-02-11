function [G,Gt] = defGGt_mod(K,sensorInf)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Operators for super-resolution
% Stanley Chan
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G  = @(x) fdown(x,K,sensorInf);
Gt = @(x) upf(x,K,sensorInf);
end

function y = fdown(x,K,sensorInf)
% tmp = imfilter(x,h,'circular');
% y = imresize(tmp, 1/K, 'nearest');
%y = downsample2(tmp,K);
y = MTF_conv_sample(x, sensorInf, K, 1);
end

function y = upf(x,K,sensorInf)
y = interpWrapper(x, K, sensorInf.upsampling);
% ratio = K;
% [r c b] = size(x);
% I1LRU = zeros(ratio.*r, ratio.*c, b);    
% I1LRU(floor(ratio/2)+1:ratio:end,floor(ratio/2)+1:ratio:end,:) = I_Interpolated;
% 
% y = interp23tap(x,K);
% %tmp = upsample2(x,K);
% tmp = imresize(x, K, 'nearest');
% y = imfilter(tmp,h,'circular');

end
