function GGt = constructGGt_mod2(K, sensorInf, rows,cols, n_band, usingHt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A modified code by Jingzhe Tao from:
% Eigen-decomposition for super-resolution
% Stanley Chan
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%hth = conv2(h2',h1);
%hth = conv2(h1, rot90(h1,2));
h1 = sensorInf.PSF_G;
if exist('n_band','var')
    h1 = h1(:,:,1:n_band);
end
if ~exist('usingHt','var')
    usingHt = 0;
end

if ~usingHt
    h2 = getInterpKernel(K, sensorInf.upsampling.interp_type, 2, sensorInf.upsampling.tap);
     hth = convn(h2,h1);
else
    for i=1:size(h1,3)
        hth(:,:,i) = conv2(h1(:,:,i), rot90(h1(:,:,i),2)); %%%%%
    end
end
    yc = ceil(size(hth,1)/2);  % mark the center coordinate
    xc = ceil(size(hth,2)/2);

    L = floor(size(hth,1)/K);  % width of the new filter 
                               %  = (1/k) with of the original filter
    
    g = zeros(L,L,size(hth,3));            % initialize new filter
    for i=-floor(L/2):floor(L/2)
        for j=-floor(L/2):floor(L/2)
            g(i+floor(L/2)+1, j+floor(L/2)+1, :) = hth(yc+K*i, xc+K*j, :);
        end
    end

    GGt = abs(fft2(g,rows/K,cols/K));
end