function PSF_G = MTF_GNyq2PSF(GNyq, tap, ratio)
if nargin < 2
   tap = 41;
end
if nargin < 3
    ratio = 4;
end
n_band = length(GNyq);

fcut = 1/ratio;
PSF_G = zeros(tap,tap,n_band);
for i = 1:n_band
    alpha = sqrt((tap*(fcut/2))^2/(-2*log(GNyq(i))));
    H = fspecial('gaussian', tap, alpha);
    Hd = H./max(H(:));
    h = fwind1(Hd,kaiser(tap));
    PSF_G(:,:,i) = real(h);
end
