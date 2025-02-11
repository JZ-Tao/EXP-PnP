function out = wrapper_median(in,sigma)
out = medfilt2(in, [sigma, sigma]);

end