function I_filtered = idfilter23tap(I,ratio)

L = 22; % (L+1) is the tap
BaseCoeff = fir1(L,1./ratio);
t = imfilter(I',BaseCoeff,'circular'); 
I_filtered = imfilter(t',BaseCoeff,'circular'); 

end