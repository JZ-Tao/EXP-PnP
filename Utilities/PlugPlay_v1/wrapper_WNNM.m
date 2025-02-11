function out = wrapper_WNNM(in_C,sigma)
nSig = NoiseLevel(in_C*255);
nSig = sigma*255;
Par   = ParSet(nSig);
out = zeros(size(in_C));
for i = 1:size(in_C,3)
    in = in_C(:,:,i);
    t = WNNM_DeNoising( in, in, Par );
    out(:,:,i) = t./255;
end
end