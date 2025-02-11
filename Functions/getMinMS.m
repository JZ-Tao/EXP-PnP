function minMS = getMinMS(I_MS)

if size(I_MS,3) == 4
    prc = 1;
    minMS = zeros(1,1,4);
    B = I_MS(:,:,1);
    G = I_MS(:,:,2);
    R = I_MS(:,:,3);
    NIR = I_MS(:,:,4);
    minMS(1,1,1) = 0.95 * prctile(B(:),prc);
    minMS(1,1,2) = 0.45 * prctile(G(:),prc);
    minMS(1,1,3) = 0.40 * prctile(R(:),prc);
    minMS(1,1,4) = 0.05 * prctile(NIR(:),prc);
else
    minMS = zeros(1,1,size(I_MS,3));
    for ii = 1 : size(I_MS, 3)
       minMS(1,1,ii) = min(min(I_MS(:,:,ii)));  
    end
end