
function I_Fus_BDSD = CS_EXP_BDSD_PC(I_MS,I_MS_LR,I_PAN,ratio,sensorInf, EXP_Opts)

I_PAN = double(I_PAN);

opts1 = optimset('display','off');

I_GT = imresize(I_MS,1/ratio);%,'nearest');
I_GT_LP = MTF_conv_sample(I_GT,sensorInf,ratio,0);
I_PAN_LR = downsampleWrapper(MTF_PAN(I_PAN,sensorInf,ratio), ratio, sensorInf.downsampling);%imresize(MTF_PAN(I_PAN,sensorInf,ratio),1/ratio,'nearest');

I_MS_en = Res_EXP(I_GT, ratio, sensorInf, EXP_Opts);
I_Fus_BDSD = zeros(size(I_MS));
gamma = zeros(size(I_MS,3)+1,size(I_MS,3));
for ii = 1 : size(I_MS,3)
    h1 = I_GT(:,:,ii);
    h2 = I_GT_LP(:,:,ii);
    H = [I_PAN_LR(:), reshape(I_GT_LP,[size(I_GT_LP,1)*size(I_GT_LP,2), size(I_GT_LP,3)])];
    A = eye(size(I_MS,3)+1);
    A(1,1) = -1;

    gamma(:,ii) = lsqlin(H,h1(:)-h2(:),A,zeros(1,size(I_MS,3)+1),[],[],[],[],[],opts1);
    I_Fus_BDSD(:,:,ii) = I_MS_en(:,:,ii) + reshape([I_PAN(:),reshape(I_MS_en,[size(I_MS,1)*size(I_MS,2), size(I_MS,3)])]*gamma(:,ii),[size(I_MS,1) size(I_MS,2)]);
end

end