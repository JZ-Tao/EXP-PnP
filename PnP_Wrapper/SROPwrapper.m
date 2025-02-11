function I_t = SROPwrapper(img_in, is_LR, method, ratio, Opts)

NormOpts.L = Opts.L;
Norm_type = Opts.normalize_type;



ss = 0;
if ss

    n_sub = floor(size(img_in,3)/2);
    SS.n_sub = n_sub;
    SS.prePCA = 0;
    basis_type = 'SVD';
    [E, iE, E_r, iE_r,  ~] = genSubspace(img_in, basis_type, [], SS);
    [iETrans, ETrans]= subspaceTransSelector(basis_type, n_sub);
    img_in_org = img_in;
    img_in = iETrans(iE, img_in_org);
    img_in_r = iETrans(iE_r, img_in_org);

end
[img_n, NormOpts] = nomalizeWrapper(img_in, 0, Norm_type, NormOpts);
if is_LR
    I_t = zeros(size(img_n,1)*ratio, size(img_n,2)*ratio, size(img_n,3));
    is_hoc = 1;
else
    I_t = zeros(size(img_n));
    is_hoc = 0;
end
done = 0;
switch(method)
    case 'bicubic'
        I_t = imresize(img_n, ratio);
 
    case 'MS-LapSRN'

        for i = 1:size(img_n,3)
            I_t(:,:,i) = wrapperMSLapSRN(img_n(:,:,i), ratio);
        end
    case 'LapSRN'

        for i = 1:size(img_n,3)
            I_t(:,:,i) = wrapperLapSRN(img_n(:,:,i), ratio);
        end

    case 'SRCNN'

        for i = 1:size(img_n,3)
            I_t(:,:,i) = wrapperSRCNN(img_n(:,:,i));
        end

    case 'FSRCNN'

        for i = 1:size(img_n,3)
            I_t(:,:,i) = wrapperFSRCNN(img_n(:,:,i),ratio);
        end
 
    case 'A-PLUS'
        I_t = A_PLUS_wrapper(img_n);
 
    case 'NEDI'
        for i = 1:size(img_n,3)
            I_t(:,:,i)=sri(img_n(:,:,i),log2(ratio));
        end
        I_t = imtranslate(I_t,[2,2]);
        is_hoc = 0;

    case 'RealESRGAN'
        
        I_t = RealESRGANWrapper(img_n, ratio);
    case 'RealESRNet'
        
        I_t = RealESRNetWrapper(img_n, ratio);        
    case 'TTST'
        I_t = TTSTWrapper(img_n, ratio);

    otherwise
        warning(['Unknown SR method in SROPWrapper:' method]);
        I_t = imresize(img_n, ratio);
  
end


I_t = nomalizeWrapper(I_t, 1, Norm_type, NormOpts);
if ss
    if is_LR
        I_t_r = interpWrapper(img_in_r, ratio, Opts.upsampling);
    else
        I_t_r = img_in_r;
    end
    I_t_s = cat(3, I_t, I_t_r);
    E = [E E_r];
    I_t = ETrans(E, I_t_s);
end

if is_hoc
    [hoc_intp1, hoc_dec1] = offset_compensation_func(1);
    t = hoc_intp1(I_t);
    I_t = hoc_dec1(t);
end