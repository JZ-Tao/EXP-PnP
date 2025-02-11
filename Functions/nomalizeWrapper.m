function [img_out, Opts] = nomalizeWrapper(img_in, is_denomalize, type, Opts)


if ~is_denomalize
    switch(type)
        case 1
            [M_Res, PS] = mapminmax(im2mat(img_in),0,1);
            img_out = mat2im(M_Res, size(img_in,1));
            Opts.PS = PS;
        case 2
            L = Opts.L;

            max_v = max(2^L-1, 1);%max(max(I_MS(:)), max(I_PAN(:)));
            img_out = img_in./max_v;
    
        case 3

            tol = [0.0001 0.9999];
            for i=1:size(img_in,3)
                b = double(img_in(:,:,i));
                t(1) =  quantile(b(:),tol(1));
                t(2)  = quantile(b(:),tol(2));       
                b(b<t(1))=t(1);
                b(b>t(2))=t(2);
    %             b = (b-t(1))/(t(2)-t(1));
                img_in(:,:,i) = reshape(b, size(img_in,1), size(img_in,2));
            end
            [M_Res, PS] = mapminmax(im2mat(img_in),0,1);
            img_out = mat2im(M_Res, size(img_in,1));
            Opts.PS = PS;
        otherwise
            img_out = img_in;
    end
    
else

    switch(type)
        case {1,3}
            img_out = mat2im(mapminmax('reverse', im2mat(img_in), Opts.PS), size(img_in,1));
        case 2
            L = Opts.L;
            max_v = max(2^L-1, 1);
            img_out = double(img_in).*max_v;
        otherwise
            img_out = img_in;
    end
end