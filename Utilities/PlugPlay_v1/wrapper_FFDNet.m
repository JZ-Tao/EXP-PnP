function out = wrapper_FFDNet(in,sigma)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% out = wrapper_FFDNet(in,sigma)
% performs FFDNet denoising
% 
% Require FFDNet package
%
% Download:
% https://github.com/cszn/FFDNet
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         L = 30/255;
%         D = diag(rand(3,1));
%         U = orth(rand(3,3));
%         imageNoiseSigma = abs(L^2*(U' * D * U));
%         assert(sum(sum(imageNoiseSigma<0)) == 0,'We assume that all the elements are non-negative.');
% 
%         % input noise Sigma
%         inputNoiseSigma =  imageNoiseSigma;
%         sigma           =  sqrt(inputNoiseSigma);
%         sigma           =  sigma(:);
%         sigma6          =  sigma([1,4,5,7,8,9]); % extract 6 parameters        
        
        
    global sigmas; % input noise level or input noise level map
    %sigmas = sigma;%sigma;
    useGPU = 1;
%     if size(in,3) == 1
%         load(fullfile('models','FFDNet_gray_200.mat')); % 'FFDNet_gray.mat'
%     elseif size(in,3) == 3
%         load(fullfile('models','FFDNet_color.mat'));
%     else
%         error('Only support image with 1 or 3 channels(ie. gray or RGB)!');
%     end
    load(fullfile('models','FFDNet_gray.mat')); % FFDNet_gray_ft_624000_G.mat FFDNet_gray
    
    
    net = vl_simplenn_tidy(net);

    % for i = 1:size(net.layers,2)
    %     net.layers{i}.precious = 1;
    % end

    if useGPU
        net = vl_simplenn_move(net, 'gpu') ;
    end
    [row, col,n_band] = size(in);
%     for i = 1:n_band
%         max_in(i) = max(max(in(:,:,i)));
%         min_in(i) = min(min(in(:,:,i)));
%         in(:,:,i) = (in(:,:,i) - min_in(i))./(max_in(i) - min_in(i));
%     end

%     [M_in, PS] = mapminmax(im2mat(in),0,1);
%     in = mat2im(M_in, size(in,1));
%     sigma = sigma.*PS.gain;
    %sigma = sigma';    
    
    %in_i = in;
    out = zeros(row,col, n_band);
    
    for i = 1:size(in,3)
       sigmas = sigma;%sigma(i);
       in_i = in(:,:,i);
%         if mod(row,2)==1
%             in_i = cat(1,in_i, in_i(end,:)) ;
%         end
%         if mod(col,2)==1
%             in_i = cat(2,in_i, in_i(:,end)) ;
%         end

        % tic;
        if useGPU
            in_i = gpuArray(in_i);
        end

        res    = vl_simplenn_FFDNet(net,in_i, [],[],'conserveMemory',true,'mode','test');
        out_i = res(end).x;

%         if mod(row,2)==1
%             out_i = out_i(1:end-1,:);
%         end
%         if mod(col,2)==1
%             out_i = out_i(:,1:end-1);
%         end

        % convert to CPU
        if useGPU
            out_i = gather(out_i);
        end
        out(:,:,i) = out_i;
    end
    %out = out_i;
    %out = mat2im(mapminmax('reverse', im2mat(out), PS), size(out,1));
%     for i = 1:n_band
%         out(:,:,i) = out(:,:,i).*(max_in(i) - min_in(i))+min_in(i);
%     end    
%     out = out;
end