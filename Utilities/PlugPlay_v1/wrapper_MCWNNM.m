function out = wrapper_MCWNNM(in,sigma)

% the parameters for denoising real images
Par.win          = 20; % Non-local patch searching window
Par.delta        = 0;   % Parameter between each iter
Par.Constant  = 2 * sqrt(2);   % Constant num for the weight vector
Par.Innerloop = 2;   % InnerLoop Num of between re-blockmatching
Par.ps            = 6;   % Patch size
Par.step         = 5;
Par.Iter          = 2;    % total iter numbers, the parameter K1 in the paper
Par.display  = true;

Par.method = 'MCWNNM_ADMM';
% the parameters for ADMM
Par.maxIter = 10;% the parameter K2 in the paper
Par.rho = 6;
Par.mu = 1;

% the parameter for estimating noise standard deviation
Par.lambda = 1.5; % We need tune this parameter for different real noisy images

i =4; % dog
Par.image = i;

Par.nim = in.*255;
Par.I = Par.nim;
Par.nlsp  =  70;                           % Initial Non-local Patch number

[h, w, ch] = size(Par.I);
for c = 1:ch
    Par.nSig0(c,1) = sigma; %NoiseEstimation(Par.nim(:, :, c), Par.ps);
end
fprintf('The noise levels are %2.2f, %2.2f, %2.2f. \n', Par.nSig0(1), Par.nSig0(2), Par.nSig0(3) );
[im_out, Par] = MCWNNM_ADMM1_NL_Denoising( Par.nim, Par.I, Par ); % WNNM denoisng function
im_out(im_out>255)=255;
im_out(im_out<0)=0;
out = im_out./255;

end