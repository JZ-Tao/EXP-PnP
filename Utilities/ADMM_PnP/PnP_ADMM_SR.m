function out = PnP_ADMM_SR(y,sensorInf,ratio,lambda,denoise,Opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% out = PlugPlayADMM_super(y,h,K,lambda,method,opts)
%
% inversion step: x=argmin_x(||Ax-y||^2+rho/2||x-(v-u)||^2)
% denoising step: v=Denoise(x+u)
%       update u: u=u+(x-v)
%
%Input:           y    -  the observed gray scale image
%                 h    -  blur kernel
%                 K    -  downsampling factor
%              lambda  -  regularization parameter
%              method  -  denoiser, e.g., 'BM3D'
%       opts.rho       -  internal parameter of ADMM {1}
%       opts.gamma     -  parameter for updating rho {1}
%       opts.maxitr    -  maximum number of iterations for ADMM {20}
%       opts.tol       -  tolerance level for residual {1e-4}   
%       ** default values of opts are given in {}. 
%
%Output:          out  -  recovered gray scale image 
%         
%
%Xiran Wang and Stanley Chan
%Copyright 2016
%Purdue University, West Lafayette, In, USA.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check inputs
% if nargin<5
%     error('not enough input, try again \n');
% elseif nargin==5
%     Opts = [];
% end

% Check defaults
if ~isfield(Opts,'rho')
    Opts.rho = 1;
end
if ~isfield(Opts,'max_itr')
    Opts.max_itr = 20;
end
if ~isfield(Opts,'tol')
    Opts.tol = 0.5e-4;%1e-4;
end
if ~isfield(Opts,'gamma')
    Opts.gamma=1;
end
if ~isfield(Opts,'print')
    Opts.print = false;
end

% set parameters
max_itr   = Opts.max_itr;
tol       = Opts.tol;
gamma     = Opts.gamma;
rho       = Opts.rho;
method = Opts.method;
%initialize variables
[rows_in,cols_in, band_in] = size(y);
rows      = rows_in*ratio;
cols      = cols_in*ratio;
N         = rows*cols*band_in;

% method = Opts.method;
% h = sensorInf.PSF_G;

GT_mode = Opts.GT_mode;
step_gamma = Opts.step_gamma;
if ~GT_mode
    step_gamma = step_gamma/(ratio^2);
    K_P = getInterpKernel(ratio, sensorInf.upsampling.interp_type, 2, sensorInf.upsampling.tap).*step_gamma;
     Gt = @(x) interpWrapper(x, ratio, sensorInf.upsampling).*step_gamma;
else
    K_P = sensorInf.PSF_G.*step_gamma;
    Gt = @(x) interpByKernel(x,ratio,K_P,sensorInf.upsampling.offset);
end
G = @(x) MTF_conv_sample(x, sensorInf, ratio, 1);
% [G,Gt]    = defGGt_mod(ratio,sensorInf);
GGt       = constructGGt_mod(ratio,sensorInf.PSF_G,K_P,rows,cols);
Gty       = Gt(y);
%v         = imresize(y,K);
v = interpWrapper(y, ratio, sensorInf.upsampling);


x         = v;
u         = zeros(size(v));
residual  = inf;
% method = 'TV';
%set function handle for denoiser
% switch method
%     case 'BM3D'
%         denoise=@wrapper_BM3D;
%     case 'TV'
%         denoise=@wrapper_TV;
%     case 'NLM'
%         denoise=@wrapper_NLM;
%     case 'RF'
%         denoise=@wrapper_RF;
%     case 'FFDNet'
%         denoise=@wrapper_FFDNet;
%     otherwise
%         error('unknown denoiser \n');
% end

% main loop
if Opts.print==true
    fprintf('Plug-and-Play ADMM --- Super Resolution \n');
    fprintf('Denoiser = %s \n\n', method);
    fprintf('itr \t ||x-xold|| \t ||v-vold|| \t ||u-uold|| \n');
end

itr = 1;

while((residual>tol) && (itr<=max_itr))
    %store x, v, u from previous iteration for psnr residual calculation
    x_old=x;
    v_old=v;
    u_old=u;
    
    %inversion step
    xtilde = v-u;
    
    rhs = Gty + rho*xtilde;
    x = real((rhs - Gt(ifft2(fft2(G(rhs))./(GGt + rho))))/rho);
    
    %denoising step
    vtilde = x+u;
%     vtilde = proj(vtilde);
    sigma  = sqrt(lambda/rho);
    v      = denoise(vtilde,sigma);
    %rho_v = 0.5; v = rho_v*v + (1-rho_v)*vtilde;
    %v = vtilde; %%%%%%%%%%%%%%%%%%%%%%%
%     if mod(itr,3)
%         figure
%         subplot(1,2,1); imshow(x(:,:,1),[]);
%         subplot(1,2,2); imshow(v(:,:,1),[]);
%         1;
%     end
    
    %update langrangian multiplier
    u      = u + (x-v);
    
    %update rho
    rho=rho*gamma;
    
    %calculate residual
    residualx = (1/sqrt(N))*(sqrt(sum(sum(sum((x-x_old).^2)))));
    residualv = (1/sqrt(N))*(sqrt(sum(sum(sum((v-v_old).^2)))));
    residualu = (1/sqrt(N))*(sqrt(sum(sum(sum((u-u_old).^2)))));
    
    residual = residualx + residualv + residualu;

    if Opts.print==true
        fprintf('%3g \t %3.5e \t %3.5e \t %3.5e \n', itr, residualx, residualv, residualu);
    end
    
    itr=itr+1;
end
out=real(x);
end