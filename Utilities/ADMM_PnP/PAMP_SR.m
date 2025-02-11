function out = PAMP_SR(y,sensorInf,ratio,lambda,denoise,Opts)
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
    Opts.tol = 1e-5;%1e-4;
end
if ~isfield(Opts,'gamma')
    Opts.gamma=1;
end
if ~isfield(Opts,'print')
    Opts.print = false;
end

% set parameters
max_itr   = Opts.max_itr;
tol       = 1e-5;
gamma     = Opts.gamma;
rho_v       = Opts.rho;
method = Opts.method;
%initialize variables
[rows_in,cols_in, band_in] = size(y);
rows      = rows_in*ratio;
cols      = cols_in*ratio;
N         = rows*cols*band_in;

% method = Opts.method;
h = sensorInf.PSF_G;

[G,Gt]    = defGGt_mod(ratio,sensorInf);
GGt       = constructGGt_mod(ratio,sensorInf,rows,cols,band_in);
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
epsilon = 1e-3;
rng(1);
% bn = zeros(size(x));
% for i = 1:size(x,3)
%     bn(:,:,i) = randn(rows, cols);
% end
bn = randn(size(x));
while((residual>tol) && (itr<=max_itr))
    %store x, v, u from previous iteration for psnr residual calculation
    x_old=x;
    v_old=v;
    u_old=u;
    
    
    %denoising step
    vtilde = x+(1/rho_v)*u;
%     vtilde = proj(vtilde);
    sigma  = sqrt(lambda/rho_v);
    v      = double(denoise(vtilde,sigma));
    if itr <= 5
        v_bn      = double(denoise(vtilde+epsilon*bn,sigma));
        dv = (v_bn - v)./(N*epsilon);
        divDv = real(bn(:)'*dv(:));
        rho = rho_v./divDv;
    end
    
    %inversion step
    xtilde = v-(1/rho).*u;
    
    rhs = Gty + rho*xtilde;
    x = (rhs - Gt(ifft2(fft2(G(rhs))./(GGt + rho))))/rho;
    if itr <= 5
        rho_v = rho/(rho+1);
    end

    
    %update langrangian multiplier
    u      = u + rho*(x-v);
    
    %update rho
    %rho=rho*gamma;
    
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