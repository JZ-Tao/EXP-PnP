function out = PnP_ADMM_deblur(y,sensorInf,lambda,denoise,opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%out = PlugPlayADMM_deblur(y,h,lambda,method,opts)
%deblurs image y by solving the ADMM:
%
%inversion step: x=argmin_x(||Ax-y||^2+rho/2||x-(v-u)||^2)
%denoising step: v=Denoise(x+u)
%      update u: u=u+(x-v)
%
%Input:           y    -  the observed gray scale image
%                 h    -  blur kernel
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
%Xiran Wang and Stanley Chan
%Copyright 2016
%Purdue University, West Lafayette, In, USA.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check inputs
if nargin<4
    error('not enough input, try again \n');
elseif nargin==4
    opts = [];
end

% Check defaults
if ~isfield(opts,'rho')
    opts.rho = 1;
end
if ~isfield(opts,'max_itr')
    opts.max_itr = 20;
end
if ~isfield(opts,'tol')
    opts.tol = 1e-4;
end
if ~isfield(opts,'gamma')
    opts.gamma=1;
end
if ~isfield(opts,'print')
    opts.print = false;
end

method = opts.method;
% set parameters
max_itr   = opts.max_itr;
tol       = opts.tol;
gamma     = opts.gamma;
rho       = opts.rho;

%initialize variables
[dim(1) dim(2) dim(3)] = size(y);          
N           = dim(1)*dim(2)*dim(3);
h = sensorInf.PSF_G;

Hty = zeros(size(y));
for ii = 1 : dim(3)
    Hty(:,:,ii) = imfilter(y(:,:,ii),real(h(:,:,ii)),'circular'); % replicate
end
% Hty = MTF(y, sensorInf, ratio, 0);
% Hty         = imfilter(y,h,'circular'); 
% eigHtH      = abs(fftn(h, dim)).^2;
eigHtH = zeros(size(y));
for ii = 1:dim(3)
    eigHtH(:,:,ii) = abs(fft2(sensorInf.PSF_G(:,:,ii), dim(1), dim(2))).^2;
end
v           = 0.5*ones(dim);
x           = v;
u           = zeros(dim);
residual    = inf;

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

if opts.print==true
    fprintf('Plug-and-Play ADMM --- Deblurring \n');
    fprintf('Denoiser = %s \n\n', method);
    fprintf('itr \t ||x-xold|| \t ||v-vold|| \t ||u-uold|| \n');
end

itr = 1;
while(residual>tol&&itr<=max_itr)
    %store x, v, u from previous iteration for psnr residual calculation
    x_old = x;
    v_old = v;
    u_old = u;
    
    %inversion step
    xtilde = v-u;
    % rhs    = fftn(Hty+rho*xtilde,dim);
    % x      = real(ifftn(rhs./(eigHtH+rho),dim));
    rhs    = fft2(Hty+rho*xtilde);
    x      = real(ifft2(rhs./(eigHtH+rho)));
%     for ii = 1:dim(3)
%         rhs    = fft2(Hty(:,:,ii)+rho*xtilde(:,:,ii),dim(1), dim(2));
%         x      = real(ifft2(rhs./(eigHtH(:,:,ii)+rho),dim(1), dim(2)));
% %         rhs    = fft2(Hty(:,:,ii)+rho*xtilde(:,:,ii),dim(1), dim(2));
% %         x      = real(ifft2(rhs./(eigHtH(:,:,ii)+rho),dim(1), dim(2)));
%     end
    %denoising step
    vtilde = x+u;
%     vtilde = proj(vtilde,[-1,1]);
    sigma  = sqrt(lambda/rho);
    v      = denoise(vtilde,sigma);
    
    %update langrangian multiplier
    u      = u + (x-v);
    
    %update rho
    rho = rho*gamma;
    
    %calculate residual
    residualx = (1/sqrt(N))*(sqrt(sum(sum(sum((x-x_old).^2)))));
    residualv = (1/sqrt(N))*(sqrt(sum(sum(sum((v-v_old).^2)))));
    residualu = (1/sqrt(N))*(sqrt(sum(sum(sum((u-u_old).^2)))));
	
    residual = residualx + residualv + residualu;

    if opts.print==true
        fprintf('%3g \t %3.5e \t %3.5e \t %3.5e \n', itr, residualx, residualv, residualu);
    end
    
    itr = itr+1;
end
out = v;
end