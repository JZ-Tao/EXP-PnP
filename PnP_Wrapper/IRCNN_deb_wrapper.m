function I_out = IRCNN_deb_wrapper(I_in, k, sigma)
if ~exist('sigma', 'var')
    sigmas      = [2, 2.55, 7.65]/255;
    sigma       = sigmas(3);
end
% (non-blind) image deblurring

% @inproceedings{zhang2017learning,
%   title={Learning Deep CNN Denoiser Prior for Image Restoration},
%   author={Zhang, Kai and Zuo, Wangmeng and Gu, Shuhang and Zhang, Lei},
%   booktitle={IEEE Conference on Computer Vision and Pattern Recognition},
%   year={2017}
% }

% If you have any question, please feel free to contact with me.
% Kai Zhang (e-mail: cskaizhang@gmail.com)

useGPU       = 1;

model_folder_full = 'C:\ELI\codes_3rd\PnP related\IRCNN-master\models\';
% folderModel  = 'models';
% taskTestCur  = 'Deblur';

% load(fullfile('kernels','Levin09.mat'));
% kernelType = 1; % 1~8
% if kernelType > 8
%     k = fspecial('gaussian', 25, 1.6);
% else
%     k = kernels{kernelType};
% end

totalIter   = 30; % default
lamda       = (sigma^2)/3; % default 3, ****** from {1 2 3 4} ******
modelSigma1 = 49; % default
modelSigma2 = 13; % ****** from {1 3 5 7 9 11 13 15} ******
modelSigmaS = logspace(log10(modelSigma1),log10(modelSigma2),totalIter);
rho         = sigma^2/((modelSigma1/255)^2);

ns          = min(25,max(ceil(modelSigmaS/2),1));
ns          = [ns(1)-1,ns];

y = double(I_in);
[w,h,c]  = size(y);
V = psf2otf(k,[w,h]);
denominator = abs(V).^2;

if c>1
    denominator = repmat(denominator,[1,1,c]);
    V = repmat(V,[1,1,c]);
end
upperleft   = conj(V).*fft2(y);

if c==1
    load([model_folder_full '\modelgray.mat']);
elseif c==3
    load([model_folder_full '\modelcolor.mat']);
end
z = single(y);
if useGPU
    z           = gpuArray(z);
    upperleft   = gpuArray(upperleft);
    denominator = gpuArray(denominator);
end
tic;
for itern = 1:totalIter
    %%% step 1
    rho = lamda*255^2/(modelSigmaS(itern)^2);
    z = real(ifft2((upperleft + rho*fft2(z))./(denominator + rho)));
    if ns(itern+1)~=ns(itern)
        [net] = loadmodel(modelSigmaS(itern),CNNdenoiser);
        net = vl_simplenn_tidy(net);
        if useGPU
            net = vl_simplenn_move(net, 'gpu');
        end
    end
    %%% step 2
    res = vl_simplenn(net, z,[],[],'conserveMemory',true,'mode','test');
    residual = res(end).x;
    z = z - residual;
end

if useGPU
    %output = im2uint8(gather(z));
    output = double(gather(z));
end
toc;

I_out = double(output);
    

