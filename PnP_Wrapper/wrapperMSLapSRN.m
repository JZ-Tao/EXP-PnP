function img_HR = wrapperMSLapSRN(img, ratio)
% -------------------------------------------------------------------------
%   Description:
%       Script to demo MS-LapSRN for one image
%
%   Citation: 
%       Fast and Accurate Image Super-Resolution with Deep Laplacian Pyramid Networks
%       Wei-Sheng Lai, Jia-Bin Huang, Narendra Ahuja, and Ming-Hsuan Yang
%       arXiv, 2017
%
%   Contact:
%       Wei-Sheng Lai
%       wlai24@ucmerced.edu
%       University of California, Merced
% -------------------------------------------------------------------------

%% parameters
%model_name  = 'MSLapSRN_D5R2';
%model_name  = 'MSLapSRN_D5R5';
model_name  = 'MSLapSRN_D5R8';

model_scale = 4; % model SR scale
%test_scale  = 4; % testing SR scale
gpu         = 1; % GPU ID

%% setup paths
%addpath(genpath('utils'));
% addpath(fullfile(pwd, 'matconvnet/matlab'));
% vl_setupnn;

%% Load pretrained multi-scale model
model_filename = fullfile('pretrained_models', sprintf('%s.mat', model_name));
fprintf('Load %s\n', model_filename);

net = load(model_filename);
net_trained = dagnn.DagNN.loadobj(net.net);

opts_filename = fullfile('pretrained_models', sprintf('%s_opts.mat', model_name));
fprintf('Load %s\n', opts_filename);

opts = load(opts_filename);
opts = opts.opts;

opts.scales = [model_scale];

%% create single-scale model
net = init_MSLapSRN_model(opts, 'test');

%% copy pretrained weights
fprintf('Copy weights to single scale model...\n');
net = copy_model_weights(net, net_trained);

if( gpu ~= 0 )
    %gpuDevice(gpu)
    net.move('gpu');
end

%% apply LapSRN
fprintf('Apply MS-LapSRN for %dx SR\n', ratio);

% %% setup
% net.mode = 'test' ;
% output_var = sprintf('x%dSR_%dx_output', model_scale, model_scale);
% output_index = net.getVarIndex(output_var);
% net.vars(output_index).precious = 1;
% 
% % extract Y
% y = single(img_U);
% 
% if( gpu )
%     y = gpuArray(y);
% end
% 
% % bicubic upsample UV
% % img_HR = imresize(img_U, test_scale);
% 
% 
% % forward
% inputs = {sprintf('x%dSR_LR', model_scale), y};
% tic;
% net.eval(inputs);
% time = toc;
% y = gather(net.vars(output_index).value);
% img_HR = double(y);
img_HR = SR_MSLapSRN(img, net, model_scale, ratio, gpu);
% if need_resize
%     img_HR = SR_MSLapSRN(img, net, model_scale, ratio, gpu);
% else
%     img_HR = SR_MSLapSRN_U(img, net, model_scale, gpu);
% end
