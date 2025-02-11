function [img_HR, time] = SR_MSLapSRN_U(y, net, model_scale, gpu)
% -------------------------------------------------------------------------
%   Description:
%       function to apply SR with MS-LapSRN
%
%   Input:
%       - img_LR        : low-resolution image
%       - net           : MS-LapSRN model
%       - model_scale   : model upsampling scale for constructing pyramid
%       - test_scale    : image upsampling scale
%       - gpu           : GPU ID
%
%   Output:
%       - img_HR: high-resolution image
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

    %% setup
    net.mode = 'test' ;
    output_var = sprintf('x%dSR_%dx_output', model_scale, model_scale);
    output_index = net.getVarIndex(output_var);
    net.vars(output_index).precious = 1;
    
    % RGB to YUV
%     if( size(img_U, 3) > 1 )
%         img_U = rgb2ycbcr(img_U);
%     end
    
    % extract Y
%     y = single(img_U(:, :, 1));
    
    if( gpu )
        y = gpuArray(single(y));
    end

    % forward
    inputs = {sprintf('x%dSR_LR', model_scale), y};
    tic;
    net.eval(inputs);
    time = toc;
    y = gather(net.vars(output_index).value);

        
    % resize if size does not match the output image
%     if( size(y, 1) ~= size(img_U, 1) )
%         y = imresize(y, [size(img_U, 1), size(img_U, 2)]);
%     end
    
    img_HR(:, :, 1) = double(y);
        
    % YUV to RGB
%     if( size(img_HR, 3) > 1 )
%         img_HR = ycbcr2rgb(img_HR);
%     end
        
    
end
