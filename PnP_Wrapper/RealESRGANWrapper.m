function img_out = RealESRGANWrapper(img_in, ratio)
current_dir = pwd;
RealESRGAN_root = [current_dir '\RealESRGAN_Temp'];
data_root = [RealESRGAN_root '\'];
exe_bin = ['C:\ELI\codes_3rd\realesrgan-ncnn-vulkan-20210801-windows\realesrgan-ncnn-vulkan.exe']
mkdir(data_root);
% 保存当前res_LR、卷积核的mat文件
% save([data_root '\img_LR' num2str(data_index) '.mat'], 'img_LR');
%method = 'realesrnet-x4plus';
method = 'realesrgan-x4plus';
img_out = zeros(size(img_in,1)*ratio, size(img_in,2)*ratio, size(img_in,3));
for i = 1:size(img_in,3)
    %imwrite(im2uint8(ImageStretchByTol(img_in(:,:,i))),[data_root '/input.png']);
    imwrite(im2uint8((img_in(:,:,i))),[data_root '/input.png']);
    command = [exe_bin ' -i "' data_root '/input.png" -o "' data_root '/output.png" -n ' method];
    status = system(command);
    t = double(imread([data_root '/output.png']));
    img_out(:,:,i) = t(:,:,1)./255;
end

%imwrite(im2uint8(ImageStretchByTol(img_in(:,:,i))),[data_root '/input.png']);
% imwrite(im2uint8((img_in)),[data_root '/input.png']);
% command = [exe_bin ' -i "' data_root '/input.png" -o "' data_root '/output.png" -n ' method];
% status = system(command);
% t = double(imread([data_root '/output.png']));
% img_out = t./255;

% img_HR = load([data_root '\img_HR' num2str(data_index) '.mat']);
