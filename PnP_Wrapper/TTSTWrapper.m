function img_out = TTSTWrapper(img_in,ratio)
current_dir = pwd;
img_out = zeros(size(img_in,1)*ratio,size(img_in,2)*ratio,size(img_in,3));
% 定义Python脚本的路径
scriptPath = 'C:\ELI\codes_3rd\TTST-main\'; % 替换为实际路径
scriptName = 'eval_4x.py';
pyCmdPath = 'C:\Users\taoji\miniconda3\envs\FunSR\python.exe';
data_root = [scriptPath 'dataset\LR\RS\'];
result_root = [scriptPath 'results\'];
cd(scriptPath);
command = ['"' pyCmdPath '"' ' ' '"' scriptPath scriptName '"'];
Res_LR = img_in;
save([data_root 'input.mat'], 'Res_LR');
[status, cmdout] = system(command);
% 检查命令执行状态
if status == 0
	disp('Python script executed successfully.');
	disp(cmdout);
else
	disp('Python script execution failed.');
	disp(num2str(status)); % 显示错误代码
end
S = load([result_root 'output.mat']);
img_out = double(S.Res_HR)./255;
% for i = 1:size(img_in,3)
%     %imwrite(im2uint8(ImageStretchByTol(img_in(:,:,i))),[data_root '/input.png']);
%     img = repmat(img_in(:,:,i), [1, 1, 3]);
%     imwrite(im2uint8((img)),[data_root 'input.png']);
% 
%     % 使用MATLAB的system函数执行命令
%     [status, cmdout] = system(command);
%     % 检查命令执行状态
%     if status == 0
%         disp('Python script executed successfully.');
%         disp(cmdout);
%     else
%         disp('Python script execution failed.');
%         disp(num2str(status)); % 显示错误代码
%     end
% 
%     t = double(imread([result_root 'RSinput.png']));
%     img_out(:,:,i) = t(:,:,1)./255;
% end

cd(current_dir);