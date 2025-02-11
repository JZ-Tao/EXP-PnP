function img_HR = wrapperSRCNN(img_U, ratio)
if ~exist('ratio', 'var')
    ratio = 4; % the magnification factor x2, x3, x4...
end
model = ['model\9-5-5(ImageNet)\x' num2str(ratio) '.mat'];

%% SRCNN
img_HR = SRCNN(model, img_U);
