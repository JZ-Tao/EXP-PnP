function img_HR = wrapperFSRCNN(img,ratio)
if ~exist('ratio', 'var')
    ratio = 4; % the magnification factor x2, x3, x4...
end
model = ['model\FSRCNN\x' num2str(ratio) '.mat'];

%% FSRCNN
img_HR = FSRCNN(model, img, ratio);
