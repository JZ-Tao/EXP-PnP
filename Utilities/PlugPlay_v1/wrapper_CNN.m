function out =  wrapper_CNN(in_C,sigma)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% wrapper file for IrCNN
% Please download IrCNN at 
% https://github.com/cszn/IRCNN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CNNSigmaSet    =  [5,10,15,25,35,50];
[~,idxmin]     =  min(abs(CNNSigmaSet-sigma*255));
CNNSigma1      =  CNNSigmaSet(min(idxmin,numel(CNNSigmaSet)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
folderModel  = 'C:\ELI\codes_3rd\PnP related\IRCNN-master\models';%'denoisers/IRCNN-master/models';
out = zeros(size(in_C));
for i = 1:size(in_C,3)
    in = in_C(:,:,i);
    [w,h,c]  = size(in);
    if c==1
        load(fullfile(folderModel,'modelgray.mat'));
    elseif c==3
        load(fullfile(folderModel,'modelcolor.mat'));
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    net            =  loadmodel(CNNSigma1,CNNdenoiser);
    net            =  vl_simplenn_tidy(net);

    res     = vl_simplenn(net,single(in),[],[],'conserveMemory',true,'mode','test');
    zhat    = double(in - res(end).x);

    out(:,:,i) = zhat;
end
end
