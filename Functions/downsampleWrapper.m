function Y = downsampleWrapper(X, ratio, Opts)
if exist('Opts', 'var')
    if ~isstruct(Opts)
        error('Incorrect use of parameter');
    end
else
    warning('empty opts');
    Opts = struct; 
end
% if ~isfield(Opts, 'using_imresize'), Opts.using_imresize = 0; end
% if ~isfield(Opts, 'offset'), Opts.offset = [2,1]; end

offset = Opts.offset;

%using_imresize = 0;
if Opts.using_imresize
    Y = imresize(X, 1/ratio, 'nearest'); % equal to Y(:,:,ii) = downsample(downsample(X(:,:,ii),ratio,2)',ratio,2)';
else    
    [row, col, band] = size(X);
    %Y = downsample(downsample(X,ratio,offset)',ratio,offset)';
    Y = zeros(floor(row/ratio), floor(col/ratio), band);
%     
    for ii = 1:size(X,3)
        Y(:,:,ii) = downsample(downsample(X(:,:,ii),ratio,offset)',ratio,offset)';
        %Y(:,:,ii) = downsample(downsample(X(:,:,ii),ratio,0)',ratio,0)';
        %Y(:,:,ii) = downsample(downsample(X(:,:,ii),ratio,1)',ratio,1)';
    end

end