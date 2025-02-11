function PAN_LP = genPAN_LP(imageHR, ratio, method, sensorInf)
PAN_LP = zeros(size(imageHR));
switch(method)
    case 'MTF'
        PAN_LP = MTF_conv_sample(imageHR,sensorInf,ratio,0);
    case 'Gauss'
        PAN_LP = LPfilterGauss(imageHR,ratio);
    case 'box'
        if ~ rem(ratio,2)
            r = ratio + 1;
        else
            r = ratio;
        end
        for i=1:size(imageHR,3)
            PAN_LP(:,:,i) = imfilter(imageHR(:,:,i),fspecial('average',[r r]),'replicate');
        end
    case {'AWT', 'WT'}
        h=[1 4 6 4 1 ]/16;
        g=[0 0 1 0 0 ]-h;
        htilde=[ 1 4 6 4 1]/16;
        gtilde=[ 0 0 1 0 0 ]+htilde;
        h=sqrt(2)*h;
        g=sqrt(2)*g;
        htilde=sqrt(2)*htilde;
        gtilde=sqrt(2)*gtilde;
        WF={h,g,htilde,gtilde};

        Levels = ceil(log2(ratio));

        for i=1:size(imageHR,3)
            WT = ndwt2_working(imageHR(:,:,i),Levels,WF);

            for ii = 2 : numel(WT.dec), WT.dec{ii} = zeros(size(WT.dec{ii})); end

            PAN_LP(:,:,i) = indwt2_working(WT,'c');
        end
end