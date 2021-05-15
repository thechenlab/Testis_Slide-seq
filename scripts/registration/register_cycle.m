% Author: Zachary Chiang, Buenrostro Lab, Harvard University
% Uses cross correlation and intensity-based registration (imregister) to 
% register one stack to another in x, y, and z, both stacks should be
% maximum intensity projected in the channel dimension and then the
% registration is applied to all channels separately

function[cycle_stackBr reg_cycle stats tform] = register_cycle(stackA, stackB, cycle, min_overlap)

    % get stack dimensions

    xlen = length(stackA(:,1,1)); ylen = length(stackA(1,:,1)); zlen = length(stackA(1,1,:));
    num_channels = size(cycle,4);
    
    % set registration parameters
    
    metric = registration.metric.MeanSquares;
    optimizer = registration.optimizer.RegularStepGradientDescent;
    optimizer.MaximumIterations = 200;
    optimizer.MaximumStepLength = 0.05;
    
    % pre-set outputs
    
    stackBr = stackB;
    cycle_stackBr = cycle;
    
	% calculate correlation prior to registration    
    
    max_corr = corr(stackA(:),stackB(:));
    pre_corr = max_corr;

    % take maximum intensity projections
        
    imgA = max(stackA,[],3);
    imgB = max(stackB,[],3);

    % apply cross correlation
    
    c = normxcorr2_general(imgB, imgA, size(imgA,1)*size(imgA,2)*min_overlap);
    [max_c, imax] = max(abs(c(:)));
    [ypeak, xpeak] = ind2sub(size(c),imax(1));
    corr_offset = [(xpeak-size(imgB,2)) (ypeak-size(imgB,1))];
    xoffset = corr_offset(1);
    yoffset = corr_offset(2);
        
    imgC = squeeze(max(stackA,[],2));
    imgD = squeeze(max(stackB,[],2));

    c2 = normxcorr2_general(imgD, imgC,size(imgC,1)*size(imgC,2)*min_overlap);
    [max_c2, imax] = max(abs(c2(:)));
    [ypeak, zpeak] = ind2sub(size(c2),imax(1));
    zoffset = zpeak-size(imgD,2);
                
    %disp(sprintf('%s: Offset is %d (x) by %d (y) by %d (z)',sec2time(toc),xoffset,yoffset,zoffset));
    
    % calculate offset bounds

    xbegin_A = max(1+xoffset, 1);
    xend_A   = min(size(stackA,2), size(stackA,2)+xoffset);
    xbegin_B = max(1-xoffset, 1);
    xend_B   = min(size(stackA,2), size(stackA,2)-xoffset);
    
    ybegin_A = max(1+yoffset, 1);
    yend_A   = min(size(stackA,1), size(stackA,1)+yoffset);
    ybegin_B = max(1-yoffset, 1);
    yend_B   = min(size(stackA,1), size(stackA,1)-yoffset);
        
    zbegin_A = max(1+zoffset, 1);
    zend_A   = min(size(stackA,3), size(stackA,3)+zoffset);
    zbegin_B = max(1-zoffset, 1);
    zend_B   = min(size(stackA,3), size(stackA,3)-zoffset);
           
    % apply offset
    
    stackBc = zeros(xlen,ylen,zlen);
    stackBc(ybegin_A:yend_A,xbegin_A:xend_A,zbegin_A:zend_A) = stackB(ybegin_B:yend_B,xbegin_B:xend_B,zbegin_B:zend_B);
    stackBc(stackBc == 0) = mode(stackB(:));

    % calculate offset correlation
    
    offset_corr = corr(stackA(:),stackBc(:));

    %disp(sprintf('%s: Offset correlation is %.05f',sec2time(toc),offset_corr));

    % only use offset if correlation is better
    
    if offset_corr > max_corr
        max_corr = offset_corr;
        stackBr = stackBc;
        
        cycle_stackBr = zeros(xlen,ylen,zlen,num_channels);
        cycle_stackBr(ybegin_A:yend_A,xbegin_A:xend_A,zbegin_A:zend_A,:) = cycle(ybegin_B:yend_B,xbegin_B:xend_B,zbegin_B:zend_B,:);
        cycle_stackBr(cycle_stackBr == 0) = mode(cycle(:));
    end
    
    % perform intensity-based registration
        
    [imreg_new, ~, tform] = imregister2(stackBc, stackA, 'affine', optimizer, metric, 'PyramidLevels', 1);
    stackBt = imwarp(stackBc, tform, 'outputView', imref3d(size(stackBc)));
    reg_corr = corr(stackA(:),stackBt(:));
        
    % apply registration to all channels if better correlation
    
    if reg_corr > max_corr
        max_corr = reg_corr;
        stackBr =  stackBt;
        for channel=1:num_channels
            cycle_stackBr(:,:,:,channel)=imwarp(cycle_stackBr(:,:,:,channel), tform, 'outputView', imref3d([size(cycle_stackBr,1) size(cycle_stackBr,2) size(cycle_stackBr,3)]));
        end
    end
        
    % output registered stack and statistics
    
    reg_cycle = stackBr;
    post_corr = max_corr;
    %disp(sprintf('%s: Final correlation is %.05f',sec2time(toc),post_corr));
    stats = sprintf('initial=%.03f, offset(%+d, %+d, %+d)=%+.03f, registration=%+.03f, final=%.03f',pre_corr,xoffset,yoffset,zoffset,offset_corr-pre_corr,reg_corr-offset_corr,post_corr);
       
end
