% Author: Zachary Chiang, Buenrostro Lab, Harvard University

function[stackBr,offsets] = simple_offset_xyz(stackA, stackB, max_dims)
    
    max_x = max_dims(1);
    max_y = max_dims(2);
    max_z = max_dims(3);
    
    stackA_resize = zeros(max_y,max_x,max_z,'uint16');
    stackB_resize = zeros(max_y,max_x,max_z,'uint16');

    stackA_resize(1:size(stackA,1),1:size(stackA,2),1:size(stackA,3)) = stackA;
    stackB_resize(1:size(stackB,1),1:size(stackB,2),1:size(stackB,3)) = stackB;
    
    stackA_resize(stackA_resize == 0) = median(stackA(:));
    stackB_resize(stackB_resize == 0) = median(stackB(:));
    
    stackA = stackA_resize;
    stackB = stackB_resize;
    
    clearvars stackA_resize stackB_resize
    
    %figure; imshowpair(max(stackA,[],3),max(stackB,[],3))

    stackBr = zeros(size(stackB),'uint16');

    imgA = squeeze(max(max(stackA,[],3),[],4));
    imgB = squeeze(max(max(stackB,[],3),[],4));
    
    imgC = squeeze(max(max(stackA,[],2),[],4));
    imgD = squeeze(max(max(stackB,[],2),[],4));
    
    c = normxcorr2_general(imgB, imgA, size(imgB,1).*size(imgB,2)*0.5);
    [max_c, imax] = max(abs(c(:)));
    [ypeak, xpeak] = ind2sub(size(c),imax(1));
    corr_offset = [(xpeak-size(imgB,2)) (ypeak-size(imgB,1))];
    xoffset = corr_offset(1);
    yoffset = corr_offset(2);
    
    c = normxcorr2_general(imgD, imgC, size(imgD,1).*size(imgD,2)*0.5);
    %figure; imagesc(c)
    [max_c, imax] = max(abs(c(:)));
    [xpeak, zpeak] = ind2sub(size(c),imax(1));
    corr_offset = [(xpeak-size(imgD,1)) (zpeak-size(imgD,2))];
    zoffset = corr_offset(2);
    
    xbegin_A = max(1+xoffset, 1);
    xend_A   = min(size(stackB,2), size(stackB,2)+xoffset);
    xbegin_B = max(1-xoffset, 1);
    xend_B   = min(size(stackB,2), size(stackB,2)-xoffset);
    
    ybegin_A = max(1+yoffset, 1);
    yend_A   = min(size(stackB,1), size(stackB,1)+yoffset);
    ybegin_B = max(1-yoffset, 1);
    yend_B   = min(size(stackB,1), size(stackB,1)-yoffset);
    
    zbegin_A = max(1+zoffset, 1);
    zend_A   = min(size(stackB,3), size(stackB,3)+zoffset);
    zbegin_B = max(1-zoffset, 1);
    zend_B   = min(size(stackB,3), size(stackB,3)-zoffset);
           
    stackBr(ybegin_A:yend_A,xbegin_A:xend_A,zbegin_A:zend_A,:) = stackB(ybegin_B:yend_B,xbegin_B:xend_B,zbegin_B:zend_B,:);
    stackBr(stackBr == 0) = median(stackB(:));
    
    offsets = [xoffset yoffset zoffset];
    
    %figure; imshowpair(max(stackA,[],3),max(stackBr,[],3))
    
end