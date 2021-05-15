function[stackBr] = simple_offset_xyz(stackB, offsets, max_dims)
    
    xoffset = offsets(1); yoffset = offsets(2); zoffset = offsets(3);

    max_x = max_dims(1);
    max_y = max_dims(2);
    max_z = max_dims(3);
    
    stackB_resize = zeros(max_y,max_x,max_z,'uint16');
    stackB_resize(1:size(stackB,1),1:size(stackB,2),1:size(stackB,3)) = stackB;
    stackB_resize(stackB_resize == 0) = median(stackB(:));
    stackB = stackB_resize;
    
    clearvars stackA_resize stackB_resize

    stackBr = zeros(size(stackB),'uint16');
    
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
    
end