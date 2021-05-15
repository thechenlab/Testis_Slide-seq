% Author: Zachary Chiang, Buenrostro Lab, Harvard University
% Uses a previous color correction transformation matrix to speed up
% current color correction between IGS cycles(see also: color_correct_initial.m)

function[corr_channel pre_corr max_corr] = color_correct_heuristic(stackA, stackB, channel, tform)

    % create intermediate and final registered stacks

    stackBr = stackB;
    stackBh = stackB;
    
    % set registration parameters
    
    metric = registration.metric.MeanSquares;
    optimizer = registration.optimizer.RegularStepGradientDescent;
    optimizer.MaximumIterations = 100;
    optimizer.MaximumStepLength = 0.05;
    
    % calculate correlation prior to registration
    
    max_corr = corr(reshape(max(stackA,[],4),[],1),reshape(max(stackB,[],4),[],1));
    pre_corr = max_corr;
    
    % perform intermediate registration using old tranformation matrix
    
    stackBh(:,:,:,channel) = imwarp(stackB(:,:,:,channel), tform, 'outputView', imref3d(size(stackB(:,:,:,channel))));
    mid_corr = corr(reshape(max(stackA,[],4),[],1),reshape(max(stackBh,[],4),[],1));
  
    % perform final registration on intermediate stack
    
    [imreg_new, ~, tform_new] = imregister2(stackBh(:,:,:,channel), max(stackA,[],4), 'affine', optimizer, metric, 'PyramidLevels', 1);
    stackBr(:,:,:,channel) = imwarp(stackBh(:,:,:,channel), tform_new, 'outputView', imref3d(size(stackBh(:,:,:,channel))));
    reg_corr = corr(reshape(max(stackA,[],4),[],1),reshape(max(stackBr,[],4),[],1));
        
    % output stack with highest correlation
    
    if reg_corr > max_corr
        max_corr = reg_corr;
        corr_channel = stackBr(:,:,:,channel);
    elseif mid_corr > max_corr
        max_corr = mid_corr;
        corr_channel = stackBh(:,:,:,channel);
    else
        max_corr = pre_corr;
        corr_channel = stackB(:,:,:,channel);
    end
       
end
