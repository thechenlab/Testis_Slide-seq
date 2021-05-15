% Author: Zachary Chiang, Buenrostro Lab, Harvard University
% Performs an initial color correction between IGS cycles and saves the 
% transformation matrix to speed up later iterations(see also: color_correct_heuristic.m)

function[corr_channel pre_corr max_corr tform_new] = color_correct_initial(stackA, stackB, channel)

    stackBr = stackB;
    
    % set registration parameters
    
    metric = registration.metric.MeanSquares;
    optimizer = registration.optimizer.RegularStepGradientDescent;
    optimizer.MaximumIterations = 300;
    optimizer.MaximumStepLength = 0.05;
    
	% calculate correlation prior to registration
    
    max_corr = corr(reshape(max(stackA,[],4),[],1),reshape(max(stackB,[],4),[],1));
    pre_corr = max_corr;
   
    % perform color correction and save tform
    
    [imreg_new, ~, tform_new] = imregister2(stackB(:,:,:,channel), max(stackA,[],4), 'affine', optimizer, metric, 'PyramidLevels', 1);
    stackBr(:,:,:,channel) = imwarp(stackB(:,:,:,channel), tform_new, 'outputView', imref3d(size(stackB(:,:,:,channel))));
    reg_corr = corr(reshape(max(stackA,[],4),[],1),reshape(max(stackBr,[],4),[],1));
        
    % output stack with highest correlation
    
    if reg_corr > max_corr
        max_corr = reg_corr;
        corr_channel = stackBr(:,:,:,channel);
    else
        max_corr = pre_corr;
        corr_channel = stackB(:,:,:,channel);
    end
       
end
