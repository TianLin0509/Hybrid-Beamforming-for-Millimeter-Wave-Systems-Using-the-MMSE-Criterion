function [stop ] = stoppingcriterion(options, info)
stop = 0;
% Target gradient norm attained
if      info.gradnorm < options.tolgradnorm
    stop = 1;
    return;
end

% Alloted iteration count exceeded
if       info.iter >= options.maxiter
    stop = 1;
    return;
end

if  abs(info.stepsize) < options.minstepsize
    stop = 1;
    return;
end

end
