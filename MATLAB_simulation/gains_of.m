function gj = gains_of(kind, shape, wc)
%GAINS_OF Channel gain vector from bandwidth wc and shape parameters.
%   kind 2 : shape = [a b cd rf bsp] -> [kp ki kd wf bsp 0]
%   kind 3/4 : shape = r = wo/wc     -> [wc wo 0 0 0 0]
switch kind
    case 2, gj = [shape(1)*wc^2, shape(2)*wc^3, shape(3)*wc, shape(4)*wc, shape(5), 0];
    otherwise, gj = [wc, shape(1)*wc, 0, 0, 0, 0];
end
end


