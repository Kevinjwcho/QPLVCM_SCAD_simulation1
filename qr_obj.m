% Objective function for quantile regression
% ========================================================================
function obj = qr_obj(x, tau)
    obj = (0.5 * abs(x) + (tau-0.5)* x);
end
% ========================================================================