function [trans_mat,det_jacob] = getCoorTransMatrix(elesize_x, elesize_y)
% getCoorTransMatrix defines the coordinate transformation matrix between
% the natural coordinates of the isoparametric element and the reference
% coordiantes.
%
%
% Created by Hao Gao (SJTU)
% Create on Feb 28, 2024
% Modified on Feb 28, 2024
% -------------------------------------------------------------------------

trans_mat = [(2/elesize_x)^2, 0, 0;
    0, (2/elesize_y)^2, 0;
    0, 0, (2/elesize_x)*(2/elesize_y)/2];

det_jacob = (2/elesize_x)*(2/elesize_y);

end