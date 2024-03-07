function [trans_mat,det_jacob] = getCoorTransMatrix(node_index, node_coord, r, s)
% getCoorTransMatrix defines the coordinate transformation matrix between
% the natural coordinates of the isoparametric element and the reference
% coordiantes.
%
%
% Created by Hao Gao (SJTU)
% Create on Feb 28, 2024
% Modified on Mar 07, 2024
% -------------------------------------------------------------------------

% trans_mat = [(2/elesize_x)^2, 0, 0;
%     0, (2/elesize_y)^2, 0;
%     0, 0, (2/elesize_x)*(2/elesize_y)];

% det_jacob = (elesize_x/2)*(elesize_y/2);

num_node = length(node_index);

% Calculate the derivatives of shape functions
dNdr = [-1/4*(1-s), 1/4*(1-s), 1/4*(1+s), -1/4*(1+s)];
dNds = [-1/4*(1-r), -1/4*(1+r), 1/4*(1+r), 1/4*(1-r)];

% Initialize derivative terms
dxdr =0; dxds =0; dydr = 0; dyds = 0;

for i_node = 1:num_node
    dxdr = dxdr+dNdr(i_node)*node_coord(node_index(i_node),2);
    dxds = dxds+dNds(i_node)*node_coord(node_index(i_node),2);
    dydr = dydr+dNdr(i_node)*node_coord(node_index(i_node),3);
    dyds = dyds+dNds(i_node)*node_coord(node_index(i_node),3);    
end

if abs(dxdr) < 1e-10
    drdx = 0;
else
    drdx = 1/dxdr;
end
if abs(dxds) < 1e-10
    dsdx = 0;
else
    dsdx = 1/dxds;
end
if abs(dydr) < 1e-10
    drdy = 0;
else
    drdy = 1/dydr;
end
if abs(dyds) < 1e-10
    dsdy = 0;
else
    dsdy = 1/dyds;
end
trans_mat = [drdx^2, dsdx^2, drdx*dsdx;
             drdy^2, dsdy^2, drdy*dsdy;
             2*drdx*drdy, 2*dsdx*dsdy, dsdx*drdy+drdx*dsdy];

jacob = [dxdr, dydr; dxds, dyds];
det_jacob = det(jacob);


end