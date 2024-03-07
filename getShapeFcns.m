 function [shapefcn,shapefcn_ddr,shapefcn_dds,shapefcn_drds] = getShapeFcns(r, s)
% getShapeFcns evaluates the shape functions and their derivatives at given
% local natural coordinates defined in the isoparametric element.
%
%
% Created by Hao Gao (SJTU)
% Create on Feb 21, 2024
% Modified on Feb 28, 2024
% -------------------------------------------------------------------------

% The shape functions are associated with the 4 DOFs of the quadrilaterial
% element. For each node, there is only ONE transverse displacement DOF.

shapefcn = zeros(1,4);
shapefcn_ddr = zeros(1,4);
shapefcn_dds = zeros(1,4);
shapefcn_drds = zeros(1,4);

% Evaluate the shape functions
shapefcn(1) = (1/2-3*r/4+r^3/4)*(1/2-3*s/4+s^3/4);
shapefcn(2) = (1/2+3*r/4-r^3/4)*(1/2-3*s/4+s^3/4);
shapefcn(4) = (1/2-3*r/4+r^3/4)*(1/2+3*s/4-s^3/4);
shapefcn(3) = (1/2+3*r/4-r^3/4)*(1/2+3*s/4-s^3/4);

% Evaluate the second derivative dd/ddxi
shapefcn_ddr(1) = (3*r/2)*(1/2-3*s/4+s^3/4);
shapefcn_ddr(2) = (-3*r/2)*(1/2-3*s/4+s^3/4);
shapefcn_ddr(4) = (3*r/2)*(1/2+3*s/4-s^3/4);
shapefcn_ddr(3) = (-3*r/2)*(1/2+3*s/4-s^3/4);

% Evaluate the second derivative dd/ddeta
shapefcn_dds(1) = (1/2-3*r/4+r^3/4)*(3*s/2);
shapefcn_dds(2) = (1/2+3*r/4-r^3/4)*(3*s/2);
shapefcn_dds(4) = (1/2-3*r/4+r^3/4)*(-3*s/2);
shapefcn_dds(3) = (1/2+3*r/4-r^3/4)*(-3*s/2);

% Evaluate the mixed derivative dd/dxideta
shapefcn_drds(1) = (3*r^2/4-3/4)*(3*s^2/4-3/4);
shapefcn_drds(2) = (3/4-3*r^2/4)*(3*s^2/4-3/4);
shapefcn_drds(4) = (3*r^2/4-3/4)*(3/4-3*s^2/4);
shapefcn_drds(3) = (3/4-3*r^2/4)*(3/4-3*s^2/4);

end