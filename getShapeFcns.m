function [shapefcn,shapefcn_ddxi,shapefcn_ddeta,shapefcn_dxideta] = getShapeFcns(xi, eta)
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
shapefcn_ddxi = zeros(1,4);
shapefcn_ddeta = zeros(1,4);
shapefcn_dxideta = zeros(1,4);

% Evaluate the shape functions
shapefcn(1) = (1-3*xi^2+2*xi^3)*(1-3*eta^2+2*eta^3);
shapefcn(2) = (3*xi^2-2*xi^3)*(1-3*eta^2+2*eta^3);
shapefcn(3) = (1-3*xi^2+2*xi^3)*(3*eta^2-2*eta^3);
shapefcn(4) = (3*xi^2-2*xi^3)*(3*eta^2-2*eta^3);

% Evaluate the second derivative dd/ddxi
shapefcn_ddxi(1) = (12*xi-6)*(1-3*eta^2+2*eta^3);
shapefcn_ddxi(2) = (6-12*xi)*(1-3*eta^2+2*eta^3);
shapefcn_ddxi(3) = (12*xi-6)*(3*eta^2-2*eta^3);
shapefcn_ddxi(4) = (6-12*xi)*(3*eta^2-2*eta^3);

% Evaluate the second derivative dd/ddeta
shapefcn_ddeta(1) = (1-3*xi^2+2*xi^3)*(12*eta-6);
shapefcn_ddeta(2) = (3*xi^2-2*xi^3)*(12*eta-6);
shapefcn_ddeta(3) = (1-3*xi^2+2*xi^3)*(6-12*eta);
shapefcn_ddeta(4) = (3*xi^2-2*xi^3)*(6-12*eta);

% Evaluate the mixed derivative dd/dxideta
shapefcn_dxideta(1) = (6*xi^2-6*xi)*(6*eta^2-6*eta);
shapefcn_dxideta(2) = (6*xi-6*xi^2)*(6*eta^2-6*eta);
shapefcn_dxideta(3) = (6*xi^2-6*xi)*(6*eta-6*eta^2);
shapefcn_dxideta(4) = (6*xi-6*xi^2)*(6*eta-6*eta^2);

end