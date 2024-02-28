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
shapefcn(1) = (1/2-3*xi/4+xi^3/4)*(1/2-3*eta/4+eta^3/4);
shapefcn(2) = (1/2+3*xi/4-xi^3/4)*(1/2-3*eta/4+eta^3/4);
shapefcn(3) = (1/2-3*xi/4+xi^3/4)*(1/2+3*eta/4-eta^3/4);
shapefcn(4) = (1/2+3*xi/4-xi^3/4)*(1/2+3*eta/4-eta^3/4);

% Evaluate the second derivative dd/ddxi
shapefcn_ddxi(1) = (3*xi/2)*(1/2-3*eta/4+eta^3/4);
shapefcn_ddxi(2) = (-3*xi/2)*(1/2-3*eta/4+eta^3/4);
shapefcn_ddxi(3) = (3*xi/2)*(1/2+3*eta/4-eta^3/4);
shapefcn_ddxi(4) = (-3*xi/2)*(1/2+3*eta/4-eta^3/4);

% Evaluate the second derivative dd/ddeta
shapefcn_ddeta(1) = (1/2-3*xi/4+xi^3/4)*(3*eta/2);
shapefcn_ddeta(2) = (1/2+3*xi/4-xi^3/4)*(3*eta/2);
shapefcn_ddeta(3) = (1/2-3*xi/4+xi^3/4)*(-3*eta/2);
shapefcn_ddeta(4) = (1/2+3*xi/4-xi^3/4)*(-3*eta/2);

% Evaluate the mixed derivative dd/dxideta
shapefcn_dxideta(1) = (3*xi^2/4-3/4)*(3*eta^2/4-3/4);
shapefcn_dxideta(2) = (3/4-3*xi^2/4)*(3*eta^2/4-3/4);
shapefcn_dxideta(3) = (3*xi^2/4-3/4)*(3/4-3*eta^2/4);
shapefcn_dxideta(4) = (3/4-3*xi^2/4)*(3/4-3*eta^2/4);

end