 function [shapefcn,shapefcn_ddr,shapefcn_dds,shapefcn_drds] = getShapeFcns16DOF(r, s, a)
% getShapeFcns evaluates the shape functions and their derivatives at given
% local natural coordinates defined in the isoparametric element.
%
%
% Created by Hao Gao (SJTU)
% Create on Jun 04, 2024
% Modified on Jun 04, 2024
% -------------------------------------------------------------------------

% The shape functions are associated with the 12 DOFs of the quadrilaterial
% element. For each node, there are a transverse displacement and two
% rotation angles.

shapefcn = zeros(1,16);
shapefcn_ddr = zeros(1,16);
shapefcn_dds = zeros(1,16);
shapefcn_drds = zeros(1,16);

% Define the 1D Hermite polynomials
Hr_01 = 1/2-3*r/4+r^3/4;
Hr_02 = 1/2+3*r/4-r^3/4;
Hr_11 = (r^3-r^2-r+1)/8*a;
Hr_12 = (r^3+r^2-r-1)/8*a;

Hs_01 = 1/2-3*s/4+s^3/4;
Hs_02 = 1/2+3*s/4-s^3/4;
Hs_11 = (s^3-s^2-s+1)/8*a;
Hs_12 = (s^3+s^2-s-1)/8*a;

Hdr_01 = 3*r^2/4-3/4; Hdr_02 = 3/4-3*r^2/4; Hdr_11 = (3*r^2-2*r-1)/8*a; Hdr_12 = (3*r^2+2*r-1)/8*a;
Hds_01 = 3*s^2/4-3/4; Hds_02 = 3/4-3*s^2/4; Hds_11 = (3*s^2-2*s-1)/8*a; Hds_12 = (3*s^2+2*s-1)/8*a;

Hddr_01 = 3*r/2; Hddr_02 =-3*r/2; Hddr_11 = (3*r-1)/4*a; Hddr_12 = (3*r+1)/4*a;
Hdds_01 = 3*s/2; Hdds_02 =-3*s/2; Hdds_11 = (3*s-1)/4*a; Hdds_12 = (3*s+1)/4*a;


% Evaluate the shape functions
shapefcn(1) = Hr_01*Hs_01;
shapefcn(2) = Hr_11*Hs_01;
shapefcn(3) = Hr_01*Hs_11;
shapefcn(4) = Hr_11*Hs_11;

shapefcn(5) = Hr_02*Hs_01;
shapefcn(6) = Hr_12*Hs_01;
shapefcn(7) = Hr_02*Hs_11;
shapefcn(8) = Hr_12*Hs_11;

shapefcn(9) = Hr_02*Hs_02;
shapefcn(10) = Hr_12*Hs_02;
shapefcn(11) = Hr_02*Hs_12;
shapefcn(12) = Hr_12*Hs_12;

shapefcn(13) = Hr_01*Hs_02;
shapefcn(14) = Hr_11*Hs_02;
shapefcn(15) = Hr_01*Hs_12;
shapefcn(16) = Hr_11*Hs_12;

% Evaluate the second derivative dd/ddr
shapefcn_ddr(1) = Hddr_01*Hs_01;
shapefcn_ddr(2) = Hddr_11*Hs_01;
shapefcn_ddr(3) = Hddr_01*Hs_11;
shapefcn_ddr(4) = Hddr_11*Hs_11;

shapefcn_ddr(5) = Hddr_02*Hs_01;
shapefcn_ddr(6) = Hddr_12*Hs_01;
shapefcn_ddr(7) = Hddr_02*Hs_11;
shapefcn_ddr(8) = Hddr_12*Hs_11;

shapefcn_ddr(9) = Hddr_02*Hs_02;
shapefcn_ddr(10) = Hddr_12*Hs_02;
shapefcn_ddr(11) = Hddr_02*Hs_12;
shapefcn_ddr(12) = Hddr_12*Hs_12;

shapefcn_ddr(13) = Hddr_01*Hs_02;
shapefcn_ddr(14) = Hddr_11*Hs_02;
shapefcn_ddr(15) = Hddr_01*Hs_12;
shapefcn_ddr(16) = Hddr_11*Hs_12;

% Evaluate the second derivative dd/dds
shapefcn_dds(1) = Hr_01*Hdds_01;
shapefcn_dds(2) = Hr_11*Hdds_01;
shapefcn_dds(3) = Hr_01*Hdds_11;
shapefcn_dds(4) = Hr_11*Hdds_11;

shapefcn_dds(5) = Hr_02*Hdds_01;
shapefcn_dds(6) = Hr_12*Hdds_01;
shapefcn_dds(7) = Hr_02*Hdds_11;
shapefcn_dds(8) = Hr_12*Hdds_11;

shapefcn_dds(9) = Hr_02*Hdds_02;
shapefcn_dds(10) = Hr_12*Hdds_02;
shapefcn_dds(11) = Hr_02*Hdds_12;
shapefcn_dds(12) = Hr_12*Hdds_12;

shapefcn_dds(13) = Hr_01*Hdds_02;
shapefcn_dds(14) = Hr_11*Hdds_02;
shapefcn_dds(15) = Hr_01*Hdds_12;
shapefcn_dds(16) = Hr_11*Hdds_12;

% Evaluate the mixed derivative dd/drds
shapefcn_drds(1) = Hdr_01*Hds_01;
shapefcn_drds(2) = Hdr_11*Hds_01;
shapefcn_drds(3) = Hdr_01*Hds_11;
shapefcn_drds(4) = Hdr_11*Hds_11;

shapefcn_drds(5) = Hdr_02*Hds_01;
shapefcn_drds(6) = Hdr_12*Hds_01;
shapefcn_drds(7) = Hdr_02*Hds_11;
shapefcn_drds(8) = Hdr_12*Hds_11;

shapefcn_drds(9) = Hdr_02*Hds_02;
shapefcn_drds(10) = Hdr_12*Hds_02;
shapefcn_drds(11) = Hdr_02*Hds_12;
shapefcn_drds(12) = Hdr_12*Hds_12;

shapefcn_drds(13) = Hdr_01*Hds_02;
shapefcn_drds(14) = Hdr_11*Hds_02;
shapefcn_drds(15) = Hdr_01*Hds_12;
shapefcn_drds(16) = Hdr_11*Hds_12;

end