function slope = calRK4Slope(state_coord_inst, inv_mass_mat_asb, stiff_mat_asb, couple_mat_asb, inv_electr_mat_inst, force_vec)
% calRK4slope calculates the slope at given time point in the 4th order
% Runge Kutta algorithm of time marching scheme.
%
% The mass, stiffness assembled matrices and the coupling matrix are
% time-invariant during the marching scheme. 
% The capacitance matrix (diagonal) is time varying and determined by the
% instant time input.
%
% Created by Hao Gao (SJTU)
% Create on Mar 04, 2024
% Modified on Mar 04, 2024
% -------------------------------------------------------------------------

num_dof_mech = size(inv_mass_mat_asb,1);
num_dof_electr = size(electr_mat_inst,1);

state_mat = [zeros(num_dof_mech), eye(num_dof_mech), zeros(num_dof_mech,num_dof_electr);
             -inv_mass_mat_asb*stiff_mat_asb, zeros(num_dof_mech), inv_mass_mat_asb*couple_mat_asb;
             zeros(num_dof_electr,num_dof_mech), -inv_electr_mat_inst*transpose(couple_mat_asb), zeros(num_dof_electr)];

slope = state_mat*state_coord_inst+force_vec;

end