function [mass_mat_asb,stiff_mat_asb,couple_mat_asb] = assembMat(mass_mat_asb, stiff_mat_asb, couple_mat_asb, node_index, piezo_index,...
    stiff_mat_shim, stiff_mat_piezo, couple_mat, mass_mat)
% assembMat assembles the mass and stiffness matrices following the order
% of the nodes.
%
%
% Created by Hao Gao (SJTU)
% Create on Feb 29, 2024
% Modified on Feb 29, 2024
% -------------------------------------------------------------------------
num_dof = length(node_index);

for i_dof = 1:num_dof
    i_dof_index = node_index(i_dof);
    for j_dof = 1:num_dof
        j_dof_index = node_index(j_dof);
        mass_mat_asb(i_dof_index,j_dof_index) = mass_mat_asb(i_dof_index,j_dof_index)+mass_mat(i_dof,j_dof);
        stiff_mat_asb(i_dof_index,j_dof_index) = stiff_mat_asb(i_dof_index,j_dof_index)+...
            stiff_mat_shim(i_dof,j_dof)+stiff_mat_piezo(i_dof,j_dof);
    end
    
    % This coupling matrix is the right-up corner stiffness matrix
    couple_mat_asb(i_dof_index,piezo_index) = couple_mat_asb(i_dof_index,piezo_index)+couple_mat(i_dof);
    
end

end