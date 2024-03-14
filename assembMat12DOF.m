function [mass_mat_asb,stiff_mat_asb,couple_mat_asb] = assembMat12DOF(mass_mat_asb, stiff_mat_asb, couple_mat_asb, node_index, piezo_index,...
    stiff_mat_shim, stiff_mat_piezo, couple_mat, mass_mat)
% assembMat assembles the mass and stiffness matrices following the order
% of the nodes.
%
%
% Created by Hao Gao (SJTU)
% Create on Mar 13, 2024
% Modified on Mar 13, 2024
% -------------------------------------------------------------------------
num_ele = size(mass_mat_asb,1)/3;
num_node_ele = length(node_index);
num_dof = num_node_ele*3;

for i_dof = 1:num_dof
    i_node = floor((i_dof-0.1)/3)+1; % Find the local node number (1-4)
    i_dof_node = i_dof-3*(i_node-1); % Find the DOF order of this node (1:w,2:theta_x,3:theta_y)
    i_node_glb = node_index(i_node); % Find the global node number 
    
    for j_dof = 1:num_dof
        j_node = floor((j_dof-0.1)/3)+1; % Find the local node number (1-4)
        j_dof_node = j_dof-3*(j_node-1); % Find the DOF order of this node (1:w,2:theta_x,3:theta_y)
        j_node_glb = node_index(j_node); % Find the global node number
        
        % Locate the assembly location 
        i_dof_glb = (i_dof_node-1)*num_ele+i_node_glb;
        j_dof_glb = (j_dof_node-1)*num_ele+j_node_glb;
         
        mass_mat_asb(i_dof_glb,j_dof_glb) = mass_mat_asb(i_dof_glb,j_dof_glb)+mass_mat(i_dof,j_dof);
        
        stiff_mat_asb(i_dof_glb,j_dof_glb) = stiff_mat_asb(i_dof_glb,j_dof_glb)+...
            stiff_mat_shim(i_dof,j_dof)+stiff_mat_piezo(i_dof,j_dof);
    end
    
    % This coupling matrix is the right-up corner stiffness matrix
    couple_mat_asb(i_node_glb,piezo_index) = couple_mat_asb(i_node_glb,piezo_index)+couple_mat(i_dof);
    
end

end