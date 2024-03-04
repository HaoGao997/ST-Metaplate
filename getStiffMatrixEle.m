function [stiff_shim_ele_mat,stiff_piezo_ele_mat,couple_ele_mat] = getStiffMatrixEle(ele_dof, ele_size, bend_shim_mat, bend_piezo_mat, couple_mat)
% getInertiMatrixEle determines the stiffness matrix of a given element.
%
% <bend_shim_mat> is a 3x3 matrix that describes the bending stiffness of the
% shim layer plate.
% <bend_piezo_mat> is a 3x3 matrix that describes the bending stiffness of
% the piezoelectric layer plate.
% <couple_mat> is a 3x1 matrix that represents the electromechanical
% coupling.
%
% Created by Hao Gao (SJTU)
% Create on Feb 29, 2024
% Modified on Feb 29, 2024
% -------------------------------------------------------------------------

% Define the Gaussian points and weights
% Here two Gaussian points are considered
gauss_pts_x = [-1/sqrt(3),1/sqrt(3)];
gauss_weight_x = [1, 1];
gauss_pts_y = [-1/sqrt(3),1/sqrt(3)];
gauss_weight_y = [1, 1];

% Initialize the element stiffness matrix and the coupling matrix
stiff_shim_ele_mat = zeros(ele_dof);
stiff_piezo_ele_mat = zeros(ele_dof);
couple_ele_mat = zeros(ele_dof,1);

for intx = 1:length(gauss_pts_x)
    for inty = 1:length(gauss_pts_y)
        
        % Extract the Gaussian point for calculation
        r_coor = gauss_pts_x(intx);
        weight_r = gauss_weight_x(intx);
        s_coor = gauss_pts_y(inty);
        weight_s = gauss_weight_y(inty);
        
        % Get the shape functions and their derivatives at the given
        % Gaussian point
        [~, ddNddr, ddNdds, ddNdrds] = getShapeFcns(r_coor, s_coor);
        
        % Define the shape function matrix
        shapefcn_mat = [ddNddr; ddNdds; 2*ddNdrds];
        
        % Get the coordinate transformation matrix
        [trans_mat, det_jacob] = getCoorTransMatrix(ele_size(1),elesize_y(2));
        
        % Sum the Gaussian quadratures
        stiff_shim_ele_mat = stiff_shim_ele_mat+weight_r*weight_s*det_jacob*...
            transpose(shapefcn_mat)*transpose(trans_mat)*bend_shim_mat*...
            trans_mat*shapefcn_mat;
        
        stiff_piezo_ele_mat = stiff_piezo_ele_mat+weight_r*weight_s*det_jacob*...
            transpose(shapefcn_mat)*transpose(trans_mat)*bend_piezo_mat*...
            trans_mat*shapefcn_mat;
        
        couple_ele_mat = couple_ele_mat+weight_r*weight_s*det_jacob*...
            transpose(couple_mat)*trans_mat*shapefcn_mat;
        
        
    end
end

end