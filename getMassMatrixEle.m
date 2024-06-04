function mass_ele_mat = getMassMatrixEle(ele_dof, node_dof, ele_size, mass_area, node_index, node_coord)
% getInertiMatrixEle determines the mass matrix of a given element.
%
%
% Created by Hao Gao (SJTU)
% Create on Feb 29, 2024
% Modified on Feb 29, 2024
% -------------------------------------------------------------------------

% Define the Gaussian points and weights
% Here two Gaussian points are considered
gauss_pts_x = [0, -sqrt(5-2*sqrt(10/7))/3, sqrt(5-2*sqrt(10/7))/3, -sqrt(5+2*sqrt(10/7))/3, sqrt(5+2*sqrt(10/7))/3];
gauss_weight_x = [128/225, (322+13*sqrt(70))/900, (322+13*sqrt(70))/900, (322-13*sqrt(70))/900, (322-13*sqrt(70))/900];

gauss_pts_y = [0, -sqrt(5-2*sqrt(10/7))/3, sqrt(5-2*sqrt(10/7))/3, -sqrt(5+2*sqrt(10/7))/3, sqrt(5+2*sqrt(10/7))/3];
gauss_weight_y = [128/225, (322+13*sqrt(70))/900, (322+13*sqrt(70))/900, (322-13*sqrt(70))/900, (322-13*sqrt(70))/900];

% Initialize the element stiffness matrix
mass_ele_mat = zeros(ele_dof);

for intx = 1:length(gauss_pts_x)
    for inty = 1:length(gauss_pts_y)
        
        % Extract the Gaussian point for calculation
        r_coor = gauss_pts_x(intx);
        weight_r = gauss_weight_x(intx);
        s_coor = gauss_pts_y(inty);
        weight_s = gauss_weight_y(inty);
        
        % Get the shape functions and their derivatives at the given
        % Gaussian point
        switch node_dof
            case 1
                [N, ~, ~, ~] = getShapeFcns(r_coor, s_coor);
            case 3
                [N, ~, ~, ~] = getShapeFcns12DOF(r_coor, s_coor, ele_size);
            case 4
                [N, ~, ~, ~] = getShapeFcns16DOF(r_coor, s_coor, ele_size);
        end
                
        % Get the coordinate transformation matrix
        [~, det_jacob] = getCoorTransMatrix(node_index, node_coord, r_coor, s_coor);
        
        % Sum the Gaussian quadratures
        mass_ele_mat = mass_ele_mat+weight_r*weight_s*det_jacob*mass_area*transpose(N)*N;
        
        
    end
end

end