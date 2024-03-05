function mass_ele_mat = getMassMatrixEle(ele_dof, ele_size, mass_area)
% getInertiMatrixEle determines the mass matrix of a given element.
%
%
% Created by Hao Gao (SJTU)
% Create on Feb 29, 2024
% Modified on Feb 29, 2024
% -------------------------------------------------------------------------

% Define the Gaussian points and weights
% Here two Gaussian points are considered
gauss_pts_x = [-sqrt(3/7-2/7*sqrt(6/5)),sqrt(3/7-2/7*sqrt(6/5)),-sqrt(3/7+2/7*sqrt(6/5)),sqrt(3/7+2/7*sqrt(6/5))];
gauss_weight_x = [(18+sqrt(30))/36, (18+sqrt(30))/36, (18-sqrt(30))/36, (18-sqrt(30))/36];

gauss_pts_y = [-sqrt(3/7-2/7*sqrt(6/5)),sqrt(3/7-2/7*sqrt(6/5)),-sqrt(3/7+2/7*sqrt(6/5)),sqrt(3/7+2/7*sqrt(6/5))];
gauss_weight_y = [(18+sqrt(30))/36, (18+sqrt(30))/36, (18-sqrt(30))/36, (18-sqrt(30))/36];

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
        [N, ~, ~, ~] = getShapeFcns(r_coor, s_coor);
                
        % Get the coordinate transformation matrix
        [~, det_jacob] = getCoorTransMatrix(ele_size(1),ele_size(2));
        
        % Sum the Gaussian quadratures
        mass_ele_mat = mass_ele_mat+weight_r*weight_s*det_jacob*mass_area*transpose(N)*N;
        
        
    end
end

end