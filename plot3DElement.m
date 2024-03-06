%% This part is a script to draw the discrete patches of elements
figure; hold all;

for i_ele = 1:MeshParam.num_ele
    
    % Get the index of the nodes of the element
    node_index = MeshParam.element(i_ele,2:5);
    
    % Get the coordinates of the nodes
    x_coord = [MeshParam.node(node_index(1),2);
              MeshParam.node(node_index(2),2);
              MeshParam.node(node_index(3),2);
              MeshParam.node(node_index(4),2)];
          
    y_coord = [MeshParam.node(node_index(1),3);
              MeshParam.node(node_index(2),3);
              MeshParam.node(node_index(3),3);
              MeshParam.node(node_index(4),3)];
    
    z_coord = V(node_index,1);
          
    patch(x_coord,y_coord,z_coord,[0.5;0.5;0.5;0.5]);
    
end

view([45,45]);


%% This part is a hacky script to verify the element stiffness and mass matrix
n_int_pts = 1001;
a = MeshParam.elesize;
int_x = linspace(-a/2,a/2,n_int_pts);
int_y = linspace(-a/2,a/2,n_int_pts);
dx = a/(n_int_pts-1);
dy = dx;

integral = zeros(4);

for i_x = 1:n_int_pts
    if (i_x-1)*(i_x-n_int_pts) == 0
        weight_x = 1/2;
    else
        weight_x = 1;
    end
    x = int_x(i_x);
    for i_y = 1:n_int_pts
        if (i_y-1)*(i_y-n_int_pts) == 0
            weight_y = 1/2;
        else
            weight_y = 1;
        end
        y = int_y(i_y);
        
        K_mat = [shapeFcn(x, a, 1, 2)*shapeFcn(y, a, 1, 0),...
                 shapeFcn(x, a, 2, 2)*shapeFcn(y, a, 1, 0),...
                 shapeFcn(x, a, 1, 2)*shapeFcn(y, a, 2, 0),...
                 shapeFcn(x, a, 2, 2)*shapeFcn(y, a, 2, 0);
                 shapeFcn(x, a, 1, 0)*shapeFcn(y, a, 1, 2),...
                 shapeFcn(x, a, 2, 0)*shapeFcn(y, a, 1, 2),...
                 shapeFcn(x, a, 1, 0)*shapeFcn(y, a, 2, 2),...
                 shapeFcn(x, a, 2, 0)*shapeFcn(y, a, 2, 2);
               2*shapeFcn(x, a, 1, 1)*shapeFcn(y, a, 1, 1),...
               2*shapeFcn(x, a, 2, 1)*shapeFcn(y, a, 1, 1),...
               2*shapeFcn(x, a, 1, 1)*shapeFcn(y, a, 2, 1),...
               2*shapeFcn(x, a, 2, 1)*shapeFcn(y, a, 2, 1)];
        
%         M_mat = [shapeFcn(x, a, 1, 0)*shapeFcn(y, a, 1, 0),...
%                  shapeFcn(x, a, 2, 0)*shapeFcn(y, a, 1, 0),...
%                  shapeFcn(x, a, 1, 0)*shapeFcn(y, a, 2, 0),...
%                  shapeFcn(x, a, 2, 0)*shapeFcn(y, a, 2, 0)];
        
        integral = integral+weight_x*weight_y*transpose(K_mat)*bend_mat_shim*K_mat*dx*dy;
        
    end
end

% local function to determine the shape functions
function n = shapeFcn(x, a, c, d)

switch c
    case 1
        switch d
            case 0
                n = 1/2-3*x/2/a+2*x^3/a^3;
            case 1
                n = 6*x^2/a^3-3/2/a;
            case 2
                n = 12*x/a^3;
        end
    case 2
        switch d
            case 0
                n = 1/2+3*x/2/a-2*x^3/a^3;
            case 1
                n = 3/2/a-6*x^2/a^3;
            case 2
                n = -12*x/a^3;
        end
end

end