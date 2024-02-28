function MeshParam = getEleNodeRel(MeshParam, node, element)
% getEleNodeRel defines the geometrical relations between the elements
% and the nodes.
% In this particular problem, we use identical rectangular quadrilaterial
% elements with 4 nodes and 12 DOF.
%
% Created by Hao Gao (SJTU)
% Create on Feb 10, 2024
% Modified on Feb 28, 2024
% -------------------------------------------------------------------------

% Metrics of the mesh
num_dof_node = 1;
num_node_ele = size(element,2)-1;

MeshParam.num_dof_node = num_dof_node;
MeshParam.num_node_ele = num_node_ele;
MeshParam.num_node = length(node);
MeshParam.num_dof = num_dof_node*MeshParam.num_node;
MeshParam.num_ele = length(element);

% Determines the element size in x and y direction 
% (This is ONLY for square element)
ele_size = zeros(MeshParam.num_ele,2);
for i_ele = 1:MeshParam.num_ele  
   pos_node_1 = node(element(i_ele,2),2:3);
   pos_node_4 = node(element(i_ele,4),2:3);
   ele_size(i_ele,:) = abs(pos_node_4-pos_node_1)*1e-2;
end

% Assign the node and element geometry information
node(:,2:3) = node(:,2:3)*1e-2;
MeshParam.node = node;
MeshParam.element = [element,ele_size];

end