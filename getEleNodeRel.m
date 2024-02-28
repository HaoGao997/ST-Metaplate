function MeshParam = getEleNodeRel(MeshParam, node, element)
% getEleNodeRel defines the geometrical relations between the elements
% and the nodes.
% In this particular problem, we use identical rectangular quadrilaterial
% elements with 4 nodes and 12 DOF.
%
% Created by Hao Gao (SJTU)
% Create on Feb 10, 2024
% Modified on Feb 10, 2024
% -------------------------------------------------------------------------

MeshParam.node = node;
MeshParam.element = element;

% Metrics of the mesh
num_dof_node = 3;
num_node_ele = 4;

MeshParam.num_dof_node = num_dof_node;
MeshParam.num_node_ele = num_node_ele;
MeshParam.num_node = length(node);
MeshParam.num_dof = num_dof_node*MeshParam.num_node;
MeshParam.num_ele = length(element);

end