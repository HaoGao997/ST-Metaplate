function MeshParam = getEleNodeRel(MeshParam, node, element, electrode, cellsize)
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
node_dof = 1;
num_node_ele = size(element,2)-1;
num_electrode = size(electrode,1);

MeshParam.node_dof = node_dof;
MeshParam.num_node_ele = num_node_ele;
MeshParam.num_node = length(node);
MeshParam.ele_dof = node_dof*num_node_ele;
MeshParam.num_dof = node_dof*MeshParam.num_node;
MeshParam.num_ele = length(element);

% Determines the element size in x and y direction
% Determine which piezoelectric patch this element belongs to
% (This is ONLY for square element)
ele_size = zeros(MeshParam.num_ele,2);
piezo_index = zeros(MeshParam.num_ele,1);
for i_ele = 1:MeshParam.num_ele
    
   pos_node_1 = 0.01*node(element(i_ele,2),2:3); % Careful for the unit!!
   pos_node_4 = 0.01*node(element(i_ele,4),2:3);
   ele_size(i_ele,:) = abs(pos_node_4-pos_node_1);
   
   center_pos = (pos_node_1+pos_node_4)/2;
   
   for i_electrode = 1:num_electrode
       if abs(center_pos(1)-electrode(i_electrode,2)) < cellsize/2
           if abs(center_pos(2)-electrode(i_electrode,3)) < cellsize/2
               piezo_index(i_ele) = electrode(i_electrode,1);
               break
           end
       end
   end
   
end

% Assign the node and element geometry information
node(:,2:3) = node(:,2:3)*1e-2;
MeshParam.node = node;
MeshParam.element = [element,ele_size,piezo_index];

end