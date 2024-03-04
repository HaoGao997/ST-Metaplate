% Transient wave propagation analysis of a spatiotemporally modulated
% piezoelectric shunt meta-plate.
%
% -------------------------------------------------------------------------
% Created by Hao Gao (SJTU)
% Create on Feb 10, 2024
% Modified on Feb 29, 2024
% -------------------------------------------------------------------------
clear; clc;

% Input the parameters of the shim plate
PlateParam.density = 2700;       % Density [kg/m^3]
PlateParam.modulus = 69e9;       % Elastic modulus [Pa]
PlateParam.poisson = 0.33;       % Poisson's ratio
PlateParam.thickness = 0.002;    % Thickness [m]
PlateParam.cellsize = 0.05;      % Size [m]
num_cell_edge = 120;
PlateParam.num_cell = num_cell_edge^2;  % Total number of cells
PlateParam.length_plate = PlateParam.cellsize*numcells_edge;

% Input the parameters of the piezoelectric layer
PiezoParam.density = 7500;       % Density [kg/m^3]
PiezoParam.modulus = 126e9;      % Elastic modulus [Pa]
PiezoParam.permitv = 17.3e-9;    % Permittivity [F/m]
PiezoParam.piezo_strain = -23.4; % Piezoelectric strain constant [C/m^2]
PiezoParam.num_patche = PlateParam.num_cell;

% Define the mesh size (regular square element)
MeshParam.elesize = 0.01;

% Load the mesh information (element/node numbering and coordinates)
load('PlateMeshData.mat');
load('ElectrodeData.mat');
MeshParam = getEleNodeRel(MeshParam, node, element, electrode, PlateParam.cellsize);

% Loop on the element assemble the mass and stiffness matrices
% Initialize the global mass and stiffness matrices
mass_mat_asb = zeros(MeshParam.num_dof);
stiff_mat_asb = zeros(MeshParam.num_dof);
couple_mat_asb = zeros(MeshParam.num_dof,PiezoParam.num_piezo);

for i_ele = 1:MeshParam.num_ele
    
    % Calculate the material properties of the shim plate and piezo
    % layers
    bend_mat_shim = PlateParam.thickness^3/12*PlateParam.elas_mat;
    bend_mat_piezo = (4*PiezoParam.thickness^2+...
        6*PlateParam.thickness*PiezoParam.thickness+...
        3*PlateParam.thickness^2)/12*PiezoParam.elas_mat;
    
    mass_area_shim = PlateParam.thickness*PlateParam.density;
    mass_area_piezo = PiezoParam.thickness*PiezoParam.density;
    
    % Get the index of the associated piezoelectric patch
    piezo_index = MeshParam.element(i_ele,8);
    
    % Numerical integration of the mass and stiffness matrices for the element
    [stiff_mat_shim,stiff_mat_piezo,couple_mat] =...
        getStiffMatrixEle(MeshParam.ele_dof, MeshParam.element(i_ele,6:7), bend_mat_shim, bend_mat_piezo, couple_param_mat);
    
    mass_mat_shim = getMassMatrixEle(MeshParam.ele_dof, MeshParam.element(i_ele,6:7), mass_area_shim);
    mass_mat_piezo = getMassMatrixEle(MeshParam.ele_dof, MeshParam.element(i_ele,6:7), mass_area_piezo);
    mass_mat = mass_mat_shim+mass_mat_piezo;
    
    [mass_mat_asb,stiff_mat_asb,couple_mat_asb] = assembMat(mass_mat_asb, stiff_mat_asb, couple_mat_asb, node_index, piezo_index,...
        stiff_mat_shim, stiff_mat_piezo, couple_mat, mass_mat);
     
end

% Constraint the excited node
ex_node_index = 10; %(TBD!!!)
stiff_vec_ex = stiff_mat_asb(:,ex_node_index);
stiff_vec_ex(ex_node_index) = [];
% Removed the excitaed node from mass and stiffness matrices
mass_mat_asb(ex_node_index,:) =[];
mass_mat_asb(:,ex_node_index) =[];
stiff_mat_asb(ex_node_index,:) =[];
stiff_mat_asb(:,ex_node_index) =[];
couple_mat_asb(ex_node_index,:) = [];

% Time domain march scheme (4th order Runge Kutta)
time_step = 1e-5;
time_sta = 0;
time_end = 20;
num_time_pts = round((time_end-time_sta)/time_step)+1;
time_pts = linspace(time_sta,time_end,num_time_pts);


% Finds the total degress of freedom of the state vector
num_states = 2*MeshParam.num_dof+PiezoParam.num_patch;
% Initialize the state coordinates
state_coord = zeros(num_states, num_time_pts);

% Calculates the inverse matrix of the mass matrix
inv_mass_mat_asb = b*inv(a);

for i_time = 1:num_time_pts
    
end






