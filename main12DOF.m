% Transient wave propagation analysis of a spatiotemporally modulated
% piezoelectric shunt meta-plate.
% 
% The 4 node 12 dof quadrilaterial element is used in the formulation.
%
% -------------------------------------------------------------------------
% Created by Hao Gao (SJTU)
% Create on Mar 13, 2024
% Modified on Mar 13, 2024
% -------------------------------------------------------------------------
clear; clc;

% Load the mesh information (element/node numbering and coordinates)
load('BeamMeshData.mat');
load('BeamElectrodeData.mat');

% Input the parameters of the shim plate
PlateParam.density = 2700;       % Density [kg/m^3]
PlateParam.modulus = 69e9;       % Elastic modulus [Pa]
PlateParam.poisson = 0.33;       % Poisson's ratio
PlateParam.thickness = 0.002;    % Thickness [m]
PlateParam.cellsize = 0.05;      % Size [m]
dim = 1;
switch dim
    case 1
        num_cell_edge = electrode(end,1);
        PlateParam.num_cell = num_cell_edge;  % Total number of cells
    case 2
        num_cell_edge = sqrt(electrode(end,1));
        PlateParam.num_cell = num_cell_edge^2;  % Total number of cells
end

PlateParam.length_plate = PlateParam.cellsize*num_cell_edge;
% Define the elasticity matrix of the shim plate
PlateParam.elas_mat = PlateParam.modulus/(1-PlateParam.poisson^2)*...
    [1, PlateParam.poisson, 0;
     PlateParam.poisson, 1, 0;
     0, 0, (1-PlateParam.poisson)/2];

% Input the parameters of the piezoelectric layer
PiezoParam.density = 7500;       % Density [kg/m^3]
PiezoParam.modulus = 126e9;      % Elastic modulus [Pa]
PiezoParam.poisson = 0.3;        % Poisson's ratio
PiezoParam.permitv = 17.3e-9;    % Permittivity [F/m]
PiezoParam.piezo_strain = -23.4; % Piezoelectric strain constant [C/m^2]
PiezoParam.thickness = 0.0005;   % Thickness [m]
PiezoParam.num_patch = PlateParam.num_cell;
PiezoParam.capacit_patch = 2*PiezoParam.permitv/PiezoParam.thickness;
% Define the elasticity matrix of the piezoelectric layer
PiezoParam.elas_mat = PiezoParam.modulus/(1-PiezoParam.poisson^2)*...
    [1, PiezoParam.poisson, 0;
     PiezoParam.poisson, 1, 0;
     0, 0, (1-PiezoParam.poisson)/2];
% Define the graded phase
phase_unit = [0,pi*2/3,4*pi/3]; % Phase of a unit
PiezoParam.phase_vec = repmat(phase_unit, 1, num_cell_edge/3);
% PiezoParam.phase_vec = zeros(1,PiezoParam.num_patch);  % TBD!!!!
% Define the modulation shunte circuit
PiezoParam.fre_mod = 10;        % Modulation frequency [rad/s]
PiezoParam.capacit_shunt = 5e-5;   % Modulation capacitance [F]

% Define the mesh parameters
MeshParam = getEleNodeRel(node, element, electrode, PlateParam.cellsize, 2);
MeshParam.elesize = MeshParam.element(1,6);

%% Loop on the element assemble the mass and stiffness matrices
% Initialize the sparse global mass and stiffness matrices
mass_mat_asb = sparse(MeshParam.num_dof,MeshParam.num_dof);
stiff_mat_asb = sparse(MeshParam.num_dof,MeshParam.num_dof);
couple_mat_asb = sparse(MeshParam.num_dof,PiezoParam.num_patch);

disp('Assemble the matrices');

for i_ele = 1:MeshParam.num_ele
    
    % Calculate the material properties of the shim plate and piezo
    % layers
    bend_mat_shim = PlateParam.thickness^3/12*PlateParam.elas_mat;
    bend_mat_piezo = 2*(4*PiezoParam.thickness^2+...
        6*PlateParam.thickness*PiezoParam.thickness+...
        3*PlateParam.thickness^2)*PiezoParam.thickness/12*PiezoParam.elas_mat; % set zero tentatively!
    
    mass_area_shim = PlateParam.thickness*PlateParam.density;
    mass_area_piezo = 2*PiezoParam.thickness*PiezoParam.density; % set zero tentatively!

    couple_param_mat = [PiezoParam.piezo_strain;PiezoParam.piezo_strain;0];
    
    % Get the index of the associated piezoelectric patch
    piezo_index = MeshParam.element(i_ele,8);
    % Get the index of nodes of the element
    node_index = MeshParam.element(i_ele,2:5);
    
    % Numerical integration of the mass and stiffness matrices for the element
    [stiff_mat_shim,stiff_mat_piezo,couple_mat] =...
        getStiffMatrixEle(MeshParam.ele_dof, MeshParam.node_dof, MeshParam.elesize, bend_mat_shim, bend_mat_piezo, couple_param_mat, node_index, MeshParam.node);
    
    mass_mat_shim = getMassMatrixEle(MeshParam.ele_dof, MeshParam.node_dof,  MeshParam.elesize, mass_area_shim, node_index, MeshParam.node);
    mass_mat_piezo = getMassMatrixEle(MeshParam.ele_dof, MeshParam.node_dof,  MeshParam.elesize, mass_area_piezo, node_index, MeshParam.node);
    mass_mat = mass_mat_shim+mass_mat_piezo;
    
    [mass_mat_asb,stiff_mat_asb,couple_mat_asb] = assembMat12DOF(mass_mat_asb, stiff_mat_asb, couple_mat_asb, node_index, piezo_index,...
        stiff_mat_shim, stiff_mat_piezo, couple_mat, mass_mat, MeshParam.num_dof);
    
    disp(['Element ',num2str(i_ele)]);
     
end

%% Constraint the excited node
ex_node_index = [362, 603, 844]; %(TBD!!!)
node_index = 1:MeshParam.num_dof;
node_index(:,ex_node_index) = [];
stiff_vec_ex = stiff_mat_asb(:,ex_node_index);
stiff_vec_ex = stiff_vec_ex(node_index,:);
% Removed the excitaed node from mass and stiffness matrices
mass_mat_asb(ex_node_index,:) =[];
mass_mat_asb(:,ex_node_index) =[];
stiff_mat_asb(ex_node_index,:) =[];
stiff_mat_asb(:,ex_node_index) =[];
couple_mat_asb(ex_node_index,:) = [];

%% Time domain march scheme (4th order Runge Kutta)
time_step = 1e-7;
time_sta = 0;
time_end = 0.001;
num_time_pts = round((time_end-time_sta)/time_step)+1;
time_pts = linspace(time_sta,time_end,num_time_pts);

% Define the displacement excitation
fre_ex = 400*2*pi;    % Excitation frequency [rad/s]

% Finds the total degress of freedom of the state vector
num_dof_mech = size(mass_mat_asb,1);
num_dof_electr = PiezoParam.num_patch;
num_states = 2*num_dof_mech+num_dof_electr;
% Initialize the state coordinates
state_coord_old = zeros(num_states,1);
state_coord_new = zeros(num_states,1);

% Define the record step size
num_pts_record = 1e2;
state_coord_record = zeros(num_states,num_pts_record);

% Calculates the inverse matrix of the mass matrix
inv_mass_mat_asb = mass_mat_asb\sparse(eye(size(mass_mat_asb)));

for i_time = 1:num_time_pts-1
    
    fid = floor((i_time-0.1)/num_pts_record);
    vecid = i_time-fid*num_pts_record;

    % Determine the slope k1 at the current time
    t_inst = time_pts(i_time);
    state_coord_inst = state_coord_old;

    % Define the time-variant capacitance matrix and its inverse
    capacit_patch_mat = sparse(diag(PiezoParam.capacit_patch*ones(1,PiezoParam.num_patch)));
    capacit_shunt_mat = sparse(diag(PiezoParam.capacit_shunt*sin(PiezoParam.fre_mod*t_inst+PiezoParam.phase_vec)));
    capacit_mat_inst = capacit_patch_mat+capacit_shunt_mat;
    inv_electr_mat_inst = capacit_mat_inst\sparse(eye(size(capacit_mat_inst)));

    % Define the constraint force as external load
    force_vec = stiff_vec_ex*[1;1;1]*sin(fre_ex*t_inst);

    % Find the slope at the given time
    k_1 = calRK4Slope(state_coord_inst, inv_mass_mat_asb, stiff_mat_asb, couple_mat_asb, inv_electr_mat_inst, force_vec);

    % Determine the slope k2 at the mid step time
    t_inst = time_pts(i_time)+time_step/2;
    state_coord_inst = state_coord_old+k_1*time_step/2;
    
    % Define the time-variant capacitance matrix and its inverse
    capacit_patch_mat = sparse(diag(PiezoParam.capacit_patch*ones(1,PiezoParam.num_patch)));
    capacit_shunt_mat = sparse(diag(PiezoParam.capacit_shunt*sin(PiezoParam.fre_mod*t_inst+PiezoParam.phase_vec)));
    capacit_mat_inst = capacit_patch_mat+capacit_shunt_mat;
    inv_electr_mat_inst = capacit_mat_inst\sparse(eye(size(capacit_mat_inst)));

    % Define the constraint force as external load
    force_vec = stiff_vec_ex*[1;1;1]*sin(fre_ex*t_inst);

    % Find the slope at the given time
    k_2 = calRK4Slope(state_coord_inst, inv_mass_mat_asb, stiff_mat_asb, couple_mat_asb, inv_electr_mat_inst, force_vec);

    % Determine the slope k3 at the mid step time
    t_inst = time_pts(i_time)+time_step/2;
    state_coord_inst = state_coord_old+k_2*time_step/2;

    % Define the time-variant capacitance matrix and its inverse
    capacit_patch_mat = sparse(diag(PiezoParam.capacit_patch*ones(1,PiezoParam.num_patch)));
    capacit_shunt_mat = sparse(diag(PiezoParam.capacit_shunt*sin(PiezoParam.fre_mod*t_inst+PiezoParam.phase_vec)));
    capacit_mat_inst = capacit_patch_mat+capacit_shunt_mat;
    inv_electr_mat_inst = capacit_mat_inst\sparse(eye(size(capacit_mat_inst)));

    % Define the constraint force as external load
    force_vec = stiff_vec_ex*[1;1;1]*sin(fre_ex*t_inst);

    % Find the slope at the given time
    k_3 = calRK4Slope(state_coord_inst, inv_mass_mat_asb, stiff_mat_asb, couple_mat_asb, inv_electr_mat_inst, force_vec);

    % Determine the slope k4 at the mid step time
    t_inst = time_pts(i_time+1);
    state_coord_inst = state_coord_old+k_3*time_step;

    % Define the time-variant capacitance matrix and its inverse
    capacit_patch_mat = sparse(diag(PiezoParam.capacit_patch*ones(1,PiezoParam.num_patch)));
    capacit_shunt_mat = sparse(diag(PiezoParam.capacit_shunt*sin(PiezoParam.fre_mod*t_inst+PiezoParam.phase_vec)));
    capacit_mat_inst = capacit_patch_mat+capacit_shunt_mat;
    inv_electr_mat_inst = capacit_mat_inst\sparse(eye(size(capacit_mat_inst)));

    % Define the constraint force as external load
    force_vec = stiff_vec_ex*[1;1;1]*sin(fre_ex*t_inst);

    % Find the slope at the given time
    k_4 = calRK4Slope(state_coord_inst, inv_mass_mat_asb, stiff_mat_asb, couple_mat_asb, inv_electr_mat_inst, force_vec);

    % March to the next step
    state_coord_new = state_coord_old+time_step/6*(k_1+2*k_2+2*k_3+k_4);
    
    % Save the state data and update the state vector
    state_coord_old = state_coord_new;
    
    state_coord_record(:,vecid) = state_coord_new;
    if vecid == num_pts_record
        fname = ['state_vec_mat_',num2str(fid)];
        save(fname,'state_coord_record');
    else
    end

end

%% Post-process
plot(time_pts(1:end), state_coord(1:num_dof_mech,1:end));




