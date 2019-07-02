clear;
clc;
load('pet.mat')
%load('pet_not_fine.mat')
%load('pet_fine.mat')
%load('pet_extra_fine.mat')
%%
addpath('./calfem/fem');
%% Properties for Al, Steel, Cu, Core
young_mod = [70, 210, 128, 500];
poss_rat = [.33, .3, .36, .45];
exp_coeff = [69e-6, 35e-6, 51e-6, 20e-6];
density = [2710, 7900, 8930, 2000];
spec_heat = [903, 460, 386, 900];
therm_cond = [238, 20, 385, 1.6];

T0 = 30+273.15; %Converted from Celsius to degrees Kelvin
Tinf = 20+273.15; %Degrees Kelvin
ac = 40;
Q = 1e5;

figure %Displays the mesh
pdemesh(p,e,t)
%%
%% Create edof matrix and plot nodes and edges
%The number of columns in t corresponds to the number of triangles in our
%mesh:
nelm=length(t(1,:)); 
edof(:,1)=1:nelm; %The first column in edof contains the element index
%Columns 2, 3 and 4 are the nodes in anti-clockwise order for the
%corresponding element (triangle):
edof(:,2:4)=t(1:3,:)'; 
coord=p'; %Column vector with all nodes (x- and y-coordinate)

%max(t(1:3,:) returns a row vector with the largest number in each column.
ndof=max(max(t(1:3,:))); 
%Gives the element coordinate matrices. Each row in Ex contains
%x-coordinates of one element:
[Ex,Ey]=coordxtr(edof,coord,(1:ndof)',3); 
eldraw2(Ex,Ey,[1,2,0]);

%% Boundaries
%Index of the nodes which x-coordinate has less distance to 0 than 1e-3.
%This gives nodes on the symmetry line, L_5 & L_6:
b_fix_y = find(abs(coord(:,2)-0)<1e-6 | abs(coord(:,2) - 0.05) < 1e-6);
%Same as above, but for y, gives lower boundary L_1:
b_fix_x = find((abs(coord(:,1)-0)<1e-6 & abs(coord(:,2) < 0.015)) | ...
               (abs(coord(:,1)-0.025)<1e-6 & abs(coord(:,2) < 0.0325)));
b_conv = find((abs(coord(:,1)-0)<1e-6 & abs(coord(:,2) > 0.015)) | ...
               (abs(coord(:,1)-0.025)<1e-6 & abs(coord(:,2) > 0.0325)));
%%
figure
plot(coord(b_fix_y,1),coord(b_fix_y,2),'rx')
hold on
plot(coord(b_fix_x,1),coord(b_fix_x,2),'bx')
plot(coord(b_conv,1),coord(b_conv,2),'gx')
legend({'fixed y', 'fixed x', 'not fixedm convection'})
%%
%% Compute global stiffness matrix (without convection terms)
b=1; %The thickness of the structure bears no relevance, set to 1 everywhere
ndep=edof(:,2:4); %Takes out only nodal dependancies from edof
K = zeros(ndof); %Create empty global stiffness matrix
for i=1:nelm %Loop over all 712 triangular elements
    k = therm_cond(t(4,i));
    D = k*eye(2); %Isotropic material
    Ke = flw2te(Ex(i,:), Ey(i,:), b, D); %Q=0 in all subdomains
    %Assemble the stiffness matrix:
    K(ndep(i,:),ndep(i,:)) = K(ndep(i,:),ndep(i,:)) + Ke; 
end
%% Convection Boundary
for k = b_conv
    elm = edof(k,2:end);
    boundries = boundary_element_type(elm) == 1;
    if sum(boundries) == 2
        if boundries(1) == 0
            NtN = [0 0 0; 0 2 1; 0 1 2];
        elseif boundries(2) == 0
            NtN = [2 0 1; 0 0 0; 1 0 2;];
        else
            NtN = [2 1 0; 1 2 0; 0 0 0];
        end
        indx = edof(k,2:end);
        alpha = 100;
        L = max([coord(elm(1),2), coord(elm(2),2), coord(elm(3),2)]) - ...
            min([coord(elm(1),2), coord(elm(2),2), coord(elm(3),2)]);
        Kb = ep*alpha*L*NtN/6;
        Fb = ep*alpha*Tinf*L*boundries'/2;
        K(indx,indx) = K(indx,indx)+Kb;
        f(indx) = f(indx)+Fb;
    end
end

