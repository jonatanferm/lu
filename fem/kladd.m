%% Material Coefficients
Esmd= 105*10^9;  Epcb= 105*10^9; Esol= 50*10^9; 
nysmd= 0.118;    nypcb= 0.136;   nysol= 0.36;
ksmd= 0.290;       kpcb= 1.059;    ksol= 66.80; 
rosmd= 1850;     ropcb= 1850;    rosol= 7265;
csmd= 950;       cpcb= 950;      csol= 210;
alphasmd= 1.2*10^-5; alphapcb= 2*10^-5; alphasol= 1.2*10^-5;

T0= 30+273.15; %Converted from Celsius to degrees Kelvin
Tinf = 20+273.15; %Degrees Kelvin
ac= 40;
qel= 9*10^3;

%% Export mesh from PDE-tool
%Imports the mesh created by PDE-tool
load('e.mat'); load('p.mat'); load('t.mat');
figure %Displays the mesh
pdemesh(p,e,t)

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
bcSym=find(abs(coord(:,1)-0)<1e-6); 
%Same as above, but for y, gives lower boundary L_1:
bcBottom=find(abs(coord(:,2)-0)<1e-6); 
%Index for nodes with x~1, i.e. right boundary L_2:
bcRight=find(abs(coord(:,1)-0.001)<1e-6);
%Index for nodes with y~0.6, i.e. upper boundary L_3:
bcTop=find(abs(coord(:,2)-0.0006)<1e-6);
%All nodes on the line y=-x+1.2:
bcConv=find(abs(coord(:,1)+coord(:,2)-0.0012)<1e-6); 
bcFlux =[];

%Loop to divide L_3 into L_3 flux and L_3 convection
for i=1:length(bcTop)
    if (coord(bcTop(i),1) < 0.0002-1e-6)
        bcFlux = [bcFlux; bcTop(i)]; %Add node number (bcTop(i) to bcFlux
    elseif (coord(bcTop(i),1) > 0.0002+1e-6)
        bcConv = [bcConv; bcTop(i)]; %Nodes that don't belong to bcFlux belong
        %to L_3 convection, except th node at x=0.2 which belongs to both
    elseif (abs(coord(bcTop(i),1)-0.0002)<1e-6)
        bcFlux = [bcFlux; bcTop(i)];
        bcConv = [bcConv; bcTop(i)];
    end
end

%Displays nodes belonging to the different boundaries
%Note that all corner nodes are double, but only shown in one color
figure
plot(coord(bcSym,1),coord(bcSym,2),'b*')
hold on
plot(coord(bcBottom,1),coord(bcBottom,2),'r*')
plot(coord(bcRight,1),coord(bcRight,2),'g*')
plot(coord(bcFlux,1),coord(bcFlux,2),'k*')
plot(coord(bcConv,1),coord(bcConv,2),'m*')
legend('L_5 and L_6 (symmetry line): insulated, u_x=0', ...
    'L_1: insulated, u_y=0', 'L_2: insulated, u_x=0', 'L_{3 flux}', ...
    'L_{3 convection} and L_4')

%% Compute global stiffness matrix (without convection terms)
b=1; %The thickness of the structure bears no relevance, set to 1 everywhere
ndep=edof(:,2:4); %Takes out only nodal dependancies from edof
K = zeros(ndof); %Create empty global stiffness matrix
for i=1:nelm %Loop over all 712 triangular elements
    if t(4,i) == 1
        %If the subdomain nr (found in 4th column of t-matrix) is 1
        k = kpcb; %Then the current subdomain is the PCB
    elseif t(4,i) == 2
        k = ksmd; %Subdomain number 2 is the SMD
    else
        k = ksol; %Subdomain number 3 is the solder
    end
    D = k*eye(2); %Isotropic material
    Ke = flw2te(Ex(i,:), Ey(i,:), b, D); %Q=0 in all subdomains
    %Assemble the stiffness matrix:
    K(ndep(i,:),ndep(i,:)) = K(ndep(i,:),ndep(i,:)) + Ke; 
end

%% Convection Boundary (using edges)
ConvEdges = [];

for i=1:length(e)
    if (e(5,i) == 4)
        ConvEdges = [ConvEdges; e(1,i) e(2,i)];
        
    elseif ( (e(5,i) == 3) && (coord(e(1,i),1) >= 0.0002-1e-6) ...
            && (coord(e(2,i),1) >= 0.0002-1e-6) )
        ConvEdges = [ConvEdges; e(1,i) e(2,i)];
    end
end


%Stiffness matrix for convection terms
startNodes = ConvEdges(:,1); %Nodes in the left column of ConvEdges
endNodes = ConvEdges(:,2); %Nodes in the right column of ConvEdges

fc = zeros(length(K),1); %Creates empty column vector for storage
%of boundary convection terms

for i= 1:length(startNodes) %length(startNodes)=length(endNodes) ...
%= # edges in ConvEdges
    %Current edge length = sqrt((startx-endx)^2 + (starty-endy)^2)
    ELength = sqrt( (coord(startNodes(i),1)-coord(endNodes(i),1))^2 ...
    + (coord(startNodes(i),2)-coord(endNodes(i),2))^2);
    Kce = ELength*(ac/6)*[2 1; 1 2];
    %Assemble the stiffness matrix:
    K(ConvEdges(i,:), ConvEdges(i,:)) = K(ConvEdges(i,:), ConvEdges(i,:)) + Kce;
    
    fce = ac*Tinf*(ELength/2)*[1; 1]; %Convection terms along current edge
    fc(ConvEdges(i,:)) = fc(ConvEdges(i,:)) + fce; %Assemble convection vector
end

%% Boundary 'force' vector, temperature
FluxEdges = [];

for i=1:length(e)
    if ( (e(5,i) == 3) && (coord(e(1,i),1) <= 0.0002+1e-6)...
            && (coord(e(2,i),1) <= 0.0002+1e-6) )
        FluxEdges = [FluxEdges; e(1,i) e(2,i)];
    end
end

%Creates empty column vector for storage of boundary 'forces':
fb = zeros(length(K),1); 
startNodes = FluxEdges(:,1); %Nodes in the left column of FluxEdges
endNodes = FluxEdges(:,2); %Nodes in the right column of FluxEdges

for i= 1:length(startNodes) %length(startNodes)=length(endNodes) ...
%= # edges in FluxEdges
    %Current edge length = sqrt((startx-endx)^2 + (starty-endy)^2)
    ELength = sqrt( (coord(startNodes(i),1)-coord(endNodes(i),1))^2...
    + (coord(startNodes(i),2)-coord(endNodes(i),2))^2);
    fbe = qel*(ELength/2)*[1; 1]; %Flux terms along current edge
    %Assemble flux vector:
    fb(FluxEdges(i,:)) = fb(FluxEdges(i,:)) + fbe; 
end

%% Transients
C = zeros(ndof); %Create empty global C-matrix
for i=1:nelm %Loop over all 712 triangular elements
    if t(4,i) == 1 
        %If the subdomain nr (found in 4th column of t-matrix) is 1:
        roc = ropcb*cpcb; %Then the current subdomain is the PCB
    elseif t(4,i) == 2
        roc = rosmd*csmd; %Subdomain number 2 is the SMD
    else
        roc = rosol*csol; %Subdomain number 3 is the solder
    end
    Ce = plantml(Ex(i,:), Ey(i,:), roc); 
    C(ndep(i,:),ndep(i,:)) = C(ndep(i,:),ndep(i,:)) + Ce; %Assemble C-matrix
end

%% Time integration
f = fc + fb; %Global force vector, temperature
%Numerical integration using implicit Euler method:
a = eulerFEMint(C, K, f, T0, 0, 100, 1000); 

%% Thermal distribution plots
ed1 = extract(edof,a(:,1));
figure
fill(Ex',Ey',ed1');
c= colorbar
label = ylabel(c,'\circC');
set(label,'Rotation',0);
title('Temperature distribution after 0.1 s');
xlabel('x (m)');
ylabel('y (m)');

ed2 = extract(edof,a(:,250));
figure
fill(Ex',Ey',ed2');
c= colorbar
label = ylabel(c,'\circC');
set(label,'Rotation',0);
title('Temperature distribution after 25 s');
xlabel('x (m)');
ylabel('y (m)');

ed3 = extract(edof,a(:,1000));
figure
fill(Ex',Ey',ed3');
c= colorbar
label = ylabel(c,'\circC');
set(label,'Rotation',0);
title('Temperature distribution after 100 s');
xlabel('x (m)');
ylabel('y (m)');

%Stationary temperature distribution
statT = solveq(K,f);
deltaT = statT-T0; %Temperature deviation from T_0 to stationary ...
%disitribution in each node
figure
fill(Ex',Ey',edstat');
c= colorbar
label = ylabel(c,'\circC');
set(label,'Rotation',0);
title('Stationary temperature distribution');
xlabel('x (m)');
ylabel('y (m)');


%% Part 2: Plane elasticity
%Stationary temperature distribution
statT = solveq(K,f);
deltaT = statT-T0; %Temperature deviation from T_0 to stationary ...
%disitribution in each node

ndofs = 2*ndof; %Degrees of freedom in solid mechanics is two per node

%Each node, i. corresponds to the degrees of freedom 2i-1 och 2i
%The new edof-matrix has dimensions [nelm,7] and is called edofs
edofs(:,1)=1:nelm;  %First column in edofs contains the elemtent indices
edofs(:,[2 4 6])=(2*t(1:3,:)-1)'; %Columns 2, 4 & 6 are 2*node index-1 ...
%for the current element
edofs(:,[3 5 7])=2*t(1:3,:)'; %Columns 3, 5 & 7 are 2*node index

%% Compute global temperature force vector:
b=1; %The thickness of the structure bears no relevance, set to 1
ep = [2, b]; %2: plane strain, b: element thickness
ndeps=edofs(:,2:7); %Takes out only nodal dependencies from edofs
fdeltaT = zeros(ndofs, 1); %Create empty global thermal strain vector
elemDeltaT = zeros(nelm,1); %Create empty vector for storing the average
%temperature deviation in each element

%Define elasticity module for the three subdomains:
Dpcb = Epcb/((nypcb+1)*(1-2*nypcb))*[1-nypcb nypcb 0;
                                     nypcb 1-nypcb 0;
                                     0 0 1/2*(1-2*nypcb)]; %Isotropic material
Dsmd = Esmd/((nysmd+1)*(1-2*nysmd))*[1-nysmd nysmd 0;nysmd 1-nysmd 0;...
    0 0 1/2*(1-2*nysmd)]; %Isotropic material
Dsol = Esol/((nysol+1)*(1-2*nysol))*[1-nysol nysol 0;nysol 1-nysol 0;...
    0 0 1/2*(1-2*nysol)]; %Isotropic material

for i=1:nelm %Loop over all 712 triangular elements
    if t(4,i) == 1 
        %If the subdomain nr (found in 4th column of t-matrix) is 1
        D = Dpcb; %Then the current subdomain is the PCB
        alpha = alphapcb;
    elseif t(4,i) == 2
        D = Dsmd; %Subdomain number 2 is the SMD
        alpha = alphasmd;
    else
        D = Dsol; %Subdomain number 3 is the solder
        alpha = alphasol;
    end 
    %Compute the average temperature deviation in the current element:
    elemDeltaT(i) = (deltaT(t(1,i)) + deltaT(t(2,i)) + deltaT(t(3,i)))/3;
    %Compute the contribution of the element to fdeltaT
    fdeltaTe = elemDeltaT(i)*plantf(Ex(i,:), Ey(i,:), ep, (alpha*D*[1; 1; 0])');
    %Assemble the global temperature strain vector:
    fdeltaT(ndeps(i,:)) = fdeltaT(ndeps(i,:)) + fdeltaTe;
end

%% Compute stiffness matrix for plane strain
%we call it Ks to differentiate it from K for temp.

Ks = zeros(ndofs);
for i=1:nelm %Loop over all 712 triangular elements
    if t(4,i) == 1 
        %If the subdomain nr (found in 4th column of t-matrix) is 1
        D = Dpcb; %Then the current subdomain is the PCB
    elseif t(4,i) == 2
        D = Dsmd; %Subdomain number 2 is the SMD
    else
        D = Dsol; %Subdomain number 3 is the solder
    end
    Ke = plante(Ex(i,:), Ey(i,:), ep, D); %Elementwise contribution
    %Assemble the stiffness matrix:
    Ks(ndeps(i,:),ndeps(i,:)) = Ks(ndeps(i,:),ndeps(i,:)) + Ke;
end

%% Boundary conditions
%Nodes where u_x=0 are in bcRight and bcSym
%On bcBottom u_y=0
bc = [2*bcRight-1; 2*bcSym-1; 2*bcBottom];
bc(:,2) = 0;
%On remaining boundaries the traction vector is 0

%%Calculate solution with solveq
[u,Q] = solveq(Ks,fdeltaT,bc);

ed = extract(edofs,u);
figure
eldraw2(Ex,Ey,[1,4,3]);
hold on
%Displays the element distortions with a factor 100 magnification:
eldisp2(Ex,Ey,ed,[1 5 3],100);
title('Deformed mesh, scale factor 100')
legend('Original mesh','Deformed mesh')


%%
es = zeros(nelm,3);
et = zeros(nelm,3);
sigmazz = zeros(nelm,1);
for i=1:nelm %Loop over all 712 triangular elements
    if t(4,i) == 1 
        %If the subdomain nr (found in 4th column of t-matrix)=1
        D = Dpcb; %Then the current subdomain is the PCB
        E = Epcb;
        ny = nypcb;
        alpha = alphapcb;
    elseif t(4,i) == 2
        D = Dsmd; %Subdomain number 2 is the SMD
        E = Esmd;
        ny = nysmd;
        alpha = alphasmd;
    else
        D = Dsol; %Subdomain number 3 is the solder
        E = Esol;
        ny = nysol;
        alpha = alphasol;
    end
    %Compute the stresses and strains of the element:
    [es(i,:), et(i,:)] = plants(Ex(i,:), Ey(i,:), ep, D, ed(i,:));
    %Take thermoelasticity into account (not done automatically by plants):
    es(i,:) = es(i,:)-(alpha*E*elemDeltaT(i))/(1-2*ny)*[1 1 0];
    %Sigmazz = ny(sigmaxx+sigmayy)-alpha*E*deltaT according to (13.42)
    sigmazz(i) = ny*(es(i,1)+es(i,2))-alpha*E*elemDeltaT(i);
end

%% Compute von Mises effective stress for each element
Seff_el = zeros(nelm, 1);
Seff_nod = zeros(length(coord),1);

for i=1:nelm %Loop over all 712 triangular elements
    %Find the effective von Mises stress for each element:
    Seff_el(i) = sqrt( es(i,1)^2 + es(i,2)^2+sigmazz(i)^2-es(i,1)*es(i,2)...
        -es(i,1)*sigmazz(i)-es(i,2)*sigmazz(i)+3*es(i,3)^2); 
end

for i=1:size(coord,1) %Loop over all nodes
[c0,c1] = find(edof(:,2:4) == i);
Seff_nod(i,1) = sum(Seff_el(c0))/size(c0,1);
end