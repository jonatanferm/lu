clear;
clc;
%load('pet.mat')
%load('pet_not_fine.mat')
load('pet_fine.mat')
%load('pet_extra_fine.mat')
%%
addpath('./calfem/fem');
%% Properties for Al, Steel, Cu, Core
young_mod = [70, 210, 128, 500]*1e9;
poss_rat = [.33, .3, .36, .45];
exp_coeff = [69e-6, 35e-6, 51e-6, 20e-6];
density = [2710, 7900, 8930, 2000];
spec_heat = [903, 460, 386, 900];
therm_cond = [238, 20, 385, 1.6];
%%
[nelm, edof, nnod, dof, ndof, coord, Ex, Ey, element_map] = extract_params(p, e, t);
%%
boundary_element_type = find_boundaries(nelm, nnod, coord);
%%
T0 = 273.15;
Tinf = T0+15;
Q = (1e5) * (element_map == 4);
[a, ~, ~, ~, ~] = heat(Tinf, nnod, nelm, ndof, edof, element_map, ...
    therm_cond, spec_heat, density, boundary_element_type, coord, Q, Ex, Ey);
eT1 = (a(edof(:,2)) + a(edof(:,3)) + a(edof(:,4)))/3 - T0;
Tdelta15_node = a-15;
Tdelta15 = eT1 - 15;
% increased current
Qincreased = (1.6^2) * Q;
[a, ~, ~, ~, ~] = heat(Tinf, nnod, nelm, ndof, edof, element_map, ...
    therm_cond, spec_heat, density, boundary_element_type, coord, Qincreased, Ex, Ey);
%%
eT2 = (a(edof(:,2)) + a(edof(:,3)) + a(edof(:,4)))/3 - T0;
%%
f = figure('visible','off');
fill(Ex',Ey',eT1','EdgeColor','none')
title('Stationary temperature distribution original current. [C]')
%colormap(gray);
colormap(hot)
colorbar;
xlabel('x-position [m]')
ylabel('y-position [m]')
axis equal
saveas(f,'temp_stat_og','png')
%%
f = figure('visible','off');
fill(Ex',Ey',eT2','EdgeColor','none')
title('Stationary temperature distribution increased current. [C]')
%colormap(gray);
colormap(hot)
colorbar;
xlabel('x-position [m]')
ylabel('y-position [m]')
axis equal
saveas(f,'temp_stat_inc','png')
%%
%{
Tinf = 25 + T0;
[~, ~, K, f, C] = heat(Tinf, nnod, nelm, ndof, edof, element_map, ...
    therm_cond, spec_heat, density, boundary_element_type, coord, Q, Ex, Ey);
d0 = Tinf * ones(ndof, 1); %Initial temp
h = 60*60;
ip = [10*60, 6*h, 0.5, [3, 0, 2*h, 4*h, 6*h]];
Tsnap = step1(K, C, d0, ip, f, []);
%%
tt = 1;
eT = (Tsnap(edof(:,2),tt) + Tsnap(edof(:,3), tt) + Tsnap(edof(:,4), tt))/3;
eT = eT-T0;
f = figure;%('visible','off');
fill(Ex',Ey',eT','EdgeColor','none')
title('Temperature distribution at 2 hours [C]')
colormap(hot);
colorbar;
xlabel('x-position [m]')
ylabel('y-position [m]')
axis equal
saveas(f,'temp_2h','png')
%%
tt = 2;
eT = (Tsnap(edof(:,2),tt) + Tsnap(edof(:,3), tt) + Tsnap(edof(:,4), tt))/3;
eT = eT-T0;
f = figure%('visible','off');
fill(Ex',Ey',eT','EdgeColor','none')
title('Temperature distribution at 4 hours [C]')
colormap(hot);
colorbar;
xlabel('x-position [m]')
ylabel('y-position [m]')
axis equal
saveas(f,'temp_4h','png')
%%
tt = 3;
eT = (Tsnap(edof(:,2),tt) + Tsnap(edof(:,3), tt) + Tsnap(edof(:,4), tt))/3;
eT = eT-T0;
f = figure%('visible','off');
fill(Ex',Ey',eT','EdgeColor','none')
title('Temperature distribution at 6 hours [C]')
colormap(hot);
colorbar;
xlabel('x-position [m]')
ylabel('y-position [m]')
axis equal
saveas(f,'temp_6h','png')
%}
%% Find matrices for elasticity formulation
%[S15, ed15, Ks1] = calc_stress(young_mod, poss_rat, exp_coeff, element_map, nnod, nelm, ndof, edof, Ex, Ey, Tdelta15, boundary_element_type);
Q = (1e5) * (element_map == 4);
edof_elast = zeros(nelm, 7);
edof_elast(:,1) = (1:nelm)';
edof_elast(:,[2 4 6]) = (2*t(1:3,:)-1)';
edof_elast(:,[3 5 7]) = 2*t(1:3,:)';
%edof_elast = [edof(:,1) edof(:,2) edof(:,2)+ndof edof(:,3) edof(:,3)+ndof edof(:,4) edof(:,4)+ndof]; 

[a, ~, ~, ~, ~] = heat(T0+25, nnod, nelm, ndof, edof, element_map, ...
    therm_cond, spec_heat, density, boundary_element_type, coord, Q, Ex, Ey);

eT1 = (a(edof(:,2)) + a(edof(:,3)) + a(edof(:,4)))/3 - T0;
Tdelta25 = eT1 - 25;
Tdelta25_node = a - 25;
Tdelta25_node(edof(:,2:end))
[Snod25, Sel25, ed25, Ks2] = calc_stress(young_mod, poss_rat, exp_coeff, element_map, nnod, nelm, ndof*2, edof_elast, Ex, Ey, Tdelta25,Tdelta25_node(edof(:,2:end)), boundary_element_type);


[Snod15, Sel15, ed15, Ks1] = calc_stress(young_mod, poss_rat, exp_coeff, element_map, nnod, nelm, ndof*2, edof_elast, Ex, Ey, Tdelta15,Tdelta15_node(edof(:,2:end)), boundary_element_type);

max(Snod25)
max(Sel25)
max(ed25)


% plot stress
% 15
f = figure('visible','off');
exd = Ex + ed15(:,1:2:end);
eyd = Ey + ed15(:,2:2:end);
fill(exd',eyd',Sel15,'EdgeColor','none')
title('Stress distribution')
%colormap(gray);
colormap(hot);
colorbar;
xlabel('x-position [m]')
ylabel('y-position [m]')
axis equal
saveas(f,'stress_15C','png')

% 25
f = figure('visible','off');
exd = Ex + ed25(:,1:2:end);
eyd = Ey + ed25(:,2:2:end);
fill(exd',eyd',Sel25,'EdgeColor','none')
title('Stress distribution')
%colormap(gray);
colormap(hot);
colorbar;
xlabel('x-position [m]')
ylabel('y-position [m]')
axis equal
saveas(f,'stress_25C','png')
%

% Calculate displaced coordinates
mag = 100; % Magnification (due to small deformations)
exd = Ex + mag*ed15(:,1:2:end);
eyd = Ey + mag*ed15(:,2:2:end);
f = figure()
colormap(jet)
fill(Ex',Ey',[0 0 0],'EdgeColor',[224,130,20]/255,'FaceAlpha',0, 'LineWidth', 1.3)
hold on
fill(exd',eyd',[0 0 0],'EdgeColor',[45,0,75]/255,'FaceAlpha',0)
axis equal
title('Displacement field for T0 = 15C [Magnitude enhancement 100]')
saveas(f,'displacement15','png')
% Calculate displaced coordinates
mag = 100; % Magnification (due to small deformations)
exd = Ex + mag*ed25(:,1:2:end);
eyd = Ey + mag*ed25(:,2:2:end);
f = figure()
colormap(jet)
fill(Ex',Ey',[0 0 0],'EdgeColor',[224,130,20]/255,'FaceAlpha',0, 'LineWidth', 1.3)
hold on
fill(exd',eyd',[0 0 0],'EdgeColor',[45,0,75]/255,'FaceAlpha',0)
axis equal
title('Displacement field for T0 = 25C [Magnitude enhancement 100]')
saveas(f,'displacement25','png')
%


function [Seff_nod, Seff_el, ed, Ks] = calc_stress(young_mod, poss_rat, exp_coeff, element_map, nnod, nelm, ndof, edof, Ex, Ey, Tdelta,Tdelta_node, boundary_element_type)
    
    thermal_stress = (young_mod(element_map).*exp_coeff(element_map) ...
        .*Tdelta')./(1-poss_rat(element_map)); %Thermal strain (13.34)
    
    size(Tdelta);
    
    Ks = spalloc(ndof, ndof, numel(edof));
    Fs = zeros(ndof, 1);
    ep = [2, 1];
    for i = 1:nelm
        indx = edof(i,2:end);
        %indx = [indx_t*2-1, indx_t*2];
        E = young_mod(element_map(i));
        v = poss_rat(element_map(i));
        alpha = exp_coeff(element_map(i));
        D = get_D(E, v);
        A = .5*det([ones(3,1) Ex(i,:)' Ey(i,:)']);
        ef=plantf(Ex(i,:), Ey(i,:), ep, (Tdelta(i)*alpha*D*[1;1;0])');
        %C = [ 1  Ex(1) Ey(1)   0          0          0  
        %    0         0        0   1   Ex(1)   Ey(1)
        %    1  Ex(2) Ey(2)   0          0          0  
        %    0         0        0   1   Ex(2)   Ey(2)
        %    1  Ex(3) Ey(3)   0          0          0  
        %    0         0        0   1   Ex(3)   Ey(3)];

        %B = [0 1 0 0 0 0
        %    0 0 0 0 0 1
        %    0 0 1 0 1 0]/C;
        
        %ef = B'*D*(1+v)*alpha*ep(2)*(1/3)*([1,1,1]*Tdelta_node(i,:)')*A*[1,1,0]';
        
        Ke = plante(Ex(i,:), Ey(i,:), ep, D);
        
        Fs(indx) = Fs(indx) + ef;
        Ks(indx, indx) = Ks(indx, indx) + Ke;
    end
    %figure;
    %spy(Ks)
    bc_x = find(boundary_element_type == 2 | boundary_element_type == 4)';
    bc_x = bc_x*2-1;
    bc_y = find(boundary_element_type == 3 | boundary_element_type == 4)';
    bc_y = bc_y*2;
    bc = [bc_x, zeros(size(bc_x)); bc_y, zeros(size(bc_y))];
    [a,Q] = solveq(Ks,Fs, bc);
    ed = extract(edof,a);
    %edx = mean(a(edof(:,2:2:end)), 2);
    %edy = mean(a(edof(:,3:2:end)), 2);
    %ed = reshape([edx edy], 1, []);
    %ed = a;
    Seff_el = zeros(nelm, 1);
    %es = zeros(nelm,3);
    %et = zeros(nelm,3);
    %sigmazz = zeros(nelm,1);
    for i = 1:nelm
        E = young_mod(element_map(i));
        v = poss_rat(element_map(i));
        alpha = exp_coeff(element_map(i));
        D = get_D(E, v);
        A=1/2*det([ones(3,1) Ex(i,:)' Ey(i,:)']);
        [sigma, ~] = plants(Ex(i,:), Ey(i,:), [2 1], D, ed(i,:));
        s_xx = sigma(1) - alpha*E*Tdelta(i)/(1-2*v);
        s_yy = sigma(2) - alpha*E*Tdelta(i)/(1-2*v);
        tauxy = sigma(3);
        s_zz = v * (s_xx + s_yy);% - alpha*E*Tdelta(i);
        sigma_eff_e = sqrt(s_xx^2 ...
                           + s_yy^2 ...
                           + s_zz^2 ... 
                           - s_xx * s_zz ...
                           - s_yy* s_zz ...
                           - s_xx*s_yy ...
                           + 3*tauxy^2);
        Seff_el(i) = sigma_eff_e;
    end
    Seff_nod = nan(nnod, 1);
    for i = 1:nnod
       [c0,~] = find(edof(:,2:2:end)==i);
       Seff_nod(i,1) = sum(Seff_el(c0))/size(c0,1); 
    end
end
%
function D = get_D(E, v)
    %D = hooke(1, E, v);
    D = (E./((v+1)*(1-2*v))) * ...
            [1-v,   v,   0;
             v,   1-v,   0;
             0,   0,     0.5*(1-2*v)];
end
%
function [boundary_element_type] = find_boundaries(nelm, nnod, coord)
    % find boundry elements
    % 0 = not boundary
    % 1 = convection, free displacement
    % 2 = isolated, no x displacement
    % 3 = isolated, no y displacement
    % 4 = isolated, no displacement
    boundary_element_type = zeros(1, nelm);
    for i=1:nnod
        x_coord = coord(i,1);
        y_coord = coord(i,2);
        if x_coord == 0 && y_coord > 0.015
            boundary_element_type(i) = 1;
        elseif x_coord == 0.025 && y_coord > 0.0325
            boundary_element_type(i) = 1;
        else
            if (x_coord == 0 || x_coord == 0.025) && ...
                    (y_coord == 0 || y_coord == 0.05) %% literally a cornercase
                boundary_element_type(i) = 4;
            elseif x_coord == 0 || x_coord == 0.025
                boundary_element_type(i) = 2;
            elseif y_coord == 0 || y_coord == 0.05
                boundary_element_type(i) = 3;
            end
        end
        if y_coord == 0.05
            boundary_element_type(i) = 3;
        end
    end
end
    
function [nelm, edof, nnod, dof, ndof, coord, Ex, Ey, element_map] = extract_params(p, e, t)
    %% extract mesh data from pdetool format
    element_map = t(4,:)'; %maps nodes to which part of the battery they belong to (core, shell, etc.)
    nelm = length(t(1,:));
    edof(:,1) = 1:nelm;
    edof(:,2:4) = t(1:3,:)';
    coord=p';
    ndof = max(max(t(1:3,:)));
    [Ex,Ey] = coordxtr(edof,coord,(1:ndof)',3);
    enod=t(1:3,:)'; % nodes of elements
    nelm=size(enod,1); % number of elements
    nnod=size(coord,1); % number of nodes
    dof=(1:nnod)'; % dof number is node number
    for ie=1:nelm
        edof(ie,:)=[ie,enod(ie,:)];
    end
end
%%
function [a, Q, K, f, C] = heat(Tinf, nnod, nelm, ndof, edof, element_map, therm_cond, spec_heat, density, boundary_element_type, coord, Q, Ex, Ey)
    %% Calculate the stationary heat distribution
    % Find K, fb, and C for heat equation
    ep = 100; %Thickness of battery, pretty arbitrary
    K = spalloc(nnod,nnod,3*nelm);
    f = zeros(ndof, 1);
    C = spalloc(nnod,nnod,3*nelm);
    
    for i = 1:nelm
        D = therm_cond(element_map(i))*eye(2);
        [Ke, Fe] = flw2te(Ex(i,:),Ey(i,:),ep,D,Q(i));
        c = ep*spec_heat(element_map(i))*density(element_map(i));
        Ce=plantml(Ex(i,:),Ey(i,:),c);
        indx = edof(i,2:end);
        %indx
        %Ke
        K(indx,indx) = K(indx,indx)+Ke;
        f(indx) = f(indx) + Fe;
        C(indx, indx) = C(indx,indx) + Ce;
    end
    % Add convection term to force vector
    for k = 1:nelm
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
    %
    [a,Q] = solveq(K,f);
end