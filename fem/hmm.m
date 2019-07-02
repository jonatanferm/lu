nelm=length(t(1,:)) 
edof(:,1)=1:nelm  
edof(:,2:4)=t(1:3,:)' 
coord=p' ;
ndof=max(max(t(1:3,:))) 
[Ex,Ey]=coordxtr(edof,coord,(1:ndof)',3);
eldraw2(Ex,Ey,[1,4,1])  


K = zeros(40);
f = K(:,1);
D = 1.059*eye(2);

for i = 1:40
    [Ke,fe] = flw2te(Ex(i,:), Ey(i,:),1,D,0);
    [K,f] = assem(edof,K,Ke,f,fe);
end
util = zeros(40);
dof = edof(:,2:4);
bc = [dof util(:,1)];
a = solveq(K,f)
%%

nelm=length(t(1,:)); 
edof(:,1)=1:nelm ; 
edof(:,2:4)=t(1:3,:)'; 
coord=p' ;
ndof=max(max(t(1:3,:))); 
[Ex,Ey]=coordxtr(edof,coord,(1:ndof)',3);
eldraw2(Ex,Ey,[1,4,1])
grid on
K = zeros(ndof);
C = zeros(ndof);
Cstiff = zeros(length(p), 3);
Kc = zeros(ndof);
alphaSOL = 1.2*10^(-5);
alphaSMD = alphaSOL;
ac = 40;
f = K(:,1);
fs = f;
fb = K(:,1);
D = 0;
x = 0;
qn = 9*10^3;  
Tinf = 20; %20 Celsius
ex = [];
ey = [];
Ex = Ex/1000;
Ey = Ey/1000;
edof2 = zeros(length(edof), 7);
edof2(:, 1) = edof(:, 1);

nodesSOL = [];
%Hitta nodes på Solderns rand
for i = 1:length(e(1,:))
    if e(6, i) == 0 && e(7, i) == 1 %Borde vara tvärtom? e verkar gå medurs
       p1 = e(1, i);
       p2 = e(2, i);
       if ~ismember(p1, nodesSOL)
           nodesSOL = [nodesSOL p1];
       end
       if ~ismember(p2, nodesSOL)
           nodesSOL = [nodesSOL p2];
       end
    end
end
%Ta reda på längden mellan noder på solderns rand
dlSOL = [];
for i = 1:length(nodesSOL)-1
   dx = p(1, nodesSOL(i+1)) - p(1, nodesSOL(i));
   dy = p(2, nodesSOL(i+1)) - p(2, nodesSOL(i));
   dl = sqrt((dx^2)+(dy^2)); 
   dlSOL = [dlSOL dl];
end
dlSOL = dlSOL/1000;

% for i = 1:length(p(1,:))
%     if p(2,i) == 0.6
%         if p(1,i) >= 0.19999 && p(1,i) <= 0.6 
%             nodesSMD1 = [nodesSMD1 i];
%         end
%     if p(1, i) <= 0.2 && p(1,i) >= 0
%         nodesSMD2 = [nodesSMD2 i];
%     end
%     end
% end

%Hitta nodes på SMD-randen, 1 med konvektion, 2 med q_n = q_el
nodesSMD1 = [];
nodesSMD2 = [];
for i = 1:length(e(1,:))
    if e(6, i) == 0 && e(7, i) == 2 
       p1 = e(1, i);
       p2 = e(2, i);
       p1coord = p(:, p1);
       p2coord = p(:, p2);
       if p1coord(2) == 0.6 && p2coord(2) == 0.6
           if p1coord(1) >= 0.19999 && p1coord(1) <= 0.6
       if ~ismember(p1, nodesSMD1)
           nodesSMD1 = [nodesSMD1 p1];
       end
           end
           if p2coord(1) >= 0.19999 && p2coord(1) <= 0.6
       if ~ismember(p2, nodesSMD1)
           nodesSMD1 = [nodesSMD1 p2];
       end
           end
           if p1coord(1) <= 0.2 && p1coord(1) >= 0
               if ~ismember(p1, nodesSMD2)
                  nodesSMD2 = [nodesSMD2 p1]; 
               end
           end
           if p2coord(1) <= 0.2 && p2coord(1) >= 0
               if ~ismember(p2, nodesSMD2)
                  nodesSMD2 = [nodesSMD2 p2]; 
               end
           end
       end
    end
end

%Ta reda på längden mellan noder på SMDs rand
dlSMD1 = [];
for i = 1:length(nodesSMD1)-1
    for j = 1:length(nodesSMD1)-1
   
   dx = p(1, nodesSMD1(i+1)) - p(1, nodesSMD1(i));
   dy = p(2, nodesSMD1(i+1)) - p(2, nodesSMD1(i));
   dl = sqrt((dx^2)+(dy^2)); 
   dlSMD1 = [dlSMD1 dl];
   
    end
end
dlSMD1 = dlSMD1/1000;
dlSMD2 = [];
for i = 1:length(nodesSMD2)-1
   dx = p(1, nodesSMD2(i+1)) - p(1, nodesSMD2(i));
   dy = p(2, nodesSMD2(i+1)) - p(2, nodesSMD2(i));
   dl = sqrt((dx^2)+(dy^2)); 
   dlSMD2 = [dlSMD2 dl];
end
  dlSMD2 = dlSMD2/1000;
%Hitta nodes på randen x = 0
 nodesMID = [];
 for i = 1:length(e(1,:))
     if (p(1, e(1,i)) == 0) && (p(1, e(2,i)) == 0)
       if (e(6, i) == 0 && e(7, i) == 2) || (e(6, i) == 0 && e(7, i) == 3)
       p1 = e(1, i);
       p2 = e(2, i);
       if ~ismember(p1, nodesMID)
           nodesMID = [nodesMID p1];
       end
       if ~ismember(p2, nodesMID)
           nodesMID = [nodesMID p2];
       end
    end
     end
 end
dlMID = [];
for i = 1:length(nodesMID)-1
   dx = p(1, nodesMID(i+1)) - p(1, nodesMID(i));
   dy = p(2, nodesMID(i+1)) - p(2, nodesMID(i));
   dl = sqrt((dx^2)+(dy^2)); 
   dlMID = [dlMID dl];
end
 dlMID = dlMID/1000;
nodesPCB = [];

%Hitta noder på PCB-boundary
for i = 1:length(e(1,:))
    if p(2,e(1,i)) == 0 && p(2,e(2,i)) == 0
            p1 = e(1, i);
            p2 = e(2, i);
            if ~ismember(p1, nodesPCB)
                nodesPCB = [nodesPCB p1];
            end
            if ~ismember(p2, nodesPCB)
                nodesPCB = [nodesPCB p2];
            end
    end
    if p(1,e(1,i)) == 1 && p(1,e(2,i)) == 1
        p1 = e(1, i);
        p2 = e(2, i);
            if ~ismember(p1, nodesPCB)
                nodesPCB = [nodesPCB p1];
            end
            if ~ismember(p2, nodesPCB)
                nodesPCB = [nodesPCB p2];
            end
    end
end

%Skapa K och C
l=1;
for i = 1:length(edof(:,1))
    if t(4, i) == 3 %PCB
        D = 1.059*eye(2); 
        x = 1850*950;
        E = 105*10^9;
        nu = 0.118;
        alpha = 2*10^(-5);
    end
    if t(4, i) == 1 %Solder
        D = 66.8*eye(2);
        x = 7265*210;   
        E = 50*10^9;
        nu = 0.36;
        alpha = 1.2*10^(-5);
    end
    if t(4, i) == 2 %SMD
        D = 0.29*eye(2);
        x = 1850*950;
        E = 105*10^9; 
        nu = 0.118;
        alpha = 1.2*10^(-5);
    end
    %SMD2 (inte konvektion)
     if (ismember(t(1, i), nodesSMD2) && ismember(t(2,i), nodesSMD2)) 
          fbSMD2 = dlSMD2(l)*(1/2)*qn*[1 1 0]';
          [K, fb] = assem(edof(i,:), K, Ke, fb, fbSMD2);
          l = l+1;
     end
     if ismember(t(1, i), nodesSMD2) && ismember(t(3,i), nodesSMD2)
         fbSMD2 = dlSMD2(l)*(1/2)*qn*[1 0 1]';
          [K, fb] = assem(edof(i,:), K, Ke, fb, fbSMD2);
          l = l+1;
     end
     if ismember(t(2, i), nodesSMD2) && ismember(t(3,i), nodesSMD2)
         fbSMD2 = dlSMD2(l)*(1/2)*qn*[0 1 1]';
          [K, fb] = assem(edof(i,:), K, Ke, fb, fbSMD2);
          l = l+1;
     end
%      x=0, y=0 -> y=0.6 (mittenranden)
%      if (ismember(t(1, i), nodesMID) && ismember(t(2,i), nodesMID)) || (ismember(t(1, i), nodesMID) && ismember(t(3,i), nodesMID)) ||(ismember(t(2, i), nodesMID) && ismember(t(3,i), nodesMID)) 
%          
%      end
    [Ke, fe] = flw2te(Ex(i, :),Ey(i,:),1,D, 0);
    [K, f] =  assem(edof(i,:),K,Ke, f, fe);
    Ce = plantml(Ex(i,:), Ey(i,:), x);
    C = assem(edof(i,:), C, Ce);
end
j = 1;
l = 1;
%Skapa Kc (konvektion)
 for k = 1:length(edof(:,1))   
     %Solder
     if (ismember(t(1, k), nodesSOL) && ismember(t(2,k), nodesSOL))
                NtNSOL = ac*dlSOL(j)*(1/6)*[2 1 0; 1 2 0; 0 0 0];
                fbSOL = ac*dlSOL(j)*(1/2)*Tinf*[1 1 0]';
                [Kc, fb] = assem(edof(k,:), Kc, NtNSOL, fb, fbSOL); 
                j = j+1;
     end
     if (ismember(t(1, k), nodesSOL) && ismember(t(3,k), nodesSOL))
                NtNSOL = ac*dlSOL(j)*(1/6)*[2 0 1; 0 0 0; 1 0 2];
                fbSOL = ac*dlSOL(j)*(1/2)*Tinf*[1 0 1]';
                [Kc, fb] = assem(edof(k,:), Kc, NtNSOL, fb, fbSOL); 
                j = j+1;
     end
     if (ismember(t(2, k), nodesSOL) && ismember(t(3,k), nodesSOL))
                NtNSOL = ac*dlSOL(j)*(1/6)*[0 0 0; 0 2 1; 0 1 2];
                fbSOL = ac*dlSOL(j)*(1/2)*Tinf*[0 1 1]';
                [Kc, fb] = assem(edof(k,:), Kc, NtNSOL, fb, fbSOL); 
                j = j+1;
     end
     %SMD1
     if (ismember(t(1, k), nodesSMD1) && ismember(t(2,k), nodesSMD1)) 
                NtNSMD1 = ac*dlSMD1(l)*(1/6)*[2 1 0; 1 2 0; 0 0 0];
                fbSMD1 = ac*dlSMD1(l)*(1/2)*Tinf*[1 1 0]';
                [Kc, fb] = assem(edof(k,:), Kc, NtNSMD1, fb, fbSMD1);
                l = l+1;
     end
     if (ismember(t(1, k), nodesSMD1) && ismember(t(3,k), nodesSMD1))
                NtNSMD1 = ac*dlSMD1(l)*(1/6)*[2 0 1; 0 0 0; 1 0 2];
                fbSMD1 = ac*dlSMD1(l)*(1/2)*Tinf*[1 0 1]';
                [Kc, fb] = assem(edof(k,:), Kc, NtNSMD1, fb, fbSMD1);
                l = l+1;
     end
     if (ismember(t(2, k), nodesSMD1) && ismember(t(3,k), nodesSMD1))
                NtNSMD1 = ac*dlSMD1(l)*(1/6)*[0 0 0; 0 2 1; 0 1 2];
                fbSMD1 = ac*dlSMD1(l)*(1/2)*Tinf*[0 1 1]';
                [Kc, fb] = assem(edof(k,:), Kc, NtNSMD1, fb, fbSMD1);
                l = l+1;
     end
 end  
  f = f + fb;
  Ktot = K+Kc;
  [statSolution, Q] = solveq(Ktot,f);
  ex = [Ex; -Ex];
  ey = [Ey; Ey];
  ed = extract(edof,statSolution);
  ed2 = [ed; ed];
  figure(1)
  fill(ex',ey',ed2')
 
 %%
 figure(2)
 N = 100;
 tf = 100000;
 t0 = 0;
 dt = (tf-t0)/N;
 T0 = 30;
 aold = T0*ones(ndof,1);
 a = [aold];
   for i = 1:N
       temp = impStep(aold, Ktot, C, f, dt);
       aold = temp;
       a = [a temp];
   end
   ex = [Ex; -Ex];
   ey = [Ey; Ey];
   edtime = extract(edof, a(:, end));
   edtime2 = [edtime; edtime];
   fill(ex', ey', edtime2');
% for i = 1:N
%     edtime = extract(edof, a(:,i));
%     edtime2 = [edtime; edtime];
%     fill(ex', ey', edtime2');
%     pause(1)
% end
%%
%Skapa edof för stress-problemet
dT = ed-T0*ones(length(ed(:,1)), 3);
ndofStress = 2*ndof;
Ks = zeros(2*ndof);
fs = zeros(2*ndof,1);
Ks = zeros(2*ndof);
edof2 = zeros(nelm, 7);
for i = 1:nelm
    edof2(i,1) = i;
   for j = 2:4
       edof2(i,2*j-2) = 2*edof(i,j)-1;
       edof2(i,2*j-1) = 2*edof(i,j);
   end
   
end
nodesPCB1 = [];
nodesPCB2 = [];
for i = 1:length(e(1,:))
    if p(2,e(1,i)) == 0 && p(2,e(2,i)) == 0
            p1 = e(1, i);
            p2 = e(2, i);
            if ~ismember(p1, nodesPCB1)
                nodesPCB1 = [nodesPCB1 p1];
            end
            if ~ismember(p2, nodesPCB1)
                nodesPCB1 = [nodesPCB1 p2];
            end
    end
    if p(1,e(1,i)) == 1 && p(1,e(2,i)) == 1
        p1 = e(1, i);
            p2 = e(2, i);
            if ~ismember(p1, nodesPCB2)
                nodesPCB2 = [nodesPCB2 p1];
            end
            if ~ismember(p2, nodesPCB2)
                nodesPCB2 = [nodesPCB2 p2];
            end
    end
end
for i = 1:length(nodesPCB1)
   PCB1free(i,1) = 2*nodesPCB1(i)-1;
   PCB1free(i,2) = 2*nodesPCB1(i);
end
for i = 1:length(nodesSOL)
   SOLfree(i,1) = 2*nodesSOL(i)-1;
   SOLfree(i,2) = 2*nodesSOL(i);
end
for i = 1:length(nodesPCB2)
   PCB2free(i,1) = 2*nodesPCB2(i)-1;
   PCB2free(i,2) = 2*nodesPCB2(i);
end
for i = 1:length(nodesSMD1)
   SMD1free(i,1) = 2*nodesSMD1(i)-1;
   SMD1free(i,2) = 2*nodesSMD1(i);
end
for i = 1:length(nodesSMD2)
   SMD2free(i,1) = 2*nodesSMD2(i)-1;
   SMD2free(i,2) = 2*nodesSMD2(i);
end
for i = 1:length(nodesMID)
   MIDfree(i,1) = 2*nodesMID(i)-1;
   MIDfree(i,2) = 2*nodesMID(i);
end
nodesPCB1 = [];
nodesPCB2 = [];
%Hitta noder på PCB-boundary
for i = 1:length(e(1,:))
    if p(2,e(1,i)) == 0 && p(2,e(2,i)) == 0
            p1 = e(1, i);
            p2 = e(2, i);
            if ~ismember(p1, nodesPCB1)
                nodesPCB1 = [nodesPCB1 p1];
            end
            if ~ismember(p2, nodesPCB1)
                nodesPCB1 = [nodesPCB1 p2];
            end
    end
    if p(1,e(1,i)) == 1 && p(1,e(2,i)) == 1
        p1 = e(1, i);
            p2 = e(2, i);
            if ~ismember(p1, nodesPCB2)
                nodesPCB2 = [nodesPCB2 p1];
            end
            if ~ismember(p2, nodesPCB2)
                nodesPCB2 = [nodesPCB2 p2];
            end
    end
end

for i = 1:nelm
    if t(4, i) == 3 %PCB
        E = 10.5*10^9;
        nu = 0.136;
        alpha = 2*10^-5;
    end
    
    if t(4,i) == 1
        E = 50*10^9;
        nu = 0.36;
        alpha = 1.2*10^-5; 
    end
    
    if t(4,i) == 2
        E = 10.5*10^9;
        nu = 0.118;
        alpha = 1.2*10^-5;
    end
    Ds = (E/((1+nu)*(1-2*nu)))*[1-nu nu 0; nu 1-nu 0; 0 0 0.5*(1-2*nu)];
    [Kse, fse] =  plante(Ex(i,:), Ey(i,:), [2 1], Ds);                                                                              
    [Ks,fs] = assem(edof2(i,:), Ks, Kse, fs, fse); 
end

% for i = 1:length(PCB2free)
%    if ismember(PCB2free(i,1) , SOLfree) && ismember(PCB2free(i,2), SOLfree)
%       corner = [PCB2free(i,1) PCB2free(i,2)];
%    end
% end

bcs = [PCB1free(:,2) zeros(length(PCB1free), 1); 
    PCB2free(:,1) zeros(length(PCB2free),1); 
    MIDfree(:,1) zeros(length(MIDfree),1)];


ep = [2 1];
estot = [];
f0 = zeros(length(Ks(:,1)),1);
for i = 1:nelm
    if t(4, i) == 3 %PCB
        E = 10.5*10^9;
        nu = 0.136;
        alpha = 2*10^-5;
    end
    if t(4,i) == 1
        E = 50*10^9;
        nu = 0.36;
        alpha = 1.2*10^-5; 
    end
    if t(4,i) == 2
        E = 10.5*10^9;
        nu = 0.118;
        alpha = 1.2*10^-5;
    end
Ds = (E/((1+nu)*(1-2*nu)))*[1-nu nu 0; nu 1-nu 0; 0 0 0.5*(1-2*nu)];
%[es, et] = plants(Ex(i,:), Ey(i,:), ep, Ds, ed4(i,:));
%estot = [estot; es]; 
dTav = mean(ed(i));
e0 = (1+nu)*alpha*dTav*[1 1 0]';
s0= Ds*e0;
ef = plantf(Ex(i,:), Ey(i,:), ep, s0');
[~, f0] = assem(edof2(i,:), Ks, Kse, f0, ef); 
end
u1 = solveq(Ks, f0, bcs);
ed5 = extract(edof2, u1);
ex = [Ex; -Ex];
ey = [Ey; Ey];
edexpanded = [ed5(:,1) ed5(:,2) ed5(:,3) ed5(:,4) ed5(:,5) ed5(:,6); 
             -ed5(:,1) ed5(:,2) -ed5(:,3) ed5(:,4) -ed5(:,5) ed5(:,6)];
 
 plotpar = [1, 4, 1];
close all
figure(3)
eldisp2(ex,ey,edexpanded,plotpar,300)

% for i=1:size(coord,1)
% [c0,c1]=find(edof(:,2:4)==i);
% Seff_nod(i,1)=sum(Seff_el(c0))/size(c0,1);
% end
%%
for i = 1:length(t)
   if ismember(21, t(:, i)) && ismember(113, t(:, i)) 
      i 
   end
end