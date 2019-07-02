tidsteg = .002;
antaltid=60;
savef=f;        %görs för att spara betydelsen av f, akoeff, etc
saveakoeff=akoeff;
savebkoeff=bkoeff;
saveanoll=anoll;
%
tid = tidsteg*(0:antaltid);
index=1:(N/2-1);
const=pi^2/L^2;
antalterm = 127;
%antalterm=input('antal termer (högst 127)')
%använd formeln från Fouriers metod för lösning av värmeledningsproblem:
takoeff = exp(-tid'*index.^2*const).*(ones(size(tid))'*akoeff);
tbkoeff = exp(-tid'*index.^2*const).*(ones(size(tid))'*bkoeff);
losn = zeros(size(tid,2),size(X,2));
for j=1:(antaltid+1)
   akoeff = takoeff(j,:);
   bkoeff = tbkoeff(j,:);
   summera;
   losnj = summ;
   losn(j,:) = losnj;
end
f=savef;
akoeff=saveakoeff;
bkoeff=savebkoeff;
anoll=saveanoll;

