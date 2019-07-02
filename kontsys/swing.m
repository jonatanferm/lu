tidsteg = .02;
antaltid=100;
savef=f;          %görs för att spara betydelsen av f, akoeff, etc
saveakoeff=akoeff;
savebkoeff=bkoeff;
saveanoll=anoll;
%
tid = tidsteg*(0:antaltid);
index=1:(N/2-1);
const=pi/L;
antalterm = 127;
%använd formeln från Fouriers metod för lösning av vågekvationen
takoeff = cos(-tid'*index*const).*(ones(size(tid))'*akoeff);
tbkoeff = cos(-tid'*index*const).*(ones(size(tid))'*bkoeff);
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



   
