hold off
inresteg=0.01;
clear summa
clear term
x=-1:inresteg:2;

term=anoll/2*ones(size(x)); %konstanta termen i Fourierserien
summa=term;
k = 0;
xdel=0:inresteg:L;
Inp='j';
%Inp=input('Är funktionen känd? (j/n) ','s');
 if Inp=='j'
     for k=1:(N/2-1),   
%     for k=1:32,   
       figure(2)
      hold off
      clf reset
      plot(x,summa); 
      axis([-1 2 -1.2 1.2])
     hold on
     plot(X,f,'r')
     title([ 'delsumma s_{' int2str(k-1) '}' ])
     hold off
     pause
     if inp~='f'     
        L=2*L;     %man måste ta hänsyn till att perioden är 2L för c och s
     end
     term=akoeff(k)*cos(2*pi*k*x/L)+bkoeff(k)*sin(2*pi*k*x/L); % uppdatera term
     summa=summa+term; % uppdatera summa
     if inp~='f'
        L=L/2;
     end
  end
else
   for k=1:(N/2-1),
      hold off
       figure(2)   %%
      plot(x,summa); 
      axis([-1 2 -2 2])
     title([ 'delsumma s_{' int2str(k-1) '}' ])
     pause
     term=akoeff(k)*cos(2*pi*k*x/L)+bkoeff(k)*sin(2*pi*k*x/L); % uppdatera term
     summa=summa+term; % uppdatera summa
  end
end
