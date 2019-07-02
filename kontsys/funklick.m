%registrera punkter från klickningar i figurfönstret
figure(1)    %ge fönstret ett nummer
clf reset
axis([0 L -1 1])  %välj axlar, [xmin xmax ymin ymax]
hold on           %se till att tidigare versioner av figuren
                  %ligger kvar när nya data läggs in

plot([0 L],[0 0],'k')
plot([0 0],[-1 1],'k')
x=[];y=[];
[xi0,yi0,button]=ginput(1);
x=[x 0];%Startpunktens x-koordinat justeras till 0
y=[y yi0];
while button==1
   [xi,yi,button]=ginput(1);
   if xi<=xi0 
      xi=xi0+0.001;
   end
   if button==1
       xi0=xi;
       x=[x xi];
       y=[y yi];
       plot(x,y,'b','Linewidth',2)
    end
 end 
 n1=size(x,2);
 x(n1)=L;
 f=interp1(x,y,X);%utfunktionsvektorn heter f, definierad i punkterna X
 %
 figure(2)
 plot(X,f)
          
