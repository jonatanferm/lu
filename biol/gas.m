[X Y]=meshgrid(0:.1:10,0:.1:10);
Z=X.*sin(4*X)+1.1*Y.*sin(2*Y);
surf(X,Y,Z)