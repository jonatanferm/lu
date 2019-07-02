M = moviein(antaltid+1);
for j=1:(antaltid+1)
   plot(X,losn(j,:))
   axis([0 L -1.1 1.1])
   M(:,j)=getframe;
end
movie(M)
