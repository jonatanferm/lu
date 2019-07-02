%beräknar summan av 'antalterm' stycken termer
summ = anoll*ones(size(X))/2;
    if inp~='f'
        L=2*L;
     end
 for k=1:antalterm
   summ = summ + akoeff(k)*cos(2*pi*k*X/L) + bkoeff(k)*sin(2*pi*k*X/L);
end
    if inp~='f'
        L=L/2;
     end