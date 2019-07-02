%Här beräknas Fourier-, sinus eller cosinuskoefficienter för en 
%funktionsvektor f, definierad i punkterna X av längd N
%
%fråga till användaren:
inp=input('Fourier-, sinus- eller cosinuskoefficienter? (f/s/c) ','s');
%om man önskar få trigonometriska Fourierkoefficienter:
if inp=='f'
fourf  = fft(f)/N; %snabb Fouriertransform av f, komplexa koefficienter
%bestämning av trigonometriska F-koeff (förutsatt att funktionen är reell):
rfourf = real(fourf);
ifourf = imag(fourf);
akoeff =  2*rfourf(2:2^(n-1));
bkoeff = -2*ifourf(2:2^(n-1));
anoll  =  2*fourf(1);
%om man önskar få koefficienter i sinusutvecklingen av f:
elseif inp=='s'
   flp=fliplr(f);        %spegla och gör en udda utvidgning
   f=[f f(N) -flp(1:N-1) ]; %f(1)=0; 
   fourf  = fft(f)/N/2;  %snabb Fouriertransform av utvidgad funktion
   %bestämning av trigonometriska Fourierkoefficienter
   rfourf = real(fourf);
   ifourf = imag(fourf);
   akoeff =  2*rfourf(2:2^n);
   bkoeff = -2*ifourf(2:2^n);
   anoll  =  2*fourf(1);
   f=f(1:N);         %återgå till halvperiod, och justera Fourierkoeff
   akoeff=akoeff(1:N/2-1);
   bkoeff=bkoeff(1:N/2-1);
   
%om man önskar få koefficienter i cosinusutvecklingen av f:   
else
   flp=fliplr(f);    %spegla och gör en jämn utvidgning 
   f=[f f(N) flp(1:N-1) ];
   fourf  = fft(f)/N/2; %snabb Fouriertransform av utvidgad funktion
   % Bestämning av trigonometriska Fourierkoefficienter
   rfourf = real(fourf);
   ifourf = imag(fourf);
   akoeff =  2*rfourf(2:2^n);
   bkoeff = -2*ifourf(2:2^n);
   anoll  =  2*fourf(1);
   f=f(1:N);           %återgå till halvperiod och justera koefficienter
   akoeff=akoeff(1:N/2-1);
   bkoeff=bkoeff(1:N/2-1);
end

hold off
%konstruera en figur med två diagram ovanför varandra, det översta med
%akoeff, det understa med bkoeff
mx=0.25;
figure(1)
subplot(211)
stem(akoeff)
title('akoeff')
axis([0 N/2-1 -mx mx])
subplot(212)
stem(bkoeff)
title('bkoeff')
axis([0 N/2-1 -mx mx])
