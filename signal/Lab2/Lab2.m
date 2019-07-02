%% 2.1
x1 = [1, 2, 1, 0, 0, 0, 0, 0, 0]';
x2 = [1, 2, 1, 0, 0, 0, 1, 2, 1]';
[c, l] = xcorr(x1, x2);
stem(l, c)
clear
%% 2.2
load man1
load man2
FT=8000;
%%
sound(man1,FT)
%%
sound(man2,FT)
%%
figure
t=(0:length(man2)-1)/FT;
plot(t,man2)
t=(0:length(man1)-1)/FT;
hold on
plot(t,man1,'r')
xlabel('time (s)')
%%
zoomman1=man1(33237:40000);
zoomman2=man2(33237:40000);
%%
sound(zoomman1,FT)
%%
sound(zoomman2,FT)
%%
figure
subplot(311)
plot(zoomman2)
subplot(312)
plot(zoomman1)
likhet=xcorr(zoomman2,zoomman1);
subplot(313)
plot(likhet)
%%
[c, l] = xcorr(zoomman2);
stem(l, c)
%%
l = 2000;
skattadman1=man2(l+1:end)-0.615*[zeros(2436-l,1);man1];
sound([skattadman1])
%% 2.3
t = (0:2*10000-1)/10000;
x = chirp(t, 300, 2, 800);
soundsc(x, 10000);
%%
D = 500;
beta = 0.4;
z = zeros(1, D);
x1 = [x, z];
x2 = [z, x];
y = x1 + beta*x2;
soundsc(y, 10000);
%%
[c, l] = xcorr(y);
stem(l, c)
%% 2.4
load gsmsig;
N = size(gsmsig, 1);
t = (0:N-1)/10000;
plot(t, gsmsig);
%%
soundsc(gsmsig, 10000);
%%
[c, l] = xcorr(gsmsig(:,1), gsmsig(:,2), 'none');
stem(l, c) %%948 stegs delay
%% 2.5
load LMEB.txt
size(LMEB)
figure
plot(LMEB)
%%
LMEB2=conv([1 -1],LMEB);
LMEB2=LMEB2(2:length(LMEB2)-1); % ta bort kantv¨arden
medel=mean(LMEB2);
LMEB2 = LMEB2 - medel;
[c, l]=xcorr(LMEB2, 'coeff');
plot(l, c)
%%
figure
subplot(211)
plot(LMEB2)
subplot(212)
plot(korrfkn)
%%
c(l == 1)
2/sqrt(size(LMEB2, 1))
%%
D1 = 5
LMEB5dmedel=filter(0.2*ones(1,D1),1,LMEB);
D2 = 25
LMEB30dmedel=filter(1/D2*ones(1,D2),1,LMEB);
%%
figure
plot(LMEB)
hold on
plot(LMEB5dmedel,'r')
plot(LMEB30dmedel,'g')
%%
Summa(1)=100;
Ready=1; % 1 n¨ar vi kan k¨opa aktier och 0 n¨ar vi redan ¨ager aktier
for dag=2:length(LMEB)
    dag
    % To buy
    if (LMEB5dmedel(dag)>LMEB30dmedel(dag))&~(LMEB5dmedel(dag-1)>LMEB30dmedel(dag-1))&Ready==1
        Summa(dag)=Summa(dag-1);
        Antal=Summa(dag)/LMEB(dag);
        Ready=0;
    elseif (LMEB5dmedel(dag)<LMEB30dmedel(dag))&~(LMEB5dmedel(dag-1)<LMEB30dmedel(dag-1))&Ready==0 % to sell
        Summa(dag)=Antal*LMEB(dag);
        Ready=1;
    else
        if Ready==1
            Summa(dag)=Summa(dag-1);
        else
            Summa(dag)=Antal*LMEB(dag);
        end;
    end;
end; %dag
figure
plot(Summa)


