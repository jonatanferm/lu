clc
clear
load('data2')
load('data3')
%load('data13')
%load('data14')
%load('data15')


c1 = [213,62,79]/255;
c2 = [50,136,189]/255;

c3 = [251,180,185;
247,104,161;
197,27,138;
122,1,119]/255;
c3 = flipud(c3);

d = data2;
d3 = data3;
%d3 = data14;
%d4 = data15;
%d2 = data2;

figure;
hold on;
plot(d(:,1), d(:,5)/3, 'LineWidth', 1.5, 'Color', c2)
plot(d(:,1), d(:,2), 'LineWidth', 1.5, 'Color', c3(1,:))
%plot(d(:,1), d2(:,2), 'LineWidth', 1.5, 'Color', c3(2,:))
%plot(d(:,1), d3(:,2), 'LineWidth', 1.5, 'Color', c3(3,:))
%plot(d(:,1), d4(:,2), 'LineWidth', 1.5, 'Color', c3(4,:))
%plot([d2(:,1); 0], [d2(:,2); 0], 'LineWidth', 1.5, 'Color', c1)
%legend({'Teoretisk', 'Uppmätt'})
%xlabel('Fjärrfält (mm)')
%ylabel('Intensitet')
hold off;