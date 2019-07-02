clear;
clc;
m = 1/1000;
al = polyshape(m.*[0, 2.5, 2.5, 22.5, 22.5, 25, 25, 0],...
        m.*[47.5, 47.5, 2.5, 2.5, 42.5, 42.5, 0, 0]) % 1 al shell
st = polyshape(m.*[0,    0,  25, 25,   10,   10], ...
        m.*[47.5, 50, 50, 42.5, 42.5, 47.5]) % 2 steel top
cu = polyshape(m.*[10,   15,   15, 12.5, 10], ...
        m.*[42.5, 42.5, 15, 12.5, 15]) % 3 copper thing
co = polyshape(m.*[15,   22.5, 22.5, 2.5, 2.5,  10,   10, 12.5, 15], ...
        m.*[42.5, 42.5, 2.5,  2.5, 47.5, 47.5, 15, 12.5, 15]) % 4 core thing
    
figure;
hold on;
plot(al, 'FaceColor', [255,255,204]/255, 'FaceAlpha', 1)
plot(st, 'FaceColor', [34,94,168]/255, 'FaceAlpha', 1)
plot(cu, 'facecolor', [65,182,196]/255, 'FaceAlpha', 1)
plot(co, 'facecolor', [161,218,180]/255, 'FaceAlpha', 1)
pbaspect([1 1 1])
legend({'Aluminim Shell', 'Steel negative terminal', 'Copper negative collector', 'Battery Core'})
xlabel('Width [m]', 'Interpreter', 'latex', 'FontSize', 15)
ylabel('Height [m]', 'Interpreter', 'latex', 'FontSize', 15)
hold off;