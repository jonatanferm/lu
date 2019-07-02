% This script is written and read by pdetool and should NOT be edited.
% There are two recommended alternatives:
% 1) Export the required variables from pdetool and create a MATLAB script
%    to perform operations on these.
% 2) Define the problem completely using a MATLAB script. See
%    http://www.mathworks.com/help/pde/examples/index.html for examples
%    of this approach.
function pdemodel
[pde_fig,ax]=pdeinit;
pdetool('appl_cb',1);
pdetool('snapon','on');
set(ax,'DataAspectRatio',[1 2.7857142857142856 1]);
set(ax,'PlotBoxAspectRatio',[35 23.333333333333336 1]);
set(ax,'XLim',[-5 30]);
set(ax,'YLim',[-5 60]);
set(ax,'XTick',[ -5,...
 -4,...
 -3,...
 -2,...
 -1,...
 0,...
 1,...
 2,...
 3,...
 4,...
 5,...
 6,...
 7,...
 8,...
 9,...
 10,...
 11,...
 12,...
 13,...
 14,...
 15,...
 16,...
 17,...
 18,...
 19,...
 20,...
 21,...
 22,...
 23,...
 24,...
 25,...
 26,...
 27,...
 28,...
 29,...
 30,...
]);
set(ax,'YTick',[ 0,...
 1,...
 2,...
 3,...
 4,...
 5,...
 6,...
 7,...
 8,...
 9,...
 10,...
 11,...
 12,...
 13,...
 14,...
 15,...
 16,...
 17,...
 18,...
 19,...
 20,...
 21,...
 22,...
 23,...
 24,...
 25,...
 26,...
 27,...
 28,...
 29,...
 30,...
 31,...
 32,...
 33,...
 34,...
 35,...
 36,...
 37,...
 38,...
 39,...
 40,...
 41,...
 42,...
 43,...
 44,...
 45,...
 46,...
 47,...
 48,...
 49,...
 50,...
 51,...
 52,...
 53,...
 54,...
 55,...
 56,...
 57,...
 58,...
 59,...
 60,...
]);
pdetool('gridon','on');

% Geometry description:
pderect([0 25 0 50],'R1');
set(findobj(get(pde_fig,'Children'),'Tag','PDEEval'),'String','R1')

% PDE coefficients:
pdeseteq(1,...
'1.0',...
'0.0',...
'10.0',...
'1.0',...
'0:10',...
'0.0',...
'0.0',...
'[0 100]')
setappdata(pde_fig,'currparam',...
['1.0 ';...
'0.0 ';...
'10.0';...
'1.0 '])

% Solve parameters:
setappdata(pde_fig,'solveparam',...
char('0','1000','10','pdeadworst',...
'0.5','longest','0','1E-4','','fixed','Inf'))

% Plotflags and user data strings:
setappdata(pde_fig,'plotflags',[1 1 1 1 1 1 1 1 0 0 0 1 1 0 0 0 0 1]);
setappdata(pde_fig,'colstring','');
setappdata(pde_fig,'arrowstring','');
setappdata(pde_fig,'deformstring','');
setappdata(pde_fig,'heightstring','');