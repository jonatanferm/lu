function fpextractdemo(action, varargin)
%	FPEXTRACTDEMO Fingerprint extracting the features on DEMO API program
%   Cheng Long Adam Wang, September 2002
%   
%   $Revision: 1.0 $  $Date: 2002.10.2  $
%   $Revision: 2.0 $  $Date: 2003.11.29 $    by Luigi Rosa
%                                            email:   luigi.rosa@tiscali.it 
%                                            mobile:  +393403463208
%                                            website: http://utenti.lycos.it/matlab
%
%   Modified:
%   - new GUI
%   - 8 Gabor filters 0 22.5 45 67.5 90 112.5 135 157.5 degree
%   - Convolution is performed in frequency domain
%   - DataBase
%   - Fingerprint matching
%   - Error management
%
% Input fingerprint should be 256 x 256 grayscale bmp image
% 8-bit grayscale image @ 500 dpi.
% If these conditions are not verified some parameters in m-functions
% should be changed in a proper way.
%  
%   Options:
%     - Centralize:    calculate the center point and the binarized image
%     - Crop:          image cropping
%     - Sectorize:     visualize the circular sectors
%     - Normalize:     normalize input image
%     - Gabor filters: visualize the Gabor filters
%     - Convolute:     compute convolution of input image and Gabor filters
%     - Features:      FingerCodes visualization
%     - FingerCode:    add the input fingerprint to Database
%     - Check:         fingerprint matching
%
%
% A crucial step in fingerprint recognition is center point determination. 
% If any error occurs while cropping image you can use the auxiliary m-file
% "vedicentro.m": it visualizes the input fingerprint and the center point
% calculated by the m-function "centralizing.m"
%
%
%   References
%
%   Cheng Long Adam Wang, researcher
%   Fingerprint Recognition System
%   http://home.kimo.com.tw/carouse9/FRS.htm
%
%   A. K. Jain, S. Prabhakar, and S. Pankanti, "A Filterbank-based Representation for 
%   Classification and Matching of Fingerprints", International Joint Conference on 
%   Neural Networks (IJCNN), pp. 3284-3285, Washington DC, July 10-16, 1999. 
%   http://www.cse.msu.edu/~prabhaka/publications.html
%
%   "Fingerprint Classification and Matching Using a Filterbank", Salil Prabhakar
%   A DISSERTATION Submitted to Michigan State University in partial fulfillment 
%   of the requirements for the degree of DOCTOR OF PHILOSOPHY, Computer 
%   Science & Engineering, 2001
%   http://biometrics.cse.msu.edu/SalilThesis.pdf
%
%   Final Report 18-551 (Spring 1999) Fingerprint Recognition Group Number 19
%   Markus Adhiwiyogo, Samuel Chong, Joseph Huang, Weechoon Teo
%   http://www.ece.cmu.edu/~ee551/Old_projects/projects/s99_19/finalreport.html
%
%   Type "fpextractdemo" on MATLAB prompt to start fingerprint processing.
%


%--------------------------------------------------------------------------
if nargin<1,
    action='InitializeFPEXTRACTDEMO';
end;

feval(action,varargin{:})
return;

%%%
%%%  Sub-function - InitializeFPEXTRACTDEMO
%%%

function InitializeFPEXTRACTDEMO()

% If fpextractdemo is already running, bring it to the foreground
h = findobj(allchild(0), 'tag', 'Extracting FingerPrint Features Demo');
if ~isempty(h)
    figure(h(1))
    return
end

screenD = get(0, 'ScreenDepth');
if screenD>8
    grayres=256;
else
    grayres=128;
end


FpextractDemoFig = figure( ...
    'Name','Extracting FingerPrint Features Demo', ...
    'NumberTitle','off', 'HandleVisibility', 'on', ...
    'tag', 'Extracting FingerPrint Features Demo', ...
    'Visible','off', 'Resize', 'off',...
    'BusyAction','Queue','Interruptible','off', ...
    'Color', [.8 .8 .8], ...
    'IntegerHandle', 'off', ...
    'Colormap', gray(grayres));

figpos = get(FpextractDemoFig, 'position');
figpos(3:4) = [1050 525];
% Adjust the size of the figure window
horizDecorations = 10;  % resize controls, etc.
vertDecorations = 45;   % title bar, etc.
screenSize = get(0,'ScreenSize');

dx = screenSize(3) - figpos(1) - figpos(3) - horizDecorations;
dy = screenSize(4) - figpos(2) - figpos(4) - vertDecorations;
if (dx < 0)
    figpos(1) = max(5,figpos(1) + dx);
end
if (dy < 0)
    figpos(2) = max(5,figpos(2) + dy);
end
set(FpextractDemoFig, 'position', figpos);

rows = figpos(4); 
cols = figpos(3);

% Colors
bgcolor = [0.45 0.45 0.45];  % Background color for frames
wdcolor = [.8 .8 .8];  % Window color
fgcolor = [1 1 1];  % For text

hs = (cols-(6*175)) / 5;        % Horizantal Spacing
vs = (rows)/8;                  % Vertical Spacing

%====================================
% Parameters for all buttons and menus

Std.Interruptible = 'off';
Std.BusyAction = 'queue';

% Defaults for image axes
Ax = Std;
Ax.Units = 'Pixels';
Ax.Parent = FpextractDemoFig;
Ax.ydir = 'reverse';
Ax.XLim = [.5 128.5];
Ax.YLim = [.5 128.5];
Ax.CLim = [0 1];
Ax.XTick = [];
Ax.YTick = [];

Img = Std;
Img.CData = [];
Img.Xdata = [1 128];
Img.Ydata = [1 128];
Img.CDataMapping = 'Scaled';
Img.Erasemode = 'none';

Ctl = Std;
Ctl.Units = 'Pixels';
Ctl.Parent = FpextractDemoFig;

Btn = Ctl;
Btn.Style = 'pushbutton';
Btn.Enable = 'off';

Edit = Ctl;
Edit.Style = 'edit';
Edit.HorizontalAlignment = 'right';
Edit.BackgroundColor = 'white';
Edit.ForegroundColor = 'black';

Menu = Ctl;
Menu.Style = 'Popupmenu';

Text = Ctl;
Text.Style = 'text';
Text.HorizontalAlignment = 'left';
Text.BackgroundColor = bgcolor;
Text.ForegroundColor = fgcolor;

%================================
% 0 degree Component 
ud.hComponent1Axes = axes(Ax, ...
    'Position', [0*vs/6 5*vs-vs/6 175 175]);
title('0 degree Component');
ud.hComponent1Image = image(Img, ...
    'Parent', ud.hComponent1Axes);
%================================
% Original FingerPrint 
ud.hOriginalAxes = axes(Ax, ...
    'Position', [cols/2-128 5*vs-vs/6-81 256 256]);
title('Original FingerPrint');
ud.hOriginalImage = image(Img, ...
    'Parent', ud.hOriginalAxes);
ud.OriginalImageIsStale = 1;

%================================
% 157.5 degree Component 
ud.hComponent8Axes = axes(Ax, ...
    'Position', [cols-175 5*vs-vs/6 175 175]);
title('157.5 degree Component');
ud.hComponent8Image = image(Img, ...
    'Parent', ud.hComponent8Axes);
%=================================
% 22.5 degree Component 
ud.hComponent2Axes = axes(Ax, ...
    'Position', [hs vs/2 175 175]);
title('22.5 degree Component');
ud.hComponent2Image = image(Img, ...
    'Parent', ud.hComponent2Axes);
%================================
% 45 degree Component 
ud.hComponent3Axes = axes(Ax, ...
    'Position', [2*hs+1*175 vs/2 175 175]);
title('45 degree Component');
ud.hComponent3Image = image(Img, ...
    'Parent', ud.hComponent3Axes);
%================================
% 67.5 degree Component 
ud.hComponent4Axes = axes(Ax, ...
    'Position', [3*hs+2*175 vs/2 175 175]);
title('67.5 degree Component');
ud.hComponent4Image = image(Img, ...
    'Parent', ud.hComponent4Axes);
%=================================
% 90 degree Component 
ud.hComponent5Axes = axes(Ax, ...
    'Position', [4*hs+3*175 vs/2 175 175]);
title('90 degree Component');
ud.hComponent5Image = image(Img, ...
    'Parent', ud.hComponent5Axes);
%=================================
% 112.5 degree Component 
ud.hComponent6Axes = axes(Ax, ...
    'Position', [5*hs+4*175 vs/2 175 175]);
title('112.5 degree Component');
ud.hComponent6Image = image(Img, ...
    'Parent', ud.hComponent6Axes);
%=================================
% 135 degree Component 
ud.hComponent7Axes = axes(Ax, ...
    'Position', [6*hs+5*175 vs/2 175 175]);
title('135 degree Component');
ud.hComponent7Image = image(Img, ...
    'Parent', ud.hComponent7Axes);
%=================================
%  The frame
ud.hControlFrame = uicontrol(Std, ...
    'Parent', FpextractDemoFig, ...
    'Style', 'Frame', ...
    'Units', 'pixels', ...
    'Position', [vs/6 5*vs-vs/6-81 200 vs+vs/8], ...
    'BackgroundColor', bgcolor);
%====================================
% Image popup menu
ud.hImgPop = uicontrol(Menu, ...
    'Position',[vs/6+vs/8 5*vs-2*vs/3+7 180 vs/16], ...
    'String','Whorl|Twin loop|Left loop|Right loop|Other image', ...
    'Callback','fpextractdemo(''LoadNewImage'')');
% Text label for Image Menu Popup
uicontrol( Text, ...
    'Position',[vs/6+vs/8 5*vs-vs/6-vs/3-2 180 vs/4], ...
    'String','Select a type of fingerprint:');
%====================================
% Extracting Step popup menu
ud.hSelectStepPop = uicontrol(Menu, ...
    'Position',[vs/6+vs/8 4*vs-7 120 vs/16], ...
    'String','Centralize|Crop|Sectorize|Normalize|Gabor filters|Convolute|Features|FingerCode|Check', ...
    'Callback','fpextractdemo(''SelectExtractingStep'')');
% Text label for Extracting Step Menu Popup
uicontrol( Text, ...
    'Position',[vs/6+vs/8 4*vs-4 90 vs/4], ...
    'String','Select step to:');
%====================================
%  Frame for Info and Close
ud.hInfoCloseFrame = uicontrol(Std, ...
    'Parent', FpextractDemoFig, ...
    'Style', 'Frame', ...
    'Units', 'pixels', ...
    'Position', [3*hs+2*175 2 vs/2+2*175 vs/2-4], ...
    'BackgroundColor', bgcolor);

%====================================
% Buttons - Info and Close
ud.hInfo=uicontrol(Btn, ...
    'Position',[3*hs+2*175+vs/2 7 vs/8+135-vs/2 vs/4], ...
    'String','Info', ...
    'Callback','helpwin fpextractdemo');

ud.hClose=uicontrol(Btn, ...
    'Position',[4*hs+3*175+vs/2 7 vs/8+135-vs/2 vs/4], ...
    'String','Close', ...
    'Callback','close(gcbf)');
%====================================
% Status bar
ud.hStatus = uicontrol(Std, ...
    'Parent', FpextractDemoFig, ...
    'Style','text', ...
    'Units','pixels', ...
    'Position',[hs vs/8 2*175-vs/8 vs/4], ...
    'Foreground', [.8 0 0], ...
    'Background',wdcolor, ...
    'Horiz','center', ...
    'Tag', 'Status', ...
    'String','Initializing fpextractdemo...');

set(FpextractDemoFig, 'UserData', ud);
set(FpextractDemoFig, 'visible','on','HandleVisibility','callback');
set([ud.hInfo ud.hClose], 'Enable', 'on');

LoadNewImage(FpextractDemoFig);
SelectExtractingStep(FpextractDemoFig);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%  Sub-Function - LoadNewImage
%%%

function LoadNewImage(DemoFig)
% Load a new image from a mat-file

if nargin<1
    DemoFig = gcbf;
end

set(DemoFig,'Pointer','watch');
ud=get(DemoFig,'Userdata');
v = get(ud.hImgPop,{'value','String'});
name = deblank(v{2}(v{1},:));
drawnow

switch name
    case 'Right loop',
        namefile='37_7.bmp';
        [img,map]=imread(namefile);
    case 'Whorl',
        namefile='19_7.bmp';
        [img,map]=imread(namefile);
    case 'Left loop',
        namefile='37_3.bmp';
        [img,map]=imread(namefile);
    case 'Twin loop',
        namefile='37_5_2.bmp';
        [img,map]=imread(namefile);
    case 'Other image',        
        [namefile,pathname]=uigetfile('*.bmp','Chose BMP GrayScale Image');
        if namefile~=0
            [img,map]=imread(strcat(pathname,namefile));
        else
            disp('   Chose a file!  ');
            [img,map]=imread('37_7.bmp');
        end
    otherwise 
        error('fpextractdemo: Unknown Image Option!');
end
% If image is N x M with  mod(N,8)~=0 or mod(M,8)~=0
% input image is resized.
imgN=size(img,1);
imgM=size(img,2);
modN=mod(imgN,8);
modM=mod(imgM,8);

%----------------------------------------
% save informations in informations.dat
if isa(img,'uint8')
    graylevmax=2^8-1;
end
if isa(img,'uint16')
    graylevmax=2^16-1;
end
if isa(img,'uint32')
    graylevmax=2^32-1;
end
save('informations.dat','graylevmax','img');
%-----------------------------------------
% resize
%-----------------------------------------
img=img(modN+1:imgN,modM+1:imgM);
%-----------------------------------------
img = double(img)/graylevmax;
set(get(ud.hOriginalAxes, 'title'), 'string', 'Original FingerPrint');
set(get(ud.hComponent1Axes, 'title'), 'string', '0 degree Component');
set(get(ud.hComponent6Axes, 'title'), 'string', '112.5 degree Component');
set(ud.hOriginalImage, 'Cdata', img);
set(DemoFig,'Pointer','arrow');
setstatus(DemoFig,'Please select a step to process...');
return;

%========================================
%%%
%%%  Sub-Function - SelectExtractingStep
%%%

function SelectExtractingStep(DemoFig)
% Load a step

if nargin<1
    DemoFig = gcbf;
end

set(DemoFig,'Pointer','watch');
ud=get(DemoFig,'Userdata');
v = get(ud.hSelectStepPop,{'value','String'});
name = deblank(v{2}(v{1},:));
drawnow

switch name
    case 'Centralize',
        Centralize(DemoFig);
    case 'Crop',
        Crop(DemoFig);
    case 'Sectorize',
        Sectorize(DemoFig);
    case 'Normalize',
        Normalize(DemoFig);
    case 'Gabor filters',
        Gaborfilter(DemoFig);
    case 'Convolute',
        Convolute(DemoFig);
    case 'Features',
        Features(DemoFig);
    case 'FingerCode',
        Fingercode(DemoFig);
    case 'Check',
        Check(DemoFig);
    otherwise 
        error('fpextractdemo: Unknown Image Option!');
end

return;

%==========================================================================
%%%
%%%  Sub-Function - Centralize
%%%

function Centralize(DemoFig)

load 'informations.dat' -mat

if nargin<1
    DemoFig = gcbf;
end
set(DemoFig,'Pointer','watch');
setstatus(DemoFig,'Centralizing..., please wait !!!');
ud=get(DemoFig,'Userdata');
fingerprint = getimage(ud.hOriginalImage);
fingerprint = fingerprint*graylevmax;

[BinarizedPrint,XofCenter,YofCenter] = centralizing(fingerprint,0);

set(get(ud.hComponent8Axes, 'title'), 'string', 'Binarized Print');
set(ud.hComponent8Image, 'Cdata', BinarizedPrint);
set(DemoFig,'Pointer','arrow');
setstatus(DemoFig,'Finished centralization');
ud.OriginalImageIsStale = 0;
set(DemoFig, 'UserData', ud);
drawnow
%==========================================================================
%%%
%%%  Sub-Function - Crop
%%%

function Crop(DemoFig)
% 
load 'informations.dat' -mat
if nargin<1
    DemoFig = gcbf;
end
set(DemoFig,'Pointer','watch');
setstatus(DemoFig,'Cropping..., please wait !!!');
ud=get(DemoFig,'Userdata');
fingerprint = getimage(ud.hOriginalImage);

fingerprint = fingerprint*graylevmax;

[BinarizedPrint,XofCenter,YofCenter]=centralizing(fingerprint,0);
[CroppedPrint]=cropping(XofCenter,YofCenter,fingerprint);

CroppedPrint = double(CroppedPrint)/graylevmax;
set(get(ud.hComponent1Axes, 'title'), 'string', 'Cropped Print');
set(ud.hComponent1Image, 'Cdata', CroppedPrint);
set(get(ud.hComponent8Axes, 'title'), 'string', 'Binarized Print');
set(ud.hComponent8Image, 'Cdata', BinarizedPrint);
set(DemoFig,'Pointer','arrow');
setstatus(DemoFig,'Finished Crop');
ud.Component1ImageIsStale = 0;
set(DemoFig, 'UserData', ud);
drawnow
%==========================================================================
%%%
%%%  Sub-Function - Sectorize
%%%

function Sectorize(DemoFig)
% 
load 'informations.dat' -mat

if nargin<1
    DemoFig = gcbf;
end
set(DemoFig,'Pointer','watch');
setstatus(DemoFig,'Sectorizing..., please wait !!!');
ud=get(DemoFig,'Userdata');
fingerprint = getimage(ud.hOriginalImage);

fingerprint = fingerprint*graylevmax;

[BinarizedPrint,XofCenter,YofCenter]=centralizing(fingerprint,0);
[CroppedPrint]=cropping(XofCenter,YofCenter,fingerprint);
for ( i=1:1:175*175)
    tmp=CroppedPrint(i);
    CroppedPrint(i)=whichsector(i);
    if (CroppedPrint(i)==36 | CroppedPrint(i)==37)
        CroppedPrint(i)=tmp/graylevmax;
    else 
        CroppedPrint(i)=CroppedPrint(i)/64;
    end
    
end

set(get(ud.hComponent1Axes, 'title'), 'string', 'SectorizedPrint');
set(ud.hComponent1Image, 'Cdata', CroppedPrint);
set(DemoFig,'Pointer','arrow');
setstatus(DemoFig,'Finished Sectorization');
ud.Component1ImageIsStale = 0;
set(DemoFig, 'UserData', ud);
drawnow
%==========================================================================
%%%
%%%  Sub-Function - Normalize
%%%

function Normalize(DemoFig)
% 
load 'informations.dat' -mat

if nargin<1
    DemoFig = gcbf;
end
set(DemoFig,'Pointer','watch');
setstatus(DemoFig,'Normalizing..., please wait !!!');
ud=get(DemoFig,'Userdata');
fingerprint = getimage(ud.hOriginalImage);

fingerprint = fingerprint*graylevmax;


[BinarizedPrint,XofCenter,YofCenter]=centralizing(fingerprint,0);
[CroppedPrint]=cropping(XofCenter,YofCenter,fingerprint);
[NormalizedPrint,vector] = sector_norm( CroppedPrint , 0 , 0);

CroppedPrint = double(CroppedPrint)/graylevmax;
NormalizedPrint = double(NormalizedPrint)/100;
set(get(ud.hComponent1Axes, 'title'), 'string', 'Cropped Print');
set(ud.hComponent1Image, 'Cdata', CroppedPrint);
set(get(ud.hComponent8Axes, 'title'), 'string', 'Normalized Print');
set(ud.hComponent8Image, 'Cdata', NormalizedPrint);
set(DemoFig,'Pointer','arrow');
setstatus(DemoFig,'Finished normalization');
ud.Component1ImageIsStale = 0;
set(DemoFig, 'UserData', ud);
drawnow
%==========================================================================
%%%
%%%  Sub-Function - Gaborfilter
%%%

function Gaborfilter(DemoFig)
% 

if nargin<1
    DemoFig = gcbf;
end
set(DemoFig,'Pointer','watch');
setstatus(DemoFig,'Gabor filter will be shown..., please wait !!!');
ud=get(DemoFig,'Userdata');

num_disk=8;

for (angle=0:1:num_disk-1)
    
    gabor=gabor2d_sub(angle,num_disk);
    gabor=gabor*128;
    switch angle<num_disk
        case (angle==0),
            set(get(ud.hComponent1Axes, 'title'), 'string', '0 degree gabor');
            set(ud.hComponent1Image, 'Cdata', gabor);
        case (angle==1),
            set(get(ud.hComponent2Axes, 'title'), 'string', '22.5 degree gabor');
            set(ud.hComponent2Image, 'Cdata', gabor);
        case (angle==2),
            set(get(ud.hComponent3Axes, 'title'), 'string', '45 degree gabor');
            set(ud.hComponent3Image, 'Cdata', gabor);
        case (angle==3),
            set(get(ud.hComponent4Axes, 'title'), 'string', '67.5 degree gabor');
            set(ud.hComponent4Image, 'Cdata', gabor);
        case (angle==4),
            set(get(ud.hComponent5Axes, 'title'), 'string', '90 degree gabor');
            set(ud.hComponent5Image, 'Cdata', gabor);
        case (angle==5),
            set(get(ud.hComponent6Axes, 'title'), 'string', '112.5 degree gabor');
            set(ud.hComponent6Image, 'Cdata', gabor);
        case (angle==6),
            set(get(ud.hComponent7Axes, 'title'), 'string', '135 degree gabor');
            set(ud.hComponent7Image, 'Cdata', gabor);
        case (angle==7),
            set(get(ud.hComponent8Axes, 'title'), 'string', '157.5 degree gabor');
            set(ud.hComponent8Image, 'Cdata', gabor);
        otherwise 
            error('Nothing !');
    end
    
end

set(DemoFig,'Pointer','arrow');
setstatus(DemoFig,'Gabor Filters were shown');
ud.OriginalImageIsStale = 0;
set(DemoFig, 'UserData', ud);
drawnow
%==========================================================================
%%%
%%%  Sub-Function - Convolute
%%%

function Convolute(DemoFig)
% 
load 'informations.dat' -mat

if nargin<1
    DemoFig = gcbf;
end
set(DemoFig,'Pointer','watch');
setstatus(DemoFig,'Convoluting with eight Gabor filters in process...');
ud=get(DemoFig,'Userdata');
fingerprint = getimage(ud.hOriginalImage);

load 'informations.dat' -mat
fingerprint = fingerprint*graylevmax;

N=175;
num_disk=8;

[BinarizedPrint,XofCenter,YofCenter]=centralizing(fingerprint,0);
[CroppedPrint]=cropping(XofCenter,YofCenter,fingerprint);
[NormalizedPrint,vector]=sector_norm(CroppedPrint,0,1);
for (angle=0:1:num_disk-1)
    
    gabor=gabor2d_sub(angle,num_disk);    
    z2=gabor;
    z1=NormalizedPrint;
    z1x=size(z1,1);
    z1y=size(z1,2);
    z2x=size(z2,1);
    z2y=size(z2,2);    
    ComponentPrint=real(ifft2(fft2(z1,z1x+z2x-1,z1y+z2y-1).*fft2(z2,z1x+z2x-1,z1y+z2y-1)));    
    px=((z2x-1)+mod((z2x-1),2))/2;
    py=((z2y-1)+mod((z2y-1),2))/2;
    ComponentPrint=ComponentPrint(px+1:px+z1x,py+1:py+z1y);
    
    
    [disk,vector]=sector_norm(ComponentPrint,1,0);
    img = double(ComponentPrint)/graylevmax;
    
    switch angle<8
        case (angle==0),
            set(get(ud.hComponent1Axes, 'title'), 'string', '0 degree Component');
            set(ud.hComponent1Image, 'Cdata', img);
        case (angle==1),
            set(get(ud.hComponent2Axes, 'title'), 'string', '22.5 degree Component');
            set(ud.hComponent2Image, 'Cdata', img);
        case (angle==2),
            set(get(ud.hComponent3Axes, 'title'), 'string', '45 degree Component');
            set(ud.hComponent3Image, 'Cdata', img);
        case (angle==3),
            set(get(ud.hComponent4Axes, 'title'), 'string', '67.5 degree Component');
            set(ud.hComponent4Image, 'Cdata', img);
        case (angle==4),
            set(get(ud.hComponent5Axes, 'title'), 'string', '90 degree Component');
            set(ud.hComponent5Image, 'Cdata', img);
        case (angle==5),
            set(get(ud.hComponent6Axes, 'title'), 'string', '112.5 degree Component');
            set(ud.hComponent6Image, 'Cdata', img);
        case (angle==6),
            set(get(ud.hComponent7Axes, 'title'), 'string', '135 degree Component');
            set(ud.hComponent7Image, 'Cdata', img);
        case (angle==7),
            set(get(ud.hComponent8Axes, 'title'), 'string', '157.5 degree Component');
            set(ud.hComponent8Image, 'Cdata', img);
        otherwise 
            error('Nothing !');
    end
    
end

set(DemoFig,'Pointer','arrow');
setstatus(DemoFig,'Finished Convolution');
set(DemoFig, 'UserData', ud);
drawnow
%==========================================================================
%%%
%%%  Sub-Function - Features
%%%

function Features(DemoFig)
% 
load 'informations.dat' -mat

if nargin<1
    DemoFig = gcbf;
end
set(DemoFig,'Pointer','watch');
setstatus(DemoFig,'Convoluting with eight Gabor filters in process...');
ud=get(DemoFig,'Userdata');
fingerprint = getimage(ud.hOriginalImage);

fingerprint = fingerprint*graylevmax;

N=175;
num_disk=8;

[BinarizedPrint,XofCenter,YofCenter]=centralizing(fingerprint,0);
[CroppedPrint]=cropping(XofCenter,YofCenter,fingerprint);
[NormalizedPrint,vector]=sector_norm(CroppedPrint,0,1);

for (angle=0:1:num_disk-1)
    
    gabor=gabor2d_sub(angle,num_disk);
    z2=gabor;
    z1=NormalizedPrint;
    z1x=size(z1,1);
    z1y=size(z1,2);
    z2x=size(z2,1);
    z2y=size(z2,2);    
    ComponentPrint=real(ifft2(fft2(z1,z1x+z2x-1,z1y+z2y-1).*fft2(z2,z1x+z2x-1,z1y+z2y-1)));    
    px=((z2x-1)+mod((z2x-1),2))/2;
    py=((z2y-1)+mod((z2y-1),2))/2;
    ComponentPrint=ComponentPrint(px+1:px+z1x,py+1:py+z1y);
    
    [disk,vector]=sector_norm(ComponentPrint,1,0);
    
    img = double(ComponentPrint)/graylevmax;
    img1 = double(disk)/51200;
    switch angle<8
        case (angle==0),
            set(get(ud.hComponent1Axes, 'title'), 'string', '0 degree Features');
            set(ud.hComponent1Image, 'Cdata', img1);
        case (angle==1),
            set(get(ud.hComponent2Axes, 'title'), 'string', '22.5 degree Features');
            set(ud.hComponent2Image, 'Cdata', img1);
        case (angle==2),
            set(get(ud.hComponent3Axes, 'title'), 'string', '45 degree Features');
            set(ud.hComponent3Image, 'Cdata', img1);
        case (angle==3),
            set(get(ud.hComponent4Axes, 'title'), 'string', '67.5 degree Features');
            set(ud.hComponent4Image, 'Cdata', img1);
        case (angle==4),
            set(get(ud.hComponent5Axes, 'title'), 'string', '90 degree Features');
            set(ud.hComponent5Image, 'Cdata', img1);
        case (angle==5),
            set(get(ud.hComponent6Axes, 'title'), 'string', '112.5 degree Features');
            set(ud.hComponent6Image, 'Cdata', img1);
        case (angle==6),
            set(get(ud.hComponent7Axes, 'title'), 'string', '135 degree Features');
            set(ud.hComponent7Image, 'Cdata', img1);
        case (angle==7),
            set(get(ud.hComponent8Axes, 'title'), 'string', '157.5 degree Features');
            set(ud.hComponent8Image, 'Cdata', img1);
        otherwise 
            error('Nothing !');
    end
    
end

set(DemoFig,'Pointer','arrow');
setstatus(DemoFig,'Features were extracted');
set(DemoFig, 'UserData', ud);
drawnow
%==========================================================================

%==========================================================================
%%%
%%%  Sub-Function - FingerCode
%%%

function Fingercode(DemoFig)
% 
load 'informations.dat' -mat

if nargin<1
    DemoFig = gcbf;
end


set(DemoFig,'Pointer','watch');
setstatus(DemoFig,'FingerCode in process...');
ud=get(DemoFig,'Userdata');
fingerprint = getimage(ud.hOriginalImage);

fingerprint = fingerprint*graylevmax;

N=175;
num_disk=8;


[BinarizedPrint,XofCenter,YofCenter]=centralizing(fingerprint,0);
[CroppedPrint]=cropping(XofCenter,YofCenter,fingerprint);
[NormalizedPrint,vector]=sector_norm(CroppedPrint,0,1);

for (angle=0:1:num_disk-1)    
    gabor=gabor2d_sub(angle,num_disk);
    z2=gabor;
    z1=NormalizedPrint;
    z1x=size(z1,1);
    z1y=size(z1,2);
    z2x=size(z2,1);
    z2y=size(z2,2);    
    ComponentPrint=real(ifft2(fft2(z1,z1x+z2x-1,z1y+z2y-1).*fft2(z2,z1x+z2x-1,z1y+z2y-1)));    
    px=((z2x-1)+mod((z2x-1),2))/2;
    py=((z2y-1)+mod((z2y-1),2))/2;
    ComponentPrint=ComponentPrint(px+1:px+z1x,py+1:py+z1y);        
    [disk,vector]=sector_norm(ComponentPrint,1,0);    
    %img = double(ComponentPrint)/graylevmax;
    %img1 = double(disk)/51200;
    finger_code1{angle+1}=vector(1:36);
end

load('informations.dat','img','-mat');
img=imrotate(img,22.5/2);
imgN=size(img,1);
imgM=size(img,2);
modN=mod(imgN,8);
modM=mod(imgM,8);
fingerprint=double(img(modN+1:imgN,modM+1:imgM));

[BinarizedPrint,XofCenter,YofCenter]=centralizing(fingerprint,0);
[CroppedPrint]=cropping(XofCenter,YofCenter,fingerprint);
[NormalizedPrint,vector]=sector_norm(CroppedPrint,0,1);

for (angle=0:1:num_disk-1)    
    gabor=gabor2d_sub(angle,num_disk);
    z2=gabor;
    z1=NormalizedPrint;
    z1x=size(z1,1);
    z1y=size(z1,2);
    z2x=size(z2,1);
    z2y=size(z2,2);    
    ComponentPrint=real(ifft2(fft2(z1,z1x+z2x-1,z1y+z2y-1).*fft2(z2,z1x+z2x-1,z1y+z2y-1)));    
    px=((z2x-1)+mod((z2x-1),2))/2;
    py=((z2y-1)+mod((z2y-1),2))/2;
    ComponentPrint=ComponentPrint(px+1:px+z1x,py+1:py+z1y);        
    [disk,vector]=sector_norm(ComponentPrint,1,0);    
    %img = double(ComponentPrint)/graylevmax;
    %img1 = double(disk)/51200;
    finger_code2{angle+1}=vector(1:36);
end
% FingerCode added to database
if (exist('fp_database.dat')==2)
    load('fp_database.dat','-mat');
    fp_number=fp_number+1;
    data{fp_number,1}=finger_code1;
    data{fp_number,2}=finger_code2;
    save('fp_database.dat','data','fp_number','-append');
else
    fp_number=1;
    data{fp_number,1}=finger_code1;
    data{fp_number,2}=finger_code2;
    save('fp_database.dat','data','fp_number');
end

message=strcat('FingerCode was succesfully added to database. Fingerprint no. ',num2str(fp_number));
msgbox(message,'FingerCode DataBase','help');

set(DemoFig,'Pointer','arrow');
setstatus(DemoFig,'FingerCode calculated');
set(DemoFig, 'UserData', ud);
drawnow
%==========================================================================

%==========================================================================
%%%
%%%  Sub-Function - Check
%%%

function Check(DemoFig)
% 
load 'informations.dat' -mat

if nargin<1
    DemoFig = gcbf;
end


set(DemoFig,'Pointer','watch');
setstatus(DemoFig,'DataBase Scanning...');
ud=get(DemoFig,'Userdata');
fingerprint = getimage(ud.hOriginalImage);

fingerprint = fingerprint*graylevmax;

N=175;
num_disk=8;


[BinarizedPrint,XofCenter,YofCenter]=centralizing(fingerprint,0);
[CroppedPrint]=cropping(XofCenter,YofCenter,fingerprint);
[NormalizedPrint,vector]=sector_norm(CroppedPrint,0,1);

for (angle=0:1:num_disk-1)    
    gabor=gabor2d_sub(angle,num_disk);
    z2=gabor;
    z1=NormalizedPrint;
    z1x=size(z1,1);
    z1y=size(z1,2);
    z2x=size(z2,1);
    z2y=size(z2,2);    
    ComponentPrint=real(ifft2(fft2(z1,z1x+z2x-1,z1y+z2y-1).*fft2(z2,z1x+z2x-1,z1y+z2y-1)));    
    px=((z2x-1)+mod((z2x-1),2))/2;
    py=((z2y-1)+mod((z2y-1),2))/2;
    ComponentPrint=ComponentPrint(px+1:px+z1x,py+1:py+z1y);        
    [disk,vector]=sector_norm(ComponentPrint,1,0);    
    %img = double(ComponentPrint)/graylevmax;
    %img1 = double(disk)/51200;
    finger_code{angle+1}=vector(1:36);
end


% FingerCode of input fngerprint has just been calculated.
% Checking with BadaBase
if (exist('fp_database.dat')==2)
    load('fp_database.dat','-mat');
    %---- alloco memoria -----------------------------------
    ruoto1=zeros(36,1);
    ruoto2=zeros(36,1);
    vettore_d1=zeros(12,1);
    vettore_d2=zeros(12,1);
    best_matching=zeros(fp_number,1);
    % start checking ---------------------------------------
    for scanning=1:fp_number
        fcode1=data{scanning,1};
        fcode2=data{scanning,2};
        for rotazione=0:1:11
            d1=0;
            d2=0;
            for disco=1:8
                f1=fcode1{disco};
                f2=fcode2{disco};
                % ora ruoto f1 ed f2  della rotazione ciclica ----------
                for old_pos=1:12
                    new_pos=mod(old_pos+rotazione,12);
                    if (new_pos==0)
                        new_pos=12;
                    end
                    ruoto1(new_pos)=f1(old_pos);
                    ruoto1(new_pos+12)=f1(old_pos+12);
                    ruoto1(new_pos+24)=f1(old_pos+24);
                    ruoto2(new_pos)=f2(old_pos);
                    ruoto2(new_pos+12)=f2(old_pos+12);
                    ruoto2(new_pos+24)=f2(old_pos+24);
                end
                %-------------------------------------------------------
                d1=d1+norm(finger_code{disco}-ruoto1);
                d2=d2+norm(finger_code{disco}-ruoto2);                
            end
            vettore_d1(rotazione+1)=d1;
            vettore_d2(rotazione+1)=d2;
        end
        [min_d1,pos_min_d1]=min(vettore_d1);
        [min_d2,pos_min_d2]=min(vettore_d2);
        if min_d1<min_d2
            minimo=min_d1;
        else
            minimo=min_d2;
        end
        best_matching(scanning)=minimo;
    end
    [distanza_minima,posizione_minimo]=min(best_matching);
    beep;
    message=strcat('The nearest fingerprint present in DataBase which matchs input fingerprint is  : ',num2str(posizione_minimo),...
                         ' with a distance of : ',num2str(distanza_minima));
    msgbox(message,'DataBase Info','help');
    %-------------------------------------------------------
    
else
    message='DataBase is empty. No check is possible.';
    msgbox(message,'FingerCode DataBase Error','warn');    
end



set(DemoFig,'Pointer','arrow');
setstatus(DemoFig,'DataBase Scanning completed');
set(DemoFig, 'UserData', ud);
drawnow
%==========================================================================