function [CroppedPrint] = cropping(XofCenter,YofCenter,CentralizedPrint)
% Modified by Luigi Rosa

N = 175;
M=size(CentralizedPrint,1);

imgN=size(CentralizedPrint,1);
imgM=size(CentralizedPrint,2);
  
if (YofCenter+30) <= M
   YofCenter = YofCenter + 20;
else
   YofCenter = M;
end

X=XofCenter-floor(N/2);
Y=YofCenter-floor(N/2);


%center point can be shown on command line
%         ----------------
%         |              |
%         |              |
%         |       .      | y=row
%         |              |
%         |              |
%         ----------------
%                x=column
%XofCenter%   column of matrix
%YofCenter%   row of matrix
  
%-------------------------------------------------------------
%          if A= 1 2 3
%                4 5 6
%                7 8 9
%             B=A(1:2,2:3)
%              = 2 3
%                5 6
%       creates B by extracting the first twos and last two 
%       columns of A
%-------------------------------------------------------------
if (YofCenter-floor(N/2)<1)||(YofCenter+floor(N/2)>imgN)||(XofCenter-floor(N/2)<1)||(XofCenter+floor(N/2)>imgM)
    message='Cropping error: when the input image is cropped an error occurs: a possible error during center point determination.';
    msgbox(message,'Cropping Error','warn');   
    CroppedPrint=zeros(175);
    return;
else
    CroppedPrint=CentralizedPrint(YofCenter-floor(N/2):YofCenter+floor(N/2),XofCenter-floor(N/2):XofCenter+floor(N/2));
end