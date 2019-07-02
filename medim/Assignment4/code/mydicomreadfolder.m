function [im,info] = mydicomreadfolder(foldername)
%[im,info] = mydicomreadfolder(foldername)
%Reads all dicom files in a folder into an image volume.
%
%-im is a three dimensional array of image data
%-info is a struct containing suitable data for voxel sizes etc.
%
%See also MYDICOMINFO, MYDICOMREAD
%
%This function is just a stub and you need to write it.

%Stub written by Einar Heiberg

%Hint:
%To get all files called in a folder, use the function 
%f = dir([foldername filesep '*.dcm'])

%Hint: Consider preallocating data for the sake of speed.

%Hint: waitbar is a good way of updating regarding progress.

%--- Initialize
im = [];
info = [];

%If called without input argument then ask for folder.
if nargin==0
  foldername = uigetdir(pwd, 'Select a folder');
end;

%Display folder name
disp(sprintf('Reading the folder %s.',foldername)); %#ok<DSPS>
first = true;
files = dir(foldername);
for k=1:length(files)
    filename = files(k).name;
    if (endsWith(filename, ".dcm"))
        [kinfo, kim] = mydicomread(files(k).folder +"/"+ filename);
        if first
            sz = size(kim);
            im = zeros(sz(1), sz(2), length(files)-2);
            info = kinfo;
            first = false;
        end
        im(:,:,k) = kim;
    end
end