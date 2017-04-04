function s=get_orient(filename)
%Function to get orientation information on a nifti/analyze file 
%
%Usage
%   s=get_orient(filename)
%       where filename is the name of the 'img' or the 'nii' file
%             s is the string specifying orientation (las or ras or rpi etc)
%
%This is meant for single volumes AXIAL slices only (not 4D) There are some
%issues with 4D that has not been sorted out yet
%

%Suresh E Joel, Sep 16, 2008

%Arguement number check
if(nargin<1), error('Too few arguements'); end;

v=spm_vol(filename);

%Make sure its only single volume
try
    size(v.mat);
catch
    error('Can use only 3D, not 4D');
end;

if(v.private.mat(1,1)<0), s='l'; else s='r'; end;
if(v.private.mat(2,2)<0), s=[s,'p']; else s=[s,'a']; end;
if(v.private.mat(3,3)<0), s=[s,'i']; else s=[s,'s']; end;