function change_orient(filename,s)
%Function to orient a nifti/analyze file so that
%matrices with differing orientations (but same dimensions) can be matched
%in orientation for calculations
%
%Usage
%   change_orient(filename,s)
%       where filename is the name of the 'img' or the 'nii' file
%             s is the string specifying orientation (las or ras or rpi etc)
%
%This is meant for single volumes only (not 4D) There are some issues with
%4D that has not been sorted out yet
%
%For native analyze formats that do not have any orientation information,
%it is assumed that the original orientation is RAS.

%Suresh E Joel, Sep 13, 2008

%Arguement number check
if(nargin<2), error('Too few arguements'); end;

v=spm_vol(filename);

%Make sure its only single volume
try
    size(v.mat);
catch
    error('Can use only 3D, not 4D');
end;

%Read the image data
i=spm_read_vols(v);

%Convert the orientation string to lowercase so its easy to compare
s=lower(s);

if( ((v.mat(1,1)<0) && s(1)=='r') || ((v.mat(1,1)>0) && s(1)=='l') )
    i=i(end:-1:1,:,:);
    v.private.mat(1,:)=-v.private.mat(1,:);
    v.private.mat0=v.private.mat;
    v.mat=v.private.mat;
end;
if( ((v.mat(2,2)<0) && s(2)=='a') || ((v.mat(2,2)>0) && s(2)=='p') )
    i=i(:,end:-1:1,:);
    v.private.mat(2,:)=-v.private.mat(2,:);
    v.private.mat0=v.private.mat;
    v.mat=v.private.mat;
end;
if( ((v.mat(3,3)<0) && s(3)=='s') || ((v.mat(1,1)>0) && s(3)=='i') )
    i=i(:,:,end:-1:1);
    v.private.mat(3,:)=-v.private.mat(3,:);
    v.private.mat0=v.private.mat;
    v.mat=v.private.mat;
end;

%Write out the file
spm_write_vol(v,i);
