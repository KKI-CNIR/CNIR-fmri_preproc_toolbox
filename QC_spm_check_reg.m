function QC_spm_check_reg(images, output)
% A visual check of image registration quality.
% FORMAT spm_check_registration
% Orthogonal views of one or more images are displayed.  Clicking in
% any image moves the centre of the orthogonal views.  Images are
% shown in orientations relative to that of the first selected image.
% The first specified image is shown at the top-left, and the last at
% the bottom right.  The fastest increment is in the left-to-right
% direction (the same as you are reading this).
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% John Ashburner
% $Id: spm_check_registration.m 755 2007-02-28 18:15:10Z john $
%June 27, 2015 modified MB Nebel: added if statement to check the operating
%system and select appropriate driver

% if numel(output)
%     flag = 1;
% else
%     flag = 0;
% end

% if nargin==0,
% 	images = spm_select([1 15],'image','Select images');
% 	if size(images,1)<1, return; end;
% 	spm_check_registration(images);
% elseif nargin > 0,
fg = spm_figure('Findwin','Graphics');
if isempty(fg),
    fg=spm_figure('Create','Graphics');
    if isempty(fg),
        error('Cant create graphics window');
    end;
else
    spm_figure('Clear','Graphics');
end
if ischar(images)
    images=spm_vol(images);
end;
spm_orthviews('Reset');
mn = length(images);
n  = round(mn^0.4);
m  = ceil(mn/n);
w  = 1/n;
h  = 1/m;
ds = (w+h)*0.02;
for ij=1:mn
    i  = 1-h*(floor((ij-1)/n)+1);
    j  = w*rem(ij-1,n);
    handle(ij) = spm_orthviews('Image', images{ij},...
        [j+ds/2 i+ds/2 w-ds h-ds]); %#ok<AGROW>
    if ij==1, spm_orthviews('Space'); end;
    spm_orthviews('AddContext',handle(ij));
    spm_orthviews('Xhairs', 'on');
end;
h = gcf;
%     disp(flag);

%     if flag == 1
% Plot points for registration`

% if(isunix)
    print(h, '-dpsc2', output, '-append');
    spm_orthviews('Reposition', [12, -29, 80])
    print(h, '-dpsc2', output, '-append');
    spm_orthviews('Reposition', [-12, -29, 80])
    print(h, '-dpsc2', output, '-append');
    spm_orthviews('Reposition', [71, -31, 2])
    print(h, '-dpsc2', output, '-append');
    spm_orthviews('Reposition', [-71, -31, 2])
    print(h, '-dpsc2', output, '-append');
    spm_orthviews('Reposition', [13, -105, 2])
    print(h, '-dpsc2', output, '-append');
    spm_orthviews('Reposition', [-13, -105, 2])
    print(h, '-dpsc2', output, '-append');
    spm_orthviews('Reposition', [10, 71, 1])
    print(h, '-dpsc2', output, '-append');
    spm_orthviews('Reposition', [-10, 71, 1])
    print(h, '-dpsc2', output, '-append');
% else
%     print(h, '-dwinc', output, '-append');
%     spm_orthviews('Reposition', [12, -29, 80])
%     print(h, '-dwinc', output, '-append');
%     spm_orthviews('Reposition', [-12, -29, 80])
%     print(h, '-dwinc', output, '-append');
%     spm_orthviews('Reposition', [71, -31, 2])
%     print(h, '-dwinc', output, '-append');
%     spm_orthviews('Reposition', [-71, -31, 2])
%     print(h, '-dwinc', output, '-append');
%     spm_orthviews('Reposition', [13, -105, 2])
%     print(h, '-dwinc', output, '-append');
%     spm_orthviews('Reposition', [-13, -105, 2])
%     print(h, '-dwinc', output, '-append');
%     spm_orthviews('Reposition', [10, 71, 1])
%     print(h, '-dwinc', output, '-append');
%     spm_orthviews('Reposition', [-10, 71, 1])
%     print(h, '-dwinc', output, '-append');
% end
% print(h, '-dpdf', output, '-append');
% spm_orthviews('Reposition', [12, -29, 80])
% print(h, '-dpdf', output, '-append');
% spm_orthviews('Reposition', [-12, -29, 80])
% print(h, '-dpdf', output, '-append');
% spm_orthviews('Reposition', [71, -31, 2])
% print(h, '-dpdf', output, '-append');
% spm_orthviews('Reposition', [-71, -31, 2])
% print(h, '-dpdf', output, '-append');
% spm_orthviews('Reposition', [13, -105, 2])
% print(h, '-dpdf', output, '-append');
% spm_orthviews('Reposition', [-13, -105, 2])
% print(h, '-dpdf', output, '-append');
% spm_orthviews('Reposition', [10, 71, 1])
% print(h, '-dpdf', output, '-append');
% spm_orthviews('Reposition', [-10, 71, 1])
% print(h, '-dpdf', output, '-append');

%     end
% else
% 	error('Incorrect Usage');
% end;


