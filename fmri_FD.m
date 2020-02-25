function FD = fmri_FD(rp_file, out_file, radius)
%Function to calculate Framewise Displacement (FD) (Power et al., 2012)from
%the six realignment parameters.  FD is calculated by summing the absolute
%value of the differenced (time_t-time_t-1) translational realignment
%parameters and the three differenced rotational parameters, which are
%converted from radians to millimeters by assuming a brain radius of 50 mm.
%
%Usage: FD = fmri_FD(rp_file, out_file)
%where  rp_file is the path/file name containing the realignment parameters
%       in columns (such as what is obtained from spm_realign or mcflirt)
%       out_file is the path/file name to be written

if(~exist('radius', 'var'))
    radius = 50;
    disp('Brain radius not defined; using default of 50 mm');
end

dat = load(rp_file);
dat = dat(:, 1:6);
order = 1;
diff_dat = abs([[0 0 0 0 0 0]; diff(dat, order, 1)]);
% 	Multiply by 50mm brain;
diff_dat(:,4:6) = diff_dat(:,4:6) * radius;
FD = sum(diff_dat, 2);

if(exist('out_file', 'var'))
writetable(table(FD), out_file)
end






