function fmri_origin(funcfile)

temp_nii = load_nii(funcfile, [], [], [], [], [], 1);
temp_nii.hdr.hist.originator(1:3) = [temp_nii.hdr.dime.dim(2:4)/2];
save_nii(temp_nii, funcfile);
