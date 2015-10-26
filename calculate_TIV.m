function calculate_TIV()
%
%        calculate_TIV
%        Copyright (C) 2012 Johns Hopkins University
%        Software by Cliff Workman
%
%        Usage: calculate_TIV
%
%        In the file selector that opens, select a single participant's GM,
%        WM, and CSF segments. This script will calculate total
%        intracranial volume (TIV) for use with VBM.
%
%        This is an adapatation of a script written by John Ashburner:
%        https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=spm;bf3307b2.0807
%        (If the link appears broken, try opening it in a system browser.)

V    = spm_vol(spm_select(Inf,'Image'));

for ii=1:1:size(V,1),
    [~,tissue_name,~] = fileparts(V(ii).fname);
    disp(tissue_name);
    clear tissue_name
end

Vols = zeros(numel(V),1);
for j=1:numel(V),
    tot = 0;
    for i=1:V(1).dim(3),
            img = spm_slice_vol(V(j),spm_matrix(...
                  [0 0 i]),V(j).dim(1:2),0);
            img = img(isfinite(img)); % <-- exclude non-finite values
            tot = tot + sum(img(:));
    end;
    voxvol = abs(det(V(j).mat))/100^3; % volume of a voxel, in litres
    Vols(j) = tot*voxvol;
end
Total = sum(Vols)
disp([Total ' litres']);
end