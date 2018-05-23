function [hs_out] = hs_detector(hsi_img,tgt_sig,mask,ems)
%
%function [hs_out] = hs_detector(hsi_img,tgt_sig,mask,ems)
%
% Hybrid Subpixel Detector
%
% ref:
% A hybrid algorithm for subpixel detection in hyperspectral imagery (inproceedings)
% Broadwater, J. and Meth, R. and Chellappa, R.
% Geoscience and Remote Sensing Symposium, 2004. IGARSS '04. Proceedings. 2004 IEEE International
% 2004 Volume 3 Pages 1601 -1604 vol.3 Month sept.
%
% inputs:
%  hsi_image - n_row x n_col x n_band hyperspectral image
%  tgt_sig - target signature (n_band x 1 - column vector)
%  mask - binary image limiting detector operation to pixels where mask is true
%         if not present or empty, no mask restrictions are used
%  ems - background endmembers
%
% outputs:
%  hs_out - detector image
%
% 8/18/2012 - Taylor C. Glenn - tcg@cise.ufl.edu
%

hs_out = img_det(@hs_det,hsi_img,tgt_sig,mask,ems);

end

function hs_data = hs_det(hsi_data,tgt_sig,ems)

[n_band,n_pix] = size(hsi_data);

params = struct();
params.sum_to_one = true;
P = unmix2(hsi_data,ems,params); %unmix2 from FUMI directory currently

% unmix data with target signature as well

targ_P = unmix2(hsi_data,[tgt_sig ems],params); 


hs_data = zeros(1,n_pix);

for i=1:n_pix
    z = hsi_data(:,i) - ems*P(i,:)';
    w = hsi_data(:,i) - [tgt_sig ems]*targ_P(i,:)';
    hs_data(i) = (z'*z) / (w'*w);
end

end
