function [hud_out] = hud_detector(hsi_img,tgt_sig,mask,ems)
%
%function [hud_out] = hud_detector(hsi_img,tgt_sig,mask,ems)
%
% Hybrid Unstructured Detector Detector
%  (from the department of redundancy department,
%   please enter your PIN number at the ATM machine)
%
%  ref:
%  Hybrid Detectors for Subpixel Targets
%  Broadwater, J. and Chellappa, R.
%  Pattern Analysis and Machine Intelligence, IEEE Transactions on
%  2007 Volume 29 Number 11 Pages 1891 -1903 Month nov.
%
% inputs:
%  hsi_image - n_row x n_col x n_band hyperspectral image
%  tgt_sig - target signature (n_band x 1 - column vector)
%  mask - binary image limiting detector operation to pixels where mask is true
%         if not present or empty, no mask restrictions are used
%  ems - background endmembers
%
% outputs:
%  hud_out - detector image
%
% 8/19/2012 - Taylor C. Glenn - tcg@cise.ufl.edu
%

hud_out = img_det(@hud_det,hsi_img,tgt_sig,mask,ems);

end

function hud_data = hud_det(hsi_data,tgt_sig,ems)

[n_band,n_pix] = size(hsi_data);

siginv = pinv(cov(hsi_data'));

params = struct();
params.sum_to_one = false;

% unmix data with target signature and background

targ_P = unmix2(hsi_data,[tgt_sig ems],params); 


hud_data = zeros(1,n_pix);

for i=1:n_pix
    s = tgt_sig*targ_P(i,1);
    x = hsi_data(:,i);
    %hud_data(i) = (x'*siginv*s)^2 / ((s'*siginv*s)*(x'*siginv*x));
    hud_data(i) = (x'*siginv*s) / (x'*siginv*x);
end

end
