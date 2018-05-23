function [hua_out] = hua_detector(hsi_img,tgt_sig,mask,ems,n_comp)
%
%function [hua_out] = hua_detector(hsi_img,tgt_sig,mask,ems,n_comp)
%
% Hybrid Unstructured Abundance Detector
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
%  hua_out - detector image
%
% 8/19/2012 - Taylor C. Glenn - tcg@cise.ufl.edu
%

hua_out = img_det(@hua_det,hsi_img,tgt_sig,mask,ems,n_comp);

end

function hua_data = hua_det(hsi_data,tgt_sig,ems,n_comp)

[n_band,n_pix] = size(hsi_data);

siginv = pinv(cov(hsi_data'));

params = struct();
params.sum_to_one = true;

% unmix data with target signature and background

targ_P = unmix2(hsi_data,[tgt_sig ems],params); 

% model abundances with and without target
P = unmix2(hsi_data,ems,params);
gmm_bg = gmdistribution.fit(P,n_comp,'Replicates',1,'Regularize',1e-6);

% compute mixture likelihood ratio of each pixel
n_em = size(ems,2);
ll_bg = log(pdf(gmm_bg,P))'/(n_em * n_comp);
ll_tgt = targ_P(:,1)' > 0.05;

hud_data = zeros(1,n_pix);

for i=1:n_pix
    s = tgt_sig*targ_P(i,1);
    x = hsi_data(:,i);
    hud_data(i) = (x'*siginv*s) / (x'*siginv*x);
end

hud_rg = max(hud_data) - min(hud_data);
bg_rg = max(ll_bg) - min(ll_bg);

hua_data = hud_data + ll_tgt*(hud_rg/3) - ll_bg*(hud_rg/(3*bg_rg));

end
