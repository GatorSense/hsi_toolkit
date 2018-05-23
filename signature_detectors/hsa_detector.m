function [hsa_out] = hsa_detector(hsi_img,tgt_sig,mask,ems,n_comp)
%
%function [hsa_out] = hsa_detector(hsi_img,tgt_sig,mask,ems,n_comp)
%
% Hybrid Structured/Abundance Detector
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
%  hsa_out - detector image
%
% 8/19/2012 - Taylor C. Glenn - tcg@cise.ufl.edu
%

hsa_out = img_det(@hsa_det,hsi_img,tgt_sig,mask,ems,n_comp);

end

function hsa_data = hsa_det(hsi_data,tgt_sig,ems,n_comp)

[n_band,n_pix] = size(hsi_data);

%siginv = pinv(cov(hsi_data'));

params = struct();
params.sum_to_one = true;

% unmix data with only background endmembers

P = unmix2(hsi_data,ems,params); %unmix2 from FUMI directory currently

% find the covariance of the background unmixing residual
res = zeros(n_band,n_pix);
for i=1:n_band
    res(i,:) = hsi_data(i,:) - ems(i,:)*P';
end

siginv = pinv(cov(res'));

% unmix data with target signature as well

targ_P = unmix2(hsi_data,[tgt_sig ems],params); 

% model background and target abundance distributions
gmm_bg = gmdistribution.fit(P,n_comp,'Replicates',1,'Regularize',1e-6);

n_em = size(ems,2);
ll_bg = log(pdf(gmm_bg,P))'; %/(n_em * n_comp);

ps = [0.05 0.95];
p_tgt = ps( (targ_P(:,1)' > 0.05) + 1);
ll_tgt = log(p_tgt);

hsd_data = zeros(1,n_pix);

for i=1:n_pix
    z = res(:,i);
    w = hsi_data(:,i) - [tgt_sig ems]*targ_P(i,:)';
    hsd_data(i) = (z'*siginv*z) / (w'*siginv*w);
end

hsd_rg = max(hsd_data) - min(hsd_data);

max_tgt = max(ll_tgt);
min_tgt = min(ll_tgt);
z_tgt = (ll_tgt - min_tgt)/(max_tgt - min_tgt);

max_bg = max(ll_bg);
min_bg = min(ll_bg);
z_bg = (ll_bg - min_bg)/(max_bg - min_bg);

hsa_data = hsd_data + z_tgt*(hsd_rg/3) - z_bg*(hsd_rg/3);

save hsa_out;
%hsa_data = hsd_data .* exp(ll_tgt - ll_bg);

end
