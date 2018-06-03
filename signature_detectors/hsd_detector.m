function [hsd_out,tgt_p] = hsd_detector(hsi_img,tgt_sig,mask,ems, siginv)
%
%function [hsd_out,tgt_p] = hsd_detector(hsi_img,tgt_sig,mask,ems)
%
% Hybrid Structured Detector 
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
%  siginv - background inverse covariance (n_band x n_band matrix) 
%
% outputs:
%  hsd_out - detector image
%  tgt_p - target proportion in unmixing
%
% 8/19/2012 - Taylor C. Glenn
% 6/2/2018 - Edited by Alina Zare
%

if ~exist('siginv','var'), siginv = []; end
addpath(fullfile('..','util'));

[hsd_out,tgt_p] = img_det(@hsd_det,hsi_img,tgt_sig,mask,ems,siginv);

end

function [hsd_data,tgt_p] = hsd_det(hsi_data,tgt_sig,ems,siginv)

[~,n_pix] = size(hsi_data);

% unmix data with only background endmembers
P = unmix(hsi_data,ems); 

% unmix data with target signature as well
targ_P = unmix(hsi_data,[tgt_sig ems]); 

if isempty(siginv)
    siginv = pinv(cov(hsi_data'));
end

hsd_data = zeros(1,n_pix);

for i=1:n_pix
    z = hsi_data(:,i) - ems*P(i,:)';
    w = hsi_data(:,i) - [tgt_sig ems]*targ_P(i,:)';
    hsd_data(i) = (z'*siginv*z) / (w'*siginv*w);
end

tgt_p = targ_P(:,1:size(tgt_sig,2))';

end

