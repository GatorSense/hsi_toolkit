function [smf_out,mu,siginv] = smf_detector(hsi_img,tgt_sig,mask,mu,siginv)
%
%function [smf_out,mu,siginv] = smf_detector(hsi_img,tgt_sig,mask,mu,siginv)
%
% Spectral Matched Filter
%
% inputs:
%  hsi_image - n_row x n_col x n_band hyperspectral image
%  tgt_sig - target signature (n_band x 1 - column vector)
%  mask - binary image limiting detector operation to pixels where mask is true
%         if not present or empty, no mask restrictions are used
%  mu - (optional) mean for filter (if not provided, computed from image)
%  siginv - (optional) inverse covariance for filter (if not provided, computed from image)
%
% outputs:
%  smf_out - detector image
%  mu - mean of background
%  siginv - inverse covariance of background
%
% 8/8/2012 - Taylor C. Glenn
% 6/2/2018 - Edited by Alina Zare
%

if ~exist('mask','var'); mask = []; end
if ~exist('mu','var'), mu = []; end
if ~exist('siginv','var'), siginv = []; end
addpath(fullfile('..','util'));

[smf_out,mu,siginv] = img_det(@smf_detector_array,hsi_img,tgt_sig,mask,mu,siginv);

end


function [smf_data,mu,siginv] = smf_detector_array(hsi_data,tgt_sig,mu,siginv)

%function [smf_data,mu,siginv] = smf_detector_array(hsi_data,tgt_sig,mu,siginv)
%
% Spectral Matched Filter for array (not image) data
%
% inputs:
%  hsi_data - n_spectra x n_band array of hyperspectral data
%  tgt_sig - target signature (n_band x 1 - column vector)
%  mu - (optional) mean for filter
%  siginv - (optional) inverse covariance for filter
%
% outputs:
%  smf_out - detector output per spectra
%  mu - mean of background
%  siginv - inverse covariance of background
%
% 8/19/2012 - Taylor C. Glenn
% 6/2/2018 - Edited by Alina Zare
%

% alternative formulation from Eismann's book, pp 653
%  subtract mean from target to reduce effect of additive model
%  on hyperspectral (non additive) data
%  also take positive square root of filter

if ~exist('mu','var') || isempty(mu)
    mu = mean(hsi_data,2);
end
if ~exist('siginv','var') || isempty(siginv)
    siginv = pinv(cov(hsi_data'));
end

s = tgt_sig - mu;
z = hsi_data - repmat(mu,[1 size(hsi_data,2)]);
f = s'*siginv / sqrt(s'*siginv*s);

smf_data = f*z;

end
