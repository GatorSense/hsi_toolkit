function [hsd_out,tgt_p] = hsd_detector(hsi_img,tgt_sig,mask,ems)
%
%function [hsd_out,tgt_p] = hsd_detector(hsi_img,tgt_sig,mask,ems)
%
% Hybrid Structured Detector Detector
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
%  hsd_out - detector image
%  tgt_p - target proportion in unmixing
%
% 8/19/2012 - Taylor C. Glenn - tcg@cise.ufl.edu
%

[hsd_out,tgt_p] = img_det(@hsd_det,hsi_img,tgt_sig,mask,ems);

end

function [hsd_data,tgt_p] = hsd_det(hsi_data,tgt_sig,ems)

[n_band,n_pix] = size(hsi_data);

% unmix data with only background endmembers
P = unmix(hsi_data,ems); 

% unmix data with target signature as well
targ_P = unmix(hsi_data,[tgt_sig ems]); 

siginv = pinv(cov(hsi_data'));
%siginv = eye(n_band);

hsd_data = zeros(1,n_pix);

for i=1:n_pix
    z = hsi_data(:,i) - ems*P(i,:)';
    w = hsi_data(:,i) - [tgt_sig ems]*targ_P(i,:)';
    hsd_data(i) = (z'*siginv*z) / (w'*siginv*w);
end

tgt_p = targ_P(:,1:size(tgt_sig,2))';

end

function [P] = unmix(data, endmembers)

%endmembers should be column vectors
X = data;

%number of endmembers
M = size(endmembers, 2);
%number of pixels
N = size(X, 2);

%set up constraint matrices
L = [];
k = [];
    
A = ones([1, M]);  % sum to one;
b = 1;

l = zeros(M,1);
u = ones(M,1);

P = zeros(N,M);

E = endmembers';
H = (2*E*E');

for i = 1:N

    Xi = X(:,i);

    F = (-2*Xi'*E')';

    P(i,:) = qpas(H, F, L, k, A, b, l, u); % from fast_spice/qpc directory
end

end
