function [acert_out,mu,siginv] = acert_max_detector(hsi_img,tgt_sigs,mask,mu,siginv)
%
%function [ace_out,mu,siginv] = ace_detector(hsi_img,tgt_sigs,mask,,mu,siginv)
%
% Adaptive Cosine/Coherence Estimator
%
% inputs:
%  hsi_image - n_row x n_col x n_band hyperspectral image
%  tgt_sig - target signature (n_band x n_sig - column vector)
%  mask - binary image limiting detector operation to pixels where mask is true
%         if not present or empty, no mask restrictions are used
%
% outputs:
%  ace_out - detector image
%  mu - mean of input data
%  siginv - inverse covariance of input data
%
% 8/8/2012 - Taylor C. Glenn - tcg@cise.ufl.edu
%

if ~exist('mask','var'), mask = []; end
if ~exist('mu','var'), mu = []; end
if ~exist('siginv','var'), siginv = []; end

[acert_out,mu,siginv] = img_det(@acert_max_det,hsi_img,tgt_sigs,mask,mu,siginv);

end

function [acert_data,mu,siginv] = acert_max_det(hsi_data,tgt_sigs,mu,siginv)

if isempty(mu)
    mu = mean(hsi_data,2);
end
if isempty(siginv)
    siginv = pinv(cov(hsi_data'));
end
n_sigs = size(tgt_sigs,2);
n_pix = size(hsi_data,2);

S = bsxfun(@minus,tgt_sigs,mu);
z = bsxfun(@minus,hsi_data,mu);

det = zeros(n_sigs,n_pix);

for i=1:n_sigs
    s = S(:,i);
    st_siginv = s'*siginv;
    st_siginv_s = s'*siginv*s;
    
    
    A = sum(st_siginv*z,1);
    B = sqrt(st_siginv_s);
    C = sqrt(sum(z.*(siginv*z),1));
    
    det(i,:) = A./(B.*C);
end

acert_data = max(det,[],1);

% ace_data = zeros(1,n_pix);
% for i=1:n_pix
%     ace_data(i) = (st_siginv*z(:,i))^2 / (st_siginv_s * z(:,i)'*siginv*z(:,i));
% end

end



