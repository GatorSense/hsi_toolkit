function [ace_out] = acert_ss_detector(hsi_img,tgt_sigs,mask,mu,siginv)
%
%function [ace_out] = acert_ss_detector(hsi_img,tgt_sigs,mask)
%
% Adaptive Cosine/Coherence Estimator - Subspace Formulation
%
% inputs:
%  hsi_image - n_row x n_col x n_band hyperspectral image
%  tgt_sigs - target signatures (n_band x M - column vector)
%  mask - binary image limiting detector operation to pixels where mask is true
%         if not present or empty, no mask restrictions are used
%
% outputs:
%  ace_out - detector image
%
% 1/17/2015 - Taylor C. Glenn - tcg@cise.ufl.edu
%

if ~exist('mask','var'); mask = []; end
if ~exist('mu','var'), mu = []; end
if ~exist('siginv','var'), siginv = []; end

acert_out = img_det(@acert_ss_det,hsi_img,tgt_sigs,mask,mu,siginv);

end

function acert_data = acert_ss_det(hsi_data,tgt_sigs,mu,siginv)
    
if isempty(mu)
    mu = mean(hsi_data,2);
end
if isempty(siginv)
    siginv = pinv(cov(hsi_data'));
end

S = bsxfun(@minus,tgt_sigs,mu);
z = bsxfun(@minus,hsi_data,mu);

G = siginv*S*pinv(S'*siginv*S)*S'*siginv;

A = sum(z.*(G*z),1);
B = sum(z.*(siginv*z),1);

acert_data = A./B;

% ace_data = zeros(1,n_pix);
% for i=1:n_pix
%     ace_data(i) = (z(:,i)'*G*z(:,i)) / (z(:,i)'*siginv*z(:,i));
% end

end



