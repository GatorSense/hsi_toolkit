function csd_out = csd_anomaly(hsi_img,n_dim_bg,n_dim_tgt,tgt_orth)
%
% function csd_out = csd_anomaly(hsi_img,n_dim_bg,n_dim_tgt,tgt_orth)
%
% Complementary Subspace Detector
%  assumes background and target are complementary subspaces
%  of PCA variance ranked space
% Ref: A. Schaum, "Joint subspace detection of hyperspectral targets," 2004 IEEE Aerospace Conference Proceedings (IEEE Cat. No.04TH8720), 2004, pp. 1824 Vol.3. doi: 10.1109/AERO.2004.1367963
%
% inputs:
%   hsi_image - n_row x n_col x n_band
%   n_dim_bg - number of leading dimensions to assign to background subspace
%   n_dim_tgt - number of dimensions to assign to target subspace
%               use empty matrix, [], to use all remaining after background assignment
%   tgt_orth - True/False, set target subspace orthogonal to background subspace
%
% 8/7/2012 - Taylor C. Glenn
% 5/5/2018 - Edited by Alina Zare

addpath(fullfile('..','util'));

[n_row,n_col,n_band] = size(hsi_img);
n_pix = n_row*n_col;

hsi_data = reshape(hsi_img,[n_pix,n_band])';

% get pca rotation, no dim reduction
[pca_data,~,vecs,vals] = pca(hsi_data,1);

% whiten the data
%   (so that later steps are equivalent to Mahalanobis distance)
z = diag(1./sqrt(vals)) * pca_data;

% figure out background and target subspaces
bg_rg = 1:n_dim_bg;

if tgt_orth
    % set target to orthogonal complement of background
    if isempty(n_dim_tgt)
        n_dim_tgt = n_band - n_dim_bg;
    end
    tgt_rg = (n_dim_bg+1):n_dim_tgt;
    
else
    % target and background overlap
    if isempty(n_dim_tgt)
        n_dim_tgt = n_band;
    end
    tgt_rg = 1:n_dim_tgt;    
end

% set background and target subspaces
B = vecs(:,bg_rg);
S = vecs(:,tgt_rg);

% run the detector
csd_data = zeros(1,n_pix);

for i=1:n_pix
    Sz = S'*z(:,i);
    Bz = B'*z(:,i);
    
    csd_data(i) = Sz'*Sz - Bz'*Bz;        
end

csd_out = reshape(csd_data,[n_row,n_col]);


