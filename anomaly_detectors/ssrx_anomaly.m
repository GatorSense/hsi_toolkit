function [ssrx_img] = ssrx_anomaly(hsi_img,n_dim_ss,guard_win,bg_win)
%
%function ssrx_img = ssrx_anomaly(hsi_img,n_dim_ss,guard_win,bg_win)
%
% Subspace Reed-Xiaoli anomaly detector
%  eliminate leading subspace as background, then
%  use local mean and covariance to determine pixel to background distance
%
% inputs:
%   hsi_image - n_row x n_col x n_band
%   n_dim_ss - number of leading dimensions to use in the background subspace
%   guard_win - guard window radius (square,symmetric about pixel of interest) 
%   bg_win - background window radius
%  
% 8/7/2012 - Taylor C. Glenn
% 5/5/2018 - Edited by Alina Zare

[n_row,n_col,n_band] = size(hsi_img);
n_pix = n_row*n_col;

hsi_data = reshape(hsi_img,[n_pix,n_band])';

% pca the data with no dim reduction
[pca_data,~,vecs,vals] = pca(hsi_data,1);

pca_img = reshape(pca_data',[n_row,n_col,n_band]);

proj = eye(n_band,n_band) - vecs(:,1:n_dim_ss)*vecs(:,1:n_dim_ss)';

% create the mask
mask_width = 1 + 2*guard_win + 2*bg_win;
half_width = guard_win + bg_win;
mask_rg = (1:mask_width)-1;

b_mask = true(mask_width,mask_width);
b_mask(bg_win+1:end-bg_win,bg_win+1:end-bg_win) = false;


% run the detector
%  (only on fully valid points)
ssrx_img = zeros(n_row,n_col);

for i=1:n_col-mask_width+1
    for j=1:n_row-mask_width+1

        b_mask_img = false(n_row, n_col);
        b_mask_img(mask_rg+j, mask_rg+i) = b_mask;
        b_mask_list = b_mask_img(:);
        
        %pull out background points        
        bg = pca_data(:,b_mask_list);
        
        %compute Mahalanobis distance
        siginv = pinv(cov(bg'));
        mu = mean(bg,2);
        z = proj*squeeze(pca_img(j+half_width,i+half_width,:)) - proj*mu;
        
        ssrx_img(j+half_width,i+half_width) = z'*siginv*z;  
    end
end
