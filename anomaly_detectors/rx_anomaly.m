function [rx_img] = rx_anomaly(hsi_img,mask,guard_win,bg_win)
%
%function rx_img = rx_anomaly(hsi_img,mask,guard_win,bg_win)
%
% Widowed Reed-Xiaoli anomaly detector
%  use local mean and covariance to determine pixel to background distance
%
% inputs:
%   hsi_image - n_row x n_col x n_band
%   mask - binary image limiting detector operation to pixels where mask is true
%         if not present or empty, no mask restrictions are used
%   guard_win - guard window radius (square,symmetric about pixel of interest) 
%   bg_win - background window radius
%  
% 8/7/2012 - Taylor C. Glenn - tcg@cise.ufl.edu
% 5/5/2018 - Edited by Alina Zare

[n_row,n_col,n_band] = size(hsi_img);
n_pix = n_row*n_col;

if ~exist('mask','var') || isempty(mask); mask = true(n_row,n_col); end    

% create the mask
mask_width = 1 + 2*guard_win + 2*bg_win;
half_width = guard_win + bg_win;
mask_rg = (1:mask_width)-1;

b_mask = true(mask_width,mask_width);
b_mask(bg_win+1:end-bg_win,bg_win+1:end-bg_win) = false;

hsi_data = reshape(hsi_img,[n_pix,n_band])';


% run the detector
%  (only on fully valid points)
rx_img = zeros(n_row,n_col);

for i=1:n_col-mask_width+1
    for j=1:n_row-mask_width+1
        
        row = j+half_width;
        col = i+half_width;
        
        if ~mask(row,col), continue; end 

        b_mask_img = false(n_row, n_col);
        b_mask_img(mask_rg+j, mask_rg+i) = b_mask;
        b_mask_list = b_mask_img(:);
        
        %pull out background points        
        bg = hsi_data(:,b_mask_list);
        
        %compute Mahalanobis distance
        siginv = pinv(cov(bg'));
        mu = mean(bg,2);
        z = squeeze(hsi_img(row,col,:)) - mu;
        
        rx_img(row,col) = z'*siginv*z;
                
    end
end
