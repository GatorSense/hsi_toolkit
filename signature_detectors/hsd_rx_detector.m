function [out] = hsd_rx_detector(hsi_img,tgt_sig,mask,guard_win,bg_win,ems,beta)
%
%function [out] = hsd_rx_detector(hsi_img,tgt_sig,mask,guard_win,bg_win,ems,beta)
%
% Hybrid Subpixel Detector with RX style local background estimation
%
% inputs:
%  hsi_image - n_row x n_col x n_band hyperspectral image
%  tgt_sig - target signature (n_band x 1 - column vector)
%  mask - binary image limiting detector operation to pixels where mask is true
%         if not present or empty, no mask restrictions are used
%  guard_win - guard window radius (square,symmetric about pixel of interest) 
%  bg_win - background window radius
%
% outputs:
%  out - detector image
%
% 1/25/2013 - Taylor C. Glenn 
% 6/3/2018 - Edited by Alina Zare
%

[n_row,n_col,n_band] = size(hsi_img);

if ~exist('mask','var') || isempty(mask), mask = true(n_row,n_col); end
if ~exist('bg_win','var'), bg_win = 4; end
if ~exist('guard_win','var'), guard_win = 2; end
if ~exist('beta','var'), beta = 0; end

hsi_data = reshape(hsi_img,n_row*n_col,n_band)';

reg = beta*eye(n_band);

% unmix data with only background endmembers
P = unmix(hsi_data,ems);

% unmix data with target signature as well
targ_P = unmix(hsi_data,[tgt_sig ems]); 

[out] = rx_det(@hsd_rx_pt,hsi_img,tgt_sig,mask,guard_win,bg_win,reg,ems,P,targ_P);

end

function [r] = hsd_rx_pt(x,ind,bg,~,args,reg,ems,P,targ_P)

if ~isempty(bg)
    sigma = cov(bg');
    siginv = pinv(sigma + reg);
else
    siginv = global_siginv;
end

z = x - ems*P(ind,:)';
w = x - [args.tgt_sigs ems]*targ_P(ind,:)';
r = (z'*siginv*z) / (w'*siginv*w);

end



