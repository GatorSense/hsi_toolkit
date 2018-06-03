function [sam_out] = sam_detector(hsi_img,tgt_sig,mask)
%
%function [sam_out] = sam_detector(hsi_img,tgt_sig,mask)
%
% Spectral Angle Mapper
%
% inputs:
%  hsi_image - n_row x n_col x n_band hyperspectral image
%  tgt_sig - target signature (n_band x 1 - column vector)
%  mask - binary image limiting detector operation to pixels where mask is true
%         if not present or empty, no mask restrictions are used
%
% outputs:
%  sam_out - detector image
%
% 8/8/2012 - Taylor C. Glenn
% 6/2/2018 - Edited by Alina Zare
%

if ~exist('mask','var'); mask = []; end
addpath(fullfile('..','util'));

sam_out = img_det(@sam_det,hsi_img,tgt_sig,mask);

end

function sam_data = sam_det(hsi_data,tgt_sig)
   
prod = tgt_sig'*hsi_data;

mag = sqrt((tgt_sig'*tgt_sig) * sum(hsi_data.^2,1));

sam_data = prod./mag;

end


