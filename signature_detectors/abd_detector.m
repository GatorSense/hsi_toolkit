function [abd_out] = abd_detector(hsi_img,tgt_sig,mask,ems)
%
%function [abd_out] = abd_detector(hsi_img,tgt_sig,mask,ems)
%
% Abundance of Target when unmixed with background endmembers
%
% inputs:
%  hsi_image - n_row x n_col x n_band hyperspectral image
%  tgt_sig - target signature (n_band x 1 - column vector)
%  mask - binary image limiting detector operation to pixels where mask is true
%         if not present or empty, no mask restrictions are used
%  ems - background endmembers
%
% outputs:
%  abd_out - detector image
%
% 8/19/2012 - Taylor C. Glenn 
% 6/2/2018 - Edited by Alina Zare
%

abd_out = img_det(@abd_det,hsi_img,tgt_sig,mask,ems);

end

function abd_data = abd_det(hsi_data,tgt_sig,ems)
addpath(fullfile('..','util'));
 
% unmix data with target signature and background
targ_P = unmix(hsi_data,double([tgt_sig ems])); 

abd_data = targ_P(:,1);

end
