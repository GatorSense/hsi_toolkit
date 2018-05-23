function [RGB_img] = get_RGB(hsi_img,wavelengths)
%
%[RGBimage] = get_RGB(img,wavelengths)
%
% creates an RGB image from a hyperspectral image
%
% inputs:
%  hsi_img - n_row x n_col x n_band hyperspectral image
%  wavelengths - 1 x n_band vector listing wavelength values for hsi_img in
%  nm
%
% outputs:
%   RGB_img - n_row x n_col x 3 RGB image
%
% 5/5/2018 - Alina Zare

RedWavelengths = [620:659];
GreenWavelengths = [550:570];
BlueWavelengths = [450:495];

RGB_img(:,:,1) = mean(get_hsi_bands(hsi_img, wavelengths, RedWavelengths),3);
RGB_img(:,:,2) = mean(get_hsi_bands(hsi_img, wavelengths, GreenWavelengths),3);
RGB_img(:,:,3) = mean(get_hsi_bands(hsi_img, wavelengths, BlueWavelengths),3);

RGB_img = ((RGB_img - min(RGB_img(:)))/(max(RGB_img(:))-min(RGB_img(:)))).^(1/1.5);
