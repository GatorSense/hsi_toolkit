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

RedWavelength = 650;
GreenWavelength = 540;
BlueWavelength = 475;

RGB_img = get_hsi_bands(hsi_img, wavelengths, [RedWavelength, GreenWavelength, BlueWavelength]);
RGB_img = sqrt((RGB_img - min(RGB_img(:)))/(max(RGB_img(:))-min(RGB_img(:))));
