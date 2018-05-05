function [RGB_img] = view_RGB(hsi_img,wavelengths)
%
%[RGBimage] = view_RGB(img,wavelengths)
%
% creates an RGB image from a hyperspectral image
%
% inputs:
%  hsi_img - n_row x n_col x n_band hyperspectral image
%  wavelengths - 1 x n_band vector listing wavelength values for hsi_img
%
% outputs:
%   RGB_img - n_row x n_col x 3 RGB image
%
% 5/5/2018 - Alina Zare

RedWavelength = 650;
GreenWavelength = 540;
BlueWavelength = 475;

RGBimage = getHSIBands(I, wavelengths, [RedWavelength, GreenWavelength, BlueWavelength]);
RGBimage = double(RGBimage)/65536;

RGB_image = sqrt((RGBimage - min(RGBimage(:)))/(max(RGBimage(:))-min(RGBimage(:))));
