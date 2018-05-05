function hsi_out = get_hsi_bands(hsi_img, wavelengths, wavelengths_to_get)
%
%hsi_out = getHSIBands(hsi_img, wavelengths, wavelengthsToGet)
%
% Return hsi stack with bands that are closest to desired wavelengths
%
% inputs:
%  hsi_img - n_row x n_col x n_band hyperspectral image
%  wavelengths - 1 x n_band vector listing wavelength values for hsi_img
%  wavengths_to_get - 1 x n_band_new vector listing desired wavelengths
%
% outputs:
%   RGB_img - n_row x n_col x n_band_new  image
%
% 5/5/2018 - Alina Zare


array_of_waves = [];
for i = 1:numel(wavelengths_to_get)
   [C,Index] = min( abs(wavelengths - wavelengths_to_get(i)) );
   array_of_waves = [array_of_waves, Index];
end

hsi_out = hsi_img(:,:,array_of_waves);
end