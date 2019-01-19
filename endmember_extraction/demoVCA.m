%{
Demo script that runs the VCA algorithm using example sub MUUFL Gulfport data
Inputs:
	hsi_img_sub - n_row x n_col x n_band hyperspectral image
	wavelengths - n_band x 1 vector listing wavelength values for hsi_img in nm
	mask_sub - n_row x n_col binary image limiting detector operation to pixels where mask is true
		       if not present or empty, no mask restrictions are used
	M - number of endmembers to compute
Outputs:
	E - n_band x M matrix of endmembers
	IdxOfE - M vector indexing the endmembers in masked_data
	Xpca - M x masked_data size matrix of data projected to M dimensions
1/17/2019 - Ronald Fick
%}

load('an_hsi_image_sub_for_demo.mat');

[x_dims, y_dims, band_dims] = size(hsi_img_sub);

mat_data = reshape(hsi_img_sub, x_dims*y_dims, band_dims);

mask_reshaped = reshape(mask_sub, x_dims*y_dims, 1);

masked_data = mat_data(mask_reshaped == 1, :);

masked_data = double(masked_data);

M = 3;

[E, IdxOfE, Xpca] = VCAbyDG(masked_data', 'Endmembers', M);

figure();
nonendmembers = setdiff(1:size(masked_data,1), IdxOfE);
scatter3(Xpca(1,nonendmembers), Xpca(2,nonendmembers), Xpca(3,nonendmembers), 5, 'b');
hold on;
scatter3(Xpca(1,IdxOfE), Xpca(2,IdxOfE), Xpca(3,IdxOfE), 40, 'r');
title('Gulfport Data Projected to 3D - Endmembers in Red');

figure();
plot(wavelengths, E);
title('Estimated Endmembers from Gulfport Data');
xlabel('Wavelength (nm)');
ylabel('Reflectance');