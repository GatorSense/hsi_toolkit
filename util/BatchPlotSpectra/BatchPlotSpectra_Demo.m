load('an_hsi_image_sub_for_demo.mat')

[x_dim, y_dim, z_dim] = size(hsi_img_sub);
data_reshaped = reshape(hsi_img_sub, x_dim*y_dim, z_dim);
mask_reshaped = reshape(mask_sub, x_dim*y_dim, 1);
masked_data = data_reshaped(mask_reshaped == 1, :);

SpecDist = BatchPlotSpectraDistribution(masked_data, wavelengths, [1,1,100], 1);