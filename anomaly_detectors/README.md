# gatorsense hsitoolkit anomaly detectors matlab version
##hsitoolkit/anomaly_detectors

Suite of anomaly detectors implemented so far: 
- beta_anomaly: global anomaly detector, fits beta distribution to each band assuming entire image is background computes negative log likelihood of each pixel in the model
- cbad_anomaly: global/cluster-based anomaly detector, Cluster Based Anomaly Detection (CBAD)
- csd_anomaly: global anomaly detector, Complementary Subspace Detector
- fcbad_anomaly: global/cluster-based anomaly detector, Fuzzy Cluster Based Anomaly Detection (FCBAD)
- gmm_anomaly: global anomaly detector, fits GMM assuming entire image is background computes negative log likelihood of each pixel in the fit model
- gmrx_anomaly: global/cluster-based anomaly detector, fits GMM assuming entire image is background assigns pixels to highest posterior probability mixture component computes pixel Mahlanobis distance to component mean
- md_anomaly: global anomaly detector, Mahalanobis Distance anomaly detector uses global image mean and covariance as background estimates
- rx_anomaly: local anomaly detector, Widowed Reed-Xiaoli anomaly detector uses local mean and covariance to determine pixel to background distance
- ssrx_anomaly: local anomaly detector, eliminate leading subspace as background, then use local mean and covariance to determine pixel to background distance

Demo: 
(1) load an_hsi_image_sub_for_demo.mat
(2) run anomaly_det_demo(hsi_img_sub, mask_sub, wavelengths)

Contact: Alina Zare, azare@ufl.edu
