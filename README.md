# SAR-RCS

SAR imaging and RCS recovery of finite-sized targets in free space

This collection of Matlab codes were used to compute the simulation results contained in the manuscript "Detecting, locating and spectrally characterizing targets with SAR" by A.D. Kim, J. Simpson and C. Tsogka.

The code "point_target_asymptotic_image.m" simulates SAR data for a point target and computes the corresponding KM image. It also produces the asymptotic KM image and computes the absolute difference between the true and asymptotic KM images. These plots correspond to Fig. 2 of the manuscript.

The code "asymptotic_KM_images.m" uses "SAR_data_MFS_3D.m" and "sphere3D.m" to simulate SAR data for dielectric spheres and computes the corresponding KM images. It also produces the asymptotic KM images and computes the absolute difference between the true and asymptotic KM images. These plots correspond to Figs. 3 and 4 of the manuscript.

The code "single_ellipsoid_image.m" uses "SAR_data_MFS_3D.m" and "ellipsoid3D_equal.m" to simulate SAR data for a dielectric ellipsoid and computes the corresponding KM image. It also uses "backscatamp3D_mfs.m" to compute the scattering amplitude for the ellipsoid in the direction of backscattering. These plots correspond to Fig. 5 of the manuscript.

The code "two_sphere_images_linearity.m" uses "SAR_data_MFS_3D.m" and "twospheres.m" to simulate SAR for dielectric spheres and computes the corresponding KM image. It also produces the modified KM image using the approximate linearity of the KM imaging function. These plots correspond to Fig. 6 of the manuscript.

The code "single_RCS_recovery.m" uses "SAR_data_MFS_3D.m", "sphere3D.m" and "ellipsoid3D_equal.m" to simulate SAR data for dielectric spheres and ellipsoids and computes the corresponding KM images. It also recovers the RCS spectrum for each target. These plots correspond to Fig. 7 of the manuscript.
