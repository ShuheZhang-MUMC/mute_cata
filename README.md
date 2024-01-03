# MUTE: Dehazing for cataractous retinal image using denoising

This is the MATLAB code for retinal image enhancement using MUTE method. We are happy that this research has been accepted and published on Medical Image Analysis: [[PAPER]](https://doi.org/10.1016/j.media.2023.102848)

USAGE:
Run main_MUTE_finished.m and select a retinal image for enhancement.

The gray value threshold in red channel for background padding should be adjusted if the retinal image has severe illumination problem. Normally, the value is set to 20. Please ensure that the background and circular field of interest is separated as the grayvalue for the background should all be zero.

Other parameters are automatically and adaptively determined.
