# MUTE: Dehazing for cataractous retinal image using denoising

This is the MATLAB code for retinal image enhancement using Double-pass fundus reflection model. We are happy that this research has been accepted and published on Medical Image Analysis: [https://doi.org/10.1016/j.sigpro.2021.108400](https://doi.org/10.1016/j.media.2023.102848)

USAGE:
Run main_MUTE_finished.m and select a retinal image for enhancement.

The gray value threshold in red channel for background padding should be adjusted if the retinal image has severe illumination problem, for example, it should be set to 2 for DiaRetdb0_image034.png. Normally, the value is set to 20.

Other parameters are automatically and adaptively determined.
