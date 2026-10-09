This directory contains FBP reconstructions based on frames of a tilt series collected from Chlamydomonas Reinhardtii. The tilt series is part of the Tomo110 dataset which was used in the tutorial for the CryoCARE denoising method. This tutorial is part of the cryoCARE_T2T GitHub repository (https://github.com/juglab/cryoCARE_T2T/tree/master). The Tomo110 dataset was downloaded from the following URL: https://download.fht.org/jug/cryoCARE/Tomo110.zip 

The files in this directory are as follows: 
	- tomo_all_frames.mrc: FBP reconstruction of a tilt series where each tilt image was obtained by summing all frames.
	- tomo_even_frames.mrc: FBP reconstruction of a tilt series where each tilt image was obtained by summing the frames with an even index.
	- tomo_odd_frames.mrc: FBP reconstruction of a tilt series where each tilt image was obtained by summing the frames with an odd index.
	- mask.mrc: A binary mask to exclude empty areas from sub-tomogram extraction. The mask was generated with IsoNet's (https://github.com/IsoNet-cryoET/IsoNet) make_mask command with the default parameters.
	- fitted_model.ckpt: A U-Net which we fitted for 1000 epochs on subtomograms extracted from tomo_even_frames.mrc and tomo_odd_frames.mrc 

All FBP reconstructions were performed as described in the CryoCARE tutorial (https://github.com/juglab/cryoCARE_T2T/tree/master/example). We cropped the FBP reconstructions to exclude empty areas.
