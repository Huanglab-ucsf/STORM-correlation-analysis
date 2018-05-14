These source codes are distributed as accompanying codes for the article "A correlation analysis framework for localization-based super-resolution microscopy" by Bo Huang group.

The provided package includes software exploiting coordinated-based correlation analysis to three different application scenarios: model-free particle alignment, diffusion analysis and colocalization analysis. All three softwares take text molecule list file saved by Insight3 as an input. We provide a C script converting tab seperated coordinate list (column order: framenumber X Y width background intensity) to Insight3 txt format. Each software is distributed with example usage as below. If you have any troubles in using this package, please do not hesitate to contact us at the email address given below.
Contact: bo.huang@ucsf.edu

STORMAlignment:
	STORMAlignment molecule list (.txt) [parameters]
	Parameters
  		-p=xxx   Image pixel size (nm) (default = 107)
  		-noEM=x  1 for EMCCD camera mode, 0 for not (default = 1)
  		-wo=x    1 for localization density weighting, 0 for not (default = 1)
  		-ts=xx   Translation step size (nm) (default = 5)
  		-tr=xxx  Translation range (nm) (default = 25)
  		-a=xxx   Rotation angle step size (degree) (default = 5)
  		-n=xxx   Maximum interation number (default = 10)
	Output:
 		*_alignedxxxxx.txt: aligned molecule list in text format 

	Example: STORMAlignment DNA-PAINT-example.txt –ts=2 –tr=10 –a=2
	Note: the computational complexity of this methods is O(N*N), where N is the total number of localizations. For an initial test, we recommend using a relative small number of particles to avoid long computation time.

DiffusionCorr:
	DiffusionCorr <molecule list (.txt)> [parameters]
	Parameters
  		-i=nnn   Image size (pixels) (default = 128)
  		-p=xxx   Image pixel size (nm) (default = 160)
  		-g=xxx   Map grid size (nm) (default = 200)
  		-os=xx   Grid oversamples image pixels (overrides -g)
  		-us=xx   Grid undersamples image pixels (overrides -g)
  		-c=xxx   XY Correlation step (nm) (default = 30)
  		-r=xxx   XY Correlation range (nm) (default = 1000)
  		-X=xxx   X dimension of correlation region (nm) (default = automatic)
  		-Y=xxx   Y dimension of correlation region (nm) (default = automatic)
  		-m=xxx   Mininum points to calculate diffusion coefficient (default = 30)
  		-ed=xx   Estimated diffusion length (stdev) for fitting (default = 200)
  		-fs=nn   Start frame to calcualte correlation (default = 1)
  		-fi=nn   Frame interval to calculate correlation (default = 1)
  		-fr=nn   Frame range to calculate correlation (default = 1)
  		-absolute   Compute absolute density (#/um^2) instead of relative (default=of)
  		-noROI   Do not correct for finite area size (default = off)
	Output:
 	*_img.spe: Diffusion map
   		frame 1: density map
   		frame 2: diffusion map (MSD per frame)
 	*_cmap.spe: Frame pair corrlation map
 	*_dif.txt: Frame pair correlation of the entire list
   		column 1: r
   		column 2: frame pair correlation
	   
	Example:DiffusionCorr 128_128_ex_8ms_vss_0.5_variable_1.80_0.20_0003_list_crop.txt -c=20 -p=183 -fi=2 -fs=1

Colocalization:
	Python source code included for calculating point-point and point-set correlation. The example usage is illustrated in colocalization_exmaple ipython notebook script.
	Python packages dependence: numpy, scipy, numba.
