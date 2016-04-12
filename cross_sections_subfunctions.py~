# -*- coding: iso-8859-1 -*-
"""
Functions to compute the mean cross-section in each bin.
"""
import numpy as np
import pdb
import matplotlib.pyplot as plt
import scipy.stats
from scipy import interpolate as interp
import cookbook


"""
***************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************
"""

def compute_band_cross_section(leftedges, rightedges, molecule):
	"""
	Objective of this code is to compute the per-molecule cross-section of a given molecule.
	Inputs:
	-left edges of the wavelength bins (nm)
	-right edges of the wavelength bins (nm)
	-molecule

	Output:
	per-molecule absorption cross-section in cm^2
	per-molecule scattering cross-section in cm^2

	Note that (at present) no pressure or temperature broadening effects are included
	"""
	import numpy as np
	import scipy.integrate
	from scipy import interpolate as interp

	n_bins=len(leftedges)
	#import data
	data=np.genfromtxt('./XCs/composite_xc_extended_'+molecule, skip_header=1, skip_footer=0)
	wav=data[:,0] #wavelengths in nm
	tot_xc=data[:,1] #total xc in cm2, rayleigh+abs
	abs_xc=data[:,2] #absorption xc in cm2
	ray_xc=data[:,3] #rayleigh scattering xc in cm2

	#form functions of cross-sections
	tot_xc_func=interp.interp1d(wav, tot_xc, kind='linear')
	abs_xc_func=interp.interp1d(wav, abs_xc, kind='linear')
	ray_xc_func=interp.interp1d(wav, ray_xc, kind='linear')
	
	#initialize variables to hold the bandpass-integrated cross-sections\
	#Reminder to self: our method assumes units of cross-sections are cm2 only (no per nm). Verify argument for this with Dimitar.
	tot_xc_band=np.zeros(n_bins)
	abs_xc_band=np.zeros(n_bins)
	ray_xc_band=np.zeros(n_bins)
	
	for ind in range(0,n_bins):
		#find average cross-sections by integrating across band and dividing by size of bandpass...
		tot_xc_band[ind]=scipy.integrate.quad(tot_xc_func, leftedges[ind], rightedges[ind])[0]/(rightedges[ind]-leftedges[ind])
		abs_xc_band[ind]=scipy.integrate.quad(abs_xc_func, leftedges[ind], rightedges[ind])[0]/(rightedges[ind]-leftedges[ind])
		ray_xc_band[ind]=scipy.integrate.quad(ray_xc_func, leftedges[ind], rightedges[ind])[0]/(rightedges[ind]-leftedges[ind])
	
	return (tot_xc_band, abs_xc_band, ray_xc_band)


def get_rugheimer_xc(leftedges, rightedges, molecule, mr_n2, mr_co2):
	"""
	Objective of this code is to load the Rugheimer molecular cross-sections. They are taken from photos.pdat. Note that this assumes preset wavelength bins.
	-left and right edges of wavelength bins in nm.
	-molecule
	-mixing ratio of n2
	-mixing ratio of co2

	Output:
	per-molecule absorption cross-section in cm^2
	per-molecule scattering cross-section in cm^2

	"""
	import numpy as np
	import scipy.integrate
	from scipy import interpolate as interp

	n_bins=len(leftedges)
	
	abs_xc_band=np.zeros(n_bins)
	ray_xc_band=np.zeros(n_bins)
	
	faruvdata=np.genfromtxt('./Raw_Data/Rugheimer_Metadata/faruvs_mod.pdat', skip_header=3, skip_footer=1)

	nearuv0=np.genfromtxt('./Raw_Data/Rugheimer_Metadata/photos.pdat',skip_header=2, skip_footer=400)
	nearuv1=np.genfromtxt('./Raw_Data/Rugheimer_Metadata/photos.pdat',skip_header=159, skip_footer=319)
	nearuv2=np.genfromtxt('./Raw_Data/Rugheimer_Metadata/photos.pdat',skip_header=397, skip_footer=70)

	if molecule=='o2':
		abs_xc_band[0:9]=faruvdata[::-1, 1]
		abs_xc_band[9:44]=nearuv1[:,1]
	if molecule=='co2':
		abs_xc_band[0:9]=faruvdata[::-1, 2]		
		abs_xc_band[9:44]=nearuv1[:,3]
	if molecule=='h2o':
		abs_xc_band[0:9]=faruvdata[::-1, 3]	
		abs_xc_band[9:44]=nearuv1[:,2]
	if molecule=='so2':
		abs_xc_band[0:9]=faruvdata[::-1, 4]	
		abs_xc_band[9:(9+68)]=nearuv2[:,2]+nearuv2[:,3]+nearuv2[:,4]	
	if molecule=='o3':
		abs_xc_band[9:(9+108)]=nearuv0[:,3]#+nearuv0[:,4]
	if molecule=='h2s':
		abs_xc_band[9:(9+68)]=nearuv2[:,5]
	if molecule=='n2':
		#no absorption from n2, only Rayleigh scattering. Take all Rayleigh to come from N2
		#Compute Rayleigh scattering according to method of SIGRAY in ltning.f and the modification in photo.f
		wavcen=0.5*(leftedges+rightedges)*1.e-3 #convert to microns
		ray_xc_band=(4.006e-28*(1.+0.0113/wavcen**2.+0.00013/wavcen**4.)/wavcen**4.)*(1.+1.5*mr_co2)/mr_n2 #scale by the mixing ratio of N2 to account for correction.
	
	tot_xc_band=abs_xc_band+ray_xc_band
	return (tot_xc_band, abs_xc_band, ray_xc_band)