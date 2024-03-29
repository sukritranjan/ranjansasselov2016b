# -*- coding: iso-8859-1 -*-
"""
This script plots our reproduction of the surface radiance from \citet{Wuttke2006} and \cite{Kerr2008}
"""
########################
###Import useful libraries
########################

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pdb
import cookbook

def cm2inch(cm): #function to convert cm to inches; useful for complying with Astrobiology size guidelines
	return cm/2.54

########################
###Specific which family of plots to generate
########################
plot_wuttke=True # Plot reproduction of Wuttke et al (2006) plot for surface radiance (Antarctica)
plot_woudc=True # Plot reproduction of direct flux for WOUDC measurements of direct flux (Toronto station, 2003/6/21, 11:54:06, SZE=20.376, O3=354 DU)
#1 means plotting happens. 0 means plotting does not happen. 

############################
#######Wuttke et al (2006) diffuse radiance measurements in Antarctica
############################

if plot_wuttke:
	importeddata=np.genfromtxt('./TwoStreamOutput/reproduce_wuttke2006.dat', skip_header=1, skip_footer=0)
	wuttke_wav=importeddata[:,2] #centers of wavelength bins, nm
	wuttke_toa_modelled=importeddata[:,3] #top-of-atmosphere intensity, erg/s/nm/cm2
	wuttke_boa_modelled=importeddata[:,4] #bottom-of-atmosphere modelled diffuse intensity, erg/s/nm/cm2
	wuttke_boa_measured=importeddata[:,5] #bottom-of-atmosphere measured diffuse intensity, erg/s/nm/cm2
	
	wuttke_boa_modelled=cookbook.movingaverage(wuttke_boa_modelled, 10)
	
	fig, (ax1, ax3)=plt.subplots(2, figsize=(cm2inch(16.5),10), sharex=True)
	ax1.set_title('Reproduction of Wuttke+2006 Diffuse Radiance Measurements')
	ax1.plot(wuttke_wav, wuttke_toa_modelled, marker='s', color='black', label='TOA Flux')
	ax1.plot(wuttke_wav, wuttke_boa_modelled, marker='s', color='red', label='Calculated Diffuse Radiance')
	ax1.plot(wuttke_wav, wuttke_boa_measured, marker='s', color='blue', label='Measured Diffuse Radiance')
	ax1.set_yscale('log')
	ax1.set_ylim([1.e-3, 1.e4])
	ax1.set_xlim([290.,900.])
	ax1.set_xlabel('nm')
	ax1.set_ylabel('erg/s/cm2/nm')
	ax1.legend(loc=0)
	#ax2.plot(wuttke_wav, np.abs(wuttke_boa_measured-wuttke_boa_modelled)/(wuttke_toa_modelled), marker='s', color='black')
	#ax2.set_xlabel('nm')
	#ax2.set_ylabel('|Obs-Model|/TOA')
	#ax2.set_yscale('log')
	##ax2.set_ylim([-0.3,.3])
	ax3.plot(wuttke_wav, np.abs(wuttke_boa_measured-wuttke_boa_modelled)/(wuttke_boa_measured), marker='s', color='black')
	ax3.set_xlabel('nm')
	ax3.set_ylabel('|Obs-Model|/Obs')
	ax3.set_yscale('linear')
	ax3.set_ylim([0,1.3])
	plt.savefig('./Plots/reproduce_wuttke_paper_alt.eps', orientation='portrait',papertype='letter', format='eps')

if plot_woudc:
	solarzenithangle=20.376*np.pi/180. #
	importeddata=np.genfromtxt('./TwoStreamOutput/reproduce_woudc.dat', skip_header=5, skip_footer=6) #skip 0 entries as measurement error.
	kerr_wav=importeddata[:,2] #centers of wavelength bins, nm
	kerr_toa_modelled=importeddata[:,3]*np.cos(solarzenithangle) #top-of-atmosphere flux, erg/s/nm/cm2
	kerr_boa_modelled=importeddata[:,4] #bottom-of-atmosphere modelled flux, erg/s/nm/cm2
	kerr_boa_measured=importeddata[:,5] #bottom-of-atmosphere measured flux, erg/s/nm/cm2

	fig1, (ax1, ax3)=plt.subplots(2, figsize=(cm2inch(16.5),10), sharex=True)
	ax1.set_title('Reproduction of WOUDC Surface Flux Measurements')
	ax1.plot(kerr_wav, kerr_toa_modelled, marker='s', color='black', label='TOA Flux')
	ax1.plot(kerr_wav, kerr_boa_modelled, marker='s', color='red', label='Calculated Surface Flux')
	ax1.plot(kerr_wav, kerr_boa_measured, marker='s', color='blue', label='Measured Surface Flux')
	ax1.set_yscale('log')
	ax1.set_ylim([1.e-3, 1.e4])
	ax1.set_xlim([292.,360.])
	ax1.set_xlabel('nm')
	ax1.set_ylabel('erg/s/cm2/nm')
	ax1.legend(loc=0)
	#ax2.plot(kerr_wav, np.abs(kerr_boa_measured-kerr_boa_modelled)/(kerr_toa_modelled), marker='s', color='black')
	#ax2.set_xlabel('nm')
	#ax2.set_ylabel('|Obs-Model|/TOA')
	#ax2.set_yscale('log')
	##ax2.set_ylim([1.e-4,2])
	ax3.plot(kerr_wav, np.abs(kerr_boa_measured-kerr_boa_modelled)/(kerr_boa_measured), marker='s', color='black')
	ax3.set_xlabel('nm')
	ax3.set_ylabel('|Obs-Model|/Obs')
	ax3.set_yscale('linear')
	ax3.set_ylim([0.,1.3])
	plt.savefig('./Plots/reproduce_woudc_paper_alt.eps', orientation='portrait',papertype='letter', format='eps')
plt.show()