# -*- coding: iso-8859-1 -*-
"""
This code demonstrates our reproduction of the Rugheimer et al results.
"""
########################
###Import useful libraries
########################

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pdb

def cm2inch(cm): #function to convert cm to inches; useful for complying with Astrobiology size guidelines
	return cm/2.54
########################
###Load in Rugheimer model.
########################
importeddata=np.genfromtxt('./LiteratureSpectra/rugheimer_earth_modern.dat', skip_header=1, skip_footer=0)
rugheimer_wav_leftedges=importeddata[:,0] #left edges of wavelength bins, nm
rugheimer_wav_rightedges=importeddata[:,1] #right edges of wavelength bins, nm
rugheimer_wav_centers=importeddata[:,2] #centers of wavelength bins, nm
rugheimer_int_toa=importeddata[:,3] #top-of-atmosphere intensity, erg/s/cm2/nm
rugheimer_int_boa=importeddata[:,4] #BOA Actinic Flux/Total Intensity, erg/s/cm2/nm. This is what we need to match.


########################
###Load in our code's reproduction (A=0.20, our cross-sections)
########################
importeddata=np.genfromtxt('./TwoStreamOutput/rugheimer_earth_modern.dat', skip_header=1, skip_footer=0)
twostr1_wav_leftedges=importeddata[:,0] #left edges of wavelength bins, nm
twostr1_wav_rightedges=importeddata[:,1] #right edges of wavelength bins, nm
twostr1_wav_centers=importeddata[:,2] #centers of wavelength bins, nm
twostr1_int_toa=importeddata[:,3] #top-of-atmosphere intensity, erg/s/cm2/nm
twostr1_flux_gnd=importeddata[:,4] #surface flux, erg/s/cm2/nm. 
twostr1_int_BOA=importeddata[:,5] #total intensity in middle of bottom layer of atmosphere
twostr1_int_gnd=importeddata[:,6] #intensity at planetary surface, =0.5*diffuse total intensity+solar intensity

########################
###Load in our code's reproduction (A=0.20, Rugheimer cross-sections)
########################
importeddata=np.genfromtxt('./TwoStreamOutput/rugheimer_earth_modern_rugheimerxcs.dat', skip_header=1, skip_footer=0)
twostr3_wav_leftedges=importeddata[:,0] #left edges of wavelength bins, nm
twostr3_wav_rightedges=importeddata[:,1] #right edges of wavelength bins, nm
twostr3_wav_centers=importeddata[:,2] #centers of wavelength bins, nm
twostr3_int_toa=importeddata[:,3] #top-of-atmosphere intensity, erg/s/cm2/nm
twostr3_flux_gnd=importeddata[:,4] #surface flux, erg/s/cm2/nm. 
twostr3_int_BOA=importeddata[:,5] #total intensity in middle of bottom layer of atmosphere
twostr3_int_gnd=importeddata[:,6] #intensity at planetary surface, =0.5*diffuse total intensity+solar intensity


########################
###Plot results
########################
fig, (ax0, ax1)=plt.subplots(2, figsize=(cm2inch(16.5),10.), sharex=True)
markersizeval=5.

ax0.plot(rugheimer_wav_centers, rugheimer_int_toa,  marker='s', markersize=markersizeval, linewidth=1, color='black', label='TOA Int.')
ax0.plot(rugheimer_wav_centers, rugheimer_int_boa,  marker='s', markersize=markersizeval, linewidth=1, color='red', label='BOA Int. \n(Rugheimer+2015)')
ax0.plot(twostr1_wav_centers, twostr1_int_BOA,  marker='s', markersize=markersizeval, linewidth=1, color='blue', label='BOA Int. (Our Model)')
ax0.plot(twostr3_wav_centers, twostr3_int_BOA,  marker='s', markersize=markersizeval, linewidth=1, color='orange', label='BOA Int. (Our Model, \nRugheimer+2013 XCs)')
ax0.set_yscale('log')
ax0.set_ylim([1.e-2, 1e4])
ax0.set_ylabel('Total Intensity \n (erg s$^{-1}$cm$^{-2}$nm$^{-1}$)')
ax0.legend(loc='lower center', fontsize=12)
ax0.set_xlim([130, 860])


#if both our codes report zero, ratio should be 1.
rugheimer_int_boa_plot_3=np.copy(rugheimer_int_boa)
rugheimer_int_boa_plot_1=np.copy(rugheimer_int_boa)
twostr1_int_BOA_plot=np.copy(twostr1_int_BOA)
twostr3_int_BOA_plot=np.copy(twostr3_int_BOA)

foo=np.where((rugheimer_int_boa==0.) & (twostr3_int_BOA==0.))
twostr3_int_BOA_plot[foo]=1.
rugheimer_int_boa_plot_3[foo]=1.

bar=np.where((rugheimer_int_boa==0.) & (twostr1_int_BOA==0.))
twostr1_int_BOA_plot[bar]=1.
rugheimer_int_boa_plot_1[bar]=1.


ax1.plot(twostr1_wav_centers, (twostr1_int_BOA-rugheimer_int_boa)/rugheimer_int_toa, marker='s', linewidth=1, color='blue', label='Our XCs')
ax1.plot(twostr3_wav_centers, (twostr3_int_BOA-rugheimer_int_boa)/rugheimer_int_toa, marker='s', linewidth=1, color='orange', label='Rugheimer+2013 XCs')
ax1.set_yscale('linear')
ax1.set_ylim([-1.1, 0.1])
ax1.set_xlabel('Wavelength (nm)')
ax1.set_ylabel('(Our Model-R+2013)/TOA Int.')
ax1.legend(loc='lower center', fontsize=12)
ax1.set_xlim([130, 860])

ax0.set_title('Reproduction of Rugheimer+2013 Model')

bunk=(twostr1_int_BOA-rugheimer_int_boa)/rugheimer_int_toa
thunk=(twostr3_int_BOA-rugheimer_int_boa)/rugheimer_int_toa
print np.max(np.abs(bunk))
print np.max(np.abs(thunk))
plt.tight_layout()
plt.savefig('./Plots/reproduce_rugheimer_modern.eps', orientation='portrait',papertype='letter', format='eps')
pdb.set_trace()
plt.show()