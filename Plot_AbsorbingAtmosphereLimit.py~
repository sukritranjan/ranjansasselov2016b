# -*- coding: iso-8859-1 -*-
"""
This code plots the absorption limit test of our code. 
"""
########################
###Import useful libraries
########################

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pdb

deg2rad=np.pi/180.

def cm2inch(cm): #function to convert cm to inches; useful for complying with Astrobiology size guidelines
	return cm/2.54

########################
###First, the A=0 case
########################
#Zenith angle of 0.
importeddata=np.genfromtxt('./TwoStreamOutput/rugheimer_earth_epoch0_w0=1e-5_a=0_z=0.dat', skip_header=1, skip_footer=0)
wav_leftedges_0_0=importeddata[:,0] #left edges of wavelength bins, nm
wav_rightedges_0_0=importeddata[:,1] #right edges of wavelength bins, nm
wav_centers_0_0=importeddata[:,2] #centers of wavelength bins, nm
flux_toa_0_0=importeddata[:,3]*np.cos(0.*deg2rad) #top-of-atmosphere flux, erg/s/cm2/nm, converted from 
flux_gnd_0_0=importeddata[:,4] #surface flux, erg/s/cm2/nm. 
flux_gnd_dir_0_0=importeddata[:,5] #Direct surface flux, erg/s/cm2/nm. 
flux_gnd_dif_0_0=importeddata[:,6] #Diffuse surface flux, erg/s/cm2/nm. 

#Zenith angle of 60
importeddata=np.genfromtxt('./TwoStreamOutput/rugheimer_earth_epoch0_w0=1e-5_a=0_z=60.dat', skip_header=1, skip_footer=0)
wav_leftedges_0_60=importeddata[:,0] #left edges of wavelength bins, nm
wav_rightedges_0_60=importeddata[:,1] #right edges of wavelength bins, nm
wav_centers_0_60=importeddata[:,2] #centers of wavelength bins, nm
flux_toa_0_60=importeddata[:,3]*np.cos(60.*deg2rad) #top-of-atmosphere flux, erg/s/cm2/nm, converted from 
flux_gnd_0_60=importeddata[:,4] #surface flux, erg/s/cm2/nm. 
flux_gnd_dir_0_60=importeddata[:,5] #Direct surface flux, erg/s/cm2/nm. 
flux_gnd_dif_0_60=importeddata[:,6] #Diffuse surface flux, erg/s/cm2/nm. 

#Zenith angle of 85
importeddata=np.genfromtxt('./TwoStreamOutput/rugheimer_earth_epoch0_w0=1e-5_a=0_z=85.dat', skip_header=1, skip_footer=0)
wav_leftedges_0_85=importeddata[:,0] #left edges of wavelength bins, nm
wav_rightedges_0_85=importeddata[:,1] #right edges of wavelength bins, nm
wav_centers_0_85=importeddata[:,2] #centers of wavelength bins, nm
flux_toa_0_85=importeddata[:,3]*np.cos(85.*deg2rad) #top-of-atmosphere flux, erg/s/cm2/nm, converted from 
flux_gnd_0_85=importeddata[:,4] #surface flux, erg/s/cm2/nm. 
flux_gnd_dir_0_85=importeddata[:,5] #Direct surface flux, erg/s/cm2/nm. 
flux_gnd_dif_0_85=importeddata[:,6] #Diffuse surface flux, erg/s/cm2/nm. 

########################
###Next, the A=0.20 case
########################
#Zenith angle of 0.
importeddata=np.genfromtxt('./TwoStreamOutput/rugheimer_earth_epoch0_w0=1e-5_a=0.2_z=0.dat', skip_header=1, skip_footer=0)
wav_leftedges_p245_0=importeddata[:,0] #left edges of wavelength bins, nm
wav_rightedges_p245_0=importeddata[:,1] #right edges of wavelength bins, nm
wav_centers_p245_0=importeddata[:,2] #centers of wavelength bins, nm
flux_toa_p245_0=importeddata[:,3]*np.cos(0.*deg2rad) #top-of-atmosphere flux, erg/s/cm2/nm, converted from 
flux_gnd_p245_0=importeddata[:,4] #surface flux, erg/s/cm2/nm. 
flux_gnd_dir_p245_0=importeddata[:,5] #Direct surface flux, erg/s/cm2/nm. 
flux_gnd_dif_p245_0=importeddata[:,6] #Diffuse surface flux, erg/s/cm2/nm. 

#Zenith angle of 60
importeddata=np.genfromtxt('./TwoStreamOutput/rugheimer_earth_epoch0_w0=1e-5_a=0.2_z=60.dat', skip_header=1, skip_footer=0)
wav_leftedges_p245_60=importeddata[:,0] #left edges of wavelength bins, nm
wav_rightedges_p245_60=importeddata[:,1] #right edges of wavelength bins, nm
wav_centers_p245_60=importeddata[:,2] #centers of wavelength bins, nm
flux_toa_p245_60=importeddata[:,3]*np.cos(60.*deg2rad) #top-of-atmosphere flux, erg/s/cm2/nm, converted from 
flux_gnd_p245_60=importeddata[:,4] #surface flux, erg/s/cm2/nm. 
flux_gnd_dir_p245_60=importeddata[:,5] #Direct surface flux, erg/s/cm2/nm. 
flux_gnd_dif_p245_60=importeddata[:,6] #Diffuse surface flux, erg/s/cm2/nm. 

#Zenith angle of 85
importeddata=np.genfromtxt('./TwoStreamOutput/rugheimer_earth_epoch0_w0=1e-5_a=0.2_z=85.dat', skip_header=1, skip_footer=0)
wav_leftedges_p245_85=importeddata[:,0] #left edges of wavelength bins, nm
wav_rightedges_p245_85=importeddata[:,1] #right edges of wavelength bins, nm
wav_centers_p245_85=importeddata[:,2] #centers of wavelength bins, nm
flux_toa_p245_85=importeddata[:,3]*np.cos(85.*deg2rad) #top-of-atmosphere flux, erg/s/cm2/nm, converted from 
flux_gnd_p245_85=importeddata[:,4] #surface flux, erg/s/cm2/nm. 
flux_gnd_dir_p245_85=importeddata[:,5] #Direct surface flux, erg/s/cm2/nm. 
flux_gnd_dif_p245_85=importeddata[:,6] #Diffuse surface flux, erg/s/cm2/nm. 

########################
###Next, the A=1.0 case
########################
#Zenith angle of 0.
importeddata=np.genfromtxt('./TwoStreamOutput/rugheimer_earth_epoch0_w0=1e-5_a=1_z=0.dat', skip_header=1, skip_footer=0)
wav_leftedges_1_0=importeddata[:,0] #left edges of wavelength bins, nm
wav_rightedges_1_0=importeddata[:,1] #right edges of wavelength bins, nm
wav_centers_1_0=importeddata[:,2] #centers of wavelength bins, nm
flux_toa_1_0=importeddata[:,3]*np.cos(0.*deg2rad) #top-of-atmosphere flux, erg/s/cm2/nm, converted from 
flux_gnd_1_0=importeddata[:,4] #surface flux, erg/s/cm2/nm. 
flux_gnd_dir_1_0=importeddata[:,5] #Direct surface flux, erg/s/cm2/nm. 
flux_gnd_dif_1_0=importeddata[:,6] #Diffuse surface flux, erg/s/cm2/nm. 

#Zenith angle of 60
importeddata=np.genfromtxt('./TwoStreamOutput/rugheimer_earth_epoch0_w0=1e-5_a=1_z=60.dat', skip_header=1, skip_footer=0)
wav_leftedges_1_60=importeddata[:,0] #left edges of wavelength bins, nm
wav_rightedges_1_60=importeddata[:,1] #right edges of wavelength bins, nm
wav_centers_1_60=importeddata[:,2] #centers of wavelength bins, nm
flux_toa_1_60=importeddata[:,3]*np.cos(60.*deg2rad) #top-of-atmosphere flux, erg/s/cm2/nm, converted from 
flux_gnd_1_60=importeddata[:,4] #surface flux, erg/s/cm2/nm. 
flux_gnd_dir_1_60=importeddata[:,5] #Direct surface flux, erg/s/cm2/nm. 
flux_gnd_dif_1_60=importeddata[:,6] #Diffuse surface flux, erg/s/cm2/nm. 

#Zenith angle of 85
importeddata=np.genfromtxt('./TwoStreamOutput/rugheimer_earth_epoch0_w0=1e-5_a=1_z=85.dat', skip_header=1, skip_footer=0)
wav_leftedges_1_85=importeddata[:,0] #left edges of wavelength bins, nm
wav_rightedges_1_85=importeddata[:,1] #right edges of wavelength bins, nm
wav_centers_1_85=importeddata[:,2] #centers of wavelength bins, nm
flux_toa_1_85=importeddata[:,3]*np.cos(85.*deg2rad) #top-of-atmosphere flux, erg/s/cm2/nm, converted from 
flux_gnd_1_85=importeddata[:,4] #surface flux, erg/s/cm2/nm. 
flux_gnd_dir_1_85=importeddata[:,5] #Direct surface flux, erg/s/cm2/nm. 
flux_gnd_dif_1_85=importeddata[:,6] #Diffuse surface flux, erg/s/cm2/nm. 

########################
###Plot results
########################
fig, ax=plt.subplots(3,2, figsize=(cm2inch(16.5),10), sharex=True)
ax00=ax[0,0]
ax01=ax[0,1]
ax10=ax[1,0]
ax11=ax[1,1]
ax20=ax[2,0]
ax21=ax[2,1]
markersizeval=5.

#First row of plots: zenith angle of 0, all different albedos
ax00.set_title(r'SZA=0$^{\circ}$')
ax00.plot(wav_centers_0_0, flux_toa_0_0,  marker='*', markersize=markersizeval, linewidth=1, color='black', label='TOA Flux')
ax00.plot(wav_centers_0_0, flux_gnd_dir_0_0,  marker='o', markersize=markersizeval, linewidth=1, color='blue', label='Direct Flux (A=0.0)')
ax00.plot(wav_centers_0_0, flux_gnd_dif_0_0,  marker='^', markersize=markersizeval, linewidth=1, linestyle='--', color='blue', label='Diffuse Flux (A=0.0)')
ax00.plot(wav_centers_p245_0, flux_gnd_dir_p245_0,  marker='o', markersize=markersizeval, linewidth=1, color='orange', label='Direct Flux (A=0.2)')
ax00.plot(wav_centers_p245_0, flux_gnd_dif_p245_0,  marker='^', markersize=markersizeval, linewidth=1, linestyle='--', color='orange', label='Diffuse Flux (A=0.2)')
ax00.plot(wav_centers_1_0, flux_gnd_dir_1_0,  marker='o', markersize=markersizeval, linewidth=1, color='red', label='Direct Flux (A=1)')
ax00.plot(wav_centers_1_0, flux_gnd_dif_1_0,  marker='^', markersize=markersizeval, linewidth=1, linestyle='--', color='red', label='Diffuse Flux (A=1)')
ax00.set_yscale('log')
ax00.set_ylim([1.e-5, 1e4])
ax00.set_ylabel(r'Flux (erg s$^{-1}$cm$^{-2}$nm$^{-1}$)', fontsize=11)

ax00.set_xlim([130, 855])

ax01.set_title(r'SZA=0$^{\circ}$')
ax01.plot(wav_centers_0_0, flux_gnd_dif_0_0/flux_toa_0_0,  marker='^', markersize=markersizeval, linewidth=1, linestyle='--', color='blue', label='A=0.0')
ax01.plot(wav_centers_p245_0, flux_gnd_dif_p245_0/flux_toa_p245_0,  marker='^', markersize=markersizeval, linewidth=1, linestyle='--', color='orange', label='A=0.245')
ax01.plot(wav_centers_1_0, flux_gnd_dif_1_0/flux_toa_1_0,  marker='^', markersize=markersizeval, linewidth=1, linestyle='--', color='red', label='A=1')
ax01.set_yscale('log')
ax01.set_ylim([1.e-13, 1e-5])
ax01.set_ylabel(r'Diffuse Surface Flux/TOA Flux', fontsize=11)
ax01.set_xlim([130, 855])


#Second row of plots: zenith angle of 60, all different albedos
ax10.set_title(r'SZA=60$^{\circ}$')
ax10.plot(wav_centers_0_60, flux_toa_0_60,  marker='*', markersize=markersizeval, linewidth=1, color='black', label='TOA Flux')
ax10.plot(wav_centers_0_60, flux_gnd_dir_0_60,  marker='o', markersize=markersizeval, linewidth=1, color='blue', label='Direct Flux\n (A=0.0)')
ax10.plot(wav_centers_0_60, flux_gnd_dif_0_60,  marker='^', markersize=markersizeval, linewidth=1, linestyle='--', color='blue', label='Diffuse Flux\n (A=0.0)')
ax10.plot(wav_centers_p245_60, flux_gnd_dir_p245_60,  marker='o', markersize=markersizeval, linewidth=1, color='orange', label='Direct Flux\n (A=0.245)')
ax10.plot(wav_centers_p245_60, flux_gnd_dif_p245_60,  marker='^', markersize=markersizeval, linewidth=1, linestyle='--', color='orange', label='Diffuse Flux\n (A=0.245)')
ax10.plot(wav_centers_1_60, flux_gnd_dir_1_60,  marker='o', markersize=markersizeval, linewidth=1, color='red', label='Direct Flux\n (A=1)')
ax10.plot(wav_centers_1_60, flux_gnd_dif_1_60,  marker='^', markersize=markersizeval, linewidth=1, linestyle='--', color='red', label='Diffuse Flux\n (A=1)')
ax10.set_yscale('log')
ax10.set_ylim([1.e-5, 1e4])
ax10.set_ylabel(r'Flux (erg s$^{-1}$cm$^{-2}$nm$^{-1}$)', fontsize=11)
ax10.set_xlim([130, 855])

ax11.set_title(r'SZA=60$^{\circ}$')
ax11.plot(wav_centers_0_60, flux_gnd_dif_0_60/flux_toa_0_60,  marker='^', markersize=markersizeval, linewidth=1, linestyle='--', color='blue', label='A=0.0')
ax11.plot(wav_centers_p245_60, flux_gnd_dif_p245_60/flux_toa_p245_60,  marker='^', markersize=markersizeval, linewidth=1, linestyle='--', color='orange', label='A=0.0')
ax11.plot(wav_centers_1_60, flux_gnd_dif_1_60/flux_toa_1_60,  marker='^', markersize=markersizeval, linewidth=1, linestyle='--', color='red', label='A=1')
ax11.set_yscale('log')
ax11.set_ylim([1.e-13, 1e-5])
ax11.set_ylabel(r'Diffuse Surface Flux/TOA Flux', fontsize=11)
ax11.set_xlim([130, 855])

#Last row of plots: zenith angle of 85, all different albedos
ax20.set_title(r'SZA=85$^{\circ}$')
ax20.plot(wav_centers_0_85, flux_toa_0_85,  marker='*', markersize=markersizeval, linewidth=1, color='black', label='TOA Flux')
ax20.plot(wav_centers_0_85, flux_gnd_dir_0_85,  marker='o', markersize=markersizeval, linewidth=1, color='blue', label='Direct Flux (A=0.0)')
ax20.plot(wav_centers_0_85, flux_gnd_dif_0_85,  marker='^', markersize=markersizeval, linewidth=1, linestyle='--', color='blue', label='Diffuse Flux (A=0.0)')
ax20.plot(wav_centers_p245_85, flux_gnd_dir_p245_85,  marker='o', markersize=markersizeval, linewidth=1, color='orange', label='Direct Flux (A=0.245)')
ax20.plot(wav_centers_p245_85, flux_gnd_dif_p245_85,  marker='^', markersize=markersizeval, linewidth=1, linestyle='--', color='orange', label='Diffuse Flux (A=0.245)')
ax20.plot(wav_centers_1_85, flux_gnd_dir_1_85,  marker='o', markersize=markersizeval, linewidth=1, color='red', label='Direct Flux (A=1)')
ax20.plot(wav_centers_1_85, flux_gnd_dif_1_85,  marker='^', markersize=markersizeval, linewidth=1, linestyle='--', color='red', label='Diffuse Flux (A=1)')
ax20.set_yscale('log')
ax20.set_ylim([1.e-5, 1e4])
ax20.set_ylabel(r'Flux (erg s$^{-1}$cm$^{-2}$nm$^{-1}$)', fontsize=11)
ax20.set_xlabel('nm')
ax20.set_xlim([130, 855])

ax21.set_title(r'SZA=85$^{\circ}$')
ax21.plot(wav_centers_0_85, flux_gnd_dif_0_85/flux_toa_0_85,  marker='^', markersize=markersizeval, linewidth=1, linestyle='--', color='blue', label='A=0.0')
ax21.plot(wav_centers_p245_85, flux_gnd_dif_p245_85/flux_toa_p245_85,  marker='^', markersize=markersizeval, linewidth=1, linestyle='--', color='orange', label='A=0.245')
ax21.plot(wav_centers_1_85, flux_gnd_dif_1_85/flux_toa_1_85,  marker='^', markersize=markersizeval, linewidth=1, linestyle='--', color='red', label='A=1')
ax21.set_yscale('log')
ax21.set_ylim([1.e-13, 1e-5])
ax21.set_ylabel(r'Diffuse Surface Flux/TOA Flux', fontsize=11)
ax21.set_xlabel('nm')
ax21.set_xlim([130, 855])

ax00.legend(bbox_to_anchor=[0, 1.22, 1., .152], loc=3, ncol=1, mode='expand', borderaxespad=0., fontsize=11)
ax01.legend(bbox_to_anchor=[0, 1.22, 1., .152], loc=3, ncol=1, mode='expand', borderaxespad=0., fontsize=11)
plt.tight_layout(rect=(0,0,1,0.8))
plt.savefig('./Plots/pure_absorption_test.eps', orientation='portrait',papertype='letter', format='eps')
plt.show()