# -*- coding: iso-8859-1 -*-
"""
This code plots the scattering limit test of our code. 
"""
########################
###Import useful libraries
########################

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pdb
import pickle

deg2rad=np.pi/180.
def cm2inch(cm): #function to convert cm to inches; useful for complying with Astrobiology size guidelines
	return cm/2.54

########################
###A=0
########################
###z=0
fnetdict=pickle.load(open('./TwoStreamOutput/rugheimer_earth_epoch0_w0=1-1e-12_a=0_z=0.p','rb'))
F_net=fnetdict['F_net'] #net flux in each layer, 0th layer is TOA, erg/s/cm2/nm
wav_leftedges=fnetdict['wav_leftedges'] #nm
wav_rightedges=fnetdict['wav_rightedges'] #nm
wav_centers=fnetdict['wav_centers'] #nm
z_lower=fnetdict['z_lower'] #cm, 0th layer is TOA
z_upper=fnetdict['z_upper'] #cm
z_center=fnetdict['z_center'] #cm
flux_toa=fnetdict['flux_toa'] #TOA "flux" (really intensity) in erg/s/cm2/nm (cgs)
solarzenithangle=fnetdict['solarzenithangle'] #radians

N_wavelengths=np.size(wav_centers) ###NOTE: We assume wavelength structure is the same for all of these (!!!)

direct_flux_toa=np.cos(solarzenithangle)*flux_toa #true TOA flux
F_net_deviation=np.zeros(np.shape(F_net))
F_net_deviation_max=np.zeros(N_wavelengths)
F_net_deviation_stddevs=np.zeros(N_wavelengths)
F_net_deviation_median=np.zeros(N_wavelengths)
for ind in range(0, N_wavelengths):
	median_val=np.median(F_net[:,ind])
	F_net_deviation_median[ind]=median_val
	F_net_deviation[:,ind]=F_net[:,ind]-median_val
	F_net_deviation_max[ind]=np.max(np.abs(F_net_deviation[:,ind]))
	F_net_deviation_stddevs[ind]=np.std(F_net[:,ind])

F_net_deviation_max_normalized_0_0=F_net_deviation_max/(direct_flux_toa)
F_net_deviation_stddevs_normalized_0_0=F_net_deviation_stddevs/(direct_flux_toa)

###z=60
fnetdict=pickle.load(open('./TwoStreamOutput/rugheimer_earth_epoch0_w0=1-1e-12_a=0_z=60.p','rb'))
F_net=fnetdict['F_net'] #net flux in each layer, 0th layer is TOA, erg/s/cm2/nm
wav_leftedges=fnetdict['wav_leftedges'] #nm
wav_rightedges=fnetdict['wav_rightedges'] #nm
wav_centers=fnetdict['wav_centers'] #nm
z_lower=fnetdict['z_lower'] #cm, 0th layer is TOA
z_upper=fnetdict['z_upper'] #cm
z_center=fnetdict['z_center'] #cm
flux_toa=fnetdict['flux_toa'] #TOA "flux" (really intensity) in erg/s/cm2/nm (cgs)
solarzenithangle=fnetdict['solarzenithangle'] #radians

direct_flux_toa=np.cos(solarzenithangle)*flux_toa #true TOA flux
F_net_deviation=np.zeros(np.shape(F_net))
F_net_deviation_max=np.zeros(N_wavelengths)
F_net_deviation_stddevs=np.zeros(N_wavelengths)
F_net_deviation_median=np.zeros(N_wavelengths)
for ind in range(0, N_wavelengths):
	median_val=np.median(F_net[:,ind])
	F_net_deviation_median[ind]=median_val
	F_net_deviation[:,ind]=F_net[:,ind]-median_val
	F_net_deviation_max[ind]=np.max(np.abs(F_net_deviation[:,ind]))
	F_net_deviation_stddevs[ind]=np.std(F_net[:,ind])

F_net_deviation_max_normalized_0_60=F_net_deviation_max/(direct_flux_toa)
F_net_deviation_stddevs_normalized_0_60=F_net_deviation_stddevs/(direct_flux_toa)

###z=85
fnetdict=pickle.load(open('./TwoStreamOutput/rugheimer_earth_epoch0_w0=1-1e-12_a=0_z=85.p','rb'))
F_net=fnetdict['F_net'] #net flux in each layer, 0th layer is TOA, erg/s/cm2/nm
wav_leftedges=fnetdict['wav_leftedges'] #nm
wav_rightedges=fnetdict['wav_rightedges'] #nm
wav_centers=fnetdict['wav_centers'] #nm
z_lower=fnetdict['z_lower'] #cm, 0th layer is TOA
z_upper=fnetdict['z_upper'] #cm
z_center=fnetdict['z_center'] #cm
flux_toa=fnetdict['flux_toa'] #TOA "flux" (really intensity) in erg/s/cm2/nm (cgs)
solarzenithangle=fnetdict['solarzenithangle'] #radians

direct_flux_toa=np.cos(solarzenithangle)*flux_toa #true TOA flux
F_net_deviation=np.zeros(np.shape(F_net))
F_net_deviation_max=np.zeros(N_wavelengths)
F_net_deviation_stddevs=np.zeros(N_wavelengths)
F_net_deviation_median=np.zeros(N_wavelengths)
for ind in range(0, N_wavelengths):
	median_val=np.median(F_net[:,ind])
	F_net_deviation_median[ind]=median_val
	F_net_deviation[:,ind]=F_net[:,ind]-median_val
	F_net_deviation_max[ind]=np.max(np.abs(F_net_deviation[:,ind]))
	F_net_deviation_stddevs[ind]=np.std(F_net[:,ind])

F_net_deviation_max_normalized_0_85=F_net_deviation_max/(direct_flux_toa)
F_net_deviation_stddevs_normalized_0_85=F_net_deviation_stddevs/(direct_flux_toa)

########################
###A=0.20
########################
###z=0
fnetdict=pickle.load(open('./TwoStreamOutput/rugheimer_earth_epoch0_w0=1-1e-12_a=0.2_z=0.p','rb'))
F_net=fnetdict['F_net'] #net flux in each layer, 0th layer is TOA, erg/s/cm2/nm
wav_leftedges=fnetdict['wav_leftedges'] #nm
wav_rightedges=fnetdict['wav_rightedges'] #nm
wav_centers=fnetdict['wav_centers'] #nm
z_lower=fnetdict['z_lower'] #cm, 0th layer is TOA
z_upper=fnetdict['z_upper'] #cm
z_center=fnetdict['z_center'] #cm
flux_toa=fnetdict['flux_toa'] #TOA "flux" (really intensity) in erg/s/cm2/nm (cgs)
solarzenithangle=fnetdict['solarzenithangle'] #radians

direct_flux_toa=np.cos(solarzenithangle)*flux_toa #true TOA flux
F_net_deviation=np.zeros(np.shape(F_net))
F_net_deviation_max=np.zeros(N_wavelengths)
F_net_deviation_stddevs=np.zeros(N_wavelengths)
F_net_deviation_median=np.zeros(N_wavelengths)
for ind in range(0, N_wavelengths):
	median_val=np.median(F_net[:,ind])
	F_net_deviation_median[ind]=median_val
	F_net_deviation[:,ind]=F_net[:,ind]-median_val
	F_net_deviation_max[ind]=np.max(np.abs(F_net_deviation[:,ind]))
	F_net_deviation_stddevs[ind]=np.std(F_net[:,ind])

F_net_deviation_max_normalized_p245_0=F_net_deviation_max/(direct_flux_toa)
F_net_deviation_stddevs_normalized_p245_0=F_net_deviation_stddevs/(direct_flux_toa)

###z=60
fnetdict=pickle.load(open('./TwoStreamOutput/rugheimer_earth_epoch0_w0=1-1e-12_a=0.2_z=60.p','rb'))
F_net=fnetdict['F_net'] #net flux in each layer, 0th layer is TOA, erg/s/cm2/nm
wav_leftedges=fnetdict['wav_leftedges'] #nm
wav_rightedges=fnetdict['wav_rightedges'] #nm
wav_centers=fnetdict['wav_centers'] #nm
z_lower=fnetdict['z_lower'] #cm, 0th layer is TOA
z_upper=fnetdict['z_upper'] #cm
z_center=fnetdict['z_center'] #cm
flux_toa=fnetdict['flux_toa'] #TOA "flux" (really intensity) in erg/s/cm2/nm (cgs)
solarzenithangle=fnetdict['solarzenithangle'] #radians

direct_flux_toa=np.cos(solarzenithangle)*flux_toa #true TOA flux
F_net_deviation=np.zeros(np.shape(F_net))
F_net_deviation_max=np.zeros(N_wavelengths)
F_net_deviation_stddevs=np.zeros(N_wavelengths)
F_net_deviation_median=np.zeros(N_wavelengths)
for ind in range(0, N_wavelengths):
	median_val=np.median(F_net[:,ind])
	F_net_deviation_median[ind]=median_val
	F_net_deviation[:,ind]=F_net[:,ind]-median_val
	F_net_deviation_max[ind]=np.max(np.abs(F_net_deviation[:,ind]))
	F_net_deviation_stddevs[ind]=np.std(F_net[:,ind])

F_net_deviation_max_normalized_p245_60=F_net_deviation_max/(direct_flux_toa)
F_net_deviation_stddevs_normalized_p245_60=F_net_deviation_stddevs/(direct_flux_toa)

###z=85
fnetdict=pickle.load(open('./TwoStreamOutput/rugheimer_earth_epoch0_w0=1-1e-12_a=0.2_z=85.p','rb'))
F_net=fnetdict['F_net'] #net flux in each layer, 0th layer is TOA, erg/s/cm2/nm
wav_leftedges=fnetdict['wav_leftedges'] #nm
wav_rightedges=fnetdict['wav_rightedges'] #nm
wav_centers=fnetdict['wav_centers'] #nm
z_lower=fnetdict['z_lower'] #cm, 0th layer is TOA
z_upper=fnetdict['z_upper'] #cm
z_center=fnetdict['z_center'] #cm
flux_toa=fnetdict['flux_toa'] #TOA "flux" (really intensity) in erg/s/cm2/nm (cgs)
solarzenithangle=fnetdict['solarzenithangle'] #radians

direct_flux_toa=np.cos(solarzenithangle)*flux_toa #true TOA flux
F_net_deviation=np.zeros(np.shape(F_net))
F_net_deviation_max=np.zeros(N_wavelengths)
F_net_deviation_stddevs=np.zeros(N_wavelengths)
F_net_deviation_median=np.zeros(N_wavelengths)
for ind in range(0, N_wavelengths):
	median_val=np.median(F_net[:,ind])
	F_net_deviation_median[ind]=median_val
	F_net_deviation[:,ind]=F_net[:,ind]-median_val
	F_net_deviation_max[ind]=np.max(np.abs(F_net_deviation[:,ind]))
	F_net_deviation_stddevs[ind]=np.std(F_net[:,ind])

F_net_deviation_max_normalized_p245_85=F_net_deviation_max/(direct_flux_toa)
F_net_deviation_stddevs_normalized_p245_85=F_net_deviation_stddevs/(direct_flux_toa)

########################
###A=1
########################
###z=0
fnetdict=pickle.load(open('./TwoStreamOutput/rugheimer_earth_epoch0_w0=1-1e-12_a=1_z=0.p','rb'))
F_net=fnetdict['F_net'] #net flux in each layer, 0th layer is TOA, erg/s/cm2/nm
wav_leftedges=fnetdict['wav_leftedges'] #nm
wav_rightedges=fnetdict['wav_rightedges'] #nm
wav_centers=fnetdict['wav_centers'] #nm
z_lower=fnetdict['z_lower'] #cm, 0th layer is TOA
z_upper=fnetdict['z_upper'] #cm
z_center=fnetdict['z_center'] #cm
flux_toa=fnetdict['flux_toa'] #TOA "flux" (really intensity) in erg/s/cm2/nm (cgs)
solarzenithangle=fnetdict['solarzenithangle'] #radians

direct_flux_toa=np.cos(solarzenithangle)*flux_toa #true TOA flux
F_net_deviation=np.zeros(np.shape(F_net))
F_net_deviation_max=np.zeros(N_wavelengths)
F_net_deviation_stddevs=np.zeros(N_wavelengths)
F_net_deviation_median=np.zeros(N_wavelengths)
for ind in range(0, N_wavelengths):
	median_val=np.median(F_net[:,ind])
	F_net_deviation_median[ind]=median_val
	F_net_deviation[:,ind]=F_net[:,ind]-median_val
	F_net_deviation_max[ind]=np.max(np.abs(F_net_deviation[:,ind]))
	F_net_deviation_stddevs[ind]=np.std(F_net[:,ind])

F_net_deviation_max_normalized_1_0=F_net_deviation_max/(direct_flux_toa)
F_net_deviation_stddevs_normalized_1_0=F_net_deviation_stddevs/(direct_flux_toa)
#pdb.set_trace()
###z=60
fnetdict=pickle.load(open('./TwoStreamOutput/rugheimer_earth_epoch0_w0=1-1e-12_a=1_z=60.p','rb'))
F_net=fnetdict['F_net'] #net flux in each layer, 0th layer is TOA, erg/s/cm2/nm
wav_leftedges=fnetdict['wav_leftedges'] #nm
wav_rightedges=fnetdict['wav_rightedges'] #nm
wav_centers=fnetdict['wav_centers'] #nm
z_lower=fnetdict['z_lower'] #cm, 0th layer is TOA
z_upper=fnetdict['z_upper'] #cm
z_center=fnetdict['z_center'] #cm
flux_toa=fnetdict['flux_toa'] #TOA "flux" (really intensity) in erg/s/cm2/nm (cgs)
solarzenithangle=fnetdict['solarzenithangle'] #radians

direct_flux_toa=np.cos(solarzenithangle)*flux_toa #true TOA flux
F_net_deviation=np.zeros(np.shape(F_net))
F_net_deviation_max=np.zeros(N_wavelengths)
F_net_deviation_stddevs=np.zeros(N_wavelengths)
F_net_deviation_median=np.zeros(N_wavelengths)
for ind in range(0, N_wavelengths):
	median_val=np.median(F_net[:,ind])
	F_net_deviation_median[ind]=median_val
	F_net_deviation[:,ind]=F_net[:,ind]-median_val
	F_net_deviation_max[ind]=np.max(np.abs(F_net_deviation[:,ind]))
	F_net_deviation_stddevs[ind]=np.std(F_net[:,ind])

F_net_deviation_max_normalized_1_60=F_net_deviation_max/(direct_flux_toa)
F_net_deviation_stddevs_normalized_1_60=F_net_deviation_stddevs/(direct_flux_toa)

###z=85
fnetdict=pickle.load(open('./TwoStreamOutput/rugheimer_earth_epoch0_w0=1-1e-12_a=1_z=85.p','rb'))
F_net=fnetdict['F_net'] #net flux in each layer, 0th layer is TOA, erg/s/cm2/nm
wav_leftedges=fnetdict['wav_leftedges'] #nm
wav_rightedges=fnetdict['wav_rightedges'] #nm
wav_centers=fnetdict['wav_centers'] #nm
z_lower=fnetdict['z_lower'] #cm, 0th layer is TOA
z_upper=fnetdict['z_upper'] #cm
z_center=fnetdict['z_center'] #cm
flux_toa=fnetdict['flux_toa'] #TOA "flux" (really intensity) in erg/s/cm2/nm (cgs)
solarzenithangle=fnetdict['solarzenithangle'] #radians

direct_flux_toa=np.cos(solarzenithangle)*flux_toa #true TOA flux
F_net_deviation=np.zeros(np.shape(F_net))
F_net_deviation_max=np.zeros(N_wavelengths)
F_net_deviation_stddevs=np.zeros(N_wavelengths)
F_net_deviation_median=np.zeros(N_wavelengths)
for ind in range(0, N_wavelengths):
	median_val=np.median(F_net[:,ind])
	F_net_deviation_median[ind]=median_val
	F_net_deviation[:,ind]=F_net[:,ind]-median_val
	F_net_deviation_max[ind]=np.max(np.abs(F_net_deviation[:,ind]))
	F_net_deviation_stddevs[ind]=np.std(F_net[:,ind])

F_net_deviation_max_normalized_1_85=F_net_deviation_max/(direct_flux_toa)
F_net_deviation_stddevs_normalized_1_85=F_net_deviation_stddevs/(direct_flux_toa)
#pdb.set_trace()
#########################
####Plot results
#########################
fig, (ax0, ax1, ax2)=plt.subplots(3, figsize=(cm2inch(16.5),10), sharex=True)
markersizeval=5.

#First row of plots: zenith angle of 0, all different albedos
ax0.set_title(r'SZA=0$^{\circ}$')
ax0.plot(wav_centers, F_net_deviation_max_normalized_0_0,  marker='o', markersize=markersizeval, linewidth=1, color='black', label='A=0')
ax0.plot(wav_centers, F_net_deviation_max_normalized_p245_0,  marker='o', markersize=markersizeval, linewidth=1, color='blue', label='A=0.2')
ax0.plot(wav_centers, F_net_deviation_max_normalized_1_0,  marker='o', markersize=markersizeval, linewidth=1, color='red', label='A=1')
ax0.set_yscale('log')
ax0.set_ylim([1.e-14, 1e-5])
ax0.set_ylabel('Maximum Fractional Columnar\nDeviation of $F_{net}$ from Median ')
#ax0.legend(loc=0, fontsize=10)
ax0.set_xlim([130, 855])

##Second row of plots: zenith angle of 60, all different albedos
ax1.set_title(r'SZA=60$^{\circ}$')
ax1.plot(wav_centers, F_net_deviation_max_normalized_0_60,  marker='o', markersize=markersizeval, linewidth=1, color='black', label='A=0')
ax1.plot(wav_centers, F_net_deviation_max_normalized_p245_60,  marker='o', markersize=markersizeval, linewidth=1, color='blue', label='A=0.2')
ax1.plot(wav_centers, F_net_deviation_max_normalized_1_60,  marker='o', markersize=markersizeval, linewidth=1, color='red', label='A=1')
ax1.set_yscale('log')
ax1.set_ylim([1.e-14, 1e-5])
ax1.set_ylabel('Maximum Fractional Columnar\nDeviation of $F_{net}$ from Median')
#ax1.legend(loc=0, fontsize=10)
ax1.set_xlim([130, 855])

#Last row of plots: zenith angle of 85, all different albedos
ax2.set_title(r'SZA=85$^{\circ}$')
ax2.plot(wav_centers, F_net_deviation_max_normalized_0_85,  marker='o', markersize=markersizeval, linewidth=1, color='black', label='A=0')
ax2.plot(wav_centers, F_net_deviation_max_normalized_p245_85,  marker='o', markersize=markersizeval, linewidth=1, color='blue', label='A=0.2')
ax2.plot(wav_centers, F_net_deviation_max_normalized_1_85,  marker='o', markersize=markersizeval, linewidth=1, color='red', label='A=1')
ax2.set_yscale('log')
ax2.set_ylim([1.e-14, 1e-5])
ax2.set_ylabel('Maximum Fractional Columnar\nDeviation of $F_{net}$ from Median')
#ax2.legend(loc=0, fontsize=10)
ax2.set_xlabel('nm')
ax2.set_xlim([130, 855])
#pdb.set_trace()

ax0.legend(bbox_to_anchor=[0, 1.22, 1., .152], loc=3, ncol=3, mode='expand', borderaxespad=0., fontsize=11)
plt.tight_layout(rect=(0,0,1,0.93))
plt.savefig('./Plots/pure_scattering_test.eps', orientation='portrait',papertype='letter', format='eps')
plt.show()