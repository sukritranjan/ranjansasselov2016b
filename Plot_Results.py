# -*- coding: iso-8859-1 -*-
"""
This code creates the Results plots for Ranjan & Sasselov 2016b.
"""
########################
###Import useful libraries
########################

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pdb
from matplotlib.pyplot import cm

def cm2inch(cm): #function to convert cm to inches; useful for complying with Astrobiology size guidelines
	return cm/2.54

########################
###Set important constants
########################

hc=1.98645e-9 #value of h*c in erg*nm, useful to convert from ergs/cm2/s/nm to photons/cm2/s/nm

########################
###Specific which family of plots to generate
########################
plot_intvsflux=False #Plot to demonstrate difference between surface radiance and surface flux.

plot_alb_zenithangle=False #plot to demonstrate the impact of albedo and zenith angle on emergent surface intensity
plot_co2_limits=False #this time, using a fixed level of N2 and variable amounts of CO2.
plot_altgas_limits=False #fixed level of N2, various levels of other gases.

plot_dosimeters_co2=True #Plots the convolution of the various action spectra with the surficial spectra for the co2 study, integrated, and normalized. Helps us see how these parameters vary. 
plot_dosimeters_h2o=True #Plots the convolution of the various action spectra with the surficial spectra for h2o
plot_dosimeters_so2=True #Plots the convolution of the various action spectra with the surficial spectra for so2
plot_dosimeters_h2s=True #Plots the convolution of the various action spectra with the surficial spectra for h2s
plot_dosimeters_ch4=False #Plots the convolution of the various action spectra with the surficial spectra for ch4
plot_dosimeters_o2=False #Plots the convolution of the various action spectra with the surficial spectra for o2
plot_dosimeters_o3=False #Plots the convolution of the various action spectra with the surficial spectra for o3


############################
#######Plot to demonstrate the difference between surface flux and surface intensity
############################
if plot_intvsflux:
	wav, toa_intensity, surface_flux, total_intensity,surface_intensity=np.genfromtxt('./TwoStreamOutput/rugheimer_earth_epoch0_a=newsnow_z=60.dat', skip_header=1, skip_footer=0, usecols=(2,3,4,5,6), unpack=True) 

	fig1, (ax1)=plt.subplots(1, figsize=(8,5), sharex=True)
	ax1.plot(wav, toa_intensity*wav/hc, linestyle='-', color='black',marker='s', label='TOA Flux')
	ax1.plot(wav, total_intensity*wav/hc, linestyle='-', color='red', marker='s',label='BOA Actinic Flux')
	ax1.plot(wav, surface_intensity*wav/hc,linestyle='-', color='purple', marker='s',label='Surface Radiance')
	ax1.plot(wav, surface_flux*wav/hc, linestyle='-', color='blue', marker='s',label='Surface Flux')
	ylimits=[1.e9, 1.e15]
	ax1.set_title('R+2015 Atmosphere, A=New Snow, SZA=60')
	ax1.legend(loc=0, ncol=1, borderaxespad=0.)
	ax1.set_yscale('log')
	ax1.set_ylabel(r'photons/s/cm$^2$/nm')
	ax1.set_ylim(ylimits)
	ax1.set_xlim([130.,500.])
	ax1.set_xlabel('nm')
	#plt.tight_layout(rect=(0,0,1,0.9))
	plt.savefig('./Plots/intvsflux.eps', orientation='portrait',papertype='letter', format='eps')
	


############################
#######Plot to demonstrate impact of albedo and zenith angle on emergent surface intensity. We use the Rugheimer prebiotic atmosphere model (mixing ratios, surface pressure, T/P profile) for this reference system. 
############################
if plot_alb_zenithangle:

	###When importing files, variables are:
	###wav_x: centers of wavelength bins, nm
	###toa_intensity_x: top-of-atmosphere intensity (incident), erg/s/nm/cm2
	###surface_flux_x: total flux incident on surface, erg/s/nm/cm2
	###surface_intensity_x: total intensity incident on the surface (direct+diffuse), erg/s/nm/cm2
	###surface_intensity_diffuse_x: diffuse total intensity incident on the surface, erg/s/nm/cm2
	###surface_intensity_direct_x: direct total intensity incident on the surface, erg/s/nm/cm2
	
	
	#####Set up info about the files to extract
	albedolist=['1', 'newsnow', 'oldsnow', '0.2', 'desert', 'ocean', 'tundra', '0'] #list of albedos we consider
	AlbedoLabels=['1', 'New Snow', 'Old Snow', '0.2', 'Desert', 'Ocean', 'Tundra', '0']#labels for the figure
	
	#####Plot Figure
	#Step 1: Initialize Figure
	fig1, (ax1, ax2, ax3)=plt.subplots(3, figsize=(cm2inch(16.5),10), sharex=True)
	colorseq1=iter(cm.rainbow(np.linspace(0,1,len(albedolist))))
	colorseq2=iter(cm.rainbow(np.linspace(0,1,len(albedolist))))
	colorseq3=iter(cm.rainbow(np.linspace(0,1,len(albedolist))))

	#Step 2: loop over all files, plot figure
	for albind in range(0, len(albedolist)):
		albedo=albedolist[albind]
		
		#Load file
		wav_z_0, toa_intensity_z_0, surface_flux_z_0,  surface_intensity_z_0, surface_intensity_diffuse_z_0, surface_intensity_direct_z_0=np.genfromtxt('./TwoStreamOutput/AlbZen/rugheimer_earth_epoch0_a='+albedo+'_z=0.dat', skip_header=1, skip_footer=0, usecols=(2,3,4,6,7,8), unpack=True) #albedo=desert, zenith angle=0 degrees
		wav_z_48p2, toa_intensity_z_48p2, surface_flux_z_48p2,  surface_intensity_z_48p2, surface_intensity_diffuse_z_48p2, surface_intensity_direct_z_48p2=np.genfromtxt('./TwoStreamOutput/AlbZen/rugheimer_earth_epoch0_a='+albedo+'_z=48.2.dat', skip_header=1, skip_footer=0, usecols=(2,3,4,6,7,8), unpack=True) #albedo=desert, zenith angle=0 degrees
		wav_z_66p5, toa_intensity_z_66p5, surface_flux_z_66p5,  surface_intensity_z_66p5, surface_intensity_diffuse_z_66p5, surface_intensity_direct_z_66p5=np.genfromtxt('./TwoStreamOutput/AlbZen/rugheimer_earth_epoch0_a='+albedo+'_z=66.5.dat', skip_header=1, skip_footer=0, usecols=(2,3,4,6,7,8), unpack=True) #albedo=desert, zenith angle=0 degrees
		
		if albind==0: #Initialize with TOA flux
			ax1.plot(wav_z_0, toa_intensity_z_0*wav_z_0/hc, marker='.', color='black', label=r'TOA Flux')
			ax2.plot(wav_z_48p2, toa_intensity_z_48p2*wav_z_48p2/hc, marker='.', color='black', label=r'TOA Flux')
			ax3.plot(wav_z_66p5, toa_intensity_z_66p5*wav_z_66p5/hc, marker='.', color='black', label=r'TOA Flux')
		
		ax1.plot(wav_z_0, surface_intensity_z_0*wav_z_0/hc, marker='.', color=next(colorseq1), label=r'A='+AlbedoLabels[albind])
		ax2.plot(wav_z_48p2, surface_intensity_z_48p2*wav_z_48p2/hc, marker='.', color=next(colorseq2), label=r'A='+AlbedoLabels[albind])
		ax3.plot(wav_z_66p5, surface_intensity_z_66p5*wav_z_66p5/hc, marker='.', color=next(colorseq3), label=r'A='+AlbedoLabels[albind])
		
	#Step 3: Clean up figure
	ylimits=[1.e9, 1.e15]
	ax1.set_title(r'z=0$^\circ$')
	ax1.legend(bbox_to_anchor=[0, 1.13, 1., .152], loc=3, ncol=3, mode='expand', borderaxespad=0., fontsize=10)
	ax2.set_title(r'z=48.2$^\circ$')
	ax3.set_title(r'z=66.5$^\circ$')
	ax1.set_yscale('log')
	ax1.set_ylabel(r'photons/s/cm$^2$/nm')
	ax1.set_ylim(ylimits)
	ax2.set_yscale('log')
	ax2.set_ylabel(r'photons/s/cm$^2$/nm')
	ax2.set_ylim(ylimits)
	ax3.set_yscale('log')
	ax3.set_ylabel(r'photons/s/cm$^2$/nm')
	ax3.set_ylim(ylimits)
	ax3.set_xlim([100.,500.])
	ax3.set_xlabel('nm')
	plt.tight_layout(rect=(0,0,1,0.9))
	plt.savefig('./Plots/paperplots_a_z_dependence.eps', orientation='portrait',papertype='letter', format='eps')

if plot_co2_limits:
	###When importing files, variables are:
	###ind 0: wav_x: centers of wavelength bins, nm
	###ind 1: toa_intensity_x: top-of-atmosphere intensity (incident), erg/s/nm/cm2
	###ind 2: surface_flux_x: total flux incident on surface, erg/s/nm/cm2
	###ind 3: surface_intensity_x: total intensity incident on the surface (direct+diffuse), erg/s/nm/cm2
	###ind 4: surface_intensity_diffuse_x: diffuse total intensity incident on the surface, erg/s/nm/cm2
	###ind 5: surface_intensity_direct_x: direct total intensity incident on the surface, erg/s/nm/cm2
	
	###############Set up info about files to extract
	N_co2_base=2.09e24 #column density of CO2 in base case (Rugheimer+2015)
	co2multiplelist=[0., 1.e-6,1.e-5, 1.e-4, 1.e-3, 0.00893, 1.e-2, 1.e-1, 0.6, 1., 1.33, 1.e1, 46.6, 1.e2, 470., 1.e3]
	co2dict={}
	isphysical=[False, False, False, False, False, True, False, False, True, True, True, False, True, False, True, False] #which of these models have a physically motivated column depth
	###############Read in base Rugheimer abundance cases
	wav_max_rugheimer, toa_intensity_max_rugheimer, surface_flux_max_rugheimer, surface_intensity_max_rugheimer, surface_intensity_diffuse_max_rugheimer, surface_intensity_direct_max_rugheimer=np.genfromtxt('./TwoStreamOutput/AlbZen/rugheimer_earth_epoch0_a=newsnow_z=0.dat', skip_header=1, skip_footer=0, usecols=(2,3,4,6,7,8), unpack=True)
	wav_min_rugheimer, toa_intensity_min_rugheimer, surface_flux_min_rugheimer, surface_intensity_min_rugheimer, surface_intensity_diffuse_min_rugheimer, surface_intensity_direct_min_rugheimer=np.genfromtxt('./TwoStreamOutput/AlbZen/rugheimer_earth_epoch0_a=tundra_z=66.5.dat', skip_header=1, skip_footer=0, usecols=(2,3,4,6,7,8), unpack=True) 
	
	##############Figure comparing surface intensity under different levels of pCO2, and different values for A and Z
	#Set up figure basics outside the loop.
	fig1, (ax1, ax2)=plt.subplots(2, figsize=(cm2inch(16.5),10), sharex=True)
	colorseq1=iter(cm.rainbow(np.linspace(0,1,len(co2multiplelist))))
	colorseq2=iter(cm.rainbow(np.linspace(0,1,len(co2multiplelist))))
	
	#Plot TOA intensities
	ax1.plot(wav_max_rugheimer, toa_intensity_max_rugheimer*wav_min_rugheimer/hc, linestyle='-', color='black', label='TOA Flux')
	ax2.plot(wav_min_rugheimer, toa_intensity_min_rugheimer*wav_min_rugheimer/hc, linestyle='-', color='black', label='TOA Flux')
	
	#In a loop, load the intensities and plot
	for ind in range(0, len(co2multiplelist)):
		multiple=co2multiplelist[ind]
		colden_co2=multiple*N_co2_base
		if isphysical[ind]: #have the physically motivated models represented differently
			linestylevar='-'
			linewidthvar=1.
		elif multiple==1.: #represent the base fiducial Rugheimer model with a different line
			linestylevar='-'
			linewidthvar=2.5
		else: #the parametric exploration 
			linestylevar='--'
			linewidthvar=1.
		
		wav_max, toa_intensity_max, surface_flux_max, surface_intensity_max, surface_intensity_diffuse_max, surface_intensity_direct_max=np.genfromtxt('./TwoStreamOutput/CO2lim/surface_intensities_co2limits_co2multiple='+str(multiple)+'_a=newsnow_z=0.dat', skip_header=1, skip_footer=0, usecols=(2,3,4,6,7,8), unpack=True) #maximum intensity for given atmosphere
		wav_min, toa_intensity_min, surface_flux_min, surface_intensity_min, surface_intensity_diffuse_min, surface_intensity_direct_min=np.genfromtxt('./TwoStreamOutput/CO2lim/surface_intensities_co2limits_co2multiple='+str(multiple)+'_a=tundra_z=66.5.dat', skip_header=1, skip_footer=0, usecols=(2,3,4,6,7,8), unpack=True) #minimum maximum intensity for given atmosphere
		co2dict[str(multiple)]=surface_intensity_max
		ax1.plot(wav_max,surface_intensity_max*wav_max/hc, linestyle=linestylevar, linewidth=linewidthvar, color=next(colorseq1), label=r'$N_{CO_{2}}=$'+'{:.2E}'.format(colden_co2)+' cm$^{-2}$')
		ax2.plot(wav_min,surface_intensity_min*wav_min/hc, linestyle=linestylevar, linewidth=linewidthvar, color=next(colorseq2), label=r'$N_{CO_{2}}=$'+'{:.2E}'.format(colden_co2)+' cm$^{-2}$')
	#print (co2dict[str(co2multiplelist[8])])/(co2dict[str(co2multiplelist[0])])
	#pdb.set_trace()
	#Set up fine detail on figure
	ylimits=[1e7, 1e15]
	ax1.set_title(r'z=0$^\circ$, A=Fresh Snow')
	ax2.set_title(r'z=66.5$^\circ$, A=Tundra')	
	ax1.legend(bbox_to_anchor=[0, 1.1, 1., .5], loc=3, ncol=2, mode='expand', borderaxespad=0., fontsize=10)
	ax1.set_yscale('log')
	ax1.set_ylabel(r'photons/s/cm$^2$/nm')
	ax1.set_ylim(ylimits)
	ax2.set_yscale('log')
	ax2.set_ylabel(r'photons/s/cm$^2$/nm')
	ax2.set_ylim(ylimits)
	ax2.set_xlim([100.,500.])
	ax2.set_xlabel('nm')
	plt.tight_layout(rect=(0,0,1,0.75))
	plt.savefig('./Plots/paperplots_co2_radiance.eps', orientation='portrait',papertype='letter', format='eps')
	
if plot_altgas_limits:
	###When importing files, variables are:
	###ind 0: wav_x: centers of wavelength bins, nm
	###ind 1: toa_intensity_x: top-of-atmosphere intensity (incident), erg/s/nm/cm2
	###ind 2: surface_flux_x: total flux incident on surface, erg/s/nm/cm2
	###ind 3: surface_intensity_x: total intensity incident on the surface (direct+diffuse), erg/s/nm/cm2
	###ind 4: surface_intensity_diffuse_x: diffuse total intensity incident on the surface, erg/s/nm/cm2
	###ind 5: surface_intensity_direct_x: direct total intensity incident on the surface, erg/s/nm/cm2
	
	#####Set up info about the files to extract ##Maximum possible natural surface radiance case (z=0, albedo=fresh snow) aka "max"
	N_tot=2.0925e25#total column density of Rugheimer+2015 model in cm**-2
	gaslist=['h2o', 'ch4', 'so2', 'o2', 'o3', 'h2s'] #list of gases we are doing this for
	gaslabellist=['H2O', 'CH4', 'SO2', 'O2', 'O3', 'H2S'] #list of nicely formated gas names for plotting
	base_abundances=np.array([4.762e-3, 1.647e-6, 3.371e-11, 2.707e-6, 9.160e-11, 6.742e-11]) #molar concentration of each of these gases in the Rugheimer model.
	#dict holding the multiples of the molar concentration we are using
	gasmultiples={}
	gasmultiples['h2o']=np.array([1.e-5, 1.e-4, 1.e-3, 1.e-2, 1.e-1, 1., 1.e1, 1.e2, 1.e3])
	gasmultiples['ch4']=np.array([1.e-2, 1.e-1, 1., 1.e1, 1.e2, 1.e3])
	gasmultiples['so2']=np.array([1., 1.e1, 1.e2, 1.e3, 1.e4, 1.e5, 1.e6, 1.e7])
	gasmultiples['o2']=np.array([1.e-5, 1.e-4, 1.e-3, 1.e-2, 1.e-1, 1., 1.e1, 1.e2, 1.e3, 1.e4, 1.e5])
	gasmultiples['o3']=np.array([1., 1.e1, 1.e2, 1.e3])
	gasmultiples['h2s']=np.array([1., 1.e1, 1.e2, 1.e3, 1.e4, 1.e5, 1.e6, 1.e7])
	
	
	#####In a loop, extract and plot the files
	for gasind in range(0, len(gaslist)):
		gas=gaslist[gasind]
		base_abundance=base_abundances[gasind]
		multiples=gasmultiples[gas]
		gaslabel=gaslabellist[gasind]
		
		#####Set up figure basics
		fig, (ax1)=plt.subplots(1, figsize=(cm2inch(16.5),7), sharex=True)
		colorseq=iter(cm.rainbow(np.linspace(0,1,len(multiples))))
		for multind in range(0, len(multiples)):
			multiple=multiples[multind]
			colden_X=base_abundance*multiple*N_tot #total column density of gas X

			datafile='./TwoStreamOutput/gaslim/surface_intensities_'+gas+'limits_'+gas+'multiple='+str(multiple)+'_a=newsnow_z=0.dat'
			wav, toa_intensity, surface_flux, surface_intensity, surface_intensity_diffuse, surface_intensity_direct=np.genfromtxt(datafile, skip_header=1, skip_footer=0, usecols=(2,3,4,6,7,8), unpack=True)
			
			if multind==0:
				ax1.plot(wav,toa_intensity*wav/hc, linestyle='-', linewidth=1, color='black', label=r'TOA Flux')

			if multiple==1.: #represent the base fiducial Rugheimer model with a different line
				linestylevar='-'
				linewidthvar=2.0
			else: #the parametric exploration 
				linestylevar='--'
				linewidthvar=1.
		
			ax1.plot(wav,surface_intensity*wav/hc, linestyle=linestylevar,linewidth=linewidthvar, color=next(colorseq), label=r'$N_{'+gaslabel+'}=$'+'{:.2E}'.format(colden_X)+' cm$^{-2}$')
		#####Finalize and save figure
		ax1.set_title(r'Varying Levels of '+gaslabel+r', (z=0$^\circ$, A=Fresh Snow)')
		ax1.legend(bbox_to_anchor=[0, 1.1, 1., .5], loc=3, ncol=2, mode='expand', borderaxespad=0., fontsize=10)
		ax1.set_yscale('log')
		ax1.set_ylabel(r'photons/s/cm$^2$/nm')
		ax1.set_ylim([1.e7, 1.e15])	
		ax1.set_xlim([100.,500.])
		ax1.set_xlabel('nm')
		plt.tight_layout(rect=(0,0,1,0.75))
		plt.savefig('./Plots/paperplots_'+gas+'_radiance.eps', orientation='portrait',papertype='letter', format='eps')


if plot_dosimeters_co2:
	###########First, import dosimeters
	SZAs, albedos, N_CO2s, rad100_165s, rad200_300s, umpgly_193s, umpgly_230s, umpgly_254s,cucn3_254s, cucn3_300s=np.genfromtxt('./Doses/co2_uv_doses.dat', skip_header=2, skip_footer=0, usecols=(0, 1, 2, 3, 4, 5, 6, 7, 8, 9), unpack=True, delimiter='	') #wavelength in nm, relative efficiency unitless
	
	max_inds=np.arange(0, len(albedos)/2) #first half of data is max radiance case (SZA=0, A=fresh snow)
	min_inds=np.arange(len(albedos)/2, len(albedos))#last half of the data is the minimum radiance case (SZA=66.5, A=tundra)
	wordsworthind=5
	kastingupperind=-4
	###########Now, plot the dosimeters vs CO2 concentration

	umpgly_193s_max_normed=umpgly_193s[max_inds]/umpgly_193s[max_inds[-7]]
	umpgly_230s_max_normed=umpgly_230s[max_inds]/umpgly_230s[max_inds[-7]]
	umpgly_254s_max_normed=umpgly_254s[max_inds]/umpgly_254s[max_inds[-7]]
	cucn3_254s_max_normed=cucn3_254s[max_inds]/cucn3_254s[max_inds[-7]]
	cucn3_300s_max_normed=cucn3_300s[max_inds]/cucn3_300s[max_inds[-7]]

	umpgly_193s_min_normed=umpgly_193s[min_inds]/umpgly_193s[min_inds[-7]]
	umpgly_230s_min_normed=umpgly_230s[min_inds]/umpgly_230s[min_inds[-7]]
	umpgly_254s_min_normed=umpgly_254s[min_inds]/umpgly_254s[min_inds[-7]]
	cucn3_254s_min_normed=cucn3_254s[min_inds]/cucn3_254s[min_inds[-7]]
	cucn3_300s_min_normed=cucn3_300s[min_inds]/cucn3_300s[min_inds[-7]]
	
	#Initialize plot basics
	fig=plt.figure(figsize=(cm2inch(16.5),7))
	gs=gridspec.GridSpec(2,2, hspace=0.40,wspace=0.35, width_ratios=[2,1], top=.77, bottom=.1, left=.1, right=.95)
	ax1=plt.subplot(gs[0])
	ax2=plt.subplot(gs[2])
	ax3=plt.subplot(gs[1])
	ax4=plt.subplot(gs[3])
	colorseq1=iter(cm.rainbow(np.linspace(0,1,5)))
	colorseq2=iter(cm.rainbow(np.linspace(0,1,5)))
	colorseq3=iter(cm.rainbow(np.linspace(0,1,5)))
	colorseq4=iter(cm.rainbow(np.linspace(0,1,5)))
	
	#Plot max case	
	ax1.plot(N_CO2s[max_inds], umpgly_193s_max_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq1), label=r'UMP Gly Bond Cleavage ($\lambda_0=193$)')
	ax1.plot(N_CO2s[max_inds], umpgly_230s_max_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq1), label=r'UMP Gly Bond Cleavage ($\lambda_0=230$)')
	ax1.plot(N_CO2s[max_inds], umpgly_254s_max_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq1), label=r'UMP Gly Bond Cleavage ($\lambda_0=254$)')
	ax1.plot(N_CO2s[max_inds], cucn3_254s_max_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq1), label=r'CuCN$_3$ Photoionization ($\lambda_0=254$)')
	ax1.plot(N_CO2s[max_inds], cucn3_300s_max_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq1), label=r'CuCN$_3$ Photoionization ($\lambda_0=300$)')
	ax1.set_title('SZA=0, Albedo=New Snow')
	ax1.axvline(N_CO2s[max_inds[-7]], color='black', linewidth=2) #Mark the Rugheimer fiducial value
	ax1.axhline(1., color='black', linewidth=1) #Mark the Rugheimer fiducial value
	ax1.axvline(N_CO2s[max_inds[wordsworthind]], color='black', linewidth=1, linestyle='--') #Wordsworth lower limit
	#ax1.axvline(N_CO2s[max_inds[kastingupperind]], color='black', linewidth=1, linestyle='--') #Kasting upper limit

	#Plot min case
	ax2.set_title('SZA=66.5, Albedo=Tundra')
	ax2.plot(N_CO2s[min_inds], umpgly_193s_min_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq2), label=r'UMP Gly Bond Cleavage ($\lambda_0=193$)')
	ax2.plot(N_CO2s[min_inds], umpgly_230s_min_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq2), label=r'UMP Gly Bond Cleavage ($\lambda_0=230$)')
	ax2.plot(N_CO2s[min_inds], umpgly_254s_min_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq2), label=r'UMP Gly Bond Cleavage ($\lambda_0=254$)')
	ax2.plot(N_CO2s[min_inds], cucn3_254s_min_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq2), label=r'CuCN$_3$ Photoionization ($\lambda_0=254$)')
	ax2.plot(N_CO2s[min_inds], cucn3_300s_min_normed,linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq2), label=r'CuCN$_3$ Photoionization ($\lambda_0=300$)')
	ax2.axvline(N_CO2s[min_inds[-7]], color='black', linewidth=2) #Mark the Rugheimer fiducial value
	ax2.axhline(1., color='black', linewidth=1) #Mark the Rugheimer fiducial value
	ax2.axvline(N_CO2s[min_inds[wordsworthind]], color='black', linewidth=1, linestyle='--') #Wordsworth lower limit
	#ax2.axvline(N_CO2s[min_inds[kastingupperind]], color='black', linewidth=1, linestyle='--') #Kasting upper limit
	
	#Plot max case	
	ax3.plot(N_CO2s[max_inds], umpgly_193s_max_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq3), label=r'UMP Gly Bond Cleavage ($\lambda_0=193$)')
	ax3.plot(N_CO2s[max_inds], umpgly_230s_max_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq3), label=r'UMP Gly Bond Cleavage ($\lambda_0=230$)')
	ax3.plot(N_CO2s[max_inds], umpgly_254s_max_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq3), label=r'UMP Gly Bond Cleavage ($\lambda_0=254$)')
	ax3.plot(N_CO2s[max_inds], cucn3_254s_max_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq3), label=r'CuCN$_3$ Photoionization ($\lambda_0=254$)')
	ax3.plot(N_CO2s[max_inds], cucn3_300s_max_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq3), label=r'CuCN$_3$ Photoionization ($\lambda_0=300$)')
	ax3.axvline(N_CO2s[max_inds[-7]], color='black', linewidth=2) #Mark the Rugheimer fiducial value
	ax3.axhline(1., color='black', linewidth=1) #Mark the Rugheimer fiducial value
	ax3.axvline(N_CO2s[max_inds[wordsworthind]], color='black', linewidth=1, linestyle='--') #Wordsworth lower limit
	#ax3.axvline(N_CO2s[max_inds[kastingupperind]], color='black', linewidth=1, linestyle='--') #Kasting upper limit
	
	#Plot min case
	ax4.plot(N_CO2s[min_inds], umpgly_193s_min_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq4), label=r'UMP Gly Bond Cleavage ($\lambda_0=193$)')
	ax4.plot(N_CO2s[min_inds], umpgly_230s_min_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq4), label=r'UMP Gly Bond Cleavage ($\lambda_0=230$)')
	ax4.plot(N_CO2s[min_inds], umpgly_254s_min_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq4), label=r'UMP Gly Bond Cleavage ($\lambda_0=254$)')
	ax4.plot(N_CO2s[min_inds], cucn3_254s_min_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq4), label=r'CuCN$_3$ Photoionization ($\lambda_0=254$)')
	ax4.plot(N_CO2s[min_inds], cucn3_300s_min_normed,linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq4), label=r'CuCN$_3$ Photoionization ($\lambda_0=300$)')
	ax4.axvline(N_CO2s[min_inds[-7]], color='black', linewidth=2) #Mark the Rugheimer fiducial value	
	ax4.axhline(1., color='black', linewidth=1)
	ax4.axvline(N_CO2s[min_inds[wordsworthind]], color='black', linewidth=1, linestyle='--') #Wordsworth lower limit
	#ax4.axvline(N_CO2s[min_inds[kastingupperind]], color='black', linewidth=1, linestyle='--') #Kasting upper limit

	#print umpgly_193s_max_normed
	#pdb.set_trace()
	
	#Finalize plot details.
	ax1.set_ylabel(r'Relative Dose Rate $D$')
	ax2.set_ylabel(r'Relative Dose Rate $D$')
	ax3.set_ylabel(r'Relative Dose Rate $D$')
	ax4.set_ylabel(r'Relative Dose Rate $D$')
	ax1.set_xlabel(r'N$_{CO2}$ (cm$^{-2}$)')
	ax2.set_xlabel(r'N$_{CO2}$ (cm$^{-2}$)')
	ax3.set_xlabel(r'N$_{CO2}$ (cm$^{-2}$)')
	ax4.set_xlabel(r'N$_{CO2}$ (cm$^{-2}$)')
	ax1.set_yscale('linear')
	ax2.set_yscale('linear')
	ax3.set_yscale('log')
	ax4.set_yscale('log')
	ax1.set_xscale('log')
	ax2.set_xscale('log')
	ax3.set_xscale('log')
	ax4.set_xscale('log')
	ax1.set_xlim([2.09e18, 2.09e24])
	ax2.set_xlim([2.09e18, 2.09e24])
	ax3.set_xlim([2.09e24, 9.85e26])
	ax4.set_xlim([2.09e24, 9.85e26])
	ax1.set_ylim([0.95,1.9])
	ax2.set_ylim([0.95,1.6])
	ax3.set_ylim([1.e-2, 1.e1])
	ax4.set_ylim([1.e-2, 1.e1])
	
	ax1.legend(bbox_to_anchor=[0, 1.2, 1.78, 0.5], loc=3, ncol=2, mode='expand', borderaxespad=0., fontsize=10)  
	plt.savefig('./Plots/paperplots_co2_uvdoses.eps', orientation='portrait',papertype='letter', format='eps')
	
	
	
	#####Now, plot stressors normalized by eustressors. In this limit, values greater than 1 mean that the relative photoreaction balance has shifted toward the eustressor.
	#Can't compare between 
	#Initialize plot basics
	fig2=plt.figure(figsize=(cm2inch(16.5),7))
	gs=gridspec.GridSpec(2,2, hspace=0.40,wspace=0.35, width_ratios=[2,1], top=.77, bottom=.1, left=.1, right=.95)
	ax1=plt.subplot(gs[0])
	ax2=plt.subplot(gs[2])
	ax3=plt.subplot(gs[1])
	ax4=plt.subplot(gs[3])
	colorseq1=iter(cm.rainbow(np.linspace(0,1,6)))
	colorseq2=iter(cm.rainbow(np.linspace(0,1,6)))
	colorseq3=iter(cm.rainbow(np.linspace(0,1,6)))
	colorseq4=iter(cm.rainbow(np.linspace(0,1,6)))
	
	#Plot max case	
	ax1.plot(N_CO2s[max_inds], umpgly_193s_max_normed/cucn3_254s_max_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq1), label=r'UMP-193/CuCN3-254')
	ax1.plot(N_CO2s[max_inds], umpgly_230s_max_normed/cucn3_254s_max_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq1), label=r'UMP-230/CuCN3-254')
	ax1.plot(N_CO2s[max_inds], umpgly_254s_max_normed/cucn3_254s_max_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq1), label=r'UMP-254/CuCN3-254')
	ax1.plot(N_CO2s[max_inds], umpgly_193s_max_normed/cucn3_300s_max_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq1), label=r'UMP-193/CuCN3-300')
	ax1.plot(N_CO2s[max_inds], umpgly_230s_max_normed/cucn3_300s_max_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq1), label=r'UMP-230/CuCN3-300')
	ax1.plot(N_CO2s[max_inds], umpgly_254s_max_normed/cucn3_300s_max_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq1), label=r'UMP-254/CuCN3-300')
	ax1.set_title('SZA=0, Albedo=New Snow')
	ax1.axvline(N_CO2s[max_inds[-7]], color='black', linewidth=2) #Mark the Rugheimer fiducial value
	ax1.axhline(1., color='black', linewidth=1) #Mark the Rugheimer fiducial value
	ax1.axvline(N_CO2s[max_inds[wordsworthind]], color='black', linewidth=1, linestyle='--') #Wordsworth lower limit
	ax1.axvline(N_CO2s[max_inds[kastingupperind]], color='black', linewidth=1, linestyle='--') #Kasting upper limit
	#Plot min case
	ax2.set_title('SZA=66.5, Albedo=Tundra')
	ax2.plot(N_CO2s[min_inds], umpgly_193s_min_normed/cucn3_254s_min_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq2), label=r'UMP-193/CuCN3-254')
	ax2.plot(N_CO2s[min_inds], umpgly_230s_min_normed/cucn3_254s_min_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq2), label=r'UMP-193/CuCN3-254')
	ax2.plot(N_CO2s[min_inds], umpgly_254s_min_normed/cucn3_254s_min_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq2), label=r'UMP-193/CuCN3-254')
	ax2.plot(N_CO2s[min_inds], umpgly_193s_min_normed/cucn3_300s_min_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq2), label=r'UMP-193/CuCN3-300')
	ax2.plot(N_CO2s[min_inds], umpgly_230s_min_normed/cucn3_300s_min_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq2), label=r'UMP-193/CuCN3-300')
	ax2.plot(N_CO2s[min_inds], umpgly_254s_min_normed/cucn3_300s_min_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq2), label=r'UMP-193/CuCN3-300')
	ax2.axvline(N_CO2s[min_inds[-7]], color='black', linewidth=2) #Mark the Rugheimer fiducial value
	ax2.axhline(1., color='black', linewidth=1) #Mark the Rugheimer fiducial value
	ax2.axvline(N_CO2s[min_inds[wordsworthind]], color='black', linewidth=1, linestyle='--') #Wordsworth lower limit
	ax2.axvline(N_CO2s[min_inds[kastingupperind]], color='black', linewidth=1, linestyle='--') #Kasting upper limit
	
	#Plot max case	
	ax3.plot(N_CO2s[max_inds], umpgly_193s_max_normed/cucn3_254s_max_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq3), label=r'UMP-193/CuCN3-254')
	ax3.plot(N_CO2s[max_inds], umpgly_230s_max_normed/cucn3_254s_max_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq3), label=r'UMP-193/CuCN3-254')
	ax3.plot(N_CO2s[max_inds], umpgly_254s_max_normed/cucn3_254s_max_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq3), label=r'UMP-193/CuCN3-254')
	ax3.plot(N_CO2s[max_inds], umpgly_193s_max_normed/cucn3_300s_max_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq3), label=r'UMP-193/CuCN3-300')
	ax3.plot(N_CO2s[max_inds], umpgly_230s_max_normed/cucn3_300s_max_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq3), label=r'UMP-193/CuCN3-300')
	ax3.plot(N_CO2s[max_inds], umpgly_254s_max_normed/cucn3_300s_max_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq3), label=r'UMP-193/CuCN3-300')
	ax3.axvline(N_CO2s[max_inds[-7]], color='black', linewidth=2) #Mark the Rugheimer fiducial value
	ax3.axhline(1., color='black', linewidth=1) #Mark the Rugheimer fiducial value
	ax3.axvline(N_CO2s[max_inds[wordsworthind]], color='black', linewidth=1, linestyle='--') #Wordsworth lower limit
	ax3.axvline(N_CO2s[max_inds[kastingupperind]], color='black', linewidth=1, linestyle='--') #Kasting upper limit
	
	#Plot min case
	ax4.plot(N_CO2s[min_inds], umpgly_193s_min_normed/cucn3_254s_min_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq4), label=r'UMP-193/CuCN3-254')
	ax4.plot(N_CO2s[min_inds], umpgly_230s_min_normed/cucn3_254s_min_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq4), label=r'UMP-193/CuCN3-254')
	ax4.plot(N_CO2s[min_inds], umpgly_254s_min_normed/cucn3_254s_min_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq4), label=r'UMP-193/CuCN3-254')
	ax4.plot(N_CO2s[min_inds], umpgly_193s_min_normed/cucn3_300s_min_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq4), label=r'UMP-193/CuCN3-300')
	ax4.plot(N_CO2s[min_inds], umpgly_230s_min_normed/cucn3_300s_min_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq4), label=r'UMP-193/CuCN3-300')
	ax4.plot(N_CO2s[min_inds], umpgly_254s_min_normed/cucn3_300s_min_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq4), label=r'UMP-193/CuCN3-300')
	ax4.axvline(N_CO2s[min_inds[-7]], color='black', linewidth=2) #Mark the Rugheimer fiducial value	
	ax4.axhline(1., color='black', linewidth=1) #Mark the Rugheimer fiducial value
	ax4.axvline(N_CO2s[min_inds[wordsworthind]], color='black', linewidth=1, linestyle='--') #Wordsworth lower limit
	ax4.axvline(N_CO2s[min_inds[kastingupperind]], color='black', linewidth=1, linestyle='--') #Kasting upper limit
	
	#Finalize plot details.
	ax1.set_ylabel(r'Relative Dose Rate $D$')
	ax2.set_ylabel(r'Relative Dose Rate $D$')
	ax3.set_ylabel(r'Relative Dose Rate $D$')
	ax4.set_ylabel(r'Relative Dose Rate $D$')
	ax1.set_xlabel(r'N$_{CO2}$ (cm$^{-2}$)')
	ax2.set_xlabel(r'N$_{CO2}$ (cm$^{-2}$)')
	ax3.set_xlabel(r'N$_{CO2}$ (cm$^{-2}$)')
	ax4.set_xlabel(r'N$_{CO2}$ (cm$^{-2}$)')
	ax1.set_yscale('linear')
	ax2.set_yscale('linear')
	ax3.set_yscale('linear')
	ax4.set_yscale('linear')
	ax1.set_xscale('log')
	ax2.set_xscale('log')
	ax3.set_xscale('log')
	ax4.set_xscale('log')
	ax1.set_xlim([2.09e18, 2.09e24])
	ax2.set_xlim([2.09e18, 2.09e24])
	ax3.set_xlim([2.09e24, 9.85e26])
	ax4.set_xlim([2.09e24, 9.85e26])
	ax1.set_ylim([0.90,1.8])
	ax2.set_ylim([0.95,1.3])
	ax3.set_ylim([0.5, 1.8])
	ax4.set_ylim([0.8, 1.2])
	
	ax1.legend(bbox_to_anchor=[0, 1.2, 1.78, 0.5], loc=3, ncol=2, mode='expand', borderaxespad=0., fontsize=10)  
	plt.savefig('./Plots/paperplots_co2_uvdoses_norm.eps', orientation='portrait',papertype='letter', format='eps')
	

if plot_dosimeters_ch4:
	###########First, import the CO2 0.1 bar base case for comparison
	SZAs, albedos, N_CO2s, rad100_165s, rad200_300s, umpgly_193s, umpgly_230s, umpgly_254s,cucn3_254s, cucn3_300s=np.genfromtxt('./Doses/co2_uv_doses.dat', skip_header=2, skip_footer=0, usecols=(0, 1, 2, 3, 4, 5, 6, 7, 8, 9), unpack=True, delimiter='	') #wavelength in nm, relative efficiency unitless
	
	basecaseind=9 #index of 0.1 bar CO2 case (NCO2=2.09e24 cm**-2), max radiance (A=new snow, sza=0) case
	
	umpgly_193s_base=umpgly_193s[basecaseind]
	umpgly_230s_base=umpgly_230s[basecaseind]
	umpgly_254s_base=umpgly_254s[basecaseind]
	cucn3_254s_base=cucn3_254s[basecaseind]
	cucn3_300s_base=cucn3_300s[basecaseind]	
	
	###########Next, set up information about CH4 file
	gas='ch4'
	gaslabel='CH4'
	minind=3 #index of minimum plausible CH4 concentration
	fiducialind=5 #index of the fiducial Rugheimer gas concentration
	maxind=8# index of maximum plausible CH4 concentration
	
	###########Next, import CH4, and normalize by the base case.

	N_gas, rad100_165s, rad200_300s, umpgly_193s, umpgly_230s, umpgly_254s,cucn3_254s, cucn3_300s=np.genfromtxt('./Doses/'+gas+'_uv_doses.dat', skip_header=2, skip_footer=0, usecols=(1, 2, 3, 4, 5, 6, 7, 8), unpack=True, delimiter=' & ') #wavelength in nm, relative efficiency unitless

	
	#values<1: the reaction proceeds slower than in the base case. Values>1: the reaction proceeds faster than in the base case.
	umpgly_193s_normed=umpgly_193s/umpgly_193s_base
	umpgly_230s_normed=umpgly_230s/umpgly_230s_base
	umpgly_254s_normed=umpgly_254s/umpgly_254s_base
	cucn3_254s_normed=cucn3_254s/cucn3_254s_base
	cucn3_300s_normed=cucn3_300s/cucn3_300s_base
	
	###########Now, plot the dosimeters vs CO2 concentration
	
	#Initialize plot basics
	fig, ax1=plt.subplots(1, figsize=(cm2inch(16.5),7))
	colorseq1=iter(cm.rainbow(np.linspace(0,1,5)))
	
	#Plot max case	
	ax1.plot(N_gas, umpgly_193s_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq1), label=r'UMP Gly Bond Cleavage ($\lambda_0=193$)')
	ax1.plot(N_gas, umpgly_230s_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq1), label=r'UMP Gly Bond Cleavage ($\lambda_0=230$)')
	ax1.plot(N_gas, umpgly_254s_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq1), label=r'UMP Gly Bond Cleavage ($\lambda_0=254$)')
	ax1.plot(N_gas, cucn3_254s_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq1), label=r'CuCN$_3$ Photoionization ($\lambda_0=254$)')
	ax1.plot(N_gas, cucn3_300s_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq1), label=r'CuCN$_3$ Photoionization ($\lambda_0=300$)')
	ax1.set_title(gaslabel+' (SZA=0, Albedo=New Snow)')
	ax1.axvline(N_gas[fiducialind], color='black', linewidth=2) #Mark the Rugheimer fiducial value
	ax1.axhline(1., color='black', linewidth=1) #Values above this: faster rxn than under 0.1 bar CO2
	
	#Finalize plot details.
	ax1.set_ylabel(r'Relative Dose Rate $D/D_{NCO2=2.09e24 cm^{-2}}$')
	ax1.set_xlabel(r'N$_{'+gaslabel+'}$ (cm$^{-2}$)')
	ax1.set_yscale('linear')
	ax1.set_xscale('log')
	ax1.set_xlim([N_gas[minind],N_gas[maxind]])
	ax1.set_ylim([0.95,1.8])

	plt.tight_layout(rect=(0,0,1., 0.7))
	ax1.legend(bbox_to_anchor=[0, 1.2, 1., 0.5], loc=3, ncol=2, mode='expand', borderaxespad=0., fontsize=10)  
	plt.savefig('./Plots/paperplots_'+gas+'_uvdoses.eps', orientation='portrait',papertype='letter', format='eps')

if plot_dosimeters_h2o:
	###########First, import the CO2 0.1 bar base case for comparison
	SZAs, albedos, N_CO2s, rad100_165s, rad200_300s, umpgly_193s, umpgly_230s, umpgly_254s,cucn3_254s, cucn3_300s=np.genfromtxt('./Doses/co2_uv_doses.dat', skip_header=2, skip_footer=0, usecols=(0, 1, 2, 3, 4, 5, 6, 7, 8, 9), unpack=True, delimiter='	') #wavelength in nm, relative efficiency unitless
	
	basecaseind=9 #index of 0.1 bar CO2 case (NCO2=2.09e24 cm**-2), max radiance (A=new snow, sza=0) case
	
	umpgly_193s_base=umpgly_193s[basecaseind]
	umpgly_230s_base=umpgly_230s[basecaseind]
	umpgly_254s_base=umpgly_254s[basecaseind]
	cucn3_254s_base=cucn3_254s[basecaseind]
	cucn3_300s_base=cucn3_300s[basecaseind]
	
	
	###########Next, set up information about H2O file
	gas='h2o'
	gaslabel='H2O'

	minind=0 #index of minimum plausible H2O concentration	
	fiducialind=5 #index of the fiducial H2O vapor concentration
	maxind=8# index of maximum plausible H2O concentration
	###########Next, import H2O, and normalize by the base case.

	N_gas, rad100_165s, rad200_300s, umpgly_193s, umpgly_230s, umpgly_254s,cucn3_254s, cucn3_300s=np.genfromtxt('./Doses/'+gas+'_uv_doses.dat', skip_header=2, skip_footer=0, usecols=(1, 2, 3, 4, 5, 6, 7, 8), unpack=True, delimiter=' & ') #wavelength in nm, relative efficiency unitless

	
	#values<1: the reaction proceeds slower than in the base case. Values>1: the reaction proceeds faster than in the base case.
	umpgly_193s_normed=umpgly_193s/umpgly_193s_base
	umpgly_230s_normed=umpgly_230s/umpgly_230s_base
	umpgly_254s_normed=umpgly_254s/umpgly_254s_base
	cucn3_254s_normed=cucn3_254s/cucn3_254s_base
	cucn3_300s_normed=cucn3_300s/cucn3_300s_base
	
	###########Now, plot the dosimeters vs CO2 concentration
	
	#Initialize plot basics
	fig, ax1=plt.subplots(1, figsize=(cm2inch(16.5),7))
	
	colorseq1=iter(cm.rainbow(np.linspace(0,1,5)))
	
	#Plot linear half	
	ax1.plot(N_gas, umpgly_193s_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq1), label=r'UMP Gly Bond Cleavage ($\lambda_0=193$)')
	ax1.plot(N_gas, umpgly_230s_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq1), label=r'UMP Gly Bond Cleavage ($\lambda_0=230$)')
	ax1.plot(N_gas, umpgly_254s_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq1), label=r'UMP Gly Bond Cleavage ($\lambda_0=254$)')
	ax1.plot(N_gas, cucn3_254s_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq1), label=r'CuCN$_3$ Photoionization ($\lambda_0=254$)')
	ax1.plot(N_gas, cucn3_300s_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq1), label=r'CuCN$_3$ Photoionization ($\lambda_0=300$)')
	ax1.set_title(gaslabel+' (SZA=0, Albedo=New Snow)')
	ax1.axvline(N_gas[fiducialind], color='black', linewidth=2) #Mark the Rugheimer fiducial value
	ax1.axhline(1., color='black', linewidth=1) #Values above this: faster rxn than under 0.1 bar CO2
	
	#Finalize plot details.
	ax1.set_ylabel(r'Relative Dose Rate $D/D_{NCO2=2.09e24 cm^{-2}}$')
	ax1.set_xlabel(r'N$_{'+gaslabel+'}$ (cm$^{-2}$)')
	ax1.set_yscale('linear')
	ax1.set_xscale('log')
	ax1.set_xlim([N_gas[minind],N_gas[maxind]])
	ax1.set_ylim([0.7,1.7])

	plt.tight_layout(rect=(0,0,1., 0.7))
	ax1.legend(bbox_to_anchor=[0, 1.2, 1., 0.5], loc=3, ncol=2, mode='expand', borderaxespad=0., fontsize=10)  
	plt.savefig('./Plots/paperplots_'+gas+'_uvdoses.eps', orientation='portrait',papertype='letter', format='eps')

if plot_dosimeters_so2:
	###########First, import the CO2 0.1 bar base case for comparison
	SZAs, albedos, N_CO2s, rad100_165s, rad200_300s, umpgly_193s, umpgly_230s, umpgly_254s,cucn3_254s, cucn3_300s=np.genfromtxt('./Doses/co2_uv_doses.dat', skip_header=2, skip_footer=0, usecols=(0, 1, 2, 3, 4, 5, 6, 7, 8, 9), unpack=True, delimiter='	') #wavelength in nm, relative efficiency unitless
	
	basecaseind=9 #index of 0.1 bar CO2 case (NCO2=2.09e24 cm**-2), max radiance (A=new snow, sza=0) case
	
	umpgly_193s_base=umpgly_193s[basecaseind]
	umpgly_230s_base=umpgly_230s[basecaseind]
	umpgly_254s_base=umpgly_254s[basecaseind]
	cucn3_254s_base=cucn3_254s[basecaseind]
	cucn3_300s_base=cucn3_300s[basecaseind]
	
	
	###########Next, set up information about H2O file
	gas='so2'
	gaslabel='SO2'

	minind=0 #index of minimum plausible H2O concentration	
	fiducialind=5 #index of the fiducial H2O vapor concentration
	maxind=12# index of maximum plausible H2O concentration
	breakind=9 #index of where to break the plots
	###########Next, import H2O, and normalize by the base case.

	N_gas, rad100_165s, rad200_300s, umpgly_193s, umpgly_230s, umpgly_254s,cucn3_254s, cucn3_300s=np.genfromtxt('./Doses/'+gas+'_uv_doses.dat', skip_header=2, skip_footer=0, usecols=(1, 2, 3, 4, 5, 6, 7, 8), unpack=True, delimiter=' & ') #wavelength in nm, relative efficiency unitless

	
	#values<1: the reaction proceeds slower than in the base case. Values>1: the reaction proceeds faster than in the base case.
	umpgly_193s_normed=umpgly_193s/umpgly_193s_base
	umpgly_230s_normed=umpgly_230s/umpgly_230s_base
	umpgly_254s_normed=umpgly_254s/umpgly_254s_base
	cucn3_254s_normed=cucn3_254s/cucn3_254s_base
	cucn3_300s_normed=cucn3_300s/cucn3_300s_base
	
	###########Now, plot the dosimeters vs CO2 concentration
	
	#Initialize plot basics
	fig2=plt.figure(figsize=(cm2inch(16.5),7))
	gs=gridspec.GridSpec(1,2, hspace=0.40,wspace=0.35, width_ratios=[breakind-minind,maxind-breakind], top=.70, bottom=.1, left=.1, right=.95)
	ax1=plt.subplot(gs[0])
	ax2=plt.subplot(gs[1])

	colorseq1=iter(cm.rainbow(np.linspace(0,1,5)))
	colorseq2=iter(cm.rainbow(np.linspace(0,1,5)))
	
	#Plot linear half	
	ax1.plot(N_gas, umpgly_193s_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq1), label=r'UMP Gly Bond Cleavage ($\lambda_0=193$)')
	ax1.plot(N_gas, umpgly_230s_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq1), label=r'UMP Gly Bond Cleavage ($\lambda_0=230$)')
	ax1.plot(N_gas, umpgly_254s_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq1), label=r'UMP Gly Bond Cleavage ($\lambda_0=254$)')
	ax1.plot(N_gas, cucn3_254s_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq1), label=r'CuCN$_3$ Photoionization ($\lambda_0=254$)')
	ax1.plot(N_gas, cucn3_300s_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq1), label=r'CuCN$_3$ Photoionization ($\lambda_0=300$)')
	ax1.set_title(gaslabel+' (SZA=0, Albedo=New Snow)')
	ax1.axvline(N_gas[fiducialind], color='black', linewidth=2) #Mark the Rugheimer fiducial value
	ax1.axhline(1., color='black', linewidth=1) #Values above this: faster rxn than under 0.1 bar CO2

	#Plot log half	
	ax2.plot(N_gas, umpgly_193s_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq2), label=r'UMP Gly Bond Cleavage ($\lambda_0=193$)')
	ax2.plot(N_gas, umpgly_230s_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq2), label=r'UMP Gly Bond Cleavage ($\lambda_0=230$)')
	ax2.plot(N_gas, umpgly_254s_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq2), label=r'UMP Gly Bond Cleavage ($\lambda_0=254$)')
	ax2.plot(N_gas, cucn3_254s_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq2), label=r'CuCN$_3$ Photoionization ($\lambda_0=254$)')
	ax2.plot(N_gas, cucn3_300s_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq2), label=r'CuCN$_3$ Photoionization ($\lambda_0=300$)')
	#ax2.set_title(gaslabel+' (SZA=0, Albedo=New Snow)')
	ax2.axvline(N_gas[fiducialind], color='black', linewidth=2) #Mark the Rugheimer fiducial value
	ax2.axhline(1., color='black', linewidth=1) #Values above this: faster rxn than under 0.1 bar CO2
	
	#Finalize plot details.
	ax1.set_ylabel(r'Relative Dose Rate $D/D_{NCO2=2.09e24 cm^{-2}}$')
	ax1.set_xlabel(r'N$_{'+gaslabel+'}$ (cm$^{-2}$)')
	ax1.set_yscale('linear')
	ax1.set_xscale('log')
	ax1.set_xlim([N_gas[minind],N_gas[breakind]])
	#ax1.set_ylim([0.4,1.8])

	ax2.set_ylabel(r'Relative Dose Rate $D/D_{NCO2=2.09e24 cm^{-2}}$')
	ax2.set_xlabel(r'N$_{'+gaslabel+'}$ (cm$^{-2}$)')
	ax2.set_yscale('log')
	ax2.set_xscale('log')
	ax2.set_xlim([N_gas[breakind],N_gas[maxind]])
	#ax2.set_ylim([0.4,1.8])
	#plt.tight_layout(rect=(0,0,1., 0.7))
	ax1.legend(bbox_to_anchor=[0, 1.2, 1.6, 0.5], loc=3, ncol=2, mode='expand', borderaxespad=0., fontsize=10)  
	plt.savefig('./Plots/paperplots_'+gas+'_uvdoses.eps', orientation='portrait',papertype='letter', format='eps')
	
	###########Now, plot the good/bad balance as a function of changing column density
	
	#Initialize plot basics
	fig3=plt.figure(figsize=(cm2inch(16.5),7))
	gs=gridspec.GridSpec(1,2, hspace=0.40,wspace=0.35, width_ratios=[breakind-minind,maxind-breakind], top=.70, bottom=.1, left=.1, right=.95)
	ax1=plt.subplot(gs[0])
	ax2=plt.subplot(gs[1])

	colorseq1=iter(cm.rainbow(np.linspace(0,1,6)))
	colorseq2=iter(cm.rainbow(np.linspace(0,1,6)))
	
	#Plot linear half	
	ax1.plot(N_gas, umpgly_193s_normed/cucn3_254s_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq1), label=r'UMP-193/CuCN3-254')
	ax1.plot(N_gas, umpgly_230s_normed/cucn3_254s_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq1), label=r'UMP-230/CuCN3-254')
	ax1.plot(N_gas, umpgly_254s_normed/cucn3_254s_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq1), label=r'UMP-254/CuCN3-254')
	ax1.plot(N_gas, umpgly_193s_normed/cucn3_300s_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq1), label=r'UMP-193/CuCN3-300')
	ax1.plot(N_gas, umpgly_230s_normed/cucn3_300s_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq1), label=r'UMP-230/CuCN3-300')
	ax1.plot(N_gas, umpgly_254s_normed/cucn3_300s_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq1), label=r'UMP-254/CuCN3-300')
	ax1.set_title(gaslabel+' (SZA=0, Albedo=New Snow)')
	ax1.axvline(N_gas[fiducialind], color='black', linewidth=2) #Mark the Rugheimer fiducial value
	ax1.axhline(1., color='black', linewidth=1) #Values above this: faster rxn than under 0.1 bar CO2

	#Plot log half	
	ax2.plot(N_gas, umpgly_193s_normed/cucn3_254s_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq2), label=r'UMP-193/CuCN3-254')
	ax2.plot(N_gas, umpgly_230s_normed/cucn3_254s_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq2), label=r'UMP-230/CuCN3-254')
	ax2.plot(N_gas, umpgly_254s_normed/cucn3_254s_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq2), label=r'UMP-254/CuCN3-254')
	ax2.plot(N_gas, umpgly_193s_normed/cucn3_300s_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq2), label=r'UMP-193/CuCN3-300')
	ax2.plot(N_gas, umpgly_230s_normed/cucn3_300s_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq2), label=r'UMP-230/CuCN3-300')
	ax2.plot(N_gas, umpgly_254s_normed/cucn3_300s_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq2), label=r'UMP-254/CuCN3-300')
	ax2.axvline(N_gas[fiducialind], color='black', linewidth=2) #Mark the Rugheimer fiducial value
	ax2.axhline(1., color='black', linewidth=1) #Values above this: faster rxn than under 0.1 bar CO2
	
	#Finalize plot details.
	ax1.set_ylabel(r'Relative Dose Rate $D^{UMP-X}/D^{CuCN3-Y}$')
	ax1.set_xlabel(r'N$_{'+gaslabel+'}$ (cm$^{-2}$)')
	ax1.set_yscale('linear')
	ax1.set_xscale('log')
	ax1.set_xlim([N_gas[minind],N_gas[breakind]])
	#ax1.set_ylim([0.4,1.8])

	ax2.set_ylabel(r'Relative Dose Rate $D^{UMP-X}/D^{CuCN3-Y}$')
	ax2.set_xlabel(r'N$_{'+gaslabel+'}$ (cm$^{-2}$)')
	ax2.set_yscale('log')
	ax2.set_xscale('log')
	ax2.set_xlim([N_gas[breakind],N_gas[maxind]])
	#ax2.set_ylim([0.4,1.8])
	#plt.tight_layout(rect=(0,0,1., 0.7))
	ax1.legend(bbox_to_anchor=[0, 1.2, 1.6, 0.5], loc=3, ncol=2, mode='expand', borderaxespad=0., fontsize=10)  
	plt.savefig('./Plots/paperplots_'+gas+'_uvdoses_balance.eps', orientation='portrait',papertype='letter', format='eps')
	
if plot_dosimeters_h2s:
	###########First, import the CO2 0.1 bar base case for comparison
	SZAs, albedos, N_CO2s, rad100_165s, rad200_300s, umpgly_193s, umpgly_230s, umpgly_254s,cucn3_254s, cucn3_300s=np.genfromtxt('./Doses/co2_uv_doses.dat', skip_header=2, skip_footer=0, usecols=(0, 1, 2, 3, 4, 5, 6, 7, 8, 9), unpack=True, delimiter='	') #wavelength in nm, relative efficiency unitless
	
	basecaseind=9 #index of 0.1 bar CO2 case (NCO2=2.09e24 cm**-2), max radiance (A=new snow, sza=0) case
	
	umpgly_193s_base=umpgly_193s[basecaseind]
	umpgly_230s_base=umpgly_230s[basecaseind]
	umpgly_254s_base=umpgly_254s[basecaseind]
	cucn3_254s_base=cucn3_254s[basecaseind]
	cucn3_300s_base=cucn3_300s[basecaseind]
	
	
	###########Next, set up information about H2O file
	gas='h2s'
	gaslabel='H2S'

	minind=0 #index of minimum plausible H2O concentration	
	fiducialind=5 #index of the fiducial H2O vapor concentration
	maxind=12# index of maximum plausible H2O concentration
	breakind=9 #index of where to break the plots
	breakind2=8 #index of where to break the second set of plots

	###########Next, import H2O, and normalize by the base case.

	N_gas, rad100_165s, rad200_300s, umpgly_193s, umpgly_230s, umpgly_254s,cucn3_254s, cucn3_300s=np.genfromtxt('./Doses/'+gas+'_uv_doses.dat', skip_header=2, skip_footer=0, usecols=(1, 2, 3, 4, 5, 6, 7, 8), unpack=True, delimiter=' & ') #wavelength in nm, relative efficiency unitless

	
	#values<1: the reaction proceeds slower than in the base case. Values>1: the reaction proceeds faster than in the base case.
	umpgly_193s_normed=umpgly_193s/umpgly_193s_base
	umpgly_230s_normed=umpgly_230s/umpgly_230s_base
	umpgly_254s_normed=umpgly_254s/umpgly_254s_base
	cucn3_254s_normed=cucn3_254s/cucn3_254s_base
	cucn3_300s_normed=cucn3_300s/cucn3_300s_base
	
	###########Now, plot the dosimeters vs CO2 concentration
	
	#Initialize plot basics
	fig2=plt.figure(figsize=(cm2inch(16.5),7))
	gs=gridspec.GridSpec(1,2, hspace=0.40,wspace=0.35, width_ratios=[breakind-minind,maxind-breakind], top=.70, bottom=.1, left=.1, right=.95)
	ax1=plt.subplot(gs[0])
	ax2=plt.subplot(gs[1])

	colorseq1=iter(cm.rainbow(np.linspace(0,1,5)))
	colorseq2=iter(cm.rainbow(np.linspace(0,1,5)))
	
	#Plot linear half	
	ax1.plot(N_gas, umpgly_193s_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq1), label=r'UMP Gly Bond Cleavage ($\lambda_0=193$)')
	ax1.plot(N_gas, umpgly_230s_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq1), label=r'UMP Gly Bond Cleavage ($\lambda_0=230$)')
	ax1.plot(N_gas, umpgly_254s_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq1), label=r'UMP Gly Bond Cleavage ($\lambda_0=254$)')
	ax1.plot(N_gas, cucn3_254s_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq1), label=r'CuCN$_3$ Photoionization ($\lambda_0=254$)')
	ax1.plot(N_gas, cucn3_300s_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq1), label=r'CuCN$_3$ Photoionization ($\lambda_0=300$)')
	ax1.set_title(gaslabel+' (SZA=0, Albedo=New Snow)')
	ax1.axvline(N_gas[fiducialind], color='black', linewidth=2) #Mark the Rugheimer fiducial value
	ax1.axhline(1., color='black', linewidth=1) #Values above this: faster rxn than under 0.1 bar CO2

	#Plot log half	
	ax2.plot(N_gas, umpgly_193s_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq2), label=r'UMP Gly Bond Cleavage ($\lambda_0=193$)')
	ax2.plot(N_gas, umpgly_230s_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq2), label=r'UMP Gly Bond Cleavage ($\lambda_0=230$)')
	ax2.plot(N_gas, umpgly_254s_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq2), label=r'UMP Gly Bond Cleavage ($\lambda_0=254$)')
	ax2.plot(N_gas, cucn3_254s_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq2), label=r'CuCN$_3$ Photoionization ($\lambda_0=254$)')
	ax2.plot(N_gas, cucn3_300s_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq2), label=r'CuCN$_3$ Photoionization ($\lambda_0=300$)')
	#ax2.set_title(gaslabel+' (SZA=0, Albedo=New Snow)')
	ax2.axvline(N_gas[fiducialind], color='black', linewidth=2) #Mark the Rugheimer fiducial value
	ax2.axhline(1., color='black', linewidth=1) #Values above this: faster rxn than under 0.1 bar CO2
	
	#Finalize plot details.
	ax1.set_ylabel(r'Relative Dose Rate $D/D_{NCO2=2.09e24 cm^{-2}}$')
	ax1.set_xlabel(r'N$_{'+gaslabel+'}$ (cm$^{-2}$)')
	ax1.set_yscale('linear')
	ax1.set_xscale('log')
	ax1.set_xlim([N_gas[minind],N_gas[breakind]])
	#ax1.set_ylim([0.4,1.8])

	ax2.set_ylabel(r'Relative Dose Rate $D/D_{NCO2=2.09e24 cm^{-2}}$')
	ax2.set_xlabel(r'N$_{'+gaslabel+'}$ (cm$^{-2}$)')
	ax2.set_yscale('log')
	ax2.set_xscale('log')
	ax2.set_xlim([N_gas[breakind],N_gas[maxind]])
	#ax2.set_ylim([0.4,1.8])
	#plt.tight_layout(rect=(0,0,1., 0.7))
	ax1.legend(bbox_to_anchor=[0, 1.2, 1.6, 0.5], loc=3, ncol=2, mode='expand', borderaxespad=0., fontsize=10)  
	plt.savefig('./Plots/paperplots_'+gas+'_uvdoses.eps', orientation='portrait',papertype='letter', format='eps')
	###########Now, plot the good/bad balance as a function of changing column density
	
	#Initialize plot basics
	fig3=plt.figure(figsize=(cm2inch(16.5),7))
	gs=gridspec.GridSpec(1,2, hspace=0.40,wspace=0.35, width_ratios=[breakind2-minind,maxind-breakind2], top=.70, bottom=.1, left=.1, right=.95)
	ax1=plt.subplot(gs[0])
	ax2=plt.subplot(gs[1])

	colorseq1=iter(cm.rainbow(np.linspace(0,1,6)))
	colorseq2=iter(cm.rainbow(np.linspace(0,1,6)))
	
	#Plot linear half	
	ax1.plot(N_gas, umpgly_193s_normed/cucn3_254s_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq1), label=r'UMP-193/CuCN3-254')
	ax1.plot(N_gas, umpgly_230s_normed/cucn3_254s_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq1), label=r'UMP-230/CuCN3-254')
	ax1.plot(N_gas, umpgly_254s_normed/cucn3_254s_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq1), label=r'UMP-254/CuCN3-254')
	ax1.plot(N_gas, umpgly_193s_normed/cucn3_300s_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq1), label=r'UMP-193/CuCN3-300')
	ax1.plot(N_gas, umpgly_230s_normed/cucn3_300s_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq1), label=r'UMP-230/CuCN3-300')
	ax1.plot(N_gas, umpgly_254s_normed/cucn3_300s_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq1), label=r'UMP-254/CuCN3-300')
	ax1.set_title(gaslabel+' (SZA=0, Albedo=New Snow)')
	ax1.axvline(N_gas[fiducialind], color='black', linewidth=2) #Mark the Rugheimer fiducial value
	ax1.axhline(1., color='black', linewidth=1) #Values above this: faster rxn than under 0.1 bar CO2

	#Plot log half	
	ax2.plot(N_gas, umpgly_193s_normed/cucn3_254s_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq2), label=r'UMP-193/CuCN3-254')
	ax2.plot(N_gas, umpgly_230s_normed/cucn3_254s_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq2), label=r'UMP-230/CuCN3-254')
	ax2.plot(N_gas, umpgly_254s_normed/cucn3_254s_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq2), label=r'UMP-254/CuCN3-254')
	ax2.plot(N_gas, umpgly_193s_normed/cucn3_300s_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq2), label=r'UMP-193/CuCN3-300')
	ax2.plot(N_gas, umpgly_230s_normed/cucn3_300s_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq2), label=r'UMP-230/CuCN3-300')
	ax2.plot(N_gas, umpgly_254s_normed/cucn3_300s_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq2), label=r'UMP-254/CuCN3-300')
	ax2.axvline(N_gas[fiducialind], color='black', linewidth=2) #Mark the Rugheimer fiducial value
	ax2.axhline(1., color='black', linewidth=1) #Values above this: faster rxn than under 0.1 bar CO2
	
	#Finalize plot details.
	ax1.set_ylabel(r'Relative Dose Rate $D^{UMP-X}/D^{CuCN3-Y}$')
	ax1.set_xlabel(r'N$_{'+gaslabel+'}$ (cm$^{-2}$)')
	ax1.set_yscale('linear')
	ax1.set_xscale('log')
	ax1.set_xlim([N_gas[minind],N_gas[breakind2]])
	ax1.set_ylim([0.,4])

	ax2.set_ylabel(r'Relative Dose Rate $D^{UMP-X}/D^{CuCN3-Y}$')
	ax2.set_xlabel(r'N$_{'+gaslabel+'}$ (cm$^{-2}$)')
	ax2.set_yscale('log')
	ax2.set_xscale('log')
	ax2.set_xlim([N_gas[breakind2],N_gas[maxind]])
	#ax2.set_ylim([0.4,1.8])
	#plt.tight_layout(rect=(0,0,1., 0.7))
	ax1.legend(bbox_to_anchor=[0, 1.2, 1.7, 0.5], loc=3, ncol=2, mode='expand', borderaxespad=0., fontsize=10)  
	plt.savefig('./Plots/paperplots_'+gas+'_uvdoses_balance.eps', orientation='portrait',papertype='letter', format='eps')

if plot_dosimeters_o2:
	###########First, import the CO2 0.1 bar base case for comparison
	SZAs, albedos, N_CO2s, rad100_165s, rad200_300s, umpgly_193s, umpgly_230s, umpgly_254s,cucn3_254s, cucn3_300s=np.genfromtxt('./Doses/co2_uv_doses.dat', skip_header=2, skip_footer=0, usecols=(0, 1, 2, 3, 4, 5, 6, 7, 8, 9), unpack=True, delimiter='	') #wavelength in nm, relative efficiency unitless
	
	basecaseind=9 #index of 0.1 bar CO2 case (NCO2=2.09e24 cm**-2), max radiance (A=new snow, sza=0) case
	
	umpgly_193s_base=umpgly_193s[basecaseind]
	umpgly_230s_base=umpgly_230s[basecaseind]
	umpgly_254s_base=umpgly_254s[basecaseind]
	cucn3_254s_base=cucn3_254s[basecaseind]
	cucn3_300s_base=cucn3_300s[basecaseind]
	
	
	###########Next, set up information about H2O file
	gas='o2'
	gaslabel='O2'

	minind=0 #index of minimum plausible H2O concentration	
	fiducialind=5 #index of the fiducial H2O vapor concentration
	maxind=10# index of maximum plausible H2O concentration
	###########Next, import H2O, and normalize by the base case.

	N_gas, rad100_165s, rad200_300s, umpgly_193s, umpgly_230s, umpgly_254s,cucn3_254s, cucn3_300s=np.genfromtxt('./Doses/'+gas+'_uv_doses.dat', skip_header=2, skip_footer=0, usecols=(1, 2, 3, 4, 5, 6, 7, 8), unpack=True, delimiter=' & ') #wavelength in nm, relative efficiency unitless

	
	#values<1: the reaction proceeds slower than in the base case. Values>1: the reaction proceeds faster than in the base case.
	umpgly_193s_normed=umpgly_193s/umpgly_193s_base
	umpgly_230s_normed=umpgly_230s/umpgly_230s_base
	umpgly_254s_normed=umpgly_254s/umpgly_254s_base
	cucn3_254s_normed=cucn3_254s/cucn3_254s_base
	cucn3_300s_normed=cucn3_300s/cucn3_300s_base
	
	###########Now, plot the dosimeters vs CO2 concentration
	
	#Initialize plot basics
	fig, ax1=plt.subplots(1, figsize=(cm2inch(16.5),7))

	colorseq1=iter(cm.rainbow(np.linspace(0,1,5)))
	
	#Plot linear half	
	ax1.plot(N_gas, umpgly_193s_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq1), label=r'UMP Gly Bond Cleavage ($\lambda_0=193$)')
	ax1.plot(N_gas, umpgly_230s_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq1), label=r'UMP Gly Bond Cleavage ($\lambda_0=230$)')
	ax1.plot(N_gas, umpgly_254s_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq1), label=r'UMP Gly Bond Cleavage ($\lambda_0=254$)')
	ax1.plot(N_gas, cucn3_254s_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq1), label=r'CuCN$_3$ Photoionization ($\lambda_0=254$)')
	ax1.plot(N_gas, cucn3_300s_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq1), label=r'CuCN$_3$ Photoionization ($\lambda_0=300$)')
	ax1.set_title(gaslabel+' (SZA=0, Albedo=New Snow)')
	ax1.axvline(N_gas[fiducialind], color='black', linewidth=2) #Mark the Rugheimer fiducial value
	ax1.axhline(1., color='black', linewidth=1) #Values above this: faster rxn than under 0.1 bar CO2

	#Finalize plot details.
	ax1.set_ylabel(r'Relative Dose Rate $D/D_{NCO2=2.09e24 cm^{-2}}$')
	ax1.set_xlabel(r'N$_{'+gaslabel+'}$ (cm$^{-2}$)')
	ax1.set_yscale('linear')
	ax1.set_xscale('log')
	ax1.set_xlim([N_gas[minind],N_gas[maxind]])
	ax1.set_ylim([0.,1.8])

	plt.tight_layout(rect=(0,0,1., 0.80))
	ax1.legend(bbox_to_anchor=[0, 1.1, 1., 0.5], loc=3, ncol=2, mode='expand', borderaxespad=0., fontsize=10)  
	plt.savefig('./Plots/paperplots_'+gas+'_uvdoses.eps', orientation='portrait',papertype='letter', format='eps')

if plot_dosimeters_o3:
	###########First, import the CO2 0.1 bar base case for comparison
	SZAs, albedos, N_CO2s, rad100_165s, rad200_300s, umpgly_193s, umpgly_230s, umpgly_254s,cucn3_254s, cucn3_300s=np.genfromtxt('./Doses/co2_uv_doses.dat', skip_header=2, skip_footer=0, usecols=(0, 1, 2, 3, 4, 5, 6, 7, 8, 9), unpack=True, delimiter='	') #wavelength in nm, relative efficiency unitless
	
	basecaseind=9 #index of 0.1 bar CO2 case (NCO2=2.09e24 cm**-2), max radiance (A=new snow, sza=0) case
	
	umpgly_193s_base=umpgly_193s[basecaseind]
	umpgly_230s_base=umpgly_230s[basecaseind]
	umpgly_254s_base=umpgly_254s[basecaseind]
	cucn3_254s_base=cucn3_254s[basecaseind]
	cucn3_300s_base=cucn3_300s[basecaseind]
	
	
	###########Next, set up information about H2O file
	gas='o3'
	gaslabel='O3'

	minind=0 #index of minimum plausible H2O concentration	
	fiducialind=5 #index of the fiducial H2O vapor concentration
	maxind=8# index of maximum plausible H2O concentration
	###########Next, import H2O, and normalize by the base case.

	N_gas, rad100_165s, rad200_300s, umpgly_193s, umpgly_230s, umpgly_254s,cucn3_254s, cucn3_300s=np.genfromtxt('./Doses/'+gas+'_uv_doses.dat', skip_header=2, skip_footer=0, usecols=(1, 2, 3, 4, 5, 6, 7, 8), unpack=True, delimiter=' & ') #wavelength in nm, relative efficiency unitless

	
	#values<1: the reaction proceeds slower than in the base case. Values>1: the reaction proceeds faster than in the base case.
	umpgly_193s_normed=umpgly_193s/umpgly_193s_base
	umpgly_230s_normed=umpgly_230s/umpgly_230s_base
	umpgly_254s_normed=umpgly_254s/umpgly_254s_base
	cucn3_254s_normed=cucn3_254s/cucn3_254s_base
	cucn3_300s_normed=cucn3_300s/cucn3_300s_base
	
	###########Now, plot the dosimeters vs CO2 concentration
	
	#Initialize plot basics
	fig, ax1=plt.subplots(1, figsize=(cm2inch(16.5),7))

	colorseq1=iter(cm.rainbow(np.linspace(0,1,5)))
	
	#Plot linear half	
	ax1.plot(N_gas, umpgly_193s_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq1), label=r'UMP Gly Bond Cleavage ($\lambda_0=193$)')
	ax1.plot(N_gas, umpgly_230s_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq1), label=r'UMP Gly Bond Cleavage ($\lambda_0=230$)')
	ax1.plot(N_gas, umpgly_254s_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq1), label=r'UMP Gly Bond Cleavage ($\lambda_0=254$)')
	ax1.plot(N_gas, cucn3_254s_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq1), label=r'CuCN$_3$ Photoionization ($\lambda_0=254$)')
	ax1.plot(N_gas, cucn3_300s_normed, linestyle='--',markersize=4, marker='s',linewidth=1, color=next(colorseq1), label=r'CuCN$_3$ Photoionization ($\lambda_0=300$)')
	ax1.set_title(gaslabel+' (SZA=0, Albedo=New Snow)')
	ax1.axvline(N_gas[fiducialind], color='black', linewidth=2) #Mark the Rugheimer fiducial value
	ax1.axhline(1., color='black', linewidth=1) #Values above this: faster rxn than under 0.1 bar CO2

	#Finalize plot details.
	ax1.set_ylabel(r'Relative Dose Rate $D/D_{NCO2=2.09e24 cm^{-2}}$')
	ax1.set_xlabel(r'N$_{'+gaslabel+'}$ (cm$^{-2}$)')
	ax1.set_yscale('log')
	ax1.set_xscale('log')
	ax1.set_xlim([N_gas[minind],N_gas[maxind]])
	ax1.set_ylim([1.8e-3,1.8])

	plt.tight_layout(rect=(0,0,1., 0.80))
	ax1.legend(bbox_to_anchor=[0, 1.1, 1., 0.5], loc=3, ncol=2, mode='expand', borderaxespad=0., fontsize=9.5)  
	plt.savefig('./Plots/paperplots_'+gas+'_uvdoses.eps', orientation='portrait',papertype='letter', format='eps')
plt.show()