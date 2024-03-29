# -*- coding: iso-8859-1 -*-
"""
Create files (from Rugheimer metadata) that give the atmospheric profile, i.e. mixing ratio, temperature and pressure as a function of altitude.
Since the Rugheimer T/P and mixing ratio files are generated from different codes, they have different abscissa, and so different files are generated for them. Interpolation is used in our code to match the two files.
"""
import numpy as np
import pdb
import matplotlib.pyplot as plt
import scipy.stats
from scipy import interpolate as interp
import cookbook

def extract_profiles_primitive_earth_rugheimer():
	"""
	Purpose of this code is to form spectra, mixing ratio files, and T/P profiles for the revised Rugheimer Epoch 0 (3.9 Ga) Earth models. This is to triangulate the sources of our differences.
	"""
	
	#####Zeroth: set value of constants, specify filenames
	import cookbook
	filename='./Raw_Data/Rugheimer_Metadata/outchem_Ep0_A0.2_Frac1.0.dat'
	bar2Ba=1.0e6 #1 bar in Ba
	k=1.3806488e-16 #Boltzmann Constant in erg/K
	
	
	#####First, form the spectra for comparison.
	importeddata=np.genfromtxt(filename, skip_header=290, skip_footer=1277)
	
	#Remove the first wavelength bin which corresponds to Lyman Alpha and which does not have a bin width that fits with its neighbors.
	rugheimer_wav_centers=importeddata[1:,1]/10. #Convert wavelengths from Angstroms to nm
	rugheimer_s=importeddata[1:,4] #ratio of 4piJ(surf)/I_0
	rugheimer_s[19]=3.16548e-128 #one element of rugheimer_s has value 3.16548e-128. Python has trouble with this and imports as a NaN. Here, we manually set its value.

	###Form wavelength bins from Rugheimer wavelength centers
	rugheimer_wav_bin_leftedges=np.zeros(len(rugheimer_wav_centers))
	rugheimer_wav_bin_rightedges=np.zeros(len(rugheimer_wav_centers))

	#First ten FUV fluxes are 5 nm (50 A) bins (email from srugheimer@gmail.com, 3/12/2015) 
	rugheimer_wav_bin_leftedges[0:9]=rugheimer_wav_centers[0:9]-2.5
	rugheimer_wav_bin_rightedges[0:9]=rugheimer_wav_centers[0:9]+2.5

	#Remainder of FUV fluxes are taken from a file that sarah sent me (srugheimer@gmail.com, 3/12/2015)
	del importeddata
	importeddata=np.genfromtxt('./Raw_Data/Rugheimer_Metadata/Active_M9_Teff2300_photo.pdat', skip_header=1, skip_footer=0)

	rugheimer_wav_bin_leftedges[9:]=importeddata[:,2]*0.1 #convert A to nm
	rugheimer_wav_bin_rightedges[9:]=importeddata[:,3]*0.1 #convert A to nm

	####Check that bins are correct:
	###print np.sum(rugheimer_wav_centers-0.5*(rugheimer_wav_bin_leftedges+rugheimer_wav_bin_rightedges)) #0 to within 1e-12 rounding error.

	###Rebin Claire et al input.
	#Import 0.01-nm resolution Claire et al 3.9 Ga Sun model.
	del importeddata
	importeddata=np.genfromtxt('./Raw_Data/Claire_Model/claire_youngsun_highres.dat', skip_header=1, skip_footer=0)
	claire_wav=importeddata[:,0] #nm, 0.01 nm resolution
	claire_fluxes=importeddata[:,1]#erg/s/cm2/nm

	#Bin Claire et al model to resolution of Rugheimer model
	claire_fluxes_rebinned=np.zeros(len(rugheimer_wav_centers))
	claire_wav_rebinned=np.zeros(len(claire_fluxes_rebinned))#This should be redundant with rugheimer_wav_centers. We include it as a check statistic that the rebinning is proceeding appropriately. 

	for ind in range(0, len(rugheimer_wav_centers)):
		min_wav=rugheimer_wav_bin_leftedges[ind]
		max_wav=rugheimer_wav_bin_rightedges[ind]
		inds=(claire_wav >= min_wav) & (claire_wav <= max_wav)
		claire_fluxes_rebinned[ind]=np.mean(claire_fluxes[inds])
		claire_wav_rebinned[ind]=np.mean(claire_wav[inds]) #check statistic.

	###print np.sum((claire_wav_rebinned-rugheimer_wav_centers)/rugheimer_wav_centers) #check statistic. Good to within 1e-5 in all cases. Any problems caused by slight misalignment from 0.01 due to rounding error. Good enough. 

	###Compute bottom-of-atmosphere actinic flux, which is what is reported in Rugheimer+2015.
	rugheimer_ground_energies=claire_fluxes_rebinned*rugheimer_s
	
	#Let's print out the results
	spectable=np.zeros([len(rugheimer_wav_bin_leftedges),5])
	spectable[:,0]=rugheimer_wav_bin_leftedges
	spectable[:,1]=rugheimer_wav_bin_rightedges
	spectable[:,2]=rugheimer_wav_centers
	spectable[:,3]=claire_fluxes_rebinned
	spectable[:,4]=rugheimer_ground_energies
	
	header='Left Bin Edge (nm)	Right Bin Edge (nm)	Bin Center (nm)		Solar Flux at Earth (erg/s/nm/cm2)		3.9 Ga BOA Intensity (erg/s/nm/cm2)\n'

	f=open('./LiteratureSpectra/rugheimer_epoch0_recomputed_A0.2.dat', 'w')
	f.write(header)
	np.savetxt(f, spectable, delimiter='		', fmt='%1.7e', newline='\n')
	f.close()
	
	###########################################################################################
	###########################################################################################
	###########################################################################################
	
	#####Second, form the mixing ratio files
	importeddata1=np.genfromtxt(filename, skip_header=779, skip_footer=873) #O2, O3, H2O
	importeddata2=np.genfromtxt(filename, skip_header=837, skip_footer=817) #CH4, SO2
	importeddata4=np.genfromtxt(filename, skip_header=958, skip_footer=704) #N2, CO2
	
	#Let's print out the results. We have established that the z values are the same, so can use a common block
	printtable=np.zeros([np.shape(importeddata1)[0],9])
	printtable[:,0]=importeddata1[:,0] #altitude in cm
	#N2 and CO2: We use the values from this block rather than block 1 because rugheimer et al force it to these values in their code, regardless of what the photochemistry code wants to do.
	printtable[:,1]=importeddata4[:,2] #N2. 
	printtable[:,2]=importeddata4[:,1] #CO2
	#The rest are normal
	printtable[:,3]=importeddata1[:,3] #H2O
	printtable[:,4]=importeddata2[:,2] #CH4
	printtable[:,5]=importeddata2[:,9] #SO2
	printtable[:,6]=importeddata1[:,2] #O2
	printtable[:,7]=importeddata1[:,8] #O3
	#printtable[:,8]# H2S; left as zeros since not included in Rugheimer model
	
	#print np.sum(printtable[:,1:],1)
	#pdb.set_trace()
	
	header0='Extracted from Rugheimer outchem_Ep0_A0.2_Frac1.0.dat\n'
	header1='Z (cm)		N2	CO2	H2O	CH4	SO2	O2	O3	H2S \n'

	f=open('./MixingRatios/rugheimer_earth_epoch0_recomputed_A0.2_mixingratios_v2.dat', 'w')
	f.write(header0)
	f.write(header1)
	np.savetxt(f, printtable, delimiter='	', fmt='%1.7e', newline='\n')
	f.close()	

	###########################################################################################
	###########################################################################################
	###########################################################################################
	#####Third, form the T/P profiles
	
	#Extract temperature and pressure profile from climate model output
	#For whatever reason the very last line of the table is doubled. We remove this. 
	importeddata=np.genfromtxt(filename, skip_header=1568, skip_footer=104)
	model_z=importeddata[:-1,0] #altitude in cm
	model_t=importeddata[:-1,1] #temperature in K
	model_n=importeddata[:-1,3] #number density in cm**-3.
	model_p=importeddata[:-1,4] #pressure, in bar (based on text in draft manuscript sent to me by Sarah Rugheimer)

	#Let's print out the results
	printtable=np.zeros([len(model_z)+1,4])
	printtable[1:,0]=model_z
	printtable[1:,1]=model_t
	printtable[1:,2]=model_n
	printtable[1:,3]=model_p
	
	#Rugheimer data file does not explicitly include t, P, n at z=0 (Surface). Our code requires z=0 data. To reconcile, we include these data manually as follows:
	printtable[0,0]=0. #z=0 case
	printtable[0,3]=1. #In the paper, p=1.0 bar at surface is specified
	printtable[0,1]=292.95 #From linear extrapolation from z=0.5 km and z=1.5 km points
	printtable[0,2]= 1.*bar2Ba/(k*292.95)#Compute number density self-consistently from temperature, pressure via Ideal Gas Law as is done elsewhere (n [cm**-3] = p [Barye]/(k*T [K])
	
	header0='Extracted from Rugheimer outchem_Ep0_A0.2_Frac1.0.dat\n'
	header1='Z (cm)	T (K)	DEN (cm**-3)	P (bar) \n'

	f=open('./TPProfiles/rugheimer_earth_epoch0_recomputed_A0.2_atmosphereprofile.dat', 'w')
	f.write(header0)
	f.write(header1)
	np.savetxt(f, printtable, delimiter='		', fmt='%1.7e', newline='\n')
	f.close()
#extract_profiles_primitive_earth_rugheimer()

def extract_profiles_modern_earth_rugheimer():
	"""
	Purpose of this code is to form spectra, mixing ratio files, and T/P profiles for the Rugheimer+2014 modern Earth surface UV models. This is a test case.
	"""
	
	#####Zeroth: set value of constants, specify filenames
	import cookbook
	filename='./Raw_Data/Rugheimer_Metadata/output_couple_Sun_100.dat'
	bar2Ba=1.0e6 #1 bar in Ba
	k=1.3806488e-16 #Boltzmann Constant in erg/K
	
	
	#####First, form the spectra for comparison. 
	#Extract spectra from Rugheimer file
	importeddata=np.genfromtxt(filename, skip_header=286, skip_footer=102)
	
	#Remove the first wavelength bin which corresponds to Lyman Alpha and which does not have a bin width that fits with its neighbors.
	spec_wav=importeddata[1:,0]*0.1 #A to nm
	spec_top=importeddata[1:,1]*1.e3 #W/m^2/nm to erg/cm^2/s/nm
	spec_gnd=importeddata[1:,2]*1.e3 #W/m^2/nm to erg/cm^2/s/nm
	#two elements of the file are not importing correctly, set them manually here	
	spec_gnd[23]=2.92059e-121*1.e3
	spec_gnd[24]=1.57780e-102 *1.e3
	
	#Next, extract the edges of the spectral bins.
	bin_left_edges=np.zeros(np.shape(spec_wav))
	bin_right_edges=np.zeros(np.shape(spec_wav))
	
	#first 9 bins are 5-nm (50 angstrom) wide bins (See faruv_sun.pdat)
	bin_left_edges[0:9]=spec_wav[0:9]-2.5
	bin_right_edges[0:9]=spec_wav[0:9]+2.5
	
	#The edges for the rest of the bins can be taken from G2V_photo.pdat:
	importeddata=np.genfromtxt('./Raw_Data/Rugheimer_Metadata/G2V_photo.pdat', skip_header=1, skip_footer=0)
	bin_left_edges[9:]=importeddata[:,2]*0.1 #convert from A to nm
	bin_right_edges[9:]=importeddata[:,3]*0.1 #convert from A to nm
	
	###let's validate our bin edges by computing the bin centers and making sure the residuals aren't too high
	##diff=(0.5*(bin_left_edges+bin_right_edges)-spec_wav)#/spec_wav
	##print diff
	##print np.max(np.abs(diff))
	###this test shows very slight offsets, at the 0.05 nm level at maximum. Should not affect results given bins are >1nm in width. 

	#Let's print out the results
	printtable=np.zeros([len(bin_left_edges),5])
	printtable[:,0]=bin_left_edges
	printtable[:,1]=bin_right_edges
	printtable[:,2]=spec_wav
	printtable[:,3]=spec_top
	printtable[:,4]=spec_gnd
	
	header='Left Bin Edge (nm)	Right Bin Edge (nm)	Bin Center (nm)		TOA Flux (erg/s/nm/cm2)		BOA Actinic Flux (erg/s/nm/cm2) \n'

	f=open('./LiteratureSpectra/rugheimer_earth_modern.dat', 'w')
	f.write(header)
	np.savetxt(f, printtable, delimiter='		', fmt='%1.7e', newline='\n')
	f.close()	
	###########################################################################################
	###########################################################################################
	###########################################################################################
	#####Second, form the mixing ratio files
	importeddata1=np.genfromtxt(filename, skip_header=78, skip_footer=323) #water, methane
	importeddata2=np.genfromtxt(filename, skip_header=182, skip_footer=222) #ozone, must derive from number density
	
	#Let's print out the results. We have established that the z values are the same, so can use a common block
	printtable=np.zeros([np.shape(importeddata1)[0],9])
	printtable[:,0]=importeddata1[:,0]*1.e5 #altitude in cm (converted from km)
	#N2 O2, and CO2: Well-mixed
	#H2O, CH4, O3: tracked through atmosphere
	#SO2: Not tracked. Assume 0. 
	printtable[:,1]=printtable[:,1]+ 0.78#N2; level tuned to assure 1 bar of surface pressure. Earth mean value given here.
	printtable[:,2]=printtable[:,2]+355.e-6 #CO2; level directly quoted in paper
	printtable[:,3]=importeddata1[:,2] #H2O
	printtable[:,4]=importeddata1[:,4] #CH4
	#printtable[:,5]=printtable[:,5] #SO2; left as zeros since not included in the model
	printtable[:,6]=printtable[:,6]+0.21 #O2; level directly quoted in paper
	printtable[:,7]=importeddata2[:,4]/importeddata2[:,2]#O3
	#printtable[:,8]=printtable[:,8]# H2S; left as zeros since not included in the model
	
	header0='Extracted from Rugheimer output_couple_Sun_100.dat\n'
	header1='Z (cm)		N2	CO2	H2O	CH4	SO2	O2	O3	H2S\n'

	f=open('./MixingRatios/rugheimer_earth_modern_mixingratios_v2.dat', 'w')
	f.write(header0)
	f.write(header1)
	np.savetxt(f, printtable, delimiter='	', fmt='%1.7e', newline='\n')
	f.close()


	###########################################################################################
	###########################################################################################
	###########################################################################################
	####Third, form the T/P profiles
	N_A=6.022e23 #Avogadro's number
	bar2Ba=1.0e6 #1 bar in Ba
	atm2bar=1.01325 #1 atm in bar
	k=83.14472/N_A #Boltzman constant in bar*cm^3/K, converted from bar*cm^3/(K*mol) (from http://www.engineeringtoolbox.com/individual-universal-gas-constant-d_588.html)
	
	#Extract temperature and pressure profile from climate model output
	importeddata=np.genfromtxt(filename, skip_header=409, skip_footer=0)
	model_z=importeddata[::-1,1]*1.e5 #altitude in cm, converted from km
	model_t=importeddata[::-1,2] #temperature in K
	model_p=importeddata[::-1,0]*atm2bar #pressure, in bar, converted from atm.
	model_n=model_p/(model_t*k) #number density in cm**-3, computed from ideal gas law.
	
	#Let's print out the results
	printtable=np.zeros([len(model_z),4])
	printtable[:,0]=model_z
	printtable[:,1]=model_t
	printtable[:,2]=model_n
	printtable[:,3]=model_p
	

	header0='Extracted from Rugheimer output_couple_Sun_100.dat\n'
	header1='Z (cm)	T (K)	DEN (cm**-3)	P (bar) \n'

	f=open('./TPProfiles/rugheimer_earth_modern_atmosphereprofile.dat', 'w')
	f.write(header0)
	f.write(header1)
	np.savetxt(f, printtable, delimiter='		', fmt='%1.7e', newline='\n')
	f.close()
#extract_profiles_modern_earth_rugheimer()

def form_profiles_wuttke():
	"""
	Purpose of this code is to form the feedstock files to replicat the Wuttke+2006 Antarctic diffuse radiance measurements
	"""
	import cookbook
	#First, form the spectral file.
	#Define spectral bins. 0.25 nm from 280-500 nm, 1 nm from 500-1000 nm. We just go to 900 since that's what our data is good to. Also we start at 292.75 because that's where our graphclicked data starts
	bin_left_edges=np.concatenate((np.arange(292.75,500.,0.25),np.arange(500., 900.,1.)))
	bin_right_edges=np.concatenate((np.arange(293.,500.25,0.25),np.arange(501., 901.,1.)))
	bin_centers=0.5*(bin_left_edges+bin_right_edges)
	
	
	#load BOA diffuse zenith flux from Wuttke+2006 (extracted via GraphClick)	
	importeddata=np.genfromtxt('./Raw_Data/UV_Surface_Measurements/wuttke.csv', skip_header=0, skip_footer=0, delimiter=',')
	dif_wav=importeddata[:,0] #nm
	dif_flux=importeddata[:,1]*2.*np.pi #mW/m2/nm/sr=erg/s/cm2/nm/sr; multiply by 2pi to convert to hemisphere-integrated total surface diffuse radiances
	dif_func=interp.interp1d(dif_wav, dif_flux, kind='linear')
	dif_flux_interp=dif_func(bin_centers)
	
	#load solar spectrum from Claire et al (2012) models, normalized to 1 au 
	importeddata=np.genfromtxt('./Raw_Data/Claire_Model/claire_modernsun_highres.dat', skip_header=1, skip_footer=0)
	claire_wav=importeddata[:,0] #nm, 0.1 nm resolution, 100-900 nm.
	claire_fluxes=importeddata[:,1]#erg/s/cm2/nm
	
	#rebin claire spectrum
	claire_fluxes_rebinned=cookbook.rebin_uneven(np.arange(99.995,900.005,0.01), np.arange(100.005, 900.015,0.01),claire_fluxes,bin_left_edges, bin_right_edges)   
	
	#Plot to make sure rebinning worked correctly
	fig, ax1=plt.subplots(1, figsize=(6,4))
	ax1.plot(claire_wav, claire_fluxes, marker='s', color='black', label='Claire Fluxes')
	ax1.plot(bin_centers, claire_fluxes_rebinned, marker='s', color='blue', label='Binned Claire Fluxes')	
	ax1.set_yscale('log')
	ax1.set_ylim([1.e-2, 1.e4])
	ax1.set_xlim([280.,900.])
	ax1.set_xlabel('nm')
	ax1.set_ylabel('erg/s/cm2/nm')
	ax1.legend(loc=0)
	plt.show()	
	
	#Let's print out the results
	spectable=np.zeros([len(bin_left_edges),5])
	spectable[:,0]=bin_left_edges
	spectable[:,1]=bin_right_edges
	spectable[:,2]=bin_centers
	spectable[:,3]=claire_fluxes_rebinned
	spectable[:,4]=dif_flux_interp
	
	header='Left Bin Edge (nm)	Right Bin Edge (nm)	Bin Center (nm)		Top of Atm Flux (erg/s/nm/cm2)		Zenith Diffuse Flux (erg/s/nm/cm2)\n'

	f=open('./LiteratureSpectra/wuttke2006.dat', 'w')
	f.write(header)
	np.savetxt(f, spectable, delimiter='		', fmt='%1.7e', newline='\n')
	f.close()
	
	###########################################################################################
	###########################################################################################
	###########################################################################################
	#####Second, form the mixing ratio files
	#####Form by replicating the Rugheimer modern Earth profile, then scaling down the H2O level and scaling up the O3 level.
	mixingratios=np.genfromtxt('./MixingRatios/rugheimer_earth_modern_mixingratios_v2.dat', skip_header=2, skip_footer=0)
	mixingratios[:,3]=mixingratios[:,3]*0.1 #scale down h2o by factor of 10
	mixingratios[:,7]=mixingratios[:,7]*1.25 #scale up ozone by factor of 1.25
	header0='Based on Rugheimer+2013 Modern Earth Model\n'
	header1='Z (cm)		N2	CO2	H2O	CH4	SO2	O2	O3	H2S\n'

	f=open('./MixingRatios/wuttke2006_mixingratios_v2.dat', 'w')
	f.write(header0)
	f.write(header1)
	np.savetxt(f, mixingratios, delimiter='	', fmt='%1.7e', newline='\n')
	f.close()
	###########################################################################################
	###########################################################################################
	###########################################################################################
	#####Finally, form TP profile
	#####Form by duplicating Rugheimer+2013 modern Earth profile
	tpprofile=np.genfromtxt('./TPProfiles/rugheimer_earth_modern_atmosphereprofile.dat', skip_header=2, skip_footer=0)
	header0='Based on Rugheimer+2013 Modern Earth Model\n'
	header1='Z (cm)	T (K)	DEN (cm**-3)	P (bar) \n'

	f=open('./TPProfiles/wuttke2006_atmosphereprofile.dat', 'w')
	f.write(header0)
	f.write(header1)
	np.savetxt(f, tpprofile, delimiter='		', fmt='%1.7e', newline='\n')
	f.close()
#form_profiles_wuttke()

def form_profiles_woudc():
	"""
	Purpose of this code is to form the feedstock files to replicate the irradiance measurements from the WOUDC website for Toronto (June 21 2003, SZA=20.376, O3=354, Brewer no. 145)
	""" 
	########First, form the spectral file.
	#load measured irradiances	
	importeddata=np.genfromtxt('./Raw_Data/UV_Surface_Measurements/woudc_toronto_2003_145_cut.dat', skip_header=1, skip_footer=0, delimiter='	')
	woudc_wav=importeddata[:,0] #nm
	woudc_flux=importeddata[:,1]*1.e3 #W/m2/nm=1000 erg/s/cm2/nm
	#woudc_func=interp.interp1d(woudc_wav, woudc_flux, kind='linear')
	#woudc_flux_interp=dif_func(bin_centers)

	#Define spectral bins.
	bin_centers=woudc_wav
	bin_left_edges=woudc_wav-0.25
	bin_right_edges=woudc_wav+0.25
	
	#load solar spectrum from Claire et al (2012) models, normalized to 1 au 
	importeddata2=np.genfromtxt('/home/sranjan/IDL/UV/YoungSun/claire_modernsun_highres.dat', skip_header=1, skip_footer=0)
	claire_wav=importeddata2[:,0] #nm, 0.1 nm resolution, 100-900 nm.
	claire_fluxes=importeddata2[:,1]#erg/s/cm2/nm

	#rebin claire spectrum
	claire_fluxes_rebinned=cookbook.rebin_uneven(np.arange(99.995,900.005,0.01), np.arange(100.005, 900.015,0.01),claire_fluxes,bin_left_edges, bin_right_edges)   
	
	#Plot to make sure rebinning worked correctly
	fig, ax1=plt.subplots(1, figsize=(6,4))
	ax1.plot(claire_wav, claire_fluxes, marker='s', color='black', label='Claire Fluxes')
	ax1.plot(bin_centers, claire_fluxes_rebinned, marker='s', color='blue', label='Binned Claire Fluxes')	
	ax1.set_yscale('log')
	ax1.set_ylim([1.e-2, 1.e4])
	ax1.set_xlim([280.,360.])
	ax1.set_xlabel('nm')
	ax1.set_ylabel('erg/s/cm2/nm')
	ax1.legend(loc=0)
	plt.show()	
	
	#Let's print out the results
	spectable=np.zeros([len(bin_left_edges),5])
	spectable[:,0]=bin_left_edges
	spectable[:,1]=bin_right_edges
	spectable[:,2]=bin_centers
	spectable[:,3]=claire_fluxes_rebinned
	spectable[:,4]=woudc_flux
	
	header='Left Bin Edge (nm)	Right Bin Edge (nm)	Bin Center (nm)		Top of Atm Flux (erg/s/nm/cm2)		Surface Flux (erg/s/nm/cm2)\n'

	f=open('./LiteratureSpectra/woudc.dat', 'w')
	f.write(header)
	np.savetxt(f, spectable, delimiter='		', fmt='%1.7e', newline='\n')
	f.close()
	
	###########################################################################################
	###########################################################################################
	###########################################################################################
	#####Second, form the mixing ratio files
	#####Form by replicating the Rugheimer modern Earth profile, then scaling down the H2O level and scaling up the O3 level.
	mixingratios=np.genfromtxt('./MixingRatios/rugheimer_earth_modern_mixingratios_v2.dat', skip_header=2, skip_footer=0)
	mixingratios[:,7]=mixingratios[:,7]*1.77 #scale up ozone by factor of 1.25
	header0='Based on Rugheimer+2013 Modern Earth Model\n'
	header1='Z (cm)		N2	CO2	H2O	CH4	SO2	O2	O3	H2S\n'

	f=open('./MixingRatios/woudc_mixingratios_v2.dat', 'w')
	f.write(header0)
	f.write(header1)
	np.savetxt(f, mixingratios, delimiter='	', fmt='%1.7e', newline='\n')
	f.close()
	###########################################################################################
	###########################################################################################
	###########################################################################################
	#####Finally, form TP profile
	#####Form by duplicating Rugheimer+2013 modern Earth profile
	tpprofile=np.genfromtxt('./TPProfiles/rugheimer_earth_modern_atmosphereprofile.dat', skip_header=2, skip_footer=0)
	header0='Based on Rugheimer+2013 Modern Earth Model\n'
	header1='Z (cm)	T (K)	DEN (cm**-3)	P (bar) \n'

	f=open('./TPProfiles/woudc_atmosphereprofile.dat', 'w')
	f.write(header0)
	f.write(header1)
	np.savetxt(f, tpprofile, delimiter='		', fmt='%1.7e', newline='\n')
	f.close()
#form_profiles_woudc()

def form_spectral_feedstock_ourwork():
	"""
	Purpose of this code is to form the spectral feedstock file to explore formally the dependence of UV surface intensity on various factors. The mixing ratio and TP profiles vary in each case though.
	"""
	import cookbook
	#Extract spectra to match and TOA intensity
	#Define spectral bins.
	bin_left_edges=np.arange(100.,500.,1.)
	bin_right_edges=np.arange(101.,501.,1.)
	bin_centers=0.5*(bin_left_edges+bin_right_edges)
	
	
	#There are no literature intensity values for this file, since at this point we are not comparing against any other datasets but are rather running our code internally. However, we can use the Rugheimer et al base case (60 degrees, 0.2) as a reference
	literature_intensities=np.zeros(np.shape(bin_centers))
	importeddata=np.genfromtxt('./TwoStreamOutput/AlbZen/rugheimer_earth_epoch0_a=0.2_z=60.dat', skip_header=1, skip_footer=0)
	basecase_wav=importeddata[:,2] #nm,
	basecase_surface_intensities=importeddata[:,6] #erg/s/cm2/nm
	
	#load solar spectrum from Claire et al (2012) models, normalized to 1 au. These are really TOA intensities. Multiply by mu_0 to get TOA fluxes. 
	importeddata=np.genfromtxt('./Raw_Data/Claire_Model/claire_youngsun_highres.dat', skip_header=1, skip_footer=0)
	claire_wav=importeddata[:,0] #nm, 0.01 nm resolution, 100-900 nm.
	claire_fluxes=importeddata[:,1]#erg/s/cm2/nm
	
	#rebin claire spectrum
	claire_fluxes_rebinned=cookbook.rebin_uneven(np.arange(99.995,900.005,0.01), np.arange(100.005, 900.015,0.01),claire_fluxes,bin_left_edges, bin_right_edges)   
	
	#Plot to make sure rebinning worked correctly
	fig, ax1=plt.subplots(1, figsize=(6,4))
	ax1.plot(claire_wav, claire_fluxes, marker='s', color='black', label='Claire Fluxes')
	ax1.plot(bin_centers, claire_fluxes_rebinned, marker='s', color='blue', label='Binned Claire Fluxes')	
	ax1.set_yscale('log')
	ax1.set_ylim([1.e-2, 1.e4])
	ax1.set_xlim([100.,500.])
	ax1.set_xlabel('nm')
	ax1.set_ylabel('erg/s/cm2/nm')
	ax1.legend(loc=0)
	plt.show()	
	
	#Let's print out the results
	spectable=np.zeros([len(bin_left_edges),5])
	spectable[:,0]=bin_left_edges
	spectable[:,1]=bin_right_edges
	spectable[:,2]=bin_centers
	spectable[:,3]=claire_fluxes_rebinned
	spectable[:,4]=basecase_surface_intensities
	
	header='Left Bin Edge (nm)	Right Bin Edge (nm)	Bin Center (nm)		Top of Atm Intensity (erg/s/nm/cm2)		3.9 Ga R+2015 Surface Intensity (erg/s/nm/cm2)\n'

	f=open('./LiteratureSpectra/general_spectral_input.dat', 'w')
	f.write(header)
	np.savetxt(f, spectable, delimiter='		', fmt='%1.7e', newline='\n')
	f.close()
#form_spectra_feedstock_ourwork()

def form_profiles_co2limtests():
	"""
	Purpose of this code is to form mixing ratio and T/P profile for our exploration of the surface environment on the 3.9 Ga Earth for a range of two-component atmospheres of CO2 and N2. N2 abundance is always fixed at 0.9 bar equivalent for consistency with Rugheimer et al (2015), while CO2 abundance varies. We derive these by reading in the values for the Rugheimer (2015) atmosphere, which is at 1 bar, and scaling it. 
	"""
	k=1.38064852e-16 #Boltzman constant in erg/K
	bar2Ba=1.0e6 #1 bar in Ba
	multiples=np.array([0., 1.e-6, 1.e-5, 1.e-4, 1.e-3, 1.e-2, 1.e-1, 1., 1.e1, 1.e2, 1.e3, 1.33, 46.6, 470., .6, 8.93e-3]) #values we will be scaling the CO2 column by 

	####################
	####Mixing ratios: 
	####################
	importeddata1=np.genfromtxt('./MixingRatios/rugheimer_earth_epoch0_recomputed_A0.2_mixingratios_v2.dat', skip_header=2, skip_footer=0) #use for z profile
	
	for ind in range(0,len(multiples)):
		#Let's print out the results.
		printtable=np.zeros([np.shape(importeddata1)[0],9])
		printtable[:,0]=importeddata1[:,0] #altitude in cm
				
		#N2 and CO2: There must be 0.9 bar n2, and the mixing ratios of CO2 and N2 must sum to zero. This is a 2x2 system, so:
		mr_n2=0.9/(0.9+0.1*multiples[ind])
		mr_co2=1.-mr_n2
		
		printtable[:,1]+=mr_n2 #N2
		printtable[:,2]+=mr_co2 #CO2
		
		#The rest are left to 0 for the purposes of these tests
		printtable[:,3]+=0 #H2O
		printtable[:,4]+=0  #CH4
		printtable[:,5]+=0  #SO2
		printtable[:,6]+=0  #O2
		printtable[:,7]+=0  #O3
		printtable[:,8]+=0  #H2S
		
		#print np.sum(printtable[:,1:],1)
		#pdb.set_trace()
		
		header0='Stipulated\n'
		header1='Z (cm)		N2	CO2	H2O	CH4	SO2	O2	O3	H2S\n'

		f=open('./MixingRatios/CO2lim/rugheimer_earth_epoch0_mixingratios_co2limits_co2multiple='+str(multiples[ind])+'.dat', 'w')
		f.write(header0)
		f.write(header1)
		np.savetxt(f, printtable, delimiter='	', fmt='%1.7e', newline='\n')
		f.close()
	####################
	####T/P Profile: Scaled from Rugheimer+2015
	####################
	#Extract temperature and pressure profile from climate model output
	importeddata=np.genfromtxt('./TPProfiles/rugheimer_earth_epoch0_recomputed_A0.2_atmosphereprofile.dat', skip_header=2, skip_footer=0)
	model_z=importeddata[:,0] #altitude in cm
	model_t=importeddata[:,1] #temperature in K
	model_n=importeddata[:,2] #number density in cm**-3.
	model_p=importeddata[:,3] #pressure, in bar (based on text in draft manuscript sent to me by Sarah Rugheimer)
	
	
	for ind in range(0, len(multiples)):
		printtable=np.zeros([len(model_z),4])
		printtable[:,0]=model_z
		printtable[:,1]=model_t
		printtable[:,2]=model_n*(0.9+0.1*multiples[ind]) #scale the densities, new total column density is 0.9 (N2) plus scaled CO2 column
		printtable[:,3]=k*printtable[:,2]*printtable[:,1]/bar2Ba #Compute the pressures from the ideal gas law, assuming identical temperature profile, and convert to bar
		
		header0='Scaled from rugheimer_earth_epoch0_recomputed_A0.2_atmosphereprofile.dat\n'
		header1='Z (cm)	T (K)	DEN (cm**-3)	P (bar) \n'
		f=open('./TPProfiles/CO2lim/rugheimer_earth_epoch0_atmosphereprofile_co2limits_co2multiple='+str(multiples[ind])+'.dat', 'w')
		f.write(header0)
		f.write(header1)
		np.savetxt(f, printtable, delimiter='		', fmt='%1.7e', newline='\n')
		f.close()

#form_profiles_co2limtests()

def form_profiles_othergases():
	"""
	Purpose of this code is to form the atmospheric profile (z, T, n, P) and the mixing ratio/molar concentration profile for our exploration of the surface environment on the 3.9 Ga Earth for a range of two-component atmospheres of N2 and X, where X is the other gases in our model (H2O, CH4, SO2, H2S, O2, O3). N2 abundance is fixed at the 0.9 bar value used in the Rugheimer+2015 model (i.e. column density of 1.883e25 cm**-2), whereas the abundance of X is allowed to vary  We derive these by reading in the values for the Rugheimer (2015) atmosphere, which is at 1 bar, and scaling them appropriately. Note that the column density of the Rugheimer atmosphere corresponds to 2.0922e25 cm**-2.
	"""
	k=1.38064852e-16 #Boltzman constant in erg/K
	bar2Ba=1.0e6 #1 bar in Ba
	g=981. #acceleration due to gravity in cm/s/s
	convfactor_Ntop=1.e-6 #(1 g)*(1 cm s**-2)*(1 cm**-2) in bar. prefactor required for the calculation of partial pressures from column densities
	amu2g=1.66054e-24 #1 amu in g
	
	####################
	####T/P Profile: Scaled from Rugheimer+2015
	####################
	#Extract temperature and pressure profile from climate model output
	importeddata=np.genfromtxt('./TPProfiles/rugheimer_earth_epoch0_recomputed_A0.2_atmosphereprofile.dat', skip_header=2, skip_footer=0)
	model_z=importeddata[:,0] #altitude in cm
	model_t=importeddata[:,1] #temperature in K
	model_n=importeddata[:,2] #number density in cm**-3.
	model_p=importeddata[:,3] #pressure, in bar (based on text in draft manuscript sent to me by Sarah Rugheimer)
	

	####################
	####Mixing ratios: 
	####################
	importeddata1=np.genfromtxt('./MixingRatios/rugheimer_earth_epoch0_recomputed_A0.2_mixingratios_v2.dat', skip_header=2, skip_footer=0) #use for z profile
	
	###########
	gaslist=['h2o', 'ch4', 'so2', 'o2', 'o3', 'h2s'] #list of gases we are doing this for
	base_abundances=np.array([4.762e-3, 1.647e-6, 3.371e-11, 2.707e-6, 9.160e-11, 6.742e-11]) #molar concentration of each of these gases in the Rugheimer model.
	m_gases=np.array([18.015, 16.042, 64.064, 31.999, 47.998, 34.081])*amu2g #mean molecular weight of each of these gases, converted from amu to g
	
	m_n2=28.013*amu2g #mean molecular weight of n2, in g (should be 4.6517e-23....and checks out!)
	N_tot=2.0925e25#total column density of Rugheimer+2015 model in cm**-2
	
	#dict holding the multiples of the molar concentration we are using
	gasmultiples={}
	gasmultiples['h2o']=np.array([1.e-5, 1.e-4, 1.e-3, 1.e-2, 1.e-1, 1., 1.e1, 1.e2, 1.e3, 1.e4, 1.e5])
	gasmultiples['ch4']=np.array([1.e-5, 1.e-4, 1.e-3, 1.e-2, 1.e-1, 1., 1.e1, 1.e2, 1.e3, 1.e4, 1.e5])
	gasmultiples['so2']=np.array([1.e-5, 1.e-4, 1.e-3, 1.e-2, 1.e-1, 1., 1.e1, 1.e2, 1.e3, 1.e4, 1.e5, 1.e6, 1.e7])
	gasmultiples['o2']=np.array([1.e-5, 1.e-4, 1.e-3, 1.e-2, 1.e-1, 1., 1.e1, 1.e2, 1.e3, 1.e4, 1.e5])
	gasmultiples['o3']=np.array([1.e-5, 1.e-4, 1.e-3, 1.e-2, 1.e-1, 1., 1.e1, 1.e2, 1.e3, 1.e4, 1.e5])
	gasmultiples['h2s']=np.array([1.e-5, 1.e-4, 1.e-3, 1.e-2, 1.e-1, 1., 1.e1, 1.e2, 1.e3, 1.e4, 1.e5, 1.e6, 1.e7, 1.e8])
	
	gasbygas_atm_params=np.zeros([1,6])#we will append to this list. table holding column densities, mean molecular masses, and gas partial pressures for each of our model atmospheres
	for gasind in range(0, len(gaslist)):
		gas=gaslist[gasind] #which gas we are dealing with right now
		multiples=gasmultiples[gas] #multiples of gas
		m_X=m_gases[gasind] #mean molecular mass of this gas
		base_abundance=base_abundances[gasind]#molar concentration of gas in Rugheimer model.
		
		for multind in range(0, len(multiples)):
			multiple=multiples[multind]
			
			abundance=base_abundance*multiple #what fraction of the Rugheimer total atmosphere column is allocated to this gas?			
			#####First, number density profile
			tnp_table=np.zeros([len(model_z),4])
			tnp_table[:,0]=model_z
			tnp_table[:,1]=model_t
			tnp_table[:,2]=model_n*(0.9+abundance) #scale the densities, new total column density is 0.9 (N2) plus scaled X column
			tnp_table[:,3]=k*tnp_table[:,2]*tnp_table[:,1]/bar2Ba #Compute the pressures from the ideal gas law, assuming identical temperature profile, and convert to bar
		
			header0='Scaled from rugheimer_earth_epoch0_recomputed_A0.2_atmosphereprofile.dat\n'
			header3='Z (cm)	T (K)	DEN (cm**-3)	P (bar) \n'
			f=open('./TPProfiles/gaslim/rugheimer_earth_epoch0_atmosphereprofile_'+gas+'limits_'+gas+'multiple='+str(multiple)+'.dat', 'w')
			f.write(header0)
			f.write(header3)
			np.savetxt(f, tnp_table, delimiter='		', fmt='%1.7e', newline='\n')
			f.close()
		
			#####Second, 'mixing ratio' profile		
			#N2 and CO2: There must be 0.9 bar-equivalent n2, and the concentrations of CO2 and N2 must sum to 1. So
			conc_n2=0.9/(0.9+abundance)
			conc_X=abundance/(0.9+abundance)
			
			conc_table=np.zeros([np.shape(importeddata1)[0],9])
			conc_table[:,0]=importeddata1[:,0] #altitude in cm			
			conc_table[:,1]+=conc_n2 #N2
			conc_table[:,gasind+3]+=conc_X #X 0=z, 1=n2, 2=co2, 3=h2o...
			#The rest are left to 0 for the purposes of these tests
			
			header0='Stipulated\n'
			header1='Z (cm)		N2	CO2	H2O	CH4	SO2	O2	O3	H2S\n'
			f=open('./MixingRatios/gaslim/rugheimer_earth_epoch0_mixingratios_'+gas+'limits_'+gas+'multiple='+str(multiple)+'.dat', 'w')
			f.write(header0)
			f.write(header1)
			np.savetxt(f, conc_table, delimiter='	', fmt='%1.7e', newline='\n')
			f.close()

			
			#####Finally, compute the mean molecular weight and surface partial pressures associated with such an atmosphere.
			N_n2=0.9*N_tot #n2 column density in cm**-2
			N_X=abundance*N_tot #X column density in cm**-2
			
			bar_m=conc_n2*m_n2+conc_X*m_X
			p_n2=bar_m*g*N_n2*convfactor_Ntop #N2 partial pressure in bar
			p_X=bar_m*g*N_X*convfactor_Ntop #X partial pressure in bar
			
			if multind==0:
				gasbygas_atm_params=np.array([multiple, N_n2, N_X, bar_m, p_n2, p_X])#we will append to this list. table holding column densities, mean molecular masses, and gas partial pressures for each of our model atmospheres
			else:
				gasbygas_atm_params=np.vstack((gasbygas_atm_params, np.array([multiple, N_n2, N_X, bar_m, p_n2, p_X])))
		f=open('./TPProfiles/gaslim/rugheimer_earth_epoch0_atmosphereprofile_'+gas+'limits_partialpressures.dat', 'w')
		f.write('multiple	N_n2 (cm-2)	N_X (cm-2)	mbar (g)	p_n2 (bar)	p_x (bar) \n')
		np.savetxt(f, gasbygas_atm_params, delimiter='	', fmt='%1.7e', newline='\n')
		f.close()
		del  gasbygas_atm_params	
form_profiles_othergases()


