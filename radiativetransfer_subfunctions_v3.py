# -*- coding: iso-8859-1 -*-
"""
This file contains subfunctions used by the main radiative transfer code
"""

import numpy as np
import pdb
from scipy import interpolate as interp
import scipy.integrate

def compute_optical_parameters(colden_species_z,xc_tot_species_wav, xc_scat_species_wav, deltaapproxflag):
	"""
	This function computes the parameters tau_n_tot_z_wav, w_0_z_wav, and g_0_z_wav. These are parameters required for the two-stream code. 

	Inputs:
		colden_species_z: column density of a given species integrated across a given atmospheric layer. First dimension indexes particles, second dimension indexes atmospheric layer. TOA=first layer, planet surface=last layer, to accord with twostream conventions. Assume units of cm**-2.
		
		xc_tot_species_wav: total cross-section (absorption+scattering) of a given species at a given wavelength. Assume units of cm**2.
		
		xc_scat_species_wav: scattering cross-section of a given species at a given wavelength. Assume units of cm**2. For now, as we are working only with gases, all scattering is assumed to be Rayleigh scattering. 
		
		deltaapprox: Flag, with value 1 or 0. If 1, it implements the delta-Eddington approximation of Joseph and Wiscombe (1976); if 0, does not. 

	Outputs:
		tau_n_tot_z_wav: total optical depth across a given atmospheric layer, as a function of z (atmospheric layer) and wavelength. Includes contributions from gas absorption and gas (Rayleigh) scattering. One day, will also include particulate scattering, but not yet implemented.
		
		tau_c_tot_z_wav: cumulative optical depth until the TOP of layer n. 
		
		w_0_z_wav: single-scattering albedo, as a function of z (atmospheric layer) and wavelength. When equal to 1, all extinction from scattering; when equal to 0, all extinction from absorption
		
		g_z_wav: asymmetry parameter for the atmosphere, as a function of z (atmospheric layer) and wavelength. When equal to zero, scattering is symmetric. For Rayleigh scattering, z=0 (need particles for asymmetry.

	Calling example:

	tau_n_tot_z_wav, tau_c_tot_z_wav, w_0_z_wav, g_z_wav=compute_optical_parameters(colden_species_z,xc_tot_species_wav, xc_scat_species_wav, deltaapprox)

	Notes:
	1. This initial version of the code assumes we are working only with gas species in the UV. Consequently, all scattering is assumed to be Rayleigh, and all scattering is assumed to be symmetric.
	2. This version of the code follows twostr.f and sets a maximum of 0.999 on w_0. This is done to avoid instabilities in the two-stream matrix inversion. Coincidentally, this also fixes any cases where the Rayleigh scattering prediction exceeds measured absorption, which we see for CO2 in some bins. 
	3. The optical depth is computed along the zenith, following the formalism of Toon et al. This enables us to account for zenith angle separately, in the two-stream code.
	4. Following the Toon formalism, tau=0 at the TOA, so the first atmospheric layer corresponds to the TOA and the last atmospheric layer corresponds to the BOA.

	"""
	
	[N_species, N_layers]=np.shape(colden_species_z) #[Number of species in optical depth calculation, number of layers in atmospheric model]
	[N_species, N_wav]=np.shape(xc_tot_species_wav) #[Number of species in optical depth calculation, number of wavelength bins in model]
	
	
	########################
	###Compute optical depths
	########################
	tau_n_tot_z_wav=np.zeros([N_layers, N_wav]) #total optical depth across each atmospheric layer as a function of wavelength.
	tau_n_scat_z_wav=np.zeros([N_layers, N_wav]) #scattering optical depth across each atmospheric layer as a function of wavelength.
	for n in range(0, N_layers):
		for wav in range(0, N_wav):
			tau_n_tot_z_wav[n,wav]=np.sum(colden_species_z[:,n]*xc_tot_species_wav[:, wav])#Sum over species
			tau_n_scat_z_wav[n,wav]=np.sum(colden_species_z[:,n]*xc_scat_species_wav[:, wav])#sum over species
		#pdb.set_trace()

	##Patel+2002 dust loading insertion
	#tau_n_scat_z_wav[-1,:]=tau_n_scat_z_wav[-1,:]+0.1 ####+.1*.62
	#tau_n_tot_z_wav[-1,:]=tau_n_tot_z_wav[-1,:]+0.1 ####+.1
	
	tau_c_tot_z_wav=np.zeros([N_layers+1, N_wav]) #Cumulative optical depth until the top of the layer n. So tau_c[0,:]=0 as it is the TOA. tau_c[N_layers,:0] corresponds to the bottom of the atmosphere, i.e. the highest optical depth possible. 
	for ind in range(0, N_layers):
		tau_c_tot_z_wav[ind+1,:]=tau_c_tot_z_wav[ind,:]+tau_n_tot_z_wav[ind, :]
	
	########################
	###Compute single-scattering albedo
	########################
	w_0_z_wav=tau_n_scat_z_wav/tau_n_tot_z_wav #single scattering albedo across each atmospheric layer as a function of wavelength
	
	w_0_z_wav[w_0_z_wav>0.999]=0.999 #The twostr.f code imposes a ceiling on w_0 of 0.999, so we have done the same.
	########################
	###Compute asymmetry parameter
	########################
	#For Rayleigh-scattering gases, g=0.
	g_z_wav=np.zeros([N_layers, N_wav]) #asymmetry parameter across each atmospheric layer as a function of wavelength. 
	
	
	########################
	###Do delta-Eddington approximation?
	########################
	if deltaapproxflag==1:
		#This approximation reparametrizes tau_n, w_0, and g, and improves accuracy of the Eddington 2-stream solutions in the case of strongly forward-scattering species.
		#For our initial work, where we are working with symmetrically-scattering gases with g=0, this approximation has no effect.
		#In the CLIMA models (twostr.f, photo_CO2.f, and siblings) these approximations are used for the quadrature case also. However, I have seen some work which suggests that they are not applicable in the IR (hemispheric mean). We will have to look into this more closely when working with the IR.
		#In the limit of g=0, this changes nothing.
		#The implementation of this approximation is based on twostr.f from the CLIMA code.
		
		f_z_wav=g_z_wav*g_z_wav #Eqn 5a of Joseph et al
		
		tau_n_tot_z_wav=tau_n_tot_z_wav*(1.-w_0_z_wav*f_z_wav) #eqn 13 of Joseph et al
		w_0_z_wav=w_0_z_wav*(1.-f_z_wav)/(1.-w_0_z_wav*f_z_wav) #Eqn 14 of Joseph et al
		g_z_wav=g_z_wav/(1.+g_z_wav) #eqn 5b of Joseph et al
		
	
	return (tau_n_tot_z_wav, tau_c_tot_z_wav, w_0_z_wav, g_z_wav)

def get_z_layers(upperlimit, step):
	"""
	This function defines the layers length-space (z) of an atmosphere for the two-stream approximation
	Inputs:
		upperlimit: upper edge of atmosphere, in cm, inclusive.
		step: width of atmospheric bin, in cm
	Outputs:
		z_lower: lower edge of atmospheric layer (cm)
		z_upper: upper edge of atmospheric layer (cm)
		z_center: center of atmospheric layer (cm)
		N_layer: number of layers
	Calling example:
		z_lower, z_center, z_upper=get_z_layers(upperlimit,step)
	Notes:
		1. upperlimit/step should be an integer.
		2. Following the Toon et al convention, the first atmospheric bin is the TOA (top of atmosphere) and the last atmospheric bin is the BOA.
		3. We assume the lower limit of the layers is the surface, because the two-stream code must go down to the surface (is a B.C.)
	"""
	z_upper=np.arange(upperlimit,0, -step)
	z_lower=np.arange(upperlimit-step, -step, -step)
	z_center=0.5*(z_upper+z_lower)
	N_layer=np.size(z_center)
	
	return (z_lower, z_center, z_upper, N_layer)


def get_atmospheric_profile(z_lower, z_upper, datafile):
	"""
	This function returns p_z, t_z, n_z, columnndensity_z, and t_c, for use for our two-stream approximation code.
	
	Inputs:
		z_lower: lower edge of atmospheric layer (cm)
		z_upper: upper edge of atmospheric layer (cm)
		datafile=name of the file you are trying to read in. 
			The code assumes the following of the datafile:
				-the first 2 lines are header metatext, and are so skipped
				-The file consists of 4 columns. In order: z_center (cm), T (K), n (cm**-3), P (bar)
				-The first entries correspond to the bottom of the atmosphere, as this is common to climate model output like the Rugheimer CLIMA output.
				 This code is most accurate when the bin edges match the datafile bin edges. 
		
	Outputs:
		n_z: density at z_center in cm**-3
		t_z: temperature at z_center in K
		p_z: pressure at z_center in Barye
		columndensity_z: column density across layer in cm**-2
		t_c: temperature at bin edges in K. t_c[n] is the temperature at the top of layer n.
	Calling form:
	n_z, t_z, p_z, columndensity_z, t_c=get_atmospheric_profile(z_lower, z_upper, datafile)
	Notes:
		1. In accordance with the Toon et al methodology we are following, the first entry of the output corresponds to the TOA and the last entry corresponds to the BOA.
	"""
	k=1.3806488e-16 #Boltzmann Constant in erg/K
	bar2Ba=1.0e6 #1 bar in Ba

	N_layers=np.size(z_lower)
	z_center=0.5*(z_lower+z_upper)
	z_widths=z_upper-z_lower
	z_edges=np.append(z_upper, z_lower[-1]) #all bin edges

	#initialize output variables
	t_c=np.zeros(N_layers+1) #temperature at the top of each bin, in K. t_c[n] is the temperature at the top of bin n. t_c[N_layers] corresponds to BOA.
	n_z=np.zeros(N_layers) #mean number density in each layer, in cm**-3
	t_z=np.zeros(N_layers) #mean temperature in each layer, in K
	p_z=np.zeros(N_layers) #mean pressure in each layer, in barye
	columndensity_z=np.zeros(N_layers) #column number density of particles in each layer, in cm**-2.
	
	importeddata=np.genfromtxt(datafile, skip_header=2)
	model_z=importeddata[:,0] #z in cm. TOA at end of file (so first value corresponds to just over the planet surface). This is reversed from the twostr. methodology, and done to accomodate common T/P file formats. 
	model_t=importeddata[:,1] #temp of band in K
	model_n=importeddata[:,2] #number density in cm**-3
	model_p=importeddata[:,3]*bar2Ba #pressure in bar, converted to barye
	####units test
	###print np.max(np.abs(model_n-model_p/(k*model_t))/model_n) #fractional difference between the n from the model, and that computed using ideal gas law. If not small, than either our gas is non-ideal OR there is some problem in the model.
	###pdb.set_trace()
	
	####It's possible for the desired abscissa, to lie outside the range of z in the model file. We deal with this by padding. We pad with constant t, linearly varying p. Note that such an approximation will obviously fail the longer the lever arm is, so input should be constructed to avoid this problem (as we presently do)
	###if model_z[-1]<z_edges[0]:
		###model_z=np.append(model_z, z_edges[0])
		###model_t=np.append(model_t, model_t[-1])
		####linear interpolation for p using as base (model_z[-1], model_p[-1]) and (model_z[-2], model_p[-2]) and extrapolating FROM (model_z[-1], model_p[-1]) TO (z_edges[0], pval)
		###pslope=(model_p[-1]-model_p[-2])/(model_z[-1]-model_z[-2])
		###pval=pslope*(z_edges[0]-model_z[-1])+model_p[-1]
		###model_p=np.append(model_p, pval)
		###print 'Warning: desired abscissa beyond upper range of atmospheric profile file. Padding used.'
	####pdb.set_trace()

	#extract functional forms of parameters
	#Here we are assuming the values from the table are samples to the intrinsic functions (as opposed to being integrated over the layer)
	model_t_func=interp.interp1d(model_z, model_t, kind='linear')
	model_p_func=interp.interp1d(model_z, model_p, kind='linear')

	#evaluate z at the bin centers of the layers
	p_z=model_p_func(z_center)
	t_z=model_t_func(z_center)
	t_c=model_t_func(z_edges)
	
	n_z=p_z/(k*t_z) #get n from ideal gas law.
	
	columndensity_z=n_z*z_widths
	#note this methodology fails if the bins are wide enough that t, p, n diverge significantly from the linear approximation. For example, if the entire atmosphere was 1 bin, this would estimate the column density as (0.5*n_0)*deltaz, whereas with an exponential atmosphere it is closer to H*n_0*(1-exp(-deltaz/H)\approx H*n_0 for delta_z large. OTOH, for delta_z small, this would estimate the column density as (n_0)*deltaz, and the true answer would be \approx deltaz*n_0 as well. So, small bin=good. Bins matching the input best of course. 
	
	return (n_z, t_z, p_z, columndensity_z, t_c)


def get_mixing_ratios(z_centers,filename, molecule):
	"""
	This function returns the mixing ratio of each of our species, per layer. It extracts it from a provided file. 
	File should have two lines of header information. First column is height in atmosphere (z) in cm, and other columns are mixing ratios in the order given below.
	"""
	importeddata=np.genfromtxt(filename, skip_header=2, skip_footer=0)
	if molecule=='n2':
		ind=1
	elif molecule=='co2':
		ind=2
	elif molecule=='h2o':
		ind=3
	elif molecule=='ch4':
		ind=4
	elif molecule=='so2':
		ind=5
	elif molecule=='o2':
		ind=6
	elif molecule=='o3':
		ind=7
	elif molecule=='h2s':
		ind=8
	else:
		print 'Invalid entry for molecule'
		return z_centers*0.0
	
	z_list=importeddata[:,0] #This is the height in atmosphere, in cm
	mr_list=importeddata[:,ind] #This is the mixing ratio of the molecule of interest
	#####pdb.set_trace()
	#It's possible for z_centers to lie outside the range of the input data file. In this case, we pad the data with the highest and lowest values of the mixing ratio, respectively
	if z_list[0]>z_centers[-1]: #If the lowest z-value in the model exceeds the lowest z-value of the desired abscissa
		z_list=np.append(0.0, z_list)
		mr_list=np.append(mr_list[0], mr_list)
		print 'Warning: desired abscissa beyond lower range of mixing ratio file. Padding used.'
	if z_list[-1]<z_centers[0]:
		z_list=np.append(z_list, z_centers[0])
		mr_list=np.append(mr_list, mr_list[-1])
		print 'Warning: desired abscissa beyond upper range of mixing ratio file. Padding used.'
	
	mr_func=interp.interp1d(z_list, mr_list, kind='linear')
	mr_evals=mr_func(z_centers)
	return mr_evals

