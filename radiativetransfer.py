# -*- coding: iso-8859-1 -*-
"""
Main code, runs the radiative transfer subroutines and writes output.
"""

########################
###Import useful libraries & define constants
########################

import numpy as np
import matplotlib.pyplot as plt
import pickle
import pdb

import cookbook
import cross_sections_subfunctions as css
import twostream_toon_func as twostr
import radiativetransfer_albedo_subfunctions as ras
import radiativetransfer_subfunctions as rts

########################
###Constants
########################

bar2Ba=1.0e6 #1 bar in Ba
amu2g=1.66054e-24 #1 amu in g

########################
###User-set parameters
########################
multiples=np.array([0., 1.e-6, 1.e-5, 1.e-4, 1.e-3, 1.e-2, 1.e-1, 1., 1.e1, 1.e2, 1.e3, 1.33, 46.6, 470., .6, 8.93e-3]) #values we will be scaling the CO2 column by 
multiple=str(multiples[10])

filename='/CO2lim/surface_intensities_co2limits_co2multiple='+multiple+'_a=newsnow_z=0' #name of file to write output, plot to

#TOA input, and comparison file.
inputspectrafile='./LiteratureSpectra/general_spectral_input.dat' #TOA Stellar Input and Reference Spectra for Comparison

#Mixing ratio file
N_species=8 #how many gas species are in our model?
mr_profilefile='./MixingRatios/CO2lim/rugheimer_earth_epoch0_mixingratios_co2limits_co2multiple='+multiple+'.dat' #Mixing ratio file

#T/P Profile.
atmoprofilefile='./TPProfiles/CO2lim/rugheimer_earth_epoch0_atmosphereprofile_co2limits_co2multiple='+multiple+'.dat' #T/P Profile File. File boundaries should match z_lower and z_upper

#RT model layers
z_upper_limit=64.e5 #Upper edge of the atmospheric layers, in cm. The bottom edge of the layers is assumed to be at 0.
z_step=1.e5 # Thickness of the atmospheric layers, in cm.

###Solar Zenith Angle
solarzenithangle=0*np.pi/180. #

###Albedo
#Uniform albedo: set flag to 'uniformalbedo' and set the subsequent albedo value
#Nonuniform albedo: set flag to 'nonuniformalbedo' and choose what mix of surface: ocean, old snow, new snow, desert, tundra. 
#set flag to 'plot' to plot albedo used. 
albedoflag='nonuniformalbedo' #value of 'nonuniformalbedo' or 'uniformalbedo'
uniformalbedo=0. #if adopting uniform albedo: what value?
frac_ocean=0. #if nonuniform albedo: what fraction of ground is ocean?
frac_tundra=0. #if nonuniform albedo: what fraction of ground is tundra?
frac_desert=0. #if nonuniform albedo: what fraction of ground is desert?
frac_oldsnow=0.#if nonuniform albedo: what fraction of ground is old snow?
frac_newsnow=1. #if nonuniform albedo: what fraction of ground is new snow?


###Temperature of Ground
temp_ground=292.5 #Temperature of the ground in K
########################
###Set output file names
########################
writefilename='./TwoStreamOutput/'+filename+'.dat' #name of file to write output to. 
plotfilename='./Plots/'+filename+'.pdf' #name of file to save plot to.

########################
###Set angular parameters relevant to the RT model
########################
mu_0=np.cos(solarzenithangle)#Cosine of solar zenith angle.
mu_1=1./np.sqrt(3.) #The Gaussian quadrature angle for the diffuse flux direction.

########################
###Load in TOA stellar input
########################
importeddata=np.genfromtxt(inputspectrafile, skip_header=1, skip_footer=0)
wav_leftedges=importeddata[:,0] #left edges of wavelength bins, nm
wav_rightedges=importeddata[:,1] #right edges of wavelength bins, nm
wav_centers=importeddata[:,2] #centers of wavelength bins, nm
intensity_toa=importeddata[:,3] #top-of-atmosphere total intensity, erg/s/cm2/nm. Multiply by mu_0 to get TOA flux.
surface_intensity_basecase=importeddata[:,4] #total surface intensity computed for the baseline Rugheimer+2015 3.9 Ga atmosphere, erg/s/cm2/nm. Included for comparison purposes.

N_wavelengths=len(wav_centers)
wav_edges=np.append(wav_leftedges, wav_rightedges[-1])
########################
###Load other key inputs
########################

###Layers of the atmosphere
#This is required for all variants of running the code
z_lower, z_center, z_upper, N_layers=rts.get_z_layers(z_upper_limit, z_step) #define layers of atmosphere.
z_edges=np.append(z_upper, z_lower[-1])


###Mixing ratios
#Set the mixing ratios of the gases. If set to a number, the gas is assumed to be well-mixed, i.e. have this mixing ratio everywhere in the column. If set to a filename, the mixing ratio is taken from that file.
mr_n2=rts.get_mixing_ratios(z_center, mr_profilefile, 'n2')
mr_co2=rts.get_mixing_ratios(z_center, mr_profilefile, 'co2')
mr_h2o=rts.get_mixing_ratios(z_center, mr_profilefile, 'h2o')
mr_ch4=rts.get_mixing_ratios(z_center, mr_profilefile, 'ch4')
mr_so2=rts.get_mixing_ratios(z_center, mr_profilefile, 'so2')
mr_o2=rts.get_mixing_ratios(z_center, mr_profilefile, 'o2')
mr_o3=rts.get_mixing_ratios(z_center, mr_profilefile, 'o3')
mr_h2s=rts.get_mixing_ratios(z_center, mr_profilefile, 'h2s')

####Albedos
albedo_dif_wav=ras.get_surface_albedo(wav_leftedges, wav_rightedges,solarzenithangle, albedoflag, uniformalbedo, frac_ocean,frac_oldsnow,frac_newsnow,frac_desert,frac_tundra, 'noplot', 'diffuse')
albedo_dir_wav=ras.get_surface_albedo(wav_leftedges, wav_rightedges,solarzenithangle, albedoflag, uniformalbedo, frac_ocean,frac_oldsnow,frac_newsnow,frac_desert,frac_tundra, 'noplot', 'direct')

########################
###Load absorption and scattering cross-sections.
########################
#Initialize variables to hold scattering cross-sections
xc_tot_species_wav=np.zeros([N_species,N_wavelengths])
xc_scat_species_wav=np.zeros([N_species,N_wavelengths])
xc_abs_species_wav=np.zeros([N_species,N_wavelengths])

#Load in band-averaged cross-sections.  
#Function called returns xc in cm2/molecule in each band (total extinction, absorption, and rayleigh scattering). tot=abs+ray.

(xc_tot_species_wav[0,:], xc_abs_species_wav[0,:], xc_scat_species_wav[0,:])=css.compute_band_cross_section(wav_leftedges, wav_rightedges, 'n2')
(xc_tot_species_wav[1,:], xc_abs_species_wav[1,:], xc_scat_species_wav[1,:])=css.compute_band_cross_section(wav_leftedges, wav_rightedges, 'co2')
(xc_tot_species_wav[2,:], xc_abs_species_wav[2,:], xc_scat_species_wav[2,:])=css.compute_band_cross_section(wav_leftedges, wav_rightedges, 'h2o')
(xc_tot_species_wav[3,:], xc_abs_species_wav[3,:], xc_scat_species_wav[3,:])=css.compute_band_cross_section(wav_leftedges, wav_rightedges, 'ch4')
(xc_tot_species_wav[4,:], xc_abs_species_wav[4,:], xc_scat_species_wav[4,:])=css.compute_band_cross_section(wav_leftedges, wav_rightedges, 'so2')
(xc_tot_species_wav[5,:], xc_abs_species_wav[5,:], xc_scat_species_wav[5,:])=css.compute_band_cross_section(wav_leftedges, wav_rightedges, 'o2')
(xc_tot_species_wav[6,:], xc_abs_species_wav[6,:], xc_scat_species_wav[6,:])=css.compute_band_cross_section(wav_leftedges, wav_rightedges, 'o3')
(xc_tot_species_wav[7,:], xc_abs_species_wav[7,:], xc_scat_species_wav[7,:])=css.compute_band_cross_section(wav_leftedges, wav_rightedges, 'h2s')

#Note that in some cases, the predicted Rayleigh scattering cross-section exceeds that measured in laboratory studies. This means that 1) the single-scattering albedo w_0>1 in these cases and 2) the absorption cross-section is negative (both unphysical). This is accounted for in compute_optical_parameters, where w_0 is set to a maximum of 0.999. The only parameters used in the radiative transfer are w_0 and tau_tot, so the negative absorption cross-sections never come into play. 


#######Try Rugheimer cross-sections and scattering instead
######(xc_tot_species_wav[0,:], xc_abs_species_wav[0,:], xc_scat_species_wav[0,:])=css.get_rugheimer_xc(wav_leftedges, wav_rightedges, 'n2',0.78, 0.)
#######(xc_tot_species_wav[0,:], xc_abs_species_wav[0,:], xc_scat_species_wav[0,:])=css.get_rugheimer_xc(wav_leftedges, wav_rightedges, 'n2',0.89, 0.1)
######(xc_tot_species_wav[1,:], xc_abs_species_wav[1,:], xc_scat_species_wav[1,:])=css.get_rugheimer_xc(wav_leftedges, wav_rightedges, 'co2',0,0)
######(xc_tot_species_wav[2,:], xc_abs_species_wav[2,:], xc_scat_species_wav[2,:])=css.get_rugheimer_xc(wav_leftedges, wav_rightedges, 'h2o',0,0)
######(xc_tot_species_wav[4,:], xc_abs_species_wav[4,:], xc_scat_species_wav[4,:])=css.get_rugheimer_xc(wav_leftedges, wav_rightedges, 'so2',0,0)
######(xc_tot_species_wav[5,:], xc_abs_species_wav[5,:], xc_scat_species_wav[5,:])=css.get_rugheimer_xc(wav_leftedges, wav_rightedges, 'o2',0,0)
######(xc_tot_species_wav[6,:], xc_abs_species_wav[6,:], xc_scat_species_wav[6,:])=css.get_rugheimer_xc(wav_leftedges, wav_rightedges, 'o3',0,0)
#######Rugheimer has no CH4 or H2S so we keep our own, but get rid of the scattering formalism 
######xc_scat_species_wav[3,:]=xc_scat_species_wav[3,:]*0.0
######xc_scat_species_wav[7,:]=xc_scat_species_wav[7,:]*0.0

########################
###Get atmospheric layers column densities
########################

#extract integrated column densities for each layer in this atmosphere
n_z, t_z, p_z, columndensity_z, t_c=rts.get_atmospheric_profile(z_lower, z_upper, atmoprofilefile)

#Compute column densities
colden_species_z=np.zeros([N_species,N_layers])
colden_species_z[0,:]=columndensity_z*mr_n2
colden_species_z[1,:]=columndensity_z*mr_co2
colden_species_z[2,:]=columndensity_z*mr_h2o
colden_species_z[3,:]=columndensity_z*mr_ch4
colden_species_z[4,:]=columndensity_z*mr_so2
colden_species_z[5,:]=columndensity_z*mr_o2
colden_species_z[6,:]=columndensity_z*mr_o3
colden_species_z[7,:]=columndensity_z*mr_h2s


########################
###Compute atmospheric optical parameters required for two-stream code: tau_n, tau_c, w0, and g in each layer as a function of wavelength.
########################

#Call subfunction to extract composite values
tau_n_tot_z_wav, tau_c_tot_z_wav, w_0_z_wav, g_z_wav=rts.compute_optical_parameters(colden_species_z,xc_tot_species_wav, xc_scat_species_wav,0)
#Reminder: Toon et al define their albedo with tau=0 at the TOA, and it is computed along the zenith direction (so solar zenith angle is accounted for separately in the code).

########Code to set scattering or absorption limits
#######w_0_z_wav=w_0_z_wav*0.0+1.e-5#pure absorption limit
#######w_0_z_wav=w_0_z_wav*0.0+1.-1.e-12 #pure scattering limit

########################
###compute the flux via the two-stream approximation
########################
F_plus_tau0=np.zeros(np.shape(tau_n_tot_z_wav)) #F_plus evaluated at tau=0 for every layer n
F_plus_taumax=np.zeros(np.shape(tau_n_tot_z_wav))#F_plus evaluated at tau=tau_n[n] for every layer n
F_minus_tau0=np.zeros(np.shape(tau_n_tot_z_wav))#F_minus evaluated at tau=0 for every layer n
F_minus_taumax=np.zeros(np.shape(tau_n_tot_z_wav))#F_minus evaluated at tau=tau_n[n] for every layer n

F_net=np.zeros(np.shape(tau_n_tot_z_wav))#Net flux at the BASE of layer n. 
AMEAN=np.zeros(np.shape(tau_n_tot_z_wav))#AMEAN, 4*pi*mean intensity at the base of layer n. 
SS=np.zeros(np.shape(intensity_toa)) #This quantity is the SS quantity from twostr.f.
surface_intensity=np.zeros(np.shape(intensity_toa)) #an  estimate of the total amount of intensity received by a point at the surface of the planet. It is equal to the direct intensity plus F_[surface]/mu_1, i.e. the downward diffuse intensity at the surface

#Core loop over wavelength:
for ind in range(0,N_wavelengths):
	wavelength=wav_centers[ind] #width doesn't matter as wav primarily matters for BB which varies in smooth way.
	solar_input=intensity_toa[ind]/np.pi #this converts the TOA flux to the F_s in Toon et al. Recall pi*F_s=solar flux (really solar intensity) in that formalism.
	w_0=w_0_z_wav[:,ind]
	
	g=g_z_wav[:,ind]
	tau_n=tau_n_tot_z_wav[:,ind]
	albedo_dif=albedo_dif_wav[ind]
	albedo_dir=albedo_dir_wav[ind]
	
	thunk1, thunk2, thunk3, thunk4, thunk5, thunk6, thunk7=twostr.twostr_func(wavelength, solar_input, solarzenithangle, albedo_dif, albedo_dir,temp_ground, w_0, g, tau_n, t_c)
	
	F_plus_tau0[:,ind]=thunk1
	F_plus_taumax[:,ind]=thunk2
	F_minus_tau0[:, ind]=thunk3
	F_minus_taumax[:,ind]=thunk4
	
	F_net[:,ind]=thunk5
	AMEAN[:,ind]=thunk6
	
	SS[ind]=np.sqrt(thunk6[-1]*thunk6[-2])
	surface_intensity[ind]=thunk7
	
	
	fnetcol=thunk5
	med_fnet=np.median(fnetcol)
	fnetdev=np.max(np.abs(fnetcol-med_fnet))/(mu_0*intensity_toa[ind])

########################
###Compute the direct fluxes throughout the atmosphere, and the surface flux.
########################

direct_flux_z_wav=mu_0*intensity_toa*np.exp(-tau_c_tot_z_wav/mu_0) #Direct flux at the boundary of each layer. First layer=TOA, last layer=surface. See: Toon et al eqn 50, and recall: flux_toa=F_s*np.pi

surface_direct_flux=direct_flux_z_wav[-1,:] #Get the direct flux at the surface
surface_diffuse_flux=F_minus_taumax[-1,:] #Get downwelling diffuse flux at the bottom layer, i.e. the surface

surface_flux=surface_diffuse_flux+surface_direct_flux


########################
###Compute the surface intensity at the base of the atmosphere, and compare to what is reported by the code.
########################
direct_intensity_z_wav=intensity_toa*np.exp(-tau_c_tot_z_wav/mu_0) #Direct intensity at the boundary of each layer, with first layer=TOA and last layer=BOA. 

surface_direct_intensity=direct_intensity_z_wav[-1,:] #direct intensity at base of atmosphere

surface_diffuse_intensity=F_minus_taumax[-1,:]/mu_1 #diffuse intensity at base of atmosphere.

surface_intensity_2=surface_direct_intensity+surface_diffuse_intensity

###Check for consistency:
surf_int_diff=(surface_intensity-surface_intensity_2)/surface_intensity
print np.max(np.abs(surf_int_diff))

########################
###Check energy conservation
########################
incoming_flux_tot=np.sum(mu_0*intensity_toa)

outgoing_flux_tot=np.sum(F_plus_tau0[0,:])

if outgoing_flux_tot <= incoming_flux_tot:
	print 'Outgoing Flux<= Incoming Flux: Consistent with Energy Conservation'
if outgoing_flux_tot > incoming_flux_tot:
	print 'Outgoing Flux > Incoming Flux: Energy Conservation Violated DANGER DANGER DANGER'

########################
###Some diagnostics, useful for the paper:
########################

columndensity_z
print 'Total column density is (cm-2):', np.sum(columndensity_z)
print 'N2 column density is (cm-2):', np.sum(colden_species_z[0,:])
print 'CO2 column density is (cm-2):', np.sum(colden_species_z[1,:])
print 'H2O column density is (cm-2):', np.sum(colden_species_z[2,:])
print 'CH4 column density is (cm-2):', np.sum(colden_species_z[3,:])
print 'SO2 column density is (cm-2):', np.sum(colden_species_z[4,:])
print 'O2 column density is (cm-2):', np.sum(colden_species_z[5,:])
print 'O3 column density is (cm-2):', np.sum(colden_species_z[6,:])
print 'H2S column density is (cm-2):', np.sum(colden_species_z[7,:])

inds=np.where(surface_intensity/intensity_toa<0.01)
indmax=np.max(inds) #here I assume monotonically decreasing
print 'index at which intensity suppressed to 0.01 x incident is:', indmax
print 'wavelength at which intensity suppressed 0.01 x incident is (nm):', wav_leftedges[indmax]

##debugging
print 'Surface Albedo (500 nm)', albedo_dif_wav[-1] #how much flux is making its way back out from the planet
print 'Planetary Albedo (500 nm)', F_plus_tau0[0,-1]/(mu_0*intensity_toa[-1]) #how much flux is making its way back out from the planet
#print F_minus_taumax[-1,:]/(mu_0*intensity_toa) #how much flux is irradiating the planet surface
#print np.max(tau_n_tot_z_wav[:,-1])
#pdb.set_trace()

############################
#######Plot results
############################

fig, (ax1, ax2)=plt.subplots(2, figsize=(8,6), sharex=True)
ax1.plot(wav_centers, intensity_toa, marker='s', color='black', label='TOA Intensity')
ax1.plot(wav_centers, surface_intensity_basecase, marker='s', color='blue', label='Surface Intensity (Rugheimer Base Case)')
ax1.plot(wav_centers, surface_intensity , marker='s', color='orange', label='Surface Intensity (This Model)')
ax1.set_yscale('log')
ax1.set_ylim([1.e-4, 1.e4])
ax1.set_xlabel('nm')
ax1.set_ylabel('erg/s/cm2/nm')
ax1.legend(loc=0)

ax2.plot(wav_centers, (surface_intensity-surface_intensity_basecase)/intensity_toa, marker='s', color='orange', label='Fractional Difference')
ax2.set_yscale('linear')
#ax2.set_ylim([-0.006, 0.006])
ax2.set_xlim([100.,500.])
ax2.set_xlabel('nm')
ax2.set_ylabel('Fractional Difference')
####plt.savefig(plotfilename, orientation='portrait',papertype='letter', format='pdf')

###############################
#######Print spectra
############################
toprint=np.zeros([np.size(wav_centers), 9])
toprint[:,0]=wav_leftedges #left edge of wavelength bin (nm)
toprint[:,1]=wav_rightedges #right edge of wavelength bin (nm)
toprint[:,2]=wav_centers #center of wavelength bin (nm)
toprint[:,3]=intensity_toa #intensity incident at top of atmosphere (erg/s/cm2/nm)
toprint[:,4]=surface_flux #flux incident at bottom of atmosphere (erg/s/cm2/nm)
toprint[:,5]=SS #total intensity (4\pi J) in middle of bottom layer of atmosphere, same as what CLIMA code reports (erg/s/cm2/nm)
toprint[:,6]=surface_intensity #total intensity incident on surface. It is equal to sum of direct intensity and diffuse downward intensity. (erg/s/cm2/nm)
toprint[:,7]=surface_diffuse_intensity #total downward diffuse intensity at surface, i.e. 2\pi*I_minus[N]. (erg/s/cm2/nm)
toprint[:,8]=surface_direct_intensity #total direct intensity incident at surface, i.e. I_0*exp(-tau_0/mu_0). (erg/s/cm2/nm)


header='Left Bin Edge (nm)	Right Bin Edge (nm)	Bin Center (nm)		Top of Atm Intensity (erg/s/nm/cm2)		Total Surface Flux (erg/s/nm/cm2)		Total Intensity at BOA (erg/s/nm/cm2)		Total Surface Intensity (erg/s/nm/cm2)		Total Surface Diffuse Intensity (erg/s/nm/cm2)		Total Surface Direct Intensity (erg/s/nm/cm2)\n'
f=open(writefilename, 'w')
f.write(header)
np.savetxt(f, toprint, delimiter='		', fmt='%1.7e', newline='\n')
f.close()

plt.show()



"""
********************************************************************************
SCRATCH CODE
********************************************************************************
"""

#############################
########Print spectra for pure absorption case test
#############################
#toprint=np.zeros([np.size(wav_centers), 7])
#toprint[:,0]=wav_leftedges #left edge of wavelength bin (nm)
#toprint[:,1]=wav_rightedges #right edge of wavelength bin (nm)
#toprint[:,2]=wav_centers #center of wavelength bin (nm)
#toprint[:,3]=intensity_toa #intensity incident at top of atmosphere (erg/s/cm2/nm)
#toprint[:,4]=surface_flux #diffuse intensity at bottom of atmosphere (erg/s/cm2/nm)
#toprint[:,5]=surface_direct_flux #reference measurements
#toprint[:,6]=surface_diffuse_flux #reference measurements
#header='Left Bin Edge (nm)	Right Bin Edge (nm)	Bin Center (nm)		Top of Atm Intensity (erg/s/nm/cm2)		Total Surface Flux (erg/s/nm/cm2)		Total Intensity at BOA (erg/s/nm/cm2)		Total Surface Flux(erg/s/nm/cm2)		Total Surface Direct Flux (erg/s/nm/cm2)		Total Surface Diffuse Flux (erg/s/nm/cm2)\n'
#f=open(writefilename, 'w')
#f.write(header)
#np.savetxt(f, toprint, delimiter='		', fmt='%1.7e', newline='\n')
#f.close()

#################################
############Save net fluxes for pure scattering case test
#################################
##Save the net flux, useful for checking the scattering limit
#writefilename='./TwoStreamOutput/'+filename+'.p'
#fnetdict={'F_net':F_net, 'wav_leftedges':wav_leftedges,'wav_rightedges':wav_rightedges, 'wav_centers':wav_centers, 'z_lower':z_lower, 'z_upper':z_upper, 'z_center':z_center, 'flux_toa':intensity_toa, 'solarzenithangle':solarzenithangle}
#pickle.dump(fnetdict, open(writefilename, 'wb'))

#F_net_deviation=np.zeros(np.shape(F_net))
#F_net_deviation_max=np.zeros(N_wavelengths)
#F_net_deviation_max_normalized=np.zeros(N_wavelengths)
#F_net_deviation_stddevs=np.zeros(N_wavelengths)
#for ind in range(0, N_wavelengths):
	#median_val=np.median(F_net[:,ind])
	
	#F_net_deviation[:,ind]=F_net[:,ind]-median_val
	#F_net_deviation_max[ind]=np.max(np.abs(F_net_deviation[:,ind]))
	#F_net_deviation_stddevs[ind]=np.std(F_net[:,ind])

#F_net_deviation_max_normalized=F_net_deviation_max/(direct_flux_z_wav[0,:])
#F_net_deviation_stddevs_normalized=F_net_deviation_stddevs/(direct_flux_z_wav[0,:])

#print np.abs(F_net_deviation_max_normalized)
#print np.min(np.abs(F_net_deviation_max_normalized))
#print np.max(np.abs(F_net_deviation_max_normalized))
##print np.max(np.abs(F_net_deviation_stddevs_normalized))

#############################
########Print spectra for Wuttke measurement:
#############################
#toprint=np.zeros([np.size(wav_centers), 6])
#toprint[:,0]=wav_leftedges #left edge of wavelength bin (nm)
#toprint[:,1]=wav_rightedges #right edge of wavelength bin (nm)
#toprint[:,2]=wav_centers #center of wavelength bin (nm)
#toprint[:,3]=intensity_toa #intensity incident at top of atmosphere (erg/s/cm2/nm)
#toprint[:,4]=surface_diffuse_intensity #diffuse intensity at bottom of atmosphere (erg/s/cm2/nm)
#toprint[:,5]=surface_intensity_basecase #reference measurements
#header='Left Bin Edge (nm)	Right Bin Edge (nm)	Bin Center (nm)		Top of Atm Intensity (erg/s/nm/cm2)		Total Surface Flux (erg/s/nm/cm2)		Total Intensity at BOA (erg/s/nm/cm2)		Total Surface Intensity (erg/s/nm/cm2)		Total Surface Diffuse Intensity (erg/s/nm/cm2)		Total Surface Direct Intensity (erg/s/nm/cm2)\n'
#f=open(writefilename, 'w')
#f.write(header)
#np.savetxt(f, toprint, delimiter='		', fmt='%1.7e', newline='\n')
#f.close()

#############################
########Print spectra for WOUDC measurement
#############################
#toprint=np.zeros([np.size(wav_centers), 6])
#toprint[:,0]=wav_leftedges #left edge of wavelength bin (nm)
#toprint[:,1]=wav_rightedges #right edge of wavelength bin (nm)
#toprint[:,2]=wav_centers #center of wavelength bin (nm)
#toprint[:,3]=intensity_toa #intensity incident at top of atmosphere (erg/s/cm2/nm)
#toprint[:,4]=surface_flux #flux incident at bottom of atmosphere (erg/s/cm2/nm)
#toprint[:,5]=surface_intensity_basecase #reference measurements
#header='Left Bin Edge (nm)	Right Bin Edge (nm)	Bin Center (nm)		Top of Atm Intensity (erg/s/nm/cm2)		Total Surface Flux (erg/s/nm/cm2)		Total Intensity at BOA (erg/s/nm/cm2)		Total Surface Intensity (erg/s/nm/cm2)		Total Surface Diffuse Intensity (erg/s/nm/cm2)		Total Surface Direct Intensity (erg/s/nm/cm2)\n'
#f=open(writefilename, 'w')
#f.write(header)
#np.savetxt(f, toprint, delimiter='		', fmt='%1.7e', newline='\n')
#f.close()