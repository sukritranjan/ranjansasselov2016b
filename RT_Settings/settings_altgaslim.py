# -*- coding: iso-8859-1 -*-
#Copy and paste into radiativetransfer_v10.py and above to reproduce given dataset. 

########################
###User-set parameters
########################
gas='so2'
multiples=np.array([1.e-5, 1.e-4, 1.e-3, 1.e-2, 1.e-1, 1., 1.e1, 1.e2, 1.e3, 1.e4, 1.e5, 1.e6, 1.e7, 1.e8])
multiple=str(multiples[0])

filename='/gaslim/surface_intensities_'+gas+'limits_'+gas+'multiple='+multiple+'_a=freshsnow_z=0' #name of file to write output, plot to

#TOA input, and comparison file.
inputspectrafile='./LiteratureSpectra/general_spectral_input/' #TOA Stellar Input and Reference Spectra for Comparison

#Mixing ratio file
N_species=8 #how many gas species are in our model?
mr_profilefile='./MixingRatios/gaslim/rugheimer_earth_epoch0_mixingratios_'+gas+'limits_'+gas+'multiple='+multiple+'.dat' #Mixing ratio file

#T/P Profile.
atmoprofilefile='./TPProfiles/gaslim/rugheimer_earth_epoch0_atmosphereprofile_'+gas+'limits_'+gas+'multiple='+multiple+'.dat' #T/P Profile File. File boundaries should match z_lower and z_upper

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
temp_ground=292.5 #Temperature of the ground in K.
