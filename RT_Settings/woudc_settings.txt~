########################
###User-set parameters
########################
filename='reproduce_woudc' #name of file to write output, plot to

#TOA input, and comparison file.
inputspectrafile='./LiteratureSpectra/woudc.dat' #TOA Stellar Input and Reference Spectra for Comparison

#Mixing ratio file
N_species=8 #how many gas species are in our model?
mr_profilefile='./MixingRatios/woudc_mixingratios_v2.dat' #Mixing ratio file

#T/P Profile.
atmoprofilefile='./TPProfiles/woudc_atmosphereprofile.dat' #T/P Profile File. File boundaries should match z_lower and z_upper

#RT model layers
z_upper_limit=60.e5  #Upper edge of the atmospheric layers, in cm. The bottom edge of the layers is assumed to be at 0.
z_step=60.e3 # Thickness of the atmospheric layers, in cm.

###Solar Zenith Angle
solarzenithangle=20.376*np.pi/180. #60 degrees (pi/3) used by Rugheimer et al in their paper.

###Albedo
#Uniform albedo: set flag to 'uniformalbedo' and set the subsequent albedo value
#Nonuniform albedo: set flag to 'nonuniformalbedo' and choose what mix of surface: ocean, old snow, new snow, desert, tundra. 
#set flag to 'plot' to plot albedo used. 
albedoflag='uniformalbedo' #value of 'nonuniformalbedo' or 'uniformalbedo'
uniformalbedo=0.04 #if adopting uniform albedo: what value?
frac_ocean=0. #if nonuniform albedo: what fraction of ground is ocean?
frac_tundra=0. #if nonuniform albedo: what fraction of ground is tundra?
frac_desert=0. #if nonuniform albedo: what fraction of ground is desert?
frac_oldsnow=0.#if nonuniform albedo: what fraction of ground is old snow?
frac_newsnow=0. #if nonuniform albedo: what fraction of ground is new snow?

###Temperature of Ground
temp_ground=288. #Temperature of the ground in K.


###REmember to uncomment the different print statement at the bottom of the file!