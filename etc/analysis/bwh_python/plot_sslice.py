# VII.3  Shell Slices
# --------------------------
# 
# **Summary:**    2-D, spherical profiles of selected output variables sampled in at discrete radii. 
# 
# **Subdirectory:**  Shell_Slices
# 
# **main_input prefix:** shellslice
# 
# **Python Class:** Shell_Slices
# 
# **Additional Namelist Variables:**  
# 
# * shellslice_levels (indicial) : indices along radial grid at which to output spherical surfaces.
# 
# * shellslice_levels_nrm (normalized) : normalized radial grid coordinates at which to output spherical surfaces.
# 
# 
# The shell-slice output type allows us to examine how the fluid properties vary on spherical surfaces.
# 
# Examining the *main_input* file, we see that the following output values have been denoted for the Shell Slices (see *rayleigh_output_variables.pdf* for mathematical formulae):
# 
# 
# | Menu Code  | Description |
# |------------|-------------|
# | 1          | Radial Velocity |
# | 2          | Theta Velocity |
# | 3          | Phi Velocity  |
# 
# 
# 
# 
# In the examples that follow, we demonstrate how to create a 2-D plot of radial velocity:
# * on a Cartesian, lat-lon grid
# * projected onto a spherical surface using Basemap   (MUST set use Basemap=True below)
# 
# 
# 
# Plotting on a lat-lon grid is straightforward and illustrated below.   The shell-slice data structure is also displayed via the help() function in the example below and contains information needed to define the spherical grid for plotting purposes.
# 
# It is worth noting the *slice_spec* keyword (described in the docstring) that can be passed to the init method.  When reading large shell slices, a user can save time and memory during the read process by specifying the slice they want to read.

# In[21]:

#####################################
#  Shell Slice
from rayleigh_diagnostics import Shell_Slices, build_file_list
import numpy
import sys
import matplotlib.pyplot as plt
from matplotlib import ticker, font_manager
# Read the data

if (len(sys.argv) > 1):		# Read the file specfied by the first command line argument
   file = '00000000'+sys.argv[1]
   file = 'Shell_Slices/'+file[-8:]
else:				# Read the last file
   files = build_file_list(0,99999999,path='Shell_Slices')
   file = files[-1]

ss = Shell_Slices(file, path = '')

ntheta = ss.ntheta
nphi = ss.nphi
costheta = ss.costheta
theta = numpy.arccos(costheta)


tindex =1 # All example quantities were output with same cadence.  Grab second time-index from all.
tindex = -1
sizetuple=(10,8)

print('Plotting iteration : %s' % ss.iters[tindex])

fig, ax = plt.subplots(ncols=1, nrows=3, figsize=sizetuple)

vr = ss.vals[:,:,0,ss.lut[1],tindex]
#img = plt.imshow(numpy.transpose(vr), extent=[0,360,-90,90])
ax[0].imshow(numpy.transpose(vr), extent=[0,360,-90,90], interpolation='none') 
ax[0].set_xlabel( 'Longitude')
ax[0].set_ylabel( 'Latitude')
ax[0].set_title(  'Radial Velocity')

vr = ss.vals[:,:,1,ss.lut[1],tindex]
ax[1].imshow(numpy.transpose(vr), extent=[0,360,-90,90], interpolation='none') 
ax[1].set_xlabel( 'Longitude')
ax[1].set_ylabel( 'Latitude')
ax[1].set_title(  'Radial Velocity')

vr = ss.vals[:,:,ss.nr-1,ss.lut[1],tindex]
ax[2].imshow(numpy.transpose(vr), extent=[0,360,-90,90], interpolation='none') 
ax[2].set_xlabel( 'Longitude')
ax[2].set_ylabel( 'Latitude')
ax[2].set_title(  'Radial Velocity')

plt.tight_layout()
savefile = 'Shell_Slices_LatLon.pdf'
print('Saving figure to: '+ savefile)
plt.savefig(savefile)


