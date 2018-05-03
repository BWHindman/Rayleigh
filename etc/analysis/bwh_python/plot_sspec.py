# VIII.  Spherical Harmonic Spectra
# ===================
# 
# **Summary:**    Spherical Harmonic Spectra sampled at discrete radii. 
# 
# **Subdirectory:**  Shell_Spectra
# 
# **main_input prefix:** shellspectra
# 
# **Python Classes:** 
# 
# * Shell_Spectra :  Complete data structure associated with Shell_Spectra outputs.
# * PowerSpectrum :  Reduced data structure -- contains power spectrum of velocity and/or magnetic fields only.
# 
# **Additional Namelist Variables:**  
# 
# * shellspectra_levels (indicial) : indices along radial grid at which to output spectra.
# 
# * shellspectra_levels_nrm (normalized) : normalized radial grid coordinates at which to output spectra.
# 
# 
# The shell-spectra output type allows us to examine the spherical harmonic decomposition of output variables at discrete radii.
# 
# Examining the *main_input* file, we see that the following output values have been denoted for the Shell Spectra (see *rayleigh_output_variables.pdf* for mathematical formulae):
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
# Spherical harmonic spectra can be read into Python using either the **Shell_Spectra** or **PowerSpectrum** classes.  
# 
# The **Shell_Spectra** class provides the full complex spectra, as a function of degree ell and azimuthal order m, for each specified output variable.   It possesses an attribute named *lpower* that contains the associated power for each variable, along with its m=0 contributions separated and removed.
# 
# The **Power_Spectrum** class can be used to read a Shell_Spectra file and quickly generate a velocity or magnetic-field power spectrum.   For this class to work correctly, your file must contain all three components of either the velocity or magnetic field.   Other variables are ignored (use Shell_Spectrum's lpower for those).
# 
# We illustrate how to use these two classes below.  As usual, we call the help() function to display the docstrings that describe the different data structures embodied by each class.

# In[23]:

from rayleigh_diagnostics import Shell_Spectra, Power_Spectrum, build_file_list
import numpy
import sys
import matplotlib.pyplot as plt
from matplotlib import ticker

nargs = len(sys.argv)
if (nargs==1):
   iter0 = 0
   iter1 = 99999999     # 100 million - 1
else:           # compile a list of files using the command line argments
   iter0 = int(sys.argv[1])
   if (nargs > 2):
      iter1 = int(sys.argv[2])
   else:
      iter1 = 99999999  # 100 million - 1

files = build_file_list(iter0,iter1,path='Shell_Spectra')
nfiles = len(files)
print(files)


# Average over time steps
rind = 0
ndata = 0
for i in range(len(files)):
    file = files[i]
    istring = file[14:len(file)]
    vpower = Power_Spectrum(istring)
    if (i==0):
       nell = vpower.lmax + 1
       nr = vpower.nr
       nspec = 3
       power = numpy.zeros((nell, nspec),dtype='float64')
       if (nfiles > 1 or nargs < 4):
          iter0 = 0
          iter1 = vpower.niter-1
       else:
          iter0 = int(sys.argv[3])-1
          iter1 = int(sys.argv[4])-1
    mn = min(iter1, vpower.niter-1) 
    print('Iterations %s through %s' % (vpower.iters[iter0], vpower.iters[mn]))
    for j0 in range(mn-iter0+1):
       j=j0+iter0
       power += vpower.power[0:nell,rind,j,:]
       ndata += 1.0
power = power*(1.0/ndata)


# Make plots
fig, ax = plt.subplots(nrows=3, figsize=(6,6))
ax[0].plot(power[:,0])
ax[0].set_xlabel(r'Degree $\ell$')
ax[0].set_xscale('log')
ax[0].set_yscale('log')
ax[0].set_title('Velocity Power (total)')


ax[1].plot(power[:,1])
ax[1].set_xlabel(r'Degree $\ell$')
ax[1].set_xscale('log')
ax[1].set_yscale('log')
ax[1].set_title('Velocity Power (m=0)')

ax[2].plot(power[:,2])
ax[2].set_xlabel(r'Degree $\ell$')
ax[2].set_xscale('log')
ax[2].set_yscale('log')
ax[2].set_title('Velocity Power ( total - {m=0} )')

plt.tight_layout()
savefile = 'Power_1D.pdf'
print('Saving figure to: '+ savefile)
plt.savefig(savefile)

fig, ax = plt.subplots()
ss = Shell_Spectra(istring)

mmax = ss.mmax
lmax = ss.lmax
tind = 0
power_spectrum = numpy.zeros((lmax+1,mmax+1),dtype='float64')

for i in range(1,4):   # i takes on values 1,2,3
    qind=ss.lut[i]
    complex_spectrum = ss.vals[:,:,rind,qind,tind]
    power_spectrum = power_spectrum+numpy.real(complex_spectrum)**2 + numpy.imag(complex_spectrum)**2

power_spectrum = numpy.transpose(power_spectrum)

tiny = 1e-6
img=ax.imshow(numpy.log10(power_spectrum+tiny), origin='lower')
ax.set_ylabel('Azimuthal Wavenumber m')
ax.set_xlabel(r'Degree $\ell$')
ax.set_title('Velocity Power Spectrum')

#colorbar ...
cbar = plt.colorbar(img) # ,shrink=0.5, aspect = 15)
cbar.set_label('Log Power')
        
tick_locator = ticker.MaxNLocator(nbins=5)
cbar.locator = tick_locator
cbar.update_ticks()
cbar.ax.tick_params()   #font size for the ticks


savefile = 'Power_2D.pdf'
print('Saving figure to: '+ savefile)
plt.savefig(savefile)
