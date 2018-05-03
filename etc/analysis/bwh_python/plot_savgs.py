# V.  Shell Averages
# ==========
# 
# **Summary:**   Spherical averages of requested output variables.  Each output variable is stored as a 1-D function of radius.
# 
# **Subdirectory:**  Shell_Avgs
# 
# **main_input prefix:** shellavg
# 
# **Python Class:** Shell_Avgs
# 
# **Additional Namelist Variables:**  
# None
# 
# The Shell-Averaged outputs are useful for examining how quantities vary as a function of radius.  They are particularly useful for examining the distribution of energy as a function of radius, or the heat flux balance established by the system.
# 
# Examining the *main_input* file, we see that the following output values have been denoted for the Shell Averages (see *rayleigh_output_variables.pdf* for mathematical formulae):
# 
# 
# | Menu Code  | Description |
# |------------|-------------|
# | 1          | Radial Velocity |
# | 2          | Theta Velocity |
# | 3          | Phi Velocity  |
# | 501        | Entropy |
# | 502        | Pressure |
# | 1433       | Volumetric Heating Flux|
# | 1455       | Enthalpy Flux|
# | 1470       | Radial Conductive Heat Flux |
# | 1923       | Radial Kinetic Energy Flux |
# | 1935       | Radial Viscous Flux |
# 
# 
# In the example that follows, we will plot the spherically-averaged velocity field as a function of radius, the mean temperature profile, and the radial heat flux.  We begin with a preamble similar to that used for the Global Averages.  Using the help function, we see that the Shell_Avgs data structure is similar to that of the G_Avgs.  There are three important differences:
# *  There is a radius attribute (necessary if we want to plot anything vs. radius)
# *  The dimensionality of the values array has changed;  radial index forms the first dimension.
# *  The second dimension of the values array has a length of 4.  In addition to the spherical mean, the 1st, 2nd and 3rd moments are stored in indices 0,1,2, and 3 respectively.

# In[8]:

from rayleigh_diagnostics import Shell_Avgs, build_file_list
import matplotlib.pyplot as plt
import numpy
import sys

# Build a list of all files ranging from iteration 0 million to 1 million
nargs = len(sys.argv)
if (nargs==1):
   iter0 = 0
   iter1 = 99999999	# 100 million - 1
else:		# compile a list of files using the command line argments
   iter0 = int(sys.argv[1])
   if (nargs > 2):
      iter1 = int(sys.argv[2])
   else:
      iter1 = 99999999	# 100 million - 1

files = build_file_list(iter0,iter1,path='Shell_Avgs')
print(files)


# 
# ***
# 
# While it can be useful to look at instaneous snapshots of Shell Averages, it's often useful to examine these outputs in a time-averaged sense.    Let's average of all 200 snapshots in the last file that was output.  We could average over data from multiple files, but since the benchmark run achieves a nearly steady state, a single file will do in this case.

# In[9]:

nfiles = len(files)


ndata = 0.0
nmom = 4
for i in range(nfiles):
    a = Shell_Avgs(filename=files[i], path='')
    if (i==0):
        nr = a.nr
        nq = a.nq
        radius = a.radius
        savg=numpy.zeros((nr,nmom,nq),dtype='float64')
        if (nfiles > 1 or nargs < 4):
           iter0 = 0
	   iter1 = a.niter-1
	else:
	   iter0 = int(sys.argv[3])-1
	   iter1 = int(sys.argv[4])-1
    mn = min(iter1,a.niter-1)
    print("Iterations %s through %s" % (a.iters[iter0], a.iters[mn]))
    for j0 in range(mn-iter0+1):
	j=j0+iter0
	savg[:,:,:] += a.vals[:,:,:,j]
	ndata += 1.0
savg = savg*(1.0/ndata)


lut = a.lut
vr = lut[1]         # Radial Velocity
vtheta = lut[2]     # Theta Velocity
vphi = lut[3]       # Phi Velocity
entropy = lut[501]  # Entropy
pressure = lut[502] # Pressure


eflux = lut[1455]  # Enthalpy Flux (radial)
cflux = lut[1470]  # Conductive Heat Flux (radial)
keflux = lut[1923] # Kinetic Energy Flux (radial)
vcflux = lut[1935] # Viscous Energy Flux (radial)
vhflux = lut[1433] # Volumetric Heating Flux


# Velocity vs. Radius
# ---------------------
# Next, we plot the mean velocity field, and its first moment, as a function of radius.   Notice that the radial and theta velocity components have a zero spherical mean.  Since we are running an incompressible model, this is a good sign!

# In[10]:

sizetuple = (10,3)
fig, ax = plt.subplots(ncols=2, figsize=sizetuple)

ax[0].plot(radius,savg[:,0,vr],label=r'$v_r$')
ax[0].plot(radius,savg[:,0,vtheta], label=r'$v_\theta$')
ax[0].plot(radius,savg[:,0,vphi], label=r'$v_\phi$')
ax[0].legend(shadow=True,loc='lower right')
ax[0].set_xlabel('Radius')
ax[0].set_ylabel('Velocity')
ax[0].set_title('Spherically-Average')

ax[1].plot(radius,savg[:,1,vr],label=r'$v_r$')
ax[1].plot(radius,savg[:,1,vtheta], label=r'$v_\theta$')
ax[1].plot(radius,savg[:,1,vphi], label=r'$v_\phi$')
ax[1].legend(shadow=True,loc='upper left')
ax[1].set_xlabel('Radius')
ax[1].set_ylabel('Velocity')
ax[1].set_title('First Spherical Moment')


plt.tight_layout()
savefile = 'velocity_variation.pdf'
print('Saving figure to: '+ savefile)
plt.savefig(savefile)



# Radial Entropy Profile
# ------------------------------
# We might also look at temperature ...

# In[11]:

fig, ax = plt.subplots(ncols=2,figsize=sizetuple)

ax[0].plot(radius,savg[:,0,entropy])
#ax[0].legend(shadow=True,loc='upper left')
ax[0].set_xlabel('Radius')
ax[0].set_ylabel('Mean Entropy')
#ax[0].set_title('Entropy Profile')

ax[1].plot(radius,numpy.sqrt(savg[:,1,entropy]))
#ax[1].legend(shadow=True,loc='upper right')
ax[1].set_xlabel('Radius')
ax[1].set_ylabel('Entropy Fluctuation')
#ax[1].set_title('Entropy Variance')

savefile = 'entropy_variation.pdf'
plt.tight_layout()
print('Saving figure to: '+ savefile)
plt.savefig(savefile)


# Heat Flux Contributions
# --------------------------
# We can also examine the balance between convective and conductive heat flux.  In this case, before plotting these quantities as a function of radius, we normalize them by the surface area of the sphere to form a luminosity.

# In[12]:

fpr=4.0*numpy.pi*radius*radius
elum = savg[:,0,eflux]*fpr
clum = savg[:,0,cflux]*fpr
klum = savg[:,0,keflux]*fpr
vlum = savg[:,0,vcflux]*fpr
hlum = savg[:,0,vhflux]*fpr
tlum = elum+clum+klum+vlum+hlum
fig, ax = plt.subplots()
ax.plot(radius,elum,label='Convection')
ax.plot(radius,clum, label='Conduction')
ax.plot(radius,klum, label='Kinetic Energy')
ax.plot(radius,vlum, label='Viscous')
ax.plot(radius,hlum, label='Radiative')
ax.plot(radius,tlum, label='Total')
ax.set_title('Flux Balance')
ax.set_ylabel(r'Energy Flux ($\times 4\pi r^2$)')
ax.set_xlabel('Radius')
ax.legend(shadow=True,loc='upper left')
plt.tight_layout()
savefile = 'energy_flux.pdf'
print('Saving figure to: '+ savefile)
plt.savefig(savefile)

