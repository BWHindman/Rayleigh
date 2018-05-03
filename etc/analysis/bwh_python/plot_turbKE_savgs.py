
#Codes from Diagnostics_TurbKE_budget

# | Menu Code   | Description 					|
# |-------------|-----------------------------------------------|
# | ----	| Turbulent KE density				|
# |-------------|-----------------------------------------------|
# | 2301	| Buoyant production of turbulent KE		|
# | ----	| Shear production of turbulent KE		|
# | 2303	| Viscous dissipation of turbulent KE		|
# |-------------|-----------------------------------------------|
# | 2304	| Pressure transport of turbulent KE 		|
# | 2305	| Viscous transport of turbulent KE		|
# | 2306	| Turbulent advective transport of turbulent KE	|
# | 2307	| Mean advective transport of turbulent KE	|
# |-------------|-----------------------------------------------|
# | 2308 (-)	| Radial flux of turb. KE due to press. trans.	|
# | 2309 (-)	| Radial flux of turb. KE due to viscous trans.	|
# | 2310 (-)	| Radial flux of turb. KE due to turb. advect.	|
# | 2311 (-)	| Radial flux of turb. KE due to mean advect.	|
# |-------------|-----------------------------------------------|
# | 2312 (-)	| Theta flux of turb. KE due to press. trans.	|
# | 2313 (-)	| Theta flux of turb. KE due to viscous trans.	|
# | 2314 (-)	| Theta flux of turb. KE due to turb. advect.	|
# | 2315 (-)	| Theta flux of turb. KE due to mean advect.	|
# |-------------|-----------------------------------------------|


#Codes from the standard outputs

# | Menu Code   | Description 					|
# |-------------|-----------------------------------------------|
# | 409		| Turbulent KE density				|
# |-------------|-----------------------------------------------|
# | 1905	| Buoyant production of turbulent KE		|
# | 1954	| Shear production of turbulent KE		|
# | ----	| Viscous dissipation of turbulent KE		|
# |-------------|-----------------------------------------------|
# | 1902	| Pressure transport of turbulent KE 		|
# | 1908 + 2303	| Viscous transport of turbulent KE		|
# | 1911 (-)	| Turbulent advective transport of turbulent KE	|
# | ----	| Mean advective transport of turbulent KE	|
# |-------------|-----------------------------------------------|
# | 1947 (-)	| Radial flux of turb. KE due to press. trans.	|
# | 1938 (-)	| Radial flux of turb. KE due to viscous trans.	|
# | 1932	| Radial flux of turb. KE due to turb. advect.	|
# | 1929	| Radial flux of turb. KE due to mean advect.	|
# |-------------|-----------------------------------------------|
# | 1948 (-)	| Theta flux of turb. KE due to press. trans.	|
# | 1939 (-)	| Theta flux of turb. KE due to viscous trans.	|
# | 1933	| Theta flux of turb. KE due to turb. advect.	|
# | 1930	| Theta flux of turb. KE due to mean advect.	|
# |-------------|-----------------------------------------------|



from rayleigh_diagnostics import Shell_Avgs, build_file_list
import matplotlib.pyplot as plt
import numpy
import sys


#################################################################################
#	Generate a list of files to average over

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
nfiles = len(files)
print(files)


#################################################################################
#	Perform the averaging

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
	   iter0 = int(sys.argv[3])-1		# Set the minimum and maximum iteration from the command line arguments
	   iter1 = int(sys.argv[4])-1
    mn = min(iter1,a.niter-1)
    print("Iterations %s through %s" % (a.iters[iter0], a.iters[mn]))
    for j0 in range(mn-iter0+1):
	j=j0+iter0
	savg[:,:,:] += a.vals[:,:,:,j]
	ndata += 1.0
savg = savg*(1.0/ndata)


#################################################################################
#	Extract the variables

lut = a.lut

# Turbulent kinetic energy diagnostics
p_buoy1 = lut[2301]	# buoyant production (BWH)
p_buoy2 = lut[1905]	# buoyant production (standard)

#p_shear1 = ???
p_shear2 = lut[1954]	# shear production (standard)

d_visc1 = lut[2303]	# viscous dissipation (BWH)
#d_visc2 = 

t_press1 = lut[2304]	# pressure transport (BWH)
t_press2 = lut[1902]	# pressure transport (standard)

t_visc1 = lut[2305]	# viscous transport (BWH)
thalf_visc2 = lut[1908]	# One piece of the viscous transport (standard): must be coupled with d_visc1 

t_tadv1 = lut[2306]	# turbulent advective transport (BWH)
tneg_tadv2 = lut[1911] 	# turbulent adventice transport (standard): need a minus sign

t_madv1 = lut[2307]	# mean advective transport (BWH)
#t_madv2 = 

fneg_press1 = lut[2308]	# pressure flux (BWH): need a minus sign
fneg_press2 = lut[1947] # pressure flux (standard): need a minus sign

fneg_visc1 = lut[2309]	# viscous flux (BWH): need a minus sign
fneg_visc2 = lut[1938]	# viscous flux (standard): need a minus sign

fneg_tadv1 = lut[2310]	# turbulent advective flux (BWH): need a minus sign
f_tadv2 = lut[1932]	# turbulent advective flux (standard)

fneg_madv1 = lut[2311]	# mean advective flux (BWH): need a minus sign
f_madv2 = lut[1929]	# mean advective flux (standard)


sizetuple = (10,4)
fig, ax = plt.subplots(ncols=2, figsize=sizetuple)

# BWH's code
shear_prod = savg[:,0,p_shear2]
buoy_prod = savg[:,0,p_buoy1]
visc_diss = -savg[:,0,d_visc1]
press_tran = savg[:,0,t_press1]
visc_tran = savg[:,0,t_visc1]
tadv_tran = savg[:,0,t_tadv1]
madv_tran = savg[:,0,t_madv1]
tot = shear_prod + buoy_prod + visc_diss + press_tran + tadv_tran + madv_tran + visc_tran

buoy_prod2 = savg[:,0,p_buoy2]
#shear_prod2 = savg[:,0,p_shear2]
press_tran2 = savg[:,0,t_press2]
visc_tran2 = savg[:,0,thalf_visc2] - visc_diss
tadv_tran2 = -savg[:,0,tneg_tadv2]
#madv_tran2 = 


ax[0].plot(radius,buoy_prod,label=r'$P_B$')
ax[0].plot(radius,shear_prod,label=r'$P_S$')
ax[0].plot(radius,visc_diss,label=r'$\Phi_T$')
ax[0].plot(radius,press_tran,label=r'$T_p$')
ax[0].plot(radius,visc_tran,label=r'$T_v$')
ax[0].plot(radius,tadv_tran,label=r'$T_T$')
ax[0].plot(radius,madv_tran,label=r'$T_M$')
#ax[0].plot(radius,tot,label=r'$Total$')

ax[0].plot(radius,buoy_prod2,linestyle='-.')
ax[0].plot(radius,press_tran2,linestyle='-.')
ax[0].plot(radius,visc_tran2,linestyle='-.')
ax[0].plot(radius,tadv_tran2,linestyle='-.')

#ax[0].legend(shadow=True,loc='lower right')
ax[0].set_xlabel('Radius')
ax[0].set_ylabel('Turbulent KE Budget')
ax[0].set_ylim([-20.0,20.0])
ax[0].set_title('Spherically-Averaged')


# Viscous Terms
ax[1].plot(radius,buoy_prod,label=r'$P_B$')
ax[1].plot(radius,shear_prod,label=r'$P_S$')
ax[1].plot(radius,visc_diss,label=r'$\Phi_T$')
ax[1].plot(radius,press_tran,label=r'$T_p$')
ax[1].plot(radius,visc_tran,label=r'$T_v$')
ax[1].plot(radius,tadv_tran,label=r'$T_T$')
ax[1].plot(radius,madv_tran,label=r'$T_M$')
#ax[1].plot(radius,tot,label=r'$Total$')

ax[1].plot(radius,buoy_prod2,linestyle='-.')
ax[1].plot(radius,press_tran2,linestyle='-.')
ax[1].plot(radius,visc_tran2,linestyle='-.')
ax[1].plot(radius,tadv_tran2,linestyle='-.')

ax[1].legend(shadow=True,loc='upper center')
ax[1].set_xlabel('Radius')
ax[1].set_ylabel('Oversized Terms')
ax[1].set_title('Spherically-Averaged')

plt.tight_layout()
savefile = 'turbKE_budget.pdf'
print('Saving figure to: '+ savefile)
plt.savefig(savefile)



# Flux Contributions
# --------------------------

fpr=4.0*numpy.pi*radius*radius
plum = -savg[:,0,fneg_press1]*fpr
vlum = -savg[:,0,fneg_visc1]*fpr
talum = -savg[:,0,fneg_tadv1]*fpr
malum = -savg[:,0,fneg_madv1]*fpr
totlum  = plum + vlum + talum + malum

plum2 = -savg[:,0,fneg_press2]*fpr
vlum2 = -savg[:,0,fneg_visc2]*fpr
talum2 = savg[:,0,f_tadv2]*fpr
malum2 = savg[:,0,f_madv2]*fpr
totlum2  = plum2 + vlum2 + talum2 + malum2

fig, ax = plt.subplots()
ax.plot(radius,plum,label='Pressure')
ax.plot(radius,plum2,linestyle = '-.')

ax.plot(radius,vlum, label='Viscous')
ax.plot(radius,vlum2,linestyle = '-.')

ax.plot(radius,talum, label='Turbulent')
ax.plot(radius,talum2,linestyle = '-.')

ax.plot(radius,malum, label='Advective')
ax.plot(radius,malum2,linestyle = '-.')

ax.plot(radius,totlum, label='Total')
ax.plot(radius,totlum2,linestyle = '-.')

ax.set_title('Flux Balance')
ax.set_ylabel(r'Energy Flux ($\times 4\pi r^2$)')
ax.set_xlabel('Radius')
ax.legend(shadow=True,loc='lower right')

plt.tight_layout()
savefile = 'turbKE_fluxes.pdf'
print('Saving figure to: '+ savefile)
plt.savefig(savefile)


