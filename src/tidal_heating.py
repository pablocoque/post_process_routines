import numpy as np
import matplotlib.pyplot as plt
import pynbody

# Load simulation data
simulation = 'follow_last'
s = pynbody.load(simulation+'/snapdir_061/snapshot_061')

# Load shock evolution data
data = np.loadtxt(simulation+'/oneclustershocks.txt')

# Extract data
time = data[:,0]
tensor = data[:,1:-4]
tidal_heating = data[:,-4]
shock_duration = data[:,-3]
Awj = data[:,-2]
time_since_last_shock = data[:,-1]

# Convert time column to Gyr
time_gyr = pynbody.analysis.cosmology.age(s,z = 1./time - 1.)

# Convert tensor data to Gyr^-2
gyr = gyr = 3.15576e16 # in s
unit_time = 3.08568e19
h = 0.6777
convert_factor = h**2/time**3 * (gyr/unit_time)**2 # to Gyr^-2

tensor *= convert_factor[:, np.newaxis]
tidal_heating *= convert_factor

# Convert shock duration to Gyr
shock_duration *= unit_time/h /gyr
time_since_last_shock *= unit_time/h /gyr

# Plot
fig, ax = plt.subplots(2,1, sharex=True, gridspec_kw={'hspace': 0}, figsize=(30,10))
for i in range(tensor.shape[1]):
    ax[0].plot(time_gyr, tensor[:,i])
ax[0].set(ylabel=r'$T_{ij}$ components [Gyr$^{-2}$]')
ax[0].set_ylim(-40000,10000)
for i,val in enumerate(tidal_heating):
    if val>0:
        ax[1].plot(time_gyr[i], val, 'x', c='k')
ax[1].set(xlabel='Time [Gyr]', ylabel=r'$I_{\rm{tid}}$ [Gyr$^{-2}$]')
ax[1].set_xlim(4,5)
ax[1].set_ylim(0,15)
plt.tight_layout()
plt.savefig(simulation+'_tidal_heating.pdf')