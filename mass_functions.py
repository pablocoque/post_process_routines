import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.interpolate import interpn
import pynbody

def density_scatter( x , y, ax = None, sort = True, bins = 20, **kwargs )   :
    """
    Scatter plot colored by 2d histogram
    """
    if ax is None :
        fig , ax = plt.subplots()
    data , x_e, y_e = np.histogram2d( x, y, bins = bins, density=False)
    z = interpn( ( 0.5*(x_e[1:] + x_e[:-1]) , 0.5*(y_e[1:]+y_e[:-1]) ) , data , np.vstack([x,y]).T , method = "splinef2d", bounds_error = False)

    #To be sure to plot all data
    z[np.where(np.isnan(z))] = 0.0

    # Sort the points by density, so that the densest points are plotted last
    if sort :
        idx = z.argsort()
        x, y, z = x[idx], y[idx], z[idx]

    h = ax.scatter( x, y, c=z, norm='log', **kwargs )

    return ax, h

# Load the simulation
simulation = 'rlx_fixed/'
s = pynbody.load(simulation+'snapdir_127/snapshot_127')

# Set the softening length from the simulation
s['eps'] = s['smooth']

# Load the halo catalog and center to the largest halo
h = s.halos()
transform = pynbody.analysis.angmom.faceon(h[0])

# Change units to physical
s.physical_units()

# Load stars that contain clusters (From main galaxy only)
mask_clusters_initial = h[1].s['InitialNumberOfClusters'] > 0
mask_clusters_final = h[1].s['NumberOfClusters'] > 0

truncation_mass = h[1].s['ICMFTruncationMass'][mask_clusters_initial]
cluster_masses = h[1].s['ClusterMass'][mask_clusters_initial].flatten()
init_cluster_masses = h[1].s['InitialClusterMass'][mask_clusters_initial].flatten()
masslostrelax = h[1].s['MassLostRelaxation'][mask_clusters_initial].flatten()
masslostshock = h[1].s['MassLostShocks'][mask_clusters_initial].flatten()

# Remove clusters that have no mass
not_empty_clusters = (init_cluster_masses>0)

cluster_masses = cluster_masses[not_empty_clusters]
init_cluster_masses = init_cluster_masses[not_empty_clusters]
masslostrelax = masslostrelax[not_empty_clusters]
masslostshock = masslostshock[not_empty_clusters]

# Simple log info
print('Maximum initial mass {:2.2e} Msun'.format(init_cluster_masses.max()))
print('Maximum present-day mass {:2.2e} Msun'.format(cluster_masses.max()))

# Extract clusters age and disruption time
clusters_disruptiontime = h[1].s['DisruptionTime'][mask_clusters_initial].flatten()
clusters_disruptiontime = clusters_disruptiontime[not_empty_clusters]
disruption_gyr = np.nan_to_num(pynbody.analysis.cosmology.age(s,z = 1./clusters_disruptiontime - 1.), nan=0.)

clusters_birthtime = []
clusters_age = []
for i, nclt in enumerate(h[1].s['InitialNumberOfClusters'][mask_clusters_initial]):
  clusters_birthtime.append(np.ones(nclt) * h[1].s['tform'].in_units('Gyr')[mask_clusters_initial][i])
  clusters_age.append(np.ones(nclt) * h[1].s['age'].in_units('Gyr')[mask_clusters_initial][i])

clusters_birthtime = np.concatenate(clusters_birthtime)
clusters_age = np.concatenate(clusters_age)

# calculate lifetime: not disrupted clusters have a lifetime of 16 Gyr
lifetime_gyr = disruption_gyr - clusters_birthtime
lifetime_gyr[np.logical_not(disruption_gyr>0)] = 16

mask_problematic = h[1].s['InitialMassFractionInClusters']>1
print('Problematic stars with clusters {:4d}'.format(mask_problematic.sum()))
print('That`s {:.3f}% of stars with clusters'.format(mask_problematic.sum()/mask_clusters_initial.sum() * 100))

# Ration of actual number of clusters
print('Ratio of survaving clusters ', (h[1].s['NumberOfClusters']*mask_clusters_final).sum()/(h[1].s['InitialNumberOfClusters']*mask_clusters_initial).sum())

# Define mass functions ranges
min_icmf = 5e3 # in solar masses
max_icmf = 1e8 # in solar masses
marray_icmf = np.logspace(np.log10(min_icmf), np.log10(max_icmf), 25)

min_gcmf = 1e2 # in solar masses
max_gcmf = 1e8 # in solar masses
marray_gcmf = np.logspace(np.log10(min_gcmf), np.log10(max_gcmf), 25)

# Calculate fraction of disrupted clusters
fraction = []
for i in range(24):
  mask = (init_cluster_masses > marray_icmf[i]) & (init_cluster_masses < marray_icmf[i+1])
  mask_disrupted = (clusters_disruptiontime[mask]>0)
  fraction.append(mask_disrupted.sum()/mask.sum())

fraction = np.array(fraction)
fraction = np.nan_to_num(fraction, nan=0.)

# And fraction of surviving mass in clusters
massfraction = []
for i in range(24):
  mask = (init_cluster_masses > marray_icmf[i]) & (init_cluster_masses < marray_icmf[i+1])
  massfraction.append(cluster_masses[mask].sum()/init_cluster_masses[mask].sum())

massfraction = np.array(massfraction)
massfraction = np.nan_to_num(massfraction, nan=0.)

# Calculate mass lost due to relaxation and shocks
mevminit = masslostrelax/init_cluster_masses
mshminit = masslostshock/init_cluster_masses

# Plot density scatter of mass lost due to different processes and density scatter of initial cluster mass vs disruption time
fig, ax = plt.subplots(1, 3, figsize=(16,5))
ax[0], cb = density_scatter(init_cluster_masses, mshminit, ax=ax[0],bins=[marray_icmf, np.linspace(0,1.1,100)], edgecolors='face', s=5)
ax[0].set(xscale='log', xlim=(5e3, 5e7), ylim=(0,1), xlabel=r'$\rm{Initial}$ $\rm{Cluster}$ $\rm{Mass}$ [$M_\odot$]', ylabel=r'$\Delta m_{\rm{sh}}/m_{\rm{init}}$')

ax[1], cb = density_scatter(init_cluster_masses, mevminit, ax=ax[1],bins=[marray_icmf, np.linspace(0,1.1,100)], edgecolors='face', s=5)
ax[1].set(xscale='log', xlim=(5e3, 5e7), ylim=(0,1), xlabel=r'$\rm{Initial}$ $\rm{Cluster}$ $\rm{Mass}$ [$M_\odot$]', ylabel=r'$\Delta m_{\rm{ev}}/m_{\rm{init}}$')

ax[2], cb = density_scatter(init_cluster_masses, lifetime_gyr, ax=ax[2], bins=[marray_icmf, np.logspace(-2,1.5,100)], edgecolors='face', s=5)
ax[2].set(xscale='log', yscale='log', xlim=(5e3, 1e7), ylim=(1e-2, 2.3e1))
ax[2].set(xlabel=r'$\rm{Initial}$ $\rm{Cluster}$ $\rm{Mass}$ [$M_\odot$]', ylabel=r'$t_{\rm{disruption}}$ [$\rm{Gyr}$]')
ax2 = ax[2].twinx()
ax2.plot((marray_icmf[:-1] + marray_icmf[1:])/2, fraction, color='r', label='Number fraction')
ax2.plot((marray_icmf[:-1] + marray_icmf[1:])/2, massfraction, color='b', label='Mass fraction')
ax2.set_ylim(0, 1)
ax2.set_ylabel('Fraction of disrupted clusters/surviving mass')
ax2.legend()
cbar = fig.colorbar(cb, ax=ax)
cbar.ax.set_ylabel('Number')
plt.tight_layout()
plt.savefig('disruption_plots.pdf')
print('Disruption plots saved')

# Plot mass functions comparison
# Load observational data
MW_gcmf = np.loadtxt('MW_GCMF.csv', delimiter=',')
M31_gcmf = np.loadtxt('M31_GCMF.csv', delimiter=',')

# Load other simulations data
empfiducial_gcmf = np.loadtxt('Fiducial_emppathfinder.csv', delimiter=',')
empschecter_gcmf = np.loadtxt('Schechter_emppathfinder.csv', delimiter=',')
emosaics_gcmf = np.loadtxt('emosaics.csv', delimiter=',')

# Plot mass function for old cluster population (t>10 Gyr)
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
old_mask = (clusters_age>10)
N, _ = np.histogram(np.log10(init_cluster_masses[old_mask]), bins=np.log10(marray_icmf))#, histtype='step', label='ICMF')
ax.plot((marray_icmf[:-1] + marray_icmf[1:])/2, N, label='ICMF')
N, _ = np.histogram(np.log10(cluster_masses[old_mask]), bins=np.log10(marray_gcmf))#, histtype='step', label='ICMF')
ax.plot((marray_gcmf[:-1] + marray_gcmf[1:])/2, N, label='GCMF')
# ax.hist(init_cluster_masses, bins=marray_icmf, histtype='step', label='ICMF')
# ax.hist(cluster_masses, bins=marray_gcmf, histtype='step', label='GCMF')
ax.plot(MW_gcmf[:,0], MW_gcmf[:,1], linestyle='dashed', color='red', label='Milky Way')
ax.plot(M31_gcmf[:,0], M31_gcmf[:,1], linestyle='dashdot', color='magenta', label='M31')
ax.plot(empfiducial_gcmf[:,0], empfiducial_gcmf[:,1], linestyle='dashed', color='cyan', label='EMP-Fiducial')
ax.plot(empschecter_gcmf[:,0], empschecter_gcmf[:,1], linestyle='dotted', color='brown', label='EMP-Schechter')
ax.plot(emosaics_gcmf[:,0], emosaics_gcmf[:,1], linestyle='dashed', color='black', label='E-MOSAICS')
ax.legend()
ax.set_xlim(min_gcmf, 6e6)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r'Cluster Mass [$M_\odot$]')
ax.set_ylabel('Mass Distribution')
plt.savefig('massfunctions.pdf')
print('Mass functions plot saved')