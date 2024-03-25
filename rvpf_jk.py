import numpy as np
import matplotlib.pyplot as plt
from cicTools import *
from scipy import spatial
import configparser
from astropy.io import ascii
import readgadget
import readfof
import redshift_space_library as RSL

config = configparser.ConfigParser()
config.read('config.ini')

write = config['PARAMS'].getboolean('write') #write files with results
plot = config['PARAMS'].getboolean('plot') #plot results for checking  

seed = int(config['PARAMS']['seed']) #random seed
lbox = float(config['PARAMS']['lbox']) #length of box
ngxs = float(config['PARAMS']['ngxs']) #dilution
zspace = config['PARAMS'].getboolean('zspace') #redshift space
zspaceAxis = config['PARAMS']['zspaceAxis'] #r-space axis
nesf = int(config['PARAMS']['nesf']) #num of test spheres
rsbin = int(config['PARAMS']['rsbin']) #num of bins of r
jk = int(config['PARAMS']['jk']) #num of bins of r
invoid = config['PARAMS'].getboolean('invoid') #redshift space
#completeRrange = config['PARAMS'].getboolean('completeRrange')
snap = int(config['PARAMS']['snap']) #snapshot number
minmass = float(config['PARAMS']['minmass']) #log10 of minimum mass
maxmass = float(config['PARAMS']['maxmass']) #log10 of maximum mass 
minradV = float(config['PARAMS']['minradV']) #minimum void radius
voidfile = str(config['PARAMS']['voidfile']) #location of voids file / which voids to use
delta = str(config['PARAMS']['delta']) #delta used in void identification
voids_zs = config['PARAMS'].getboolean('voids_zs') #read voids identified in z-space
evolDelta = config['PARAMS'].getboolean('evolDelta') #read voids identified with evolved integrated delta

rs = np.geomspace(40,4000,rsbin)

#
#-----------
# Read data from Illustris
#-----------
#
snapdir = '/home/fdavilakurban/mnt/clemente/quijote/Halos/FoF/fiducial/0/' #folder hosting the catalogue
snapnum = 4   
# get the name of the corresponding snapshot
snapshot = '/home/fdavilakurban/mnt/clemente/quijote/Snapshots/fiducial/0/snapdir_%03d/snap_%03d'%(snapnum,snapnum)

# read the redshift, boxsize, cosmology...etc in the header
header   = readgadget.header(snapshot)
BoxSize  = header.boxsize/1e3  #Mpc/h
Nall     = header.nall         #Total number of particles
Masses   = header.massarr*1e10 #Masses of the particles in Msun/h
Omega_m  = header.omega_m      #value of Omega_m
Omega_l  = header.omega_l      #value of Omega_l
h        = header.hubble       #value of h
redshift = header.redshift     #redshift of the snapshot
Hubble   = 100.0*np.sqrt(Omega_m*(1.0+redshift)**3+Omega_l)#Value of H(z) in km/s/(Mpc/h)

print('BoxSize = %.3f Mpc/h'%BoxSize)
print('Number of particles in the snapshot:',Nall)
print('Omega_m = %.3f'%Omega_m)
print('Omega_l = %.3f'%Omega_l)
print('h = %.3f'%h)
print('redshift = %.1f'%redshift)

# read the halo catalogue
FoF = readfof.FoF_catalog(snapdir, snapnum, long_ids=False,
                          swap=False, SFR=False, read_IDs=False)
# get the properties of the halos
gxs = FoF.GroupPos/1e3  #Halo positions in Mpc/h

#
#-----------
# Replicate box edges periodically
#-----------
#
print('Replicating box:')
newgxs = perrep(gxs,lbox,np.max(rs))
print(f'Num of original gxs in box: {len(gxs)}\n\
Num of gxs after replication: {len(newgxs)}')

#
#-----------
# Distant observer aproximation for z-space
#-----------
#
if zspace == True:
    H0 = .06774
    axis = zspaceAxis
    vaxis = 'v'+axis
    newgxs[axis]+=newgxs[vaxis]/H0
    newgxs[axis][np.where(newgxs[axis]<0.)[0]]+=lbox
    newgxs[axis][np.where(newgxs[axis]>lbox)[0]]-=lbox


#
#-----------
# VPF calculations
#-----------
#
pos = np.column_stack((newgxs['x'],newgxs['y'],newgxs['z']))

tree = spatial.cKDTree(pos)

chi = np.zeros(len(rs))
NXi = np.zeros(len(rs))
P0 = np.zeros(len(rs))
N_mean = np.zeros(len(rs))
xi_mean = np.zeros(len(rs))

if invoid == False:
    if jk!= 0:
        chi_std = np.zeros(len(rs))
        NXi_std = np.zeros(len(rs))
        P0_std = np.zeros(len(rs))
        N_mean_std = np.zeros(len(rs))
        xi_mean_std = np.zeros(len(rs))

        print('Calculating JK cic statistics...')
        
        for i,r in enumerate(rs):
            chi[i], NXi[i], P0[i], N_mean[i], xi_mean[i],\
                    chi_std[i], NXi_std[i], P0_std[i], N_mean_std[i], xi_mean_std[i]\
                        = cic_stats_jk(tree, nesf, r, lbox, jk)
    else:
        print('Calculating cic statistics...')
        for i,r in enumerate(rs):
            chi[i], NXi[i], P0[i], N_mean[i], xi_mean[i],\
                        = cic_stats(tree, nesf, r, lbox)

if invoid == True:

    # Comentado porque lo copie al principio cuando defino 'rs'
    # Guardado por las dudas
    # ----------------------------------------------------------
    # voids = ascii.read(voidsfile,\
    #     names=['r','x','y','z','vx','vy','vz',\
    #         'deltaint_1r','maxdeltaint_2-3r','log10Poisson','Nrecenter'])
    # voids = voids[voids['r']>=minradV]
    # print('N of voids:',len(voids))
    # voids['r'] = voids['r']*1000 #Converts kpc to Mpc
    # voids['x'] = voids['x']*1000
    # voids['y'] = voids['y']*1000
    # voids['z'] = voids['z']*1000

    if jk==3:
        chi_std = np.zeros(len(rs))
        NXi_std = np.zeros(len(rs))
        P0_std = np.zeros(len(rs))
        N_mean_std = np.zeros(len(rs))
        xi_mean_std = np.zeros(len(rs))

        print('Calculating JK invoid cic statistics...')
        
        for i,r in enumerate(rs):
            chi[i], NXi[i], P0[i], N_mean[i], xi_mean[i],\
                    chi_std[i], NXi_std[i], P0_std[i], N_mean_std[i], xi_mean_std[i]\
                        = cic_stats_invoid_jk(voids, tree, nesf, r)

    else:
        print('Calculating invoid cic statistics...')
        for i,r in enumerate(rs):
            chi[i], NXi[i], P0[i], N_mean[i], xi_mean[i],\
                        = cic_stats_invoid(voids, tree, nesf, r)


#
#-----------
# Writing file
#-----------
#
if write==True:
    print(f'Creating {namefile}')
    if jk!=0:
        np.savez(namefile,chi,chi_std,NXi,NXi_std,P0,P0_std,N_mean,N_mean_std,xi_mean,xi_mean_std,rs)
    else:
        np.savez(namefile,chi,NXi,P0,N_mean,xi_mean,rs)

print(chi, rs)