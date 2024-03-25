#%%
import numpy as np
import matplotlib.pyplot as plt
from cicTools import *
from scipy import spatial
import configparser
from astropy.io import ascii

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

if invoid==True:
    print(f"""
        ngxs = {ngxs}
        nesf = {nesf}
        zspace = {zspace}
        zspaceAxis = {zspaceAxis}
        Num of JK resamplings = {jk}^3
        invoid = {invoid}
        snap = {snap}
        minmass = {minmass}
        maxmass = {maxmass}
        minradV = {minradV}
        evolDelta = {evolDelta}
        voidfile = {voidfile}
        """)
elif invoid==False:
    print(f"""
        ngxs = {ngxs}
        nesf = {nesf}
        zspace = {zspace}
        zspaceAxis = {zspaceAxis}
        Num of JK resamplings = {jk}^3
        snap = {snap}
        minmass = {minmass}
        """)

#
#-----------
# Print voids name file
#-----------
#
if invoid==True:
    #if voids_zs==False: 
    # Comento esta linea y el elif de abajo porque estoy dejando implícito 
    # que si zspace==True entonces leo los voids identificados en zspace
    if zspace==False:
        if evolDelta==False:
            if delta=='09':
                if voidfile=='1e9': voidsfile='../data/tng300-1_voids.dat'
                elif voidfile=='1e10': voidsfile='../data/voids_1e10.dat'
                elif voidfile=='1e11': voidsfile='../data/voids_1e11.dat'
            if delta=='08':
                if voidfile=='1e9': voidsfile='../data/voids_1e9_08.dat'
                elif voidfile=='1e10': voidsfile='../data/voids_1e10_08.dat'
                elif voidfile=='1e11': voidsfile='../data/voids_1e11_08.dat'
            if delta=='07':
                if voidfile=='1e9': voidsfile='../data/voids_1e11_07.dat'
                elif voidfile=='1e10': voidsfile='../data/voids_1e10_07.dat'
                elif voidfile=='1e11': voidsfile='../data/voids_1e11_07.dat'
        if evolDelta==True:
            if voidfile!='1e11': raise Exception('Voids not identified with evolved delta for this "voidfile" value')
            voidsfile = f'../data/voids_1e11_snap{snap}.dat'

    #elif voids_zs==True:
    elif zspace==True:
        if evolDelta==False:
            if delta=='09':
                if voidfile=='1e9': voidsfile='../data/voids_zs_1e9_09.dat'
                elif voidfile=='1e10': voidsfile='../data/voids_zs_1e10_09.dat'
                elif voidfile=='1e11': voidsfile='../data/voids_zs_1e11_09.dat'
            if delta=='08':
                if voidfile=='1e9': voidsfile='../data/voids_zs_1e9_08.dat'
                elif voidfile=='1e10': voidsfile='../data/voids_zs_1e10_08.dat'
                elif voidfile=='1e11': voidsfile='../data/voids_zs_1e11_08.dat'
            if delta=='07':
                if voidfile=='1e9': voidsfile='../data/voids_zs_1e9_09.dat'
                elif voidfile=='1e10': voidsfile='../data/voids_zs_1e10_07.dat'
                elif voidfile=='1e11': voidsfile='../data/voids_zs_1e11_07.dat'
        if evolDelta==True:
            if voidfile!='1e11': raise Exception('Voids not identified with evolved delta for this "voidfile" value')
            voidsfile = f'../data/voids_zs_1e11_snap{snap}.dat'
    print('void file location:', voidsfile)

#
#-----------
# Namefile
#-----------
#
if write==True:
    if ngxs!=0:
        namefile = f'../data/dilut{ngxs}_nesf{nesf}'
    else:
        namefile = f'../data/allgxs_nesf{nesf}'

    if zspace==True: 
        namefile += f'_redshift{zspaceAxis}'

    if invoid == True:
        namefile+= '_invoid'

    #if completeRrange == True:
    #    namefile+='_allR'

    if jk!=0:
        namefile += '_jk'

    if snap!=99:
        namefile+=f'_snap{snap}'

    if minmass==0.:
        namefile+=f'_minMass1e10'
    elif minmass==1.:
        namefile+=f'_minMass1e11'
    elif minmass==2.:
        namefile+=f'_minMass1e12'
    elif minmass==-1.5:
        namefile+=f'_minMass1e8.5'
    elif minmass==-2.:
        namefile+=f'_minMass1e8'
    elif minmass not in [-2,-1.5,-1.,0.,1.,2.]:
        raise Exception('Invalid minmass value')
    
    if maxmass!=3.:
        if maxmass==2.:
            namefile+='_maxMass1e12'
        elif maxmass==1.:
            namefile+='_maxMass1e11'
        elif maxmass==0.:
            namefile+='_maxMass1e10'
        elif minmass not in [0.,1.,2.]:
            raise Exception('Invalid "maxmass" value')

    if invoid==True:
        if evolDelta==False:
            if voidfile=='1e9': namefile+='_v1e9'
            elif voidfile=='1e10': namefile+='_v1e10'
            elif voidfile=='1e11': namefile+='_v1e11'
            #if voids_zs==True: namefile+='zs'
            namefile += f'_minradV{minradV}'
        if evolDelta==True:
            if voidfile=='1e11': 
                namefile+='_v1e11EvolDelta'
            else: 
                raise Exception('Voids not identified with evolved delta for this "voidfile" value')
            namefile += f'_minradV{minradV}'

    if delta!='09':
        namefile += f'_d{delta}'

    namefile += '.npz'

    print('Filename to be created:',namefile)




#%%
"""
Voy probando rangos de radio para calcular chi
Radio mínimo: tal que en el eje x (xi*Nmean) me de alrededor de 0.1
Radio maximo: tal que la chi no me de inf*
Estos depende del tamaño de la muestra (ngxs)

*tambien sucede que si rmax es muy grande P0 es muy chico, y ln(P0)
crece asintoticamente

rs range for ngxs=10000: np.geomspace(1500,16000,x)
rs range for ngxs=100000: np.geomspace(500,9000,x)
rs range for ngxs=1000000: np.geomspace(200,5000,x)
"""
# if ngxs==0: rs = np.geomspace(40,4000,rsbin) 
# elif ngxs==10000000: rs = np.geomspace(30,5000,rsbin) 
# elif ngxs==1000000: rs = np.geomspace(190,5800,rsbin)
# elif ngxs==100000: rs = np.geomspace(800,9100,rsbin)
# elif ngxs==10000: rs = np.geomspace(2000,17100,rsbin)
# elif ngxs==1000: rs = np.geomspace(7000,27800,rsbin)

# if invoid==True:
#     if ngxs==0: rs = np.geomspace(250,2800,10) 
#     elif ngxs==10000000: rs = np.geomspace(300,3000,10) 
#     elif ngxs==1000000: rs = np.geomspace(700,3500,10)

#
#-----------
# Determine probing range
#-----------
#
# if completeRrange==False: 
#     rs = np.geomspace(250,2500,rsbin) #Dejo esto para que tome algún valor en caso que invoid==False
#     if invoid==True:
#         if minradV==7.:
#             rs = np.geomspace(250,2500,rsbin) 
#         elif minradV==9.:
#             rs = np.geomspace(1500,4500,rsbin)

# if completeRrange==True: rs = np.geomspace(40,4000,rsbin)
if invoid==True:

    #Read Voids
    # if 'snap' in voidsfile: #If 'snap' is in 'voidsfile' it is reading voids id'd with evolved delta
    #     voids = ascii.read(voidsfile,\
    #         names=['r','x','y','z','vx','vy','vz',\
    #             'deltaint_1r','maxdeltaint_2-3r'])

    # else:
    #     voids = ascii.read(voidsfile,\
    #         names=['r','x','y','z','vx','vy','vz',\
    #             'deltaint_1r','maxdeltaint_2-3r','log10Poisson','Nrecenter'])

    try:
        voids = ascii.read(voidsfile,\
            names=['r','x','y','z','vx','vy','vz',\
                'deltaint_1r','maxdeltaint_2-3r','log10Poisson','Nrecenter'])
    except:
        voids = ascii.read(voidsfile,\
            names=['r','x','y','z','vx','vy','vz',\
                'deltaint_1r','maxdeltaint_2-3r'])


    voids = voids[voids['r']>=minradV]
    print('N of voids:',len(voids))
    voids['r'] = voids['r']*1000 #Converts Mpc to kpc
    voids['x'] = voids['x']*1000
    voids['y'] = voids['y']*1000
    voids['z'] = voids['z']*1000

    # r_max:= maximum probing radius
    # Rv:= void radius
    # N_esf,v:= minimum number of probing spheres per void
    # r_max = Rv/(N_esf,v)**(1/3)

    n_invoid = int(nesf/len(voids))
    r_max = minradV*1000/(n_invoid/.64)**(1./3)
    rs = np.geomspace(200.,r_max,rsbin)


    # minrad = round(np.min(voids['r']))
    # rs = np.geomspace(minrad/40.,minrad/4.,rsbin)
    # print('rs:',rs)
else:
    rs = np.geomspace(40,4000,rsbin)

# rs = np.geomspace(1300,13000,rsbin)
print('rs:',rs)

#
#-----------
# Read data from Illustris
#-----------
#
gxs = readTNG(snap=snap,minmass=minmass)
if ngxs!=0:
    np.random.seed(seed)
    ids = np.random.choice(len(gxs),size=int(len(gxs)*ngxs))
    gxs = gxs[ids]

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

#%%
if plot==True:
    x = np.geomspace(1E-2,1E3,50)
    c='k'
    #chi = -np.log(P0)/N_mean
    #NE = N_mean*xi_mean

    plt.plot(x,np.log(1+x)/x,label='Negative Binomial',c=c)
    #a=.3
    #plt.plot(x,(1/((1-a)*(x/a)))*((1+x/a)**(1-a)-1),label='Generalized Hierarhichal',c=c,ls='--')
    #plt.plot(x,(1-np.e**(-x))/x,label='Minimal')
    plt.plot(x,(np.sqrt(1+2*x)-1)/x,label='Thermodynamical',c=c,ls='-.')
    #plt.plot(x[:-15],1-x[:-15]/2,label='Gauss',c=c)
    # Q=1
    # plt.plot(x,1-(np.euler_gamma+np.log(4*Q*x))/(8*Q),label='BBGKY',c=c,ls=':')

    plt.plot(NXi,chi,lw=2)
    plt.xscale('log')
    plt.legend(loc=3)
    plt.show()
# %%
