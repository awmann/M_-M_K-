import numpy as np
import pkg_resources
from astropy.io import fits
##figure out where the big fits files are in this installation
datapath = pkg_resources.resource_filename('mk_mass','resources')

def posterior(K,dist,ek=0.0,edist=0.0,feh=None,efeh=None,oned=False,silent=False,post=None):
    '''
    Calculate a mass posterior given an input K magnitude and distance
    Based on The M_{K_S}-Mass relation from Mann et al. (2018b).

    http://adsabs.harvard.edu/abs/PAPER

    See function tester below for some example uses

    PARAMETERS:
        K: float,int
            2MASS Ks-band (magnitudes)
        dist: float,int
            distance (parsecs)
        ek: float,int
            error on K (magnitudes). If not provided assumed to be 0
        edist: float,int
            error on distance (parcsecs). If not provided assumed to be 0


        feh:  float,optional default None
            Iron abundance. If provided the code will use the
            Mk-Mass-[Fe/H] relation instead
        efeh: float,optional default None
            error on iron abundance (dex). If not provided, but feh
            is provied, assumed to be 0
        post: array_like default None
            MCMC posterior. Will read it in if not provided. This is
            useful if you need to run a lot of values (so it does not
            need to be read in many times).


        oned:boolean,optional default False
            Returns a simple 1D error instead of a posterior

        silent:boolean default False
            surpress the many (annoying) warnings

    Returns
        mass:
            if /oned is set this is a 2-element array containing the
            median and standard deviation of the mass posterior. Else
            this contains the posterior on mass of 1600000
            elements (the length of the trimmed posterior if you input one yourself).

    EXAMPLE (see tester below for more examples):
        Typical case
        k = 8.782
        ek = 0.02
        dist = 14.55
        edist = 0.13
        mass = mk_mass.posterior(k,dist,ek,edist)
        print,median(mass),stdev(mass)
    If you have a posterior on K and distance instead of 1D errors:
        post = mrdfits('',/silent)
        mass = ()
        for i in range(len(l)): mass += (mk_mass.posterior(k[i],dist[i],0.0,0.0,oned=True)[0],)

    MODIFICATION HISTORY:
        May 10 2018: Ported from scratch code by A. Mann
        May 14 2018: Added testing modules. A. Mann
        May 14 2018: Python version made. A Rizzuto


    If you use this code, please cite the relevant paper:

    '''

    ##Check inputs for correct types first
    try:
        K     = float(K)
        dist  = float(dist)
        ek    = float(ek)
        edist = float(edist)
    except:
        if silent == False: print ('Inputs K, dist,eK, edist; all should be integers or floats')
        return -1

    ##SETUP and WARNINGS
    if (ek < 0.00000001):
        if silent == False: print ('Warning, assuming no error on K, mass errors underestimated!')
        ek = 0.0
    if edist < 0.000000001:
        if silent == False: print ('Warning, assuming no error on K, mass errors underestimated!')
        edist  = 0.0
    if (feh != None) & (efeh == None):
        if silent == False: print ('Warning, assume no error on [Fe/H]. mass errors underestimated!')
        efeh = 0.0
        try:
            feh = float(feh)
        except:
            if silent == False: print ('Warning, feh should be a scalar float or integer!')
    if efeh != None:
        try:
            efeh = float(efeh)
        except:
            if silent == False: print ('Warning, efeh should be a scalar float or integer!' )

    ##Unpack the relation
    if post == None:
        if feh != None:
            post = fits.open(datapath +'/Mk-M_8_feh_trim.fits')[0].data
        else: post = fits.open(datapath +'/Mk-M_7_trim.fits')[0].data
    ntot = post.shape[0]
    a0    = post[:,0]
    a1    = post[:,1]
    a2    = post[:,2]
    a3    = post[:,3]
    a4    = post[:,4]
    a5    = post[:,5]
    if feh == None:
        sige   = np.exp(post[:,6])
        f   = a0*0.0
        feh = 0.0
    else:
        e = 0.0*a0
        f = post[:,6]
        sige = np.exp(post[:,7])
    ##Compute the posterior
    m        = 0.0
    kmag     = np.random.normal(K,ek,a0.shape[0])
    distance = np.random.normal(dist,edist,a0.shape[0])
    mk       = kmag + 5 - 5*np.log10(distance)
    zp       = 7.5
    m        = (10.0**(a0+a1*(mk-zp)+a2*(mk-zp)**2.0+a3*(mk-zp)**3.0+a4*(mk-zp)**4.0+a5*(mk-zp)**5.0))*(1.0+feh*f)
    m       += np.median(sige)*m*np.random.normal(size=a0.shape[0])


    if oned == True: m = [np.mean(m),np.std(m)]

    if np.isnan(m).any() == True: print ('Warning, some outputs are NaNs!!?!!')
    return m

def tester():
    '''
    This is a small code to help test the output of the mk_mass code
    It just runs a few well-known systems (mostly planet hosts) and
    checkes how the masses compare to semi-independent determinations
    (e.g., using models, transit-fit density, etc.)
    '''
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    mpl.rcParams['lines.linewidth']   =3
    mpl.rcParams['axes.linewidth']    = 2
    mpl.rcParams['xtick.major.width'] =2
    mpl.rcParams['ytick.major.width'] =2
    mpl.rcParams['ytick.labelsize'] = 10
    mpl.rcParams['xtick.labelsize'] = 10
    mpl.rcParams['axes.labelsize'] = 12
    mpl.rcParams['legend.numpoints'] = 1
    mpl.rcParams['axes.labelweight']='bold'
    mpl.rcParams['mathtext.fontset']='stix'
    mpl.rcParams['font.weight'] = 'semibold'


    ##TRAPPIST-1
    k,ek  = 10.296,0.023
    dist  = 12.23989539
    edist = 0.018910230
    mass  = posterior(k,dist,ek,edist)
    print ('Trappist-1:')
    print ('Our mass: ' + str(np.median(mass)) + '+/-' + str(np.std(mass)))
    print ('Van Grootel et al: 0.089+/-0.006')

    mVG = np.random.normal(0.089,0.006,mass.shape[0])
    fig,ax = plt.subplots()
    __,thebins,__ = ax.hist(mass,bins=400,color='r',label='Our Mass',alpha=0.5)
    ax.hist(mVG,bins=thebins,color='b',label='Van Grootel et al.',alpha=0.5)
    ax.legend()
    ax.set_xlabel(r'Mass ($M_\odot$)')
    ax.set_title('Trappist-1')

    print ('Diff = '+str(np.median(mass-mVG)))
    print ('Sig = '+str(np.mean(mass-mVG)/np.std(mass-mVG)))
    plt.show()

    ##GJ1214
    k,ek  = 8.782,0.02
    dist  = 14.55
    edist = 0.13
    feh,efeh=0.3,0.1
    mass      = posterior(k,dist,ek,edist)
    mass_feh  = posterior(k,dist,ek,edist,feh,efeh)
    print ('GJ1214:')
    print ('Our mass: ' + str(np.median(mass)) + '+/-' + str(np.std(mass)))
    print ('Our mass with [Fe/H]: ' + str(np.median(mass_feh)) + '+/-' + str(np.std(mass_feh)))
    print ('Anglada-Escude: +0.176+/-0.009')

    mVG = np.random.normal(0.176,0.009,mass.shape[0])
    fig2,ax2 = plt.subplots()
    __,thebins,__ = ax2.hist(mass,bins=400,color='r',label='Our Mass',alpha=0.5)
    ax2.hist(mass_feh,bins=thebins,color='g',label='Our Mass with [Fe/H]',alpha=0.5)
    ax2.hist(mVG,bins=thebins,color='b',label='Anglada Escude et al',alpha=0.5)
    ax2.legend()
    ax2.set_xlabel(r'Mass ($M_\odot$)')
    ax2.set_title('GJ1214')
    print ('Diff = '+str(np.median(mass-mVG)))
    print ('Sig = '+str(np.mean(mass-mVG)/np.std(mass-mVG)))
    print ('Diff(feh) = '+str(np.median(mass_feh-mVG)))
    print ('Sig(feh) = '+str(np.mean(mass_feh-mVG)/np.std(mass_feh-mVG)))
    plt.show()

    return -1
