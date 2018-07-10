#!/usr/bin/env python
################################################################################
# Plot FOV of HAWC in Equatorial Coordinates for a given DateTime or MJD
# compared to some region of interest.
# healpy and numpy packages need to be installed.
#
# Adapted from plotMollweide.py
################################################################################

__version__ = "$Id: $"

try:
    import os

    import argparse
    import numpy as np
    import healpy as hp
    from math import log
    import re

    import matplotlib as mpl

    if not os.getenv("DISPLAY"):
        mpl.use('PDF')

    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates

    from datetime import datetime, timedelta

    from hawc import hawcnest, data_structures, astro_service, detector_service
    from hawc.hawcnest import HAWCUnits as U
    from hawc.data_structures import *
    from HAWCNest import HAWCNest

except ImportError as e:
    print(e)
    raise SystemExit

## Try importing TeVCat
try:
    import TeVCat
    haveTeVCat = True
except ImportError as e:
    haveTeVCat = False
    print(e)

## Try importing the Fermi catalogs
try:
    import FGLCatalog
    import FHLCatalog
    haveFermi = True
except ImportError as e:
    haveFermi = False
    print(e)

## Connect to the database ##
useDB = False
try:
    import mysql.connector as mysql
    useDB = True
except ImportError as e:
    print(e)

## Functions

def gpsTimeToDatetime(sec,ns):
    """Convert GPS seconds and nanoseconds to a datetime object."""
    GPS_EPOCH = datetime(1980, 1, 6)
    return GPS_EPOCH + timedelta(seconds=sec, microseconds=1e-3*ns)

def setupAxes(fig, coord="C"):
    """Create the color bar for a HEALPix map image."""
    for ax in fig.get_axes():
        if type(ax) is hp.projaxes.HpxMollweideAxes:
            if coord == "C":
                ax.annotate("0$^\circ$", xy=(1.8, 0.625), size="large")
                ax.annotate("360$^\circ$", xy=(-1.95, 0.625), size="large")
                ax.annotate("Equatorial", xy=(1.65, -0.85), size="large")
            elif coord == "N":
                ax.annotate("-180$^\circ$", xy=(1.8, 0.625), size="large")
                ax.annotate("180$^\circ$", xy=(-1.95, 0.625), size="large")
                ax.annotate("Equatorial", xy=(1.65, -0.85), size="large")
            elif coord == "E":
                ax.annotate("0$^\circ$", xy=(1.8, 0.625), size="large")
                ax.annotate("360$^\circ$", xy=(-1.95, 0.625), size="large")
                ax.annotate("Eclipitc", xy=(1.65, -0.85), size="large")
            elif coord == "G":
                ax.annotate("-180$^\circ$", xy=(1.65, 0.625), size="large")
                ax.annotate("180$^\circ$", xy=(-1.9, 0.625), size="large")
                ax.annotate("Galactic", xy=(1.65, -0.85), size="large")

def plotOtherFOV(lat, lon, alt, clr, mjd, j2000, astr, crds, upgoing=False):
    localeOther = LatLonAlt( lat*U.degree, lon*U.degree, alt*U.meter) #
    equOther = EquPoint()
    z = -1. if upgoing else 1.  #upward going events z=-1
    astr.Loc2Equ(mjd, localeOther, Vector(0.,0.,z), equOther, toJ2000=j2000) 
    otherVec = hp.ang2vec(U.halfpi - equOther.dec, equOther.ra)
    for zen in np.arange(15,91.,15):
        hp.projplot(
                     [ np.arange(0,360, 1.), np.ones(360)*(90-zen)],
                     lonlat=True, color=clr, coord=crds,
                     linestyle="-", linewidth=1.,
                     rot= (equOther.ra/U.degree, equOther.dec/U.degree,0)
                   )

def unixTimeToDatetime(sec,ns):
    """Convert UNIX seconds and nanoseconds to a datetime object."""
    UNIX_EPOCH = datetime(1970, 1, 1)
    return UNIX_EPOCH + timedelta(seconds=sec, microseconds=1e-3*ns)

def main():
    #Important
    p = argparse.ArgumentParser(description="Plot maps in Mollweide projection")
    p.add_argument("fitsfile", nargs="?", help="Input HEALPix FITS file")
    p.add_argument("-c", "--coords", dest="coords", default="C",
                   help="C=equatorial, G=galactic, E=ecliptic, N=negative RA values [-180,180]")
    p.add_argument("--location", dest="location", type=float, nargs=2, default=[83.63,22.01],
                   help="location to mark on map as RA, Dec pair [degrees]")
    p.add_argument("--datetime", dest="datetime", default=None,
                   help="datetime of transient as 'YYYY MM DD HH MM SS'")
    p.add_argument("--mjd", dest="mjd", default=None, type=float,
                   help="Modified Julian Date")
    p.add_argument("--errorCirc", dest="errorCirc", type=float, default=None,
                   help="Error circle around location [degrees]")
    p.add_argument("-o", "--output", dest="output", default=None,
                   help="Output image file (.png preferred)")
    p.add_argument("-T", "--tmax", dest="tmax", default=6., type=float,
                   help="+/- time range [hours]")
    p.add_argument("--j2000", default=False, action="store_true",
                   help="location is given in J2000 coordinates and FOV is plotted in J2000")
    p.add_argument("--antares", default=False, action="store_true",
                   help="plot antares FOV")
    p.add_argument("--fact", default=False, action="store_true",
                   help="plot FACT FOV")
    p.add_argument("--hess", default=False, action="store_true",
                   help="plot HESS FOV")
    p.add_argument("--icecube", default=False, action="store_true",
                   help="plot icecube FOV")
    p.add_argument("--magic", default=False, action="store_true",
                   help="plot MAGIC FOV")
    p.add_argument("--veritas", default=False, action="store_true",
                   help="plot VERITAS FOV")
    

    # Less Important
    try:
        pwdfile = "%s/etc/pwd.cfg" % os.environ['HOME']
    except Exception, e:
        print "Env var $HOME not found."
        pwdfile = None
    p.add_argument("--pwd",dest="pwd", default=pwdfile,
                   help = "location of file with DB password")
    p.add_argument("-z", "--maxZen", dest="maxZen", default=45., type=float,
                   help="Maximum zenith angle for FOV map [degrees]")
    p.add_argument("-Z", "--zenStep", dest="zenStep", default=2.5, type=float,
                   help="Step size for zenith band in FOV map [degrees]")
    p.add_argument("-S", "--size", dest="xsize", type=int, default=1600,
                   help="Map x-axis pixel number; increase for high resolution")
#    p.add_argument("-T", "--title", dest="title", default="",
#                   help="Plot title")
    p.add_argument("--fermicat", dest="fermicat", default=None,
                   help="Fermi xFGL catalog FITS file name")
    p.add_argument("--fermicat-labels", dest="fermicatLabels",
                  action="store_true", default=False,
                  help="Draw Fermi sources with labels.")
    p.add_argument("--download-tevcat", dest="dtevcat", action="store_true",
                   default=False,
                   help="Download data from tevcat.uchicago.edu and exit")
    p.add_argument("--tevcat", dest="tevcat", default=None,
                   help="Draw TeVCat sources using TeVCat ASCII file")
    p.add_argument("--tevcat-labels", dest="tevcatLabels", default=None,
                   help="Draw TeVCat sources with labels w/ TeVCat ASCII file")
    p.add_argument("--hotspots", dest="hotspots", default=None,
                   help="Hotspot coordinates in an ASCII file")
    p.add_argument("--preliminary", action="store_true", dest="preliminary",
                   default=False,
                   help="Add a watermark indicating the plot as preliminary")
    p.add_argument("--alpha", action="store_true", dest="alpha",
                   default=0.25,
                   help="Alpha used for drawn lines.")
    p.add_argument("--logz", action="store_true", dest="logz",
                   help="Set color scale in log.")

    # Not Important
    p.add_argument("-m", "--min", dest="min", type=float, default=None,
                   help="Plot minimum value")
    p.add_argument("-M", "--max", dest="max", type=float, default=None,
                   help="Plot maximum value")
    p.add_argument("-n", "--ncolors", dest="ncolors", type=int, default=256,
                   help="Number of colors to use in the color map")

    args = p.parse_args()

    # Get Local Sidereal Angle (return other stats)
    ###
    hawcnest.SetLoggingLevel(5, False)
    nest = HAWCNest()
    nest.Service("StdAstroService", "astro")
    nest.Service("ConfigDirDetectorService", "det", configDir=os.environ["CONFIG_HAWC"])
    nest.Configure()
    nest.Finish()
    astro = astro_service.GetService("astro")

    if args.mjd and args.datetime:
        print "Only specify MJD or Datetime. Not both!"
        raise SystemExit
        
    if args.mjd:
        mjd = ModifiedJulianDate(args.mjd*U.day)
        gps = mjd.timestamp
        utc = UTCDateTime(gps)
    elif args.datetime:
        dt  =  [ int(i) for i in args.datetime.split()  ]
        utc = UTCDateTime(dt[0], dt[1], dt[2], dt[3], dt[4], dt[5] )
        gps = utc.timestamp
        mjd = ModifiedJulianDate(utc)
    else:
        utc = GetCurrentTime()
        gps = utc.timestamp
        mjd = ModifiedJulianDate(utc)

    detector = detector_service.GetService("det").GetDetector(gps)
    locale = detector.lat_lon_alt
    NS = "N" if locale.latitude >= 0. else "S"
    EW = "E" if locale.longitude >= 0.  else "W"
    
    axis = Vector()
    equ = EquPoint( args.location[0]*U.degree, args.location[1]*U.degree)
    astro.Equ2Loc(mjd, locale, equ, axis, fromJ2000=args.j2000)

    equCurr = EquPoint()
    astro.Loc2Equ(mjd, locale, axis, equCurr)

    gmst = astro.GetGMST(mjd)
    lst  = gmst + locale.longitude
    if lst < 0:
        lst += U.twopi
    if lst > U.twopi:
        lst -= U.twopi
    lha   = lst - equCurr.ra
    if ( lha < -U.pi ):
      lha += U.twopi;
    if ( lha > U.pi ):
      lha -= U.twopi;
    status = "rising" if lha < 0. else "setting"
    
    raHr  = equ.ra*12/np.pi
    raMin = ( raHr - int(raHr) ) * 60
    raSec = ( raMin - int(raMin) ) * 60

    dec    = equ.dec/U.degree
    decMin = ( dec - int(dec) ) * 60
    decSec = ( decMin - int(decMin) ) * 60

    msg = "".ljust(60,'=')
    msg+= "\n" + ("Source at Zenith of %.02f deg and %s" % (axis.theta/U.degree, status)).center(60,' ')
    msg+= "\n" + "".ljust(60,'=')

    msg+= "\n  Detector coordinates:"
    msg+= "\n\t- Latitude ".ljust(25,'_') + " %.02f %s" % (abs(locale.latitude/U.degree), NS)
    msg+= "\n\t- Longitude ".ljust(25,'_') + " %.02f %s" % (abs(locale.longitude/U.degree), EW)
    msg+= "\n\t- Altitude ".ljust(25,'_') + " %d m.a.s.l." % (locale.altitude)

    msg+= "\n\n  Trigger Time:"
    msg+= "\n\t- Mod. Julian Day ".ljust(25,'_')   + " {:.4f} UT".format(mjd.get_date()/U.day)
    msg+= "\n\t- UTC Date Time ".ljust(25,'_')     + " %s" % utc
    msg+= "\n\t- GM Sidereal Time ".ljust(25,'_')  + " {:6.2f} deg = {:6.2f} hrs".format(gmst/U.degree, gmst*12./np.pi)
    msg+= "\n\t- Loc Sidereal Time ".ljust(25,'_') + " {:6.2f} deg = {:6.2f} hrs".format(lst/U.degree, lst*12./np.pi)

    if args.j2000:
        msg+= "\n\n  Source position (J2000):"
    else:
        msg+= "\n\n  Source position (current epoch):"
    msg+= "\n\t- Right Ascension ".ljust(25,'_') + " {:7.2f} deg = {:3d}hrs {:2d}m {:2d}s".format(equ.ra/U.degree, int(raHr), int(raMin), int(raSec))
    msg+= "\n\t- Declination ".ljust(25,'_') + " {:7.2f} deg = {:3d}deg {:2d}m {:2d}s".format(equ.dec/U.degree, int(dec), int(decMin), int(decSec))
    
    msg+= "\n\t- Local HA ".ljust(25,'_') + " {:7.2f} deg = {:6.2f} hrs".format(lha/U.degree, lha*12./np.pi)
    msg+= "\n\t- Local Zenith ".ljust(25,'_') + " {:7.2f} deg".format(axis.theta/U.degree)
    msg+= "\n\t- Local Azimuth ".ljust(25,'_') + " {:7.2f} deg".format(axis.phi/U.degree)
    msg+= "\n" + "".ljust(60,'=')
    msg+= "\n" + "".ljust(60,'=')

    print msg
    logfile = msg
    
    ###


    # Download TeVCat
    if args.dtevcat:
        if haveTeVCat:
            print("Fetching data from tevcat.uchicago.edu.")
            tevcat = TeVCat.TeVCat()
            return None
        else:
            print("Sorry, AERIE TeVCat python module is not available.")
            raise SystemExit

    # Start normal processing

    znorm = None
    if (args.logz):
        znorm = 'log'

    # Set up the figure frame and coordinate rotation angle
    if args.coords == "C":
        fig = plt.figure(1, figsize=(8.5,9.5), dpi=100, tight_layout=True)
        fig.add_subplot(211)
        rotation = 180
        coords = "C"
    elif args.coords == "N":
        fig = plt.figure(1, figsize=(8.5,9.5), dpi=100, tight_layout=True)
        fig.add_subplot(211)
        rotation = 0
        coords = "C"
    elif args.coords == "G":
        fig = plt.figure(1, figsize=(8.5,9.5), dpi=100, tight_layout=True)
        fig.add_subplot(211)
        coords = "CG"
        rotation = 0.
    else:
        fig = plt.figure(1, figsize=(8.5,9.5), dpi=100, tight_layout=True)
        fig.add_subplot(211)
        coords = "CE"
        rotation = 180.

    # Set up FOV map
    FOVmap = np.zeros(512*512*12)

    cosDep  = np.cos( np.arange(0,args.maxZen,args.zenStep)*U.degree )
    sinDep  = np.sin( np.arange(0,args.maxZen,args.zenStep)*U.degree )
    linDep  = np.arange(0,args.maxZen,args.zenStep)*U.degree 

    # Add HAWC FOV Map
    equHAWC = EquPoint()
    astro.Loc2Equ(mjd, locale, Vector(0.,0.,1.), equHAWC, toJ2000=args.j2000)
    hawcVec = hp.ang2vec(U.halfpi - equHAWC.dec, equHAWC.ra)

    for z in linDep:
        pix = hp.query_disc(512, hawcVec, z)
        FOVmap[pix] += 1.
    FOVmap[np.where( FOVmap[:] <  1.)] = hp.UNSEEN

    # Set up the map minimum and maximum values
    notMask = FOVmap[FOVmap > hp.UNSEEN]
    dMin, dMax = (notMask[np.argmin(notMask)], notMask[np.argmax(notMask)])
    if args.min is not None:
        dMin = args.min
        FOVmap[(FOVmap<dMin) & (FOVmap > hp.UNSEEN)] = dMin
    if args.max is not None:
        dMax = args.max
        FOVmap[(FOVmap>dMax) & (FOVmap > hp.UNSEEN)] = dMax

    mpl.rc("font", family="serif", size=14)


    hp.mollview(FOVmap, fig=1,
                        hold=True,
                        xsize=args.xsize,
                        coord=coords,
                        min=dMin,
                        max=dMax,
                        rot=rotation,
                        #title= r"Trigger at $(\alpha=%.02f^\circ, \delta = %.02f^\circ)$ on %s [UTC]" % (args.location[0], args.location[1],utc),
                        title= "%s [UTC]" % (utc),
                        cbar=False,
                        notext=True,
                        norm=znorm)

    hp.graticule(local=True);

    # Add ANTARES FOV Ring
    if args.antares:
        plotOtherFOV( 42.+48./60., 6.+10./60., -2475, '#ffd966', mjd, args.j2000, astro, coords, upgoing=True)
    # Add FACT/MAGIC FOV Ring
    if args.fact or args.magic:
        plotOtherFOV( 28.757 , -17.887, 2200, 'darkgreen', mjd, args.j2000, astro, coords, upgoing=False)
    # Add HESS FOV Ring
    if args.hess:
        plotOtherFOV(  -23-16.28/60., 16+30./60., 1800, '#cc0000', mjd, args.j2000, astro, coords, upgoing=False)
    # Add IceCube FOV Ring
    if args.icecube:
        plotOtherFOV(  -89-59.5/60., -63, 750, 'lightblue', mjd, args.j2000, astro, coords, upgoing=True)
    # Add VERITASFOV Ring
    if args.veritas:
        plotOtherFOV(  31+40./60., -110-57./60., 1270, '#002699', mjd, args.j2000, astro, coords, upgoing=False)

    if args.preliminary:
        hp.projtext(95*U.degree, 240*U.degree,
                    "PRELIMINARY",
                    coord=coords,
                    color=textcolor,
                    alpha=0.5,
                    rotation=0,
                    fontdict={"family":"sans-serif", "weight":"bold", "size":42})

    # If TeVCat data are available, plot them
    if args.tevcat or args.tevcatLabels:
        if haveTeVCat:
            try:
                if args.tevcat:
                    tevcat = TeVCat.TeVCat(args.tevcat)
                elif args.tevcatLabels:
                    tevcat = TeVCat.TeVCat(args.tevcatLabels)
            except IOError as e:
                print(e)
                print("Downloading data from tevcat.uchicago.edu")
                tevcat = TeVCat.TeVCat()
            except:
                print("Why caught here?")
                print("Downloading data from tevcat.uchicago.edu")
                tevcat = TeVCat.TeVCat()

            for cId in (1,2):
                catalog = tevcat.GetCatalog(cId)
                ra = catalog.GetRA()
                dec = catalog.GetDec()
                assoc = catalog.GetCanonicalName()

                hp.projscatter(90*U.degree - dec, ra, coord=coords,
                               color=textcolor,
                               facecolors="none",
                               marker="s")
                if args.tevcatLabels:
                    for r, d, s in zip(ra, dec, assoc):
                        hp.projtext(90*U.degree - d, r, s, coord=coords,
                                    color=textcolor,
                                    rotation=10,
                                    fontdict={"size":14})
                # Print significance for TeVCat sources above 3 sigma
                sThr = 3.
                if cId == 1:
                    print("TeVCat sources above %.2f sigma:" % sThr)
                for r, d, s in zip(ra, dec, assoc):
                    valueTeVCat = hp.get_interp_val(FOVmap, 90*U.degree - d, r)
                    if valueTeVCat > sThr:
                        print(" - %s: %.2f sigma" % (s, valueTeVCat))
        else:
            print("Sorry, TeVCat could not be loaded.")

    # If Fermi data are available, plot them
    if args.fermicat:
        if haveFermi:
            fcat = None
            try:
                fcat = FGLCatalog.FGLCatalog(args.fermicat)
                aflag = fcat.GetAnalysisFlags()
                acut = (aflag == 0)                     # cut analysis errors

                flux = fcat.GetFlux1000()               # 1-100 GeV flux
                dflx = fcat.GetFlux1000Error()          # flux uncertainty
                fcut = dflx/flux < 0.5                  # cut poorly measured srcs

                cuts = np.logical_and(acut, fcut)
                print('Using FGL')
            except:
                try:
                    fcat = FHLCatalog.FHLCatalog(args.fermicat)
                    cuts = (fcat.GetFlux() > 0.)        # Dummy cut
                    print('Using FHL')
                except:
                    print('Fermi catalog could not be loaded!')
        if fcat != None:
            # Don't show TeV associations if plotting from TeVCat
            if args.tevcat or args.tevcatLabels:
                tcfg = fcat.GetTeVCatFlag()
                tcut = (tcfg == "N") | (tcfg == "C")
                cuts = np.logical_and(cuts, tcut)
            ra = fcat.GetRA()[cuts]
            dec = fcat.GetDec()[cuts]
            assoc = fcat.GetSourceName()[cuts]
            catnms = fcat.GetCatalogName()[cuts]
            for i in xrange(len(assoc)):
                if assoc[i] == '':
                    assoc[i] = catnms[i]

            hp.projscatter(90*U.degree - dec, ra, coord=coords,
                           color=textcolor,
                           facecolors="none",
                           marker="o")
            if args.fermicatLabels:
              for r, d, s in zip(ra, dec, assoc):
                hp.projtext(90*U.degree - d, r, s, coord=coords,
                            color=textcolor,
                            rotation=10,
                            fontdict={"size":14})
        else:
            print("Sorry, the Fermi xFGL catalog could not be loaded.")

    # If a hotspot list is given, plot it
    if args.hotspots:
        fhot = open(args.hotspots, "r")
        ra = []
        dec = []
        for line in fhot:
            if line.startswith("#"):
                continue
            r, d = [float(t)*U.degree for t in line.strip().split()[1:3]]
            ra.append(r)
            dec.append(d)
        ra = np.array(ra)
        dec = np.array(dec)
        hp.projscatter(90*U.degree - dec, ra, coord=coords,
                       facecolors="none",
                       marker="o",
                       s=40)

    # Adjust plot so that empty southern hemisphere is cut out
    ax = fig.gca()
    # Hide not used space of image
    imgp = ax.get_images()[0]
    imgp.get_cmap().set_under('w',alpha=0.)
    imgp.get_cmap().set_over('w',alpha=0.)
    imgp.get_cmap().set_bad('w',alpha=0.)
    setupAxes(fig,coord=args.coords)

    if (args.location):
        (ra, dec) = args.location
        hp.projscatter((90 - dec)*U.degree, ra*U.degree, coord=coords,
                        color='indigo',
                        facecolors='mediumorchid',
                        marker="*",
                        s=200)
#        hp.projtext((90 - dec)*U.degree + 25*U.degree, ra*U.degree + 35*U.degree, 
#                   r"$\alpha=%.02f^\circ$"  % (args.location[0]) + "\n" + "$\delta = %.02f^\circ$" % (args.location[1]) , 
#                   coord=coords,
#                   color='black',
#                   rotation=0,
#                   fontdict={"family":"sans-serif", "weight":"bold", "size":12})
        descr = "HAWC Field of View\n(Moving Left)" if coords == "C" else "HAWC Field of View" 
        hp.projtext( U.halfpi - locale.latitude + 55*U.degree, lst,
                   descr,
                   coord=coords,
                   color='black',
                   rotation=0,
                   fontdict={"family":"sans-serif", "weight":"bold", "size":12})
        # plot error circle
        if args.errorCirc:
                hp.projplot( 
                             [ np.arange(0,360, 1.), np.ones(360)*(90-args.errorCirc)], 
                             lonlat=True, color='mediumorchid', coord=coords,
                             linestyle="-", linewidth=1.5,
                             rot= (ra,dec,0)
                           )
              
#              for axM in fig.get_axes():
#                  if type(axM) is hp.projaxes.HpxMollweideAxes:
#                      axM.fill_between( 
#                                np.arange(0,360, 1.),np.ones(360)*(90-args.errorCirc), np.ones(360)*(90-args.errorCirc)*.001, 
#                               #lonlat=True, color='mediumorchid', coord="C",
#                               #linestyle="-", linewidth=1,
#                               #rot= (ra,dec,0)
#                                  )
      
    time = []
    zen  = []
    for delhr in np.arange(-1.0*args.tmax, args.tmax+0.1, 0.1):

        mjdIt = ModifiedJulianDate( mjd.get_date() + delhr*U.hour )
        utcIt = mjdIt.datetime

        axisIt = Vector()
        equIt = EquPoint( args.location[0]*U.degree, args.location[1]*U.degree)
        astro.Equ2Loc(mjdIt, locale, equIt, axisIt, fromJ2000=args.j2000)

        time.append( datetime.strptime(str(utcIt), "%Y-%m-%d %H:%M:%S") ) 
        zen.append( axisIt.theta/U.degree )

    zen = np.array(zen)
    ax2 = fig.add_subplot(212)
    ax2.plot(time,np.ma.masked_where(zen>90,zen),"b-",linewidth=2.5, label= r"$\alpha=%.02f^\circ$"  % (args.location[0]) + "\n" + "$\delta = %.02f^\circ$" % (args.location[1]) )
    ax2.fill_between(time, args.maxZen*np.ones(len(time)), 90*np.ones(len(time)), color='gray', alpha=.5) # gray out high zenith
    ax2.text(time[0], args.maxZen+4.5, "Outside optimal viewing angle", color='black',alpha=.6)
    ax2.set_xlabel("Hour [UTC]")
    ax2.set_ylabel( "Zenith Angle [degrees]", color='blue' )
    ax2.set_ylim( [110,0] )
    ax2.set_yticks(np.arange(90,-.1,-10), minor=False)
    ax2.set_yticks([105], minor=True)
    ax2.set_yticklabels( ['ON'], color='green',minor=True)
    hours = mdates.HourLocator( interval= int(args.tmax/8.) + 1) # every hour
    months = mdates.MonthLocator() # every month
    days = mdates.DayLocator()   # every day
    hourFmt = mdates.DateFormatter('%H:%M')
    ax2.xaxis.set_major_locator(hours)
    ax2.xaxis.set_major_formatter(hourFmt)
    plt.xticks(rotation=60)
    ax2.grid()
    trigtime=[ datetime.strptime(str(utc), "%Y-%m-%d %H:%M:%S") for i in range(0,90) ]
    ax2.plot( trigtime, np.arange(90,0,-1), linestyle="--", linewidth =3.5, color='mediumorchid', label="Trigger Time")
    ax2.plot( [time[0], time[-1]], [90,90], 'k-')
    ax2.text(time[0],98,'Detector Status', color='green')

    leg = ax2.legend(numpoints=1, loc="upper left")
    leg.get_frame().set_alpha(.7)

######
    if useDB:
        host = 'hawcmon.umd.edu'
        user = 'hawc_user'
        dbname = 'athena'
        table = 'measurements_expstatus_1'
        try:
            f = open(args.pwd, "r")
            #f = open("/var/www/html/hawcview/pwd.cfg", "r")
            password = f.readline().strip()
        except:
            print '**Create file at $HOME/etc/pwd.cfg (or at a custom location using the pwd argument) which contains the database password to avoid this prompt***'
            password = raw_input('Enter database password: ')
        try:
            connection = mysql.connect(host=host,user=user,passwd=password,db=dbname)
            cursor = connection.cursor()
        except Exception,e:
            print 'Failed to connect to %s database: %s' % (dbname, e)
            raise Exception
        try:
            mjd0 = ModifiedJulianDate(mjd.get_date() - args.tmax*U.hour)
            mjd1 = ModifiedJulianDate(mjd.get_date() + args.tmax*U.hour)
            utc0 = UTCDateTime(mjd0.timestamp)
            utc1 = UTCDateTime(mjd1.timestamp)
            uni0 = utc0.unixsecond
            uni1 = utc1.unixsecond
            
            # run a select query
            query  = "select runNumber, measurementTime, runStatus"
            query += " from %s" % table
            query += " WHERE runStatus LIKE '\\'RUNNING\\'' AND measurementTime > %d AND measurementTime < %d" % (uni0,uni1)
            query += " ORDER BY runNumber, measurementTime;"
            cursor.execute(query)
            results = cursor.fetchall()
        except Exception, e:
            results = []
            print "Failed to query table %s : %s" % (table, e)

       
        if len(results) > 0: 
            sTimes = []
            runs = []
            t0 = results[0][1]
            t1 = results[0][1]
            prevRun = results[0][0] 
            run    = results[0][0]
            for r in results:
                run = r[0]
                unixT = r[1]
            
                if run != prevRun:
                   #runs.append(run)
                   d0 = unixTimeToDatetime(t0,0)
                   d1 = unixTimeToDatetime(t1,0)
                   sTimes.append( [d0,d1] ) 
                   t0=unixT
            
                t1 = unixT
                prevRun = run
            
            #runs.append(run)
            d0 = unixTimeToDatetime(t0,0)
            d1 = unixTimeToDatetime(t1,0)
            sTimes.append( [d0,d1]) 

            d0prev=time[0]
            d1prev=time[0]
            
            at2 = ax2.twinx()
            for d0,d1 in sTimes:

                at2.fill_between( [d0,d1], [0,0], [3.6,3.6],facecolor='g', edgecolor='None')
                d0prev = d0
                d1prev = d1

            at2.set_ylim( [0, 40] )
            at2.set_xlim( [ time[0], time[-1] ] )
            at2.set_yticks( [] )
        else:
            sTimes = []
            print "No statuses available"


        dbname = "QMDB"
        table = 'hawconline'
        results = []
        try:
            connection = mysql.connect(host=host,user=user,passwd=password,db=dbname)
            cursor = connection.cursor()
        except Exception,e:
            print "Failed to connect to %s database: %s" % (dbname, e)
            raise Exception
        try:
            # run a select query
            query  = "select Run_number,SubRun_number, EventCount, UTCTimeStart_datetime, UTCTimeEnd_datetime"
            query += " from %s" % table
            query += " WHERE GPSTimeStart_sec>0 AND UTCTimeStart_datetime > \"%s\" AND UTCTimeEnd_datetime < \"%s\"" % (utc0,utc1)
            query += " ORDER BY Run_number, SubRun_number;"
            cursor.execute(query)
            results = cursor.fetchall()
        except Exception,e:
            print "Failed to query table %s: %s" % (table, e)


        if len(results) > 0:
            tcens = []
            rate_kHz  = []
            for r in results:
                counts = int( r[2] )
                d0     =  r[3] 
                d1     =  r[4] 
                if d0 != d1:
                  rate_kHz.append( counts/(d1-d0).total_seconds()/1000. )
                  tcens.append(    d0 + (d1-d0)/2 )

            if len(sTimes) == 0:
                at2 = ax2.twinx()
                at2.set_ylim( [0, 40] )
                at2.set_xlim( [ time[0], time[-1] ] )
            at2.plot(tcens, rate_kHz,'r.',mec="None",alpha=0.7)
            at2.set_ylabel("Detector Rate [kHz]", color='red')
            at2.set_yticks( np.arange(0,40.1,10) )
        else:
            print "\nNo rates available\n"
    
    if args.output:
       fig.savefig(args.output, dpi=400)
       f = open("%s.log" % os.path.splitext(args.output)[0], "w")
       f.write(logfile)
       f.close()
    else:
        plt.show()

if __name__ == "__main__":
    main()

