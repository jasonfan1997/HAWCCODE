#!/usr/bin/env python
################################################################################
# Plot HEALPix map in FITS column format using a Mollweide projection.  The
# healpy and numpy packages need to be installed.
#
# The plot has built in masks and rescaling optimized for the HAWC detector
# exposure.
################################################################################

__version__ = "$Id: plotMollweide.py 41205 2017-10-19 21:36:26Z sybenzvi $"

try:
    import argparse
    import numpy as np
    import healpy as hp
    from math import log

    import os

    import matplotlib as mpl

    if not os.getenv("DISPLAY"):
        mpl.use('AGG')

    import matplotlib.pyplot as plt
    import re
    
    import MapPalette
except ImportError as e:
    print(e)
    raise SystemExit

# Try importing TeVCat
try:
    import TeVCat
    haveTeVCat = True
except ImportError as e:
    haveTeVCat = False
    print(e)

# Try importing the Fermi catalogs
try:
    import FGLCatalog
    import FHLCatalog
    haveFermi = True
except ImportError as e:
    haveFermi = False
    print(e)

# Try importing precession functions
try:
    import precess
    canPrecess = True
except ImportError as e:
    canPrecess = False
    print(e)
    
# Try importing the Gumble distribution
try:
    from scipy.stats import gumbel_r
    haveGumbel = True
except ImportError as e:
    haveGumbel = False
    print(e)

# Try importing various ephem
try:
    import ephem
    haveEphem = True
except ImportError as e:
    haveEphem = False
    print(e)

degree = np.pi/180.0

def gal2equ(l, b):
    """Convert Galactic coordinates to equatorial coordinates"""
    ra, dec = [0, 0]
    if haveEphem:
        gal = ephem.Galactic(l*degree, b*degree)
        ra, dec = [x/degree for x in gal.to_radec()]
    return ra, dec

def maskHAWC(skymap):
    """Set pixels outside the HAWC mask region to UNSEEN"""
    npix = skymap.size
    nside = hp.npix2nside(npix)

    pxls = np.arange(skymap.size)
    nZeroPix = pxls[(skymap != 0)]
    pxTh, pxPh = hp.pix2ang(nside,pxls)
    
    # Mask outside FOV
    cut = (pxTh > pxTh[nZeroPix].max()) | (pxTh < pxTh[nZeroPix].min()) \
          | (skymap == 1e20)
    skymap[cut] = hp.UNSEEN

    return skymap

def setupColorbar(fig, title=None, ticks=None, coord="C"):
    """Create the color bar for a HEALPix map image."""
    for ax in fig.get_axes():
        if type(ax) is hp.projaxes.HpxMollweideAxes:
            cb = fig.colorbar(ax.get_images()[0], ax=ax,
                              orientation="horizontal",
                              shrink=0.8,
                              aspect=40,
                              pad=0.05,
                              fraction=0.1,
                              ticks=ticks,
                              format=mpl.ticker.FormatStrFormatter("%g"))
            for label in cb.ax.get_xticklabels():
                label.set_fontsize("large")
            cb.set_label(title, size="large")
            if coord == "C" or coord == "E":
                ax.annotate("0$^\circ$", xy=(1.8, 0.625), size="large")
                ax.annotate("360$^\circ$", xy=(-1.95, 0.625), size="large")
            else:
                ax.annotate("-180$^\circ$", xy=(1.65, 0.625), size="large")
                ax.annotate("180$^\circ$", xy=(-1.9, 0.625), size="large")


def main():
    p = argparse.ArgumentParser(description="Plot maps in Mollweide projection")
    p.add_argument("fitsfile", nargs="?", help="Input HEALPix FITS file")
    p.add_argument("-c", "--coords", default="C",
                   help="C=equatorial, G=galactic, E=ecliptic")
    p.add_argument("-C", "--column", dest="fitscol", type=int, default=0,
                   help="FITS file column number of HEALPix map")
    p.add_argument("-D", "--coldiff", default=-1, type=int,
                   help="FITS column in which healpix map to be subtracted "
                   "(if any) resides")
    p.add_argument("--file-diff", dest="diff", default="", type=str,
                   help="FITS File to take the difference from (optional) otherwise uses positional argument")
    p.add_argument("--decMin", type=float, default=None,
                   help="Plot minimum declination value")
    p.add_argument("--decMax", type=float, default=None,
                   help="Plot maximum declination value")
    p.add_argument("--mjd", default=None, type=float,
                   help="MJD of the map. Will be plotted J2000")
    p.add_argument("--norm-sqdeg", dest="sqdeg", action="store_true",
                   help="Normalize values to square degree, according to nSides "
                   "(makes sense only certain dinds of maps)")
    p.add_argument("-L", "--label", default=None,
                   help="Color bar label")
    p.add_argument("-m", "--min", type=float, default=None,
                   help="Plot minimum value")
    p.add_argument("-M", "--max", type=float, default=None,
                   help="Plot maximum value")
    p.add_argument("-n", "--ncolors", type=int, default=256,
                   help="Number of colors to use in the color map")
    p.add_argument("-o", "--output", default=None,
                   help="Output image file (.png preferred)")
    p.add_argument("-s", "--scale", type=float, default=1.0,
                   help="Scale the map after input by some value")
    p.add_argument("-S", "--size", dest="xsize", type=int, default=1600,
                   help="Map x-axis pixel number; increase for high resolution")
    p.add_argument("-t", "--ticks", nargs="+", type=float,
                   help="List of ticks for color bar, space-separated")
    p.add_argument("-T", "--title", default="",
                   help="Plot title")
    p.add_argument("-x", "--threshold", dest="thresh", default=None, type=float,
                   help="Apply a (lower) threshold to the map")
    p.add_argument("-X", "--abthresh", dest="athresh", default=None, type=float,
                   help="Apply an absolute (lower/upper) threshold to the map")

    p.add_argument("--gamma", action="store_true", 
                   help="Use Fermi/HESS/VERITAS gamma-ray color scale")
    p.add_argument("--milagro", action="store_true", dest="milagro",
                   help="Use Milagro color scale")

    p.add_argument("--fermicat", default=None,
                   help="Fermi xFGL catalog FITS file name")
    p.add_argument("--fermicat-labels", dest="fermicatLabels",
                  action="store_true",
                  help="Draw Fermi sources with labels.")
    p.add_argument("--download-tevcat", dest="dtevcat", action="store_true",
                   help="Download data from tevcat.uchicago.edu and exit")
    p.add_argument("--tevcat", default=None,
                   help="Draw TeVCat sources using TeVCat ASCII file")
    p.add_argument("--tevcat-labels", dest="tevcatLabels", action="store_true",
                   help="Draw TeVCat sources with labels w/ TeVCat ASCII file")

    # Highlight hotspots listed in an ASCII file
    p.add_argument("--hotspots", default=None,
                   help="Hotspot coordinates in an ASCII file( Format:   source name, ra, dec   )")
    p.add_argument("--hotspot-labels", dest="hotspotLabels",
                   action="store_true", 
                   help="Draw hotspots sources with labels.")

    p.add_argument("--glines", action="store_true",
                   help="Draw galactic lines of latitude (celestial only)")
    p.add_argument("--gridlines", action="store_true",
                   help="Draw gridlines")
    p.add_argument("--gplane", action="store_true",
                   help="Draw the galactic plane and galactic center")
    p.add_argument("--nomask", action="store_false", dest="mask",
                   default=True,
                   help="Mask part of the sky")
    p.add_argument("--preliminary", action="store_true", 
                   help="Add a watermark indicating the plot as preliminary")
    p.add_argument("--alpha", type=float, default=0.7,
                   help="Alpha used for drawn Galactic lines.")
    p.add_argument("--logz", action="store_true",
                   help="Set color scale in log.")

    # Plot P-Value instead of sigma, expects location and scale of
    # a right-skewed Gumble distribution
    p.add_argument("--gumbel", type=float, nargs=2,
                   help="Plot P-Value instead of significance. "
                   "Arguments: location and scale of right-skewed "
                   "Gumble distribution.")

    args = p.parse_args()

    #Sanity checks
    if (args.mjd is not None) and (not canPrecess):
        print("Missing necessary packages to make precession")
        raise SystemExit
    
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
    if args.fitsfile == None:
        print("Please specify an FITS file to plot.")
        raise SystemExit

    znorm = None
    if (args.logz):
        znorm = 'log'

    skymap, skymapHeader = hp.read_map(args.fitsfile, args.fitscol, h=True)
    nside1 = hp.get_nside(skymap)
    # remove naughty values
    skymap[np.isnan(skymap)] = hp.UNSEEN

    if args.decMax:
        imin = hp.ang2pix(nside1,90*degree-args.decMax*degree,0)
        skymap[0:imin] = hp.UNSEEN
    if args.decMin:
        imax = hp.ang2pix(nside1,90*degree-args.decMin*degree,0)
        skymap[imax:] = hp.UNSEEN

    if args.coldiff > -1:
        if os.path.isfile(args.diff):
            print("Taking difference with {0}".format(args.diff))
            skymap2 = hp.read_map(args.diff, args.coldiff)
        else:
            print("No extra file provided, using same file as input")
            skymap2 = hp.read_map(args.fitsfile, args.coldiff)

        print("Subtracting column {0} from skymap...".format(args.coldiff))
        skymap -= skymap2

    if (args.gumbel):
        assert haveGumbel
        gumbel_location, gumbel_scale = args.gumbel
        gumbel = gumbel_r(loc=gumbel_location, scale=gumbel_scale)
        skymap = gumbel.logsf(skymap)/log(10)

        def inf_suppressor(x):
            return x if np.isfinite(x) else 0.

        inf_suppressor = np.vectorize(inf_suppressor)
        skymap = inf_suppressor(skymap)

    if args.sqdeg:
        # Normalize value to square degree
        pixsizestr = 4*np.pi / (12*nside1**2)
        str2sqdeg = (180/np.pi)**2
        pixsizesqdeg = pixsizestr * str2sqdeg
        skymap /= pixsizesqdeg

    # I did not find header handler, thats all I came up with...
    toFind = 'TTYPE' + str(args.fitscol+1)
    skymapName = dict(skymapHeader)[toFind]
    skymap[skymap != hp.UNSEEN] *= args.scale

    # Find FOV
    pxls = np.arange(skymap.size)
    nZeroPix = pxls[(skymap != 0)]
    pxTh, pxPh = hp.pix2ang(nside1,pxls)
    
    # Mask outside FOV
    if args.mask:
        skymap = maskHAWC(skymap)

    print("Trials: %d" % skymap[skymap != hp.UNSEEN].size)
    print("Counts: %g" % skymap[skymap != hp.UNSEEN].sum())

    # Set up the map minimum and maximum values
    notMask = skymap[skymap > hp.UNSEEN]
    dMin, dMax = (notMask[np.argmin(notMask)], notMask[np.argmax(notMask)])
    print(dMin, dMax)
    if args.min is not None:
        dMin = args.min
        skymap[(skymap<dMin) & (skymap > hp.UNSEEN)] = dMin
    if args.max is not None:
        dMax = args.max
        skymap[(skymap>dMax) & (skymap > hp.UNSEEN)] = dMax

    mpl.rc("font", family="serif", size=14)

    # Set up the color map for plotting
    textcolor, colmap = MapPalette.setupDefaultColormap(args.ncolors)

    # Use the Fermi/HESS/VERITAS purply-red-yellow color map
    if args.gamma:
        textcolor, colmap = MapPalette.setupGammaColormap(args.ncolors)

    # Use the Milagro color map
    if args.milagro:
        dMin = args.min if args.min != None else -5
        dMax = args.max if args.max != None else 15
        thresh = args.thresh if args.thresh != None else 2.
        textcolor, colmap = \
            MapPalette.setupMilagroColormap(dMin, dMax, thresh, args.ncolors)
    # Use a thresholded grayscale map with colors for extreme values
    else:
        if args.thresh != None:
            textcolor, colmap = \
                MapPalette.setupThresholdColormap(dMin, dMax, args.thresh,
                                                  args.ncolors)
        elif args.athresh != None:
            textcolor, colmap = \
                MapPalette.setupAbsThresholdColormap(dMin, dMax, args.athresh,
                                                     args.ncolors)

    # Set up the figure frame and coordinate rotation angle
    if args.coords == "C":
        fig = plt.figure(1, figsize=((11.75,6.5)), dpi=100)
        rotation = 180.
        coords = "C"
    elif args.coords == "G":
        fig = plt.figure(1, figsize=(9.5,6), dpi=100)
        coords = "CG"
        rotation = 0.
    else:
        fig = plt.figure(1, figsize=(9.5,6), dpi=100)
        coords = "CE"
        rotation = 180.

    if args.mjd is None:
        rotMap = rotation
    else:
        rotMap=precess.mjd2J2000ang(mjd=args.mjd,coord=coords,rot=rotation)
        rotPrec=hp.rotator.Rotator(rot=precess.mjd2J2000ang(mjd=args.mjd))
        
    hp.mollview(skymap, fig=1,
                        xsize=args.xsize,
                        coord=coords,
                        min=dMin,
                        max=dMax,
                        rot=rotMap,
                        title=args.title.replace("\\n","\n"),
                        cmap=colmap,
                        cbar=False,
                        notext=True,
                        norm=znorm)
    if args.gridlines:
      hp.graticule()

    if args.preliminary:
        hp.projtext(95*degree, 240*degree,
                    "PRELIMINARY",
                    coord=coords,
                    color=textcolor,
                    alpha=0.5,
                    rotation=0,
                    fontdict={"family":"sans-serif", "weight":"bold", "size":42})

    # If TeVCat data are available, plot them
    if args.tevcat:
        if haveTeVCat:
            try:
                tevcat = TeVCat.TeVCat(args.tevcat)
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

                if args.mask:
                    cut = np.logical_and(dec<64*degree, dec>-26*degree)
                    dec = dec[cut]
                    ra = ra[cut]
                    assoc = assoc[cut]

                if args.mjd is None:
                    theta,phi=(90*degree - dec,ra)
                else:
                    theta,phi=rotPrec.I(90*degree - dec,ra)
                    
                hp.projscatter(theta, phi, coord=coords,
                               color=textcolor,
                               facecolors="none",
                               marker="s")
                if args.tevcatLabels:
                    for th, ph, s in zip(theta, phi, assoc):
                        hp.projtext(th, ph, s, coord=coords,
                                    color=textcolor,
                                    rotation=10,
                                    fontdict={"size":14})
                # Print significance for TeVCat sources above 3 sigma
                sThr = 3.
                if cId == 1:
                    print("TeVCat sources above %.2f sigma:" % sThr)
                for r, d, s in zip(ra, dec, assoc):
                    valueTeVCat = hp.get_interp_val(skymap, 90*degree - d, r)
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

            if args.mjd is None:
                theta,phi=(90*degree - dec,ra)
            else:
                theta,phi=rotPrec.I(90*degree - dec,ra)
            
            for i in xrange(len(assoc)):
                if assoc[i] == '':
                    assoc[i] = catnms[i]

            if args.mask:
                cut = np.logical_and(dec<64*degree, dec>-26*degree)
                dec = dec[cut]
                ra = ra[cut]
                assoc = assoc[cut]
            hp.projscatter(theta, phi, coord=coords,
                           color=textcolor,
                           facecolors="none",
                           marker="o")
            if args.fermicatLabels:
              for th, ph, s in zip(theta, phi, assoc):
                hp.projtext(th, ph, s, coord=coords,
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
        assoc = []
        for line in fhot:
            if line.startswith("#"):
                continue
            r, d = [float(t)*degree for t in line.strip().split(",")[1:3]]
            ra.append(r)
            dec.append(d)
            assoc.append(line.strip().split(",")[0].strip())
        ra = np.array(ra)
        dec = np.array(dec)
        assoc = np.array(assoc)

        if args.mjd is None:
            theta,phi=(90*degree - dec,ra)
        else:
            theta,phi=rotPrec.I(90*degree - dec,ra)
            
        if args.mask:
            cut = np.logical_and(dec<64*degree, dec>-26*degree)
            dec = dec[cut]
            ra = ra[cut]
            assoc = assoc[cut]
            theta = theta[cut]
            phi = phi[cut]
        hp.projscatter(theta, phi, coord=coords,
                       color='darkorchid',
                       facecolors="none",
                       marker="o",
                       s=40)
        if args.hotspotLabels:
            for th, ph, r, d, s in zip(theta, phi, ra, dec, assoc):
                print(r, d, s)
                hp.projtext(th, ph, s+'  .', coord=coords,
                            color='darkorchid',
                            rotation=4,
                            fontdict={'family': 'sans-serif',
                                      'size': 10,
                                      'weight': 'bold'})

    # Adjust plot so that empty southern hemisphere is cut out
    ax = fig.gca()
    # Hide not used space of image
    imgp = ax.get_images()[0]
    imgp.get_cmap().set_under('w',alpha=0.)
    imgp.get_cmap().set_over('w',alpha=0.)
    imgp.get_cmap().set_bad('w',alpha=0.)

    if args.coords == "C":
        ax.set_ylim([-0.5,1])

    # Draw the galactic plane, latitude lines, and the galactic center
    # (only if underlying grid is celestial)
    if args.glines and args.coords == "C":
        for b in np.arange(-90., 90.1, 30):
            raGP = []
            decGP = []
            if b == 0:
                ls = "-"
                lw = 2
            else:
                ls = "--"
                lw = 1
            for l in np.arange(0., 360., 1.):
                r, d = gal2equ(l, b)
                raGP.append(r)
                decGP.append(d)

            if args.mjd is not None:    
                raGP,decGP = rotPrec.I(raGP,decGP,lonlat=True)
                
            hp.projplot([raGP, decGP], lonlat=True, color=textcolor, coord="C",
                        alpha=args.alpha,linestyle=ls, linewidth=lw)

        galCenter=gal2equ(0., 0.);
        if args.mjd is not None:
            galCenter=rotPrec.I(galCenter,lonlat=True)
            
        hp.projscatter(gal2equ(0., 0.), lonlat=True, color=textcolor,
                       alpha=args.alpha, s=50, marker='o',coord="C")

        # Label the galactic latitude lines
        marks = [-60, -30, 0, 30, 0, -30, -60]
        xposn = [-1.56,-1.20,-0.79, 0.10, 0.60, 0.97, 1.33]
        for m, x in zip(marks, xposn):
            ax.annotate("%d$^\circ$" % m, xy=(x,-0.475), size="small",
                        fontweight="bold", color=textcolor)
    # Draw lines above/below galactic plane, and the galactic center
    # (only if underlying grid is celestial)
    elif args.gplane and args.coords == "C":
        #for b in np.arange(-90., 90.1, 30):
        for b in [-5, 5]:
            raGP = []
            decGP = []
            if b == 0:
                ls = "-"
                lw = 2
            else:
                ls = "--"
                lw = 1
            for l in np.arange(0., 360., 1.):
                r, d = gal2equ(l, b)
                raGP.append(r)
                decGP.append(d)

            if args.mjd is not None:    
                raGP,decGP = rotPrec.I(raGP,decGP,lonlat=True)

            hp.projplot([raGP, decGP], lonlat=True, color=textcolor, coord="C",
                        linestyle=ls, linewidth=lw, alpha=args.alpha)

        galCenter=gal2equ(0., 0.);
        if args.mjd is not None:
            galCenter=rotPrec.I(galCenter,lonlat=True)
            
        hp.projscatter(galCenter, lonlat=True, color=textcolor,
                       s=50, marker='o', alpha=args.alpha,
                       coord="C")

        # Label the galactic latitude lines
        marks = [-5, 5, 5, -5]
        xposn = [-1.05, -0.7, 0.4, 0.7]
        for m, x in zip(marks, xposn):
            ax.annotate("%d$^\circ$" % m, xy=(x,-0.475), size="small",
                        fontweight="bold", color=textcolor)

    # Set up the color bar and tick axis
    if args.label:
        paletteLabel = args.label
    else:
        if re.match("significance", skymapName):
            paletteLabel = r"significance [$\sigma$]"
        else:
            paletteLabel = skymapName
    if args.gumbel:
        paletteLabel = "log10(P-Value)"

    # Clean up automatic tick definitions
    cticks = args.ticks
    if args.ticks == None:
        clrmin = np.floor(dMin)
        clrmax = np.ceil(dMax)
        ctspc = np.round((clrmax - clrmin)/10.)
        if ctspc == 0.:
            ctspc = 1.
        if clrmin < 0 and clrmax > 0:
            cticks = -np.arange(0,-clrmin,ctspc)
            cticks = np.sort(np.unique(np.append(cticks, np.arange(0, clrmax, ctspc))))
        else:
            cticks = np.arange(clrmin,clrmax+ctspc,ctspc)

    setupColorbar(fig, title=paletteLabel, ticks=cticks, coord=args.coords)

    if args.output:
        fig.savefig(args.output, dpi=400)
        print("File %s created" % args.output)
    else:
        plt.show()

if __name__ == "__main__":
    main()

