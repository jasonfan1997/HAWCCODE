#!/usr/bin/env python
################################################################################
# Plot Histogram and Gaussian Fit of a FITS Healpix Map
#
# Author: Dan Fiorino (dan.fiorino@icecube.wisc.edu) and Segev BenZvi
################################################################################

__version__ = "$Id: makeSignificanceHistogram.py 43395 2018-06-03 02:22:55Z henrike $"

try:
    import argparse
    import healpy as hp
    import numpy  as np
    import scipy.stats
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    from math import sqrt, pi, erf

    degree = pi/180.
except ImportError,e:
    print e
    raise SystemExit

def equ2gal(ra, dec):
    """Convert equatorial (J2000.0) coordinates to galactic.  This is taken
    from the eqgal.f program in SLALIB.

    Inputs expected in degrees.
    """
    rmtx = np.matrix([[-0.054875539726,  0.494109453312, -0.867666135858],
                     [-0.873437108010, -0.444829589425, -0.198076386122],
                     [-0.483834985808,  0.746982251810,  0.455983795705]])
    cosr = np.cos(ra*degree)
    sinr = np.sin(ra*degree)
    cosd = np.cos(dec*degree)
    sind = np.sin(dec*degree)
    evec = np.matrix([[cosr*cosd], [sinr*cosd], [sind]])
    gvec = rmtx.transpose() * evec

    x, y, z =  (gvec.item(0), gvec.item(1), gvec.item(2))
    r = np.sqrt(x*x + y*y)
    l = 0.
    b = 0.
    if r != 0.:
        l = np.arctan2(y, x) / degree
        if l < 0:
            l += 360.
    if z != 0:
        b = np.arctan2(z, r) / degree
    return (l, b)

def maskFOV(skymap):
    """Set pixels outside the FOV mask region to UNSEEN"""
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


if __name__ == "__main__":
    # Set up input
    p = argparse.ArgumentParser(description="Make significance map histogram")
    p.add_argument("fitsfile", nargs="?", help="Input HEALPix FITS file")
    p.add_argument("-C", "--column", dest="col", default=0, type=int,
                   help="FITS column in which healpix map resides")
    p.add_argument("-D", "--coldiff", dest="coldif", default=-1, type=int,
                   help="FITS column in which healpix map to be subtracted "
                   "(if any) resides")
    p.add_argument("-m", "--minDec", dest="minDec", type=float, default=-90.,
                   help="Minimum sky map declination")
    p.add_argument("-M", "--maxDec", dest="maxDec", type=float, default=90.,
                   help="Maximum sky map declination")
    p.add_argument("--minRA", dest="minRA", type=float, default = None,
                   help="Minimum sky map right ascension")
    p.add_argument("--maxRA", dest="maxRA", type=float, default= None,
                   help="Maximum sky map right ascension")
    p.add_argument("-o", "--output", dest="output",
                   help="image output name")
    p.add_argument("-s", "--minSig", dest="minSig", type=float, default=None,
                   help="min significance for histogram (exclusive)")
    p.add_argument("-S", "--maxSig", dest="maxSig", type=float, default=None,
                   help="max significance for histogram (exclusive)" )
    p.add_argument("-T", "--title",  dest="title", default="1D Significance Histogram",
                   help="Map title for legend")
    p.add_argument("-X", "--xtitle",  dest="xtitle", default=None,
                   help="X-axis title")
    p.add_argument("--batch", action="store_true", dest="batch",
                   default=False, help="only display width and mean of fit")
    p.add_argument("--excludeArea", action="store_true", dest="excludeArea",
                   default=False,
                   help="Exclude a circular area from hist")
    p.add_argument("--plotErrors", action="store_true", dest="plotErrors",
                   default=False,
                   help="Plot statistical error bars")
    p.add_argument("-b", "--minb", dest="minb", type=float,
                   help="Minimum galactic latitude")
    p.add_argument("-B", "--maxb", dest="maxb", type=float,
                   help="Maximum galactic latitude")
    p.add_argument("-l", "--minl", dest="minl", type=float,
                   help="Minimum galactic longitude")
    p.add_argument("-L", "--maxl", dest="maxl", type=float,
                   help="Maximum galactic longitude")
    p.add_argument("-r", "--range", dest="range", default="I",
                   help="I=inside (b,B), O=outside of (b,B)")
    p.add_argument("--multipole",dest="multipole",action="store_true",
                   default=False, help="significance from multipole fit")
    p.add_argument("--includeArea", action="store_true", dest="includeArea",
                   default=False,
                   help="Only include a circular area in hist")
    p.add_argument("--templateROI", dest="templateROI", help="Fits file with ROI")
    p.add_argument("--ROIthreshold", dest="ROIthreshold", help="Threshold for ROI", default=0.5, type=float)
    p.add_argument("--ra", dest="targetRA", type=float, default=83.63,
                   help="Target RA for --inclueArea or --excludeArea")
    p.add_argument("--dec", dest="targetDec", type=float, default=22.01,
                   help="Target Dec for --inclueArea or --excludeArea")
    p.add_argument("--binsize", dest="binsize", type=float, default=20.,
                   help="Circular bin radius [deg] (use with --includeArea)")
    p.add_argument("--nomask", action="store_false", dest="mask",
                   default=True,
                   help="Mask part of the sky")
    # Reinterpret sigma as Test Statistics
    p.add_argument('--teststat', dest='teststat', default=False,
                   action='store_true',
                   help='Display test statistics')
    p.add_argument('--realNorm', dest='realNorm', default=False,
                   action='store_true',
                   help='Obtain the norm of the expected distribution from the integral of the data histogram, instead of from the fit')
    p.add_argument('--noFit', dest='plotFit', default=True,
                   action='store_false',
                   help='Do not show fitted curve and parameters')
    
    args = p.parse_args()
    if not args.plotFit and not args.realNorm:
      print "Must use realNorm if no fit is performed."
      args.realNorm = True
    args.title = args.title.replace('%%MATHMODE%%','$') # Mathmode breaks ROOT
    # Imports containing ROOT objects must be after the parser, because
    # TApplication messes up with sys.args.
    try:
        from ROOT import TH1D
    except ImportError,e:
        print e
        raise SystemExit

    # Start normal processing
    if args.fitsfile == None:
        print("Please specify an FITS file to plot.")
        raise SystemExit

    if args.multipole:
        sigmap,sigmapHeader = hp.read_map(args.fitsfile,3,h=True)
    else:
        sigmap,sigmapHeader = hp.read_map(args.fitsfile,args.col,h=True)
        if args.coldif > -1:
            sigmap2 = hp.read_map(args.fitsfile, args.coldif)
            sigmap -= sigmap2

    if args.teststat:
        def correctSigma(sigma,oldDOF=1,newDOF=2):
            #if sigma > 30:
            #    return sigma-0.15
            return np.copysign(scipy.sqrt(scipy.stats.chi2.isf(scipy.stats.chi2.sf(np.power(sigma,2),newDOF),oldDOF)),sigma)
        vCorrectSigma = np.vectorize(correctSigma)
        sigmap = correctSigma(sigmap)

    # Mask outside FOV
    if args.mask:
      sigmap = maskFOV(sigmap)
    npix   = len(sigmap)
    nside  = hp.npix2nside( npix )
    startpix = hp.ang2pix(nside, pi/2. - args.maxDec*pi/180., 0.)
    endpix   = hp.ang2pix(nside, pi/2. - args.minDec*pi/180., 0.)

    minSig = args.minSig
    if minSig == None:
        minSig = np.floor(sigmap[ (sigmap > hp.UNSEEN)].min())
    else:
        minSig = np.nextafter(minSig,np.inf) #Make it exclusive

    maxSig = args.maxSig
    if maxSig == None:
      maxSig = np.ceil(sigmap[ (sigmap > hp.UNSEEN)].max())


    if (startpix > endpix):
        print "error: switch your bounds"

    # Fill root histogram and array
    nbins   = 500
    hist    = TH1D("hist", ";S;counts", nbins, minSig, maxSig)
    siglist = []

    if args.excludeArea:
        target = hp.ang2pix(nside, np.pi/2.-args.targetDec*degree, args.targetRA*degree)
        sink   = hp.query_disc(nside, hp.pix2vec(nside, target), args.binsize*degree)

    source = None
    if args.includeArea:
        # Save significances with specified deg of the Target RA/Dec
        target = hp.ang2pix(nside,np.pi/2.-args.targetDec*degree, args.targetRA*degree)
        source = hp.query_disc(nside, hp.pix2vec(nside, target), args.binsize*degree)
    
    if args.templateROI:
        roimap = hp.fitsfunc.read_map(args.templateROI )
        print "Reading ROI from %s (threshold %.2f)" %( args.templateROI, args.ROIthreshold  )
        roi_npix   = len(roimap)
        roi_nside  = hp.npix2nside( roi_npix )
        if not roi_npix == npix:
          print "ROI map and significance map do not have the same NSIDE, will try to up/degrade the ROI map."
          roimap = hp.pixelfunc.ud_grade( roimap, nside_out = nside )
        source = np.where( roimap > args.ROIthreshold )[0]
        print np.sum(source)
    
    if source is not None:
        max = 0.
        maxPix =0
        for bin in source:
            siglist.append(sigmap[bin])
            hist.Fill(sigmap[bin])
            if sigmap[bin] > max:
                max = sigmap[bin]
                maxPix = bin
        if args.includeArea:
          print "Max significance within %.2f deg of target at (%.2f, %.2f) is %f" % (args.binsize, args.targetRA, args.targetDec, max)
          print "Significance at target = %f" % sigmap[target]

    else:
        # Loop over whole sky
        for bin in range(startpix, endpix):
            if args.excludeArea and bin in sink:
              continue
              
            if any( t is not None for t in [args.minb, args.maxb, args.minl, args.maxl, args.minRA, args.maxRA]):
              dec, ra = hp.pix2ang(nside, bin)
              ra  = ra/degree
              dec = 90. - dec/degree
              l, b = equ2gal(ra, dec)
              
              if (args.minb is not None and args.maxb is not None):
                if args.range == "I" and not (b>args.minb and b<args.maxb):
                  continue
                if args.range == "O" and not( b<args.minb or b>args.maxb ):
                  continue
 
              if (args.minl is not None and args.maxl is not None):
                if args.range == "I" and not (l>args.minl and l<args.maxb):
                  continue
                if args.range == "O" and not( l<args.minl or l>args.maxb ):
                  continue
 
              if (args.minRA is not None and args.maxRA is not None):
                if args.range == "I" and not (ra>args.minRA and ra<args.maxRA):
                  continue
                if args.range == "O" and not( ra<args.minRA or ra>args.maxRA ):
                  continue
 
            siglist.append(sigmap[bin])
            hist.Fill(sigmap[bin])

    print hist.GetEntries()
    # Get gaussian fit
    if args.plotFit:
      hist.Fit("gaus", "0","",-2.,2.) #only fit +/- 2sigma
      #hist.Fit("gaus", "0")
      fit      = hist.GetFunction("gaus")
      norm     = fit.GetParameter(0)
      mean     = fit.GetParameter(1)
      meanerr  = fit.GetParError(1)
      sigma    = fit.GetParameter(2)
      sigmaerr = fit.GetParError(2)

    if args.realNorm:
        realNorm = hist.Integral()
        print "Integral = %f" % realNorm
        if args.minSig != None:
            if args.maxSig != None:
                realNorm = 2 * realNorm / ( erf(args.maxSig / sqrt(2))  - erf(args.minSig / sqrt(2)) )
            else:
                realNorm = 2 * realNorm / ( 1  - erf(args.minSig / sqrt(2)) )

            print "Integral (corrected for integration boundaries) = %f" % realNorm
        elif args.maxSig != None:
            realNorm = 2 * realNorm / ( 1 + erf(args.maxSig / sqrt(2)) )
            print "Integral (with integration limits correction) = %f" % realNorm

        realNorm = hist.GetBinWidth(1) * realNorm
            
    if args.batch:
      print "mean %f \n width %f" % (mean,sigma)
    else:
        # Plot
        # Plot gaussian fit
      x        = np.arange(minSig, maxSig,
                             (maxSig-minSig)/nbins,
                             dtype=float)
      if args.realNorm:
            y = (realNorm/sqrt(2*np.pi)) * np.exp(-0.5*((x))**2)
      else:
            y = norm * np.exp(-0.5*((x))**2)

      expPlot, = plt.plot(x, y, "r--", linewidth=1,label="Expectation")

      if args.plotFit:
        y2        = norm * np.exp(-0.5*((x-mean)/sigma)**2)
        fitPlot, = plt.plot(x, y2, "g-.", linewidth=1,label="Fit")

        # Plot histogram
      A   = plt.hist(siglist, nbins,
                   (minSig, maxSig),
                   histtype='step',
                   edgecolor='black',label="Data")

        # Plot counting statistic errorbars
      if args.plotErrors:
            data     = A[0]
            stderrs  = []
            for i in range(0,len(data)):
                stderrs.append(sqrt(data[i]))
            plt.errorbar(x, y, fmt=None, yerr=stderrs)

      # Title & axis options
      if args.xtitle:
            xtitle = args.xtitle
      else:
            toFind = 'TTYPE' + str(args.col+1)
            xtitle = dict(sigmapHeader)[toFind]
        
      mpl.rc("font", family="serif", size=14)
      plt.title(args.title.replace("\\n","\n"))
      ax = plt.gca()
      ax.set_xlabel(xtitle)
      ax.set_xlim([minSig, maxSig])
      ax.set_ylabel("Number of Pixels")
      if nside == 64:
              ax.set_ylim([8e-1, 2e3])
      #else:
      #    ax.set_ylim([8e-1, 1.5e5])
      ax.set_yscale("symlog")
      xlms = ax.get_xlim()
      ylms = ax.get_ylim()

      # Legend options
      mpl.rcParams['legend.loc'] = 'best'
        
      dummyHandle = Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none', linewidth=0)

      if args.plotFit:
            lg = ax.legend([A[2][0],expPlot,fitPlot,dummyHandle],["Data",
                                                                  "Expectation",
                                                                  "Fit",
                                                                  "mean = %.3f $\pm$ %.3f \nwidth = %.3f $\pm$ %.3f"%(mean,meanerr,sigma,sigmaerr)])
            fitParam = lg.get_texts()[3]
            fitParam.set_color('green')
            fitParam._fontproperties = lg.get_texts()[0]._fontproperties.copy()
            fitParam.set_size('medium')
      else:
            lg = ax.legend([A[2][0],expPlot,dummyHandle],["Data","Expectation"])
            
      lg.get_frame().set_linewidth(0)

      if args.output:
            plt.savefig(args.output, dpi=300)
      else:
            plt.show()
