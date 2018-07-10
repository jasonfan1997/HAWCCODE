#!/usr/bin/env python
################################################################################
# Plot HEALPix Map in column of FITS file on a Mercator Projection.
# Requires healpy to be installed.
################################################################################

__version__ = "$Id: plotGalPlane.py 28513 2015-12-17 19:41:07Z sybenzvi $"

try:
    import argparse
    import numpy as np
    import healpy as hp
    import os
    from math import log
    import pickle

    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import re

    import MapPalette

    degree = np.pi / 180.
except ImportError,e:
    print e
    raise SystemExit

# Try importing TeVCat
try:
    import TeVCat
    haveTeVCat = True
except ImportError,e:
    haveTeVCat = False
    print e

# Try importing the Fermi catalogs
try:
    import FGLCatalog
    import FHLCatalog
    haveFermi = True
except ImportError,e:
    haveFermi = False
    print e

# Try importing the Gumble distribution
try:
    from scipy.stats import gumbel_r
    haveGumbel = True
except ImportError,e:
    haveGumbel = False
    print e

def main():
    p = argparse.ArgumentParser(description="Plot map with Mercator projection")
    p.add_argument("fitsfile", nargs="?", help="Input HEALPix FITS file")
    p.add_argument("-a", "--abthresh", dest="athresh", default=None, type=float,
                   help="Apply an absolute (lower/upper) threshold to the map")
    p.add_argument("-c", "--coords", dest="coords", default="C",
                   help="C=equatorial, G=galactic")
    p.add_argument("-C", "--column", dest="col", default=0, type=int,
                   help="FITS column in which healpix map resides")
    p.add_argument("-D", "--coldiff", dest="coldif", default=-1, type=int,
                   help="FITS column in which healpix map to be subtracted "
                   "(if any) resides")
    p.add_argument("-l", "--threshold", dest="thresh", type=float, default=None,
                   help="Apply a lower threshold to the plot.")
    p.add_argument("--norm-sqdeg", dest="sqdeg", action="store_true",
                   default=False,
                   help="Normalize values to square degree, according to nSides "
                   "(makes sense only certain dinds of maps)")
    p.add_argument("--milagro", action="store_true", dest="milagro", default=False,
                   help="Use Milagro color scale")
    p.add_argument("--contours", dest="contours", type=float, nargs=3, 
                   help="plot 3 contours, e.g. --contours 1 2 3 [sigma]")
    p.add_argument("--interpolate", action="store_true", dest="interp",
                   default=False,
                   help="Use bilinear interpolation to smear HEALPix pixels")
    p.add_argument("--gamma", action="store_true", dest="gamma", default=False,
                   help="Use GeV/TeV gamma-ray color scale")
    p.add_argument("-L", "--label", dest="label",
                   default=None, help="Color bar label")
    p.add_argument("--moon", action="store_true", default=False,
                   dest="isMoon", help="is in moon centered coordinates")
    p.add_argument("-m", "--min", dest="min", type=float, default=None,
                   help="Plot minimum value")
    p.add_argument("-M", "--max", dest="max", type=float, default=None,
                   help="Plot maximum value")
    p.add_argument("-n", "--ncolors", dest="ncolors", type=int, default=256,
                   help="Number of contours to use in the colormap")
    p.add_argument("-o", "--output", dest="output", default=None,
                   help="Output file name (optional)")
    p.add_argument("--preliminary", action="store_true", dest="preliminary",
                   default=False,
                   help="Add a watermark indicating the plot as preliminary")
    p.add_argument("-s","--scale",dest="scale",type=float,default=1.,
                   help="scale up or down values in map")
    p.add_argument("--sun", action="store_true", default=False,
                   dest="isSun", help="is in Sun centered coordinates")
                   
    p.add_argument("--dpar", dest="dpar", type=float,default=5.,
                  help="Interval in degrees between parallels")
    p.add_argument("--dmer", dest="dmer", type=float, default=10.,
                  help="Interval in degrees between meridians.")
    p.add_argument("--nogrid", action="store_true", default=False,
                   dest="nogrid", help="Do not plot gridlines.")
    p.add_argument("--alpha", action="store_true", dest="alpha",
                   default=0.25,
                   help="Alpha used for drawn lines.")
    p.add_argument("--interpolation", dest="interpolation", \
                  default=False, action='store_true',
                  help="Uses bilinear interpolation in data map.")
    p.add_argument("--xsize", dest="xsize", type=int, default=1000,
                  help="Number of X pixels, Y will be scaled acordingly.")
    p.add_argument("-t", "--ticks", dest="ticks", nargs="+", type=float,
                   help="List of ticks for color bar, space-separated")
    p.add_argument("-T","--title",dest="title",
                   help="title of map")

    # Mutually exclusive option to specify plot xy range OR center + WxH
    argRange = p.add_mutually_exclusive_group(required=False)
    argRange.add_argument("--xyrange", dest="xyrange", type=float, nargs=4,
                          help="Xmin Xmax Ymin Ymax plot range [degree]")
    argRange.add_argument("--origin", dest="origin", type=float, nargs=4,
                          help="Central (x,y) coords + width + height [degree]")

    # Plot P-Value instead of sigma, expects location and scale of
    # a right-skewed Gumble distribution
    p.add_argument("--gumbel", dest="gumbel", type=float, nargs=2,
                   help="Plot P-Value instead of significance. "
                   "Arguments: location and scale of right-skewed "
                   "Gumble distribution.")

    # Download/plotting options for Fermi catalog sources
    p.add_argument("--fermicat", dest="fermicat", default=None,
                   help="Fermi xFGL catalog FITS file name")
    p.add_argument("--fermicat-labels", dest="fermicatLabels",
                  action="store_true", default=False,
                  help="Draw Fermi sources with labels.")
  

    # Download/plotting options for TeVCat sources
    p.add_argument("--download-tevcat", dest="dtevcat", action="store_true",
                   default=False,
                   help="Download data from tevcat.uchicago.edu and exit")
    p.add_argument("--tevcat", dest="tevcat", default=None,
                   help="Draw TeVCat sources using TeVCat ASCII file")
    p.add_argument("--tevcat-labels", dest="tevcatLabels",
                   action="store_true", default=False,
                   help="Draw TeVCat sources with labels.")

    # Highlight hotspots listed in an ASCII file
    p.add_argument("--hotspots", dest="hotspots", default=None,
                   help="Hotspot coordinates in an ASCII file")

    args = p.parse_args()

    # Download TeVCat
    if args.dtevcat:
      if haveTeVCat:
        print "Fetching data from tevcat.uchicago.edu."
        tevcat = TeVCat.TeVCat()
        return None
      else:
        print "Sorry, AERIE TeVCat python module is not available."
        raise SystemExit


    # Start normal processing
    if args.fitsfile == None:
      print "Please specify an FITS file to plot."
      raise SystemExit

    # Read in the skymap and mask out empty pixels
    skymap, skymapHeader = hp.read_map(args.fitsfile, args.col, h=True)
    # remove naughty values
    skymap[np.isnan(skymap)] = hp.UNSEEN
    skymap *= args.scale
    nside1 = hp.get_nside(skymap)
    npix   = hp.nside2npix(nside1)
    
    if args.coldif > -1:
        skymap2 = hp.read_map(args.fitsfile, args.coldif)
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
    toFind = 'TTYPE' + str(args.col+1)
    skymapName = dict(skymapHeader)[toFind]

    # Find FOV
    pxls = np.arange(skymap.size)
    nZeroPix = pxls[(skymap != 0)]
    pxTh, pxPh = hp.pix2ang(nside1,pxls)

    # Mask outside FOV
    values = np.ma.masked_where((pxTh > pxTh[nZeroPix].max())
                                | (pxTh < pxTh[nZeroPix].min()), skymap)
                                
    # Plot the skymap as an image, setting up the color palette on the way
    mpl.rc("font", size=14, family="serif")

    coords = ["C","G"]
    rotation = 0.
    figsize = (9,9)
    fig   = plt.figure(num=1,figsize=figsize)

    # Defin longitudinal ranges to plot
    #lonrngs=[[-10,50],[50,110],[110,170]]
    lonrngs=[[-180,-90], [-10,80], [80,170]]
    # Get extrema
    imgpix = np.array([],dtype=int)
    for cxm,cxx in lonrngs:
      angverts = [[cxm,-10],[cxx,-10],\
                  [cxx,10],[cxm,10]]
      vertlist = [ ]
      cRot = hp.Rotator(coord=coords[::-1])
      for x,y in angverts:
        ctht,cph = np.deg2rad((y,x))
        ctht = 0.5*np.pi  - ctht
        if cph < 0:
          cph = 2.*np.pi + cph
        vertlist.append(cRot(hp.ang2vec(ctht,cph)))
      # Get pixels in image
      imgpix = np.concatenate((imgpix,hp.query_polygon(nside1,
                                vertlist,
                                inclusive=True)))
    imgpix = np.unique(imgpix)
    # Set or Find min/max value of map
    dMin, dMax = (values[imgpix].min(), values[imgpix].max())
    print 'Min: %g Max: %g'%(dMin, dMax)

    if args.min is not None:
      dMin = args.min
      values[(skymap<dMin) & (values > hp.UNSEEN)] = dMin
    if args.max is not None:
      dMax = args.max
      values[(skymap>dMax) & (values > hp.UNSEEN)] = dMax

    textcolor, colormap = MapPalette.setupDefaultColormap(args.ncolors)


    # Use the Fermi/HESS/VERITAS purply-red-yellow color map
    if args.gamma:
        textcolor, colormap = MapPalette.setupGammaColormap(args.ncolors)

    # Use the Milagro color map
    if args.milagro:
      dMin = -5
      dMax = 15
      dMin = args.min if args.min != None else -5
      dMax = args.max if args.max != None else 15
      thresh = args.thresh if args.thresh != None else 2.
      textcolor, colormap = \
        MapPalette.setupMilagroColormap(dMin, dMax, thresh, args.ncolors)
      print "Milagro",dMin,dMax,thresh
    # Use a thresholded grayscale map with colors for extreme values
    else:
      if args.thresh != None:
        textcolor, colormap = \
            MapPalette.setupThresholdColormap(dMin, dMax, args.thresh,
                                              args.ncolors)
      elif args.athresh != None:
        textcolor, colormap = \
            MapPalette.setupAbsThresholdColormap(dMin, dMax, args.athresh,
                                                 args.ncolors)

    spltidx = 1
    ymin = -10
    ymax = 10
    xlabel = r"$l$ [$^\circ$]"
    ylabel = r"$b$ [$^\circ$]"
    for lonra in lonrngs:
      xmin, xmax = lonra
      faspect = abs(xmax - xmin)/abs(ymax-ymin)
      if args.interpolation:
        cRot = hp.Rotator(coord=coords[::-1],deg=False)
        phi   = np.linspace(np.deg2rad(xmax), np.deg2rad(xmin), args.xsize)
        theta = 0.5*np.pi  - np.linspace(np.deg2rad(ymin), np.deg2rad(ymax),
                                         args.xsize/faspect)
        Phi, Theta = np.meshgrid(phi, theta)
        rTheta,rPhi = cRot(Theta.reshape(phi.size*theta.size),\
                          Phi.reshape(phi.size*theta.size))
        rotimg = hp.get_interp_val(values,rTheta.reshape(Phi.shape),\
                                  rPhi.reshape(Theta.shape))
      else:
        tfig   = plt.figure(num=2,figsize=figsize)
        rotimg = hp.cartview(values, fig=2,coord=coords,title="",\
                             cmap=colormap, cbar=False,\
                             lonra=[xmin,xmax],latra=[ymin,ymax],rot=rotation,
                             notext=True, xsize=args.xsize,
                             min=dMin, max=dMax,return_projected_map=True)
        plt.close(tfig)
      # Enforce limits image limits
      rotimg[(rotimg > dMax) & (rotimg != 1e20)] = dMax
      rotimg[rotimg < dMin] = dMin
      ax = fig.add_subplot(3,1,spltidx)
      if spltidx==1:
        # X axis label
        ax.set_xlabel(xlabel, color='k')
      if spltidx==2:
        # Y axis label
        ax.set_ylabel(ylabel, color='k')
      if spltidx==3:
        # Title
        if args.title != None:
          ax.set_title(args.title,color='k')

      spltidx += 1
      imgp = ax.imshow(rotimg,extent=[xmax, xmin, ymax, ymin],\
                       vmin=dMin,vmax=dMax,cmap=colormap)
      xlms = ax.get_xlim()
      ylms = ax.get_ylim()
      imgp.get_cmap().set_under('w',alpha=0.)
      imgp.get_cmap().set_over('w',alpha=0.)
      xts = np.arange(np.floor(xmin), np.ceil(xmax+args.dmer), args.dmer)[1:-1]
      yts = np.arange(np.floor(ymin), np.ceil(ymax+args.dpar), args.dpar)[1:-1]
      ax.xaxis.set_ticks(xts)
      ax.yaxis.set_ticks(yts)
      if args.nogrid == False:
        ax.grid(color=textcolor)

      # If TeVCat data are available, plot them
      if args.tevcat or args.tevcatLabels:
        if haveTeVCat:
          try:
            if args.tevcat:
              tevcat = TeVCat.TeVCat(args.tevcat)
            elif args.tevcatLabels:
              tevcat = TeVCat.TeVCat(args.tevcatLabels)
          except IOError,e:
            print e
            print "Downloading data from tevcat.uchicago.edu"
            tevcat = TeVCat.TeVCat()
          except:
            print "Why caught here?"
            print "Downloading data from tevcat.uchicago.edu"
            tevcat = TeVCat.TeVCat()
          cRot = hp.Rotator(coord=["C",coords[-1]])
          for cId in (1,2):
            catalog = tevcat.GetCatalog(cId)
            ra = catalog.GetRA()
            dec = catalog.GetDec()
            assoc = catalog.GetCanonicalName()
            y, x = cRot(np.pi*0.5 - dec,ra)
            x = np.rad2deg(x)
            y = 90. - np.rad2deg(y)
            cut = (x > xmin) & (x < xmax) & (y>ymin) & (y<ymax)
            x = x[cut]
            y = y[cut]
            assoc = assoc[cut]
            ax.scatter(x,y,
                    color=textcolor,
                    facecolors="none",
                    marker="s")
            if args.tevcatLabels:
              for r, d, s in zip(x, y, assoc):
                ax.text(r,d, s,
                           color=textcolor,
                           rotation=45,
                           fontdict={"size":6})
        else:
          print "Sorry, TeVCat could not be loaded."

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
                    cuts = (fcat.GetFlux() > 0.) # Dummy cut
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
          
          cRot = hp.Rotator(coord=["C",coords[-1]])
          y, x = cRot(np.pi*0.5 - dec,ra)
          x = np.rad2deg(x)
          y = 90. - np.rad2deg(y)
          cut = (x > xmin) & (x < xmax) & (y>ymin) & (y<ymax)
          x = x[cut]
          y = y[cut]
          assoc = assoc[cut]
          catnms = fcat.GetCatalogName()[cuts]
          for i in xrange(len(assoc)):
              if assoc[i] == '':
                  assoc[i] = catnms[i]

          ax.scatter(x,y,
                     color=textcolor,
                     facecolors="none",
                     marker="s")
          if args.fermicatLabels:
            for r, d, s in zip(x, y, assoc):
              ax.text(r,d, s,
                      color=textcolor,
                      rotation=90,
                      fontdict={"size":14})
        else:
          print "Sorry, the Fermi xFGL catalog could not be loaded."

      # If a hotspot list is given, plot it
      if args.hotspots:
        fhot = open(args.hotspots, "r")
        ra = []
        dec = []
        for line in fhot:
          if line.startswith("#"):
            continue
          r, d = [float(t) for t in line.strip().split()[1:3]]
          ra.append(r)
          dec.append(d)
        ra = np.deg2rad(ra)
        dec = np.deg2rad(dec)
        cRot = hp.Rotator(coord=["C",coords[-1]])
        y, x = cRot(np.pi*0.5 - dec,ra)
        x = np.rad2deg(x)
        y = 90. - np.rad2deg(y)
        cut = (x > xmin) & (x < xmax) & (y>ymin) & (y<ymax)
        x = x[cut]
        y = y[cut]
        ax.scatter(x,y,
                   color=textcolor,
                   facecolors="none",
                   marker="o")
      ax.set_xlim(xlms)
      ax.set_ylim(ylms)

    
    # Set up the color bar
    nax = 0
    fhgt = 0.25
    for ax in fig.get_axes():
      ax.set_position((0.05,0.15+nax*(fhgt+0.033),0.9,fhgt))
      nax += 1

    cticks = args.ticks
    if args.ticks == None:
      clrmin = np.floor(dMin)
      clrmax = np.ceil(dMax)
      ctspc = np.round((clrmax - clrmin)/10.)
      if ctspc == 0.:
        ctspc = 1.
      if clrmin < 0 and clrmax > 0:
        cticks = -np.arange(0,-clrmin,ctspc)
        cticks = np.sort(np.unique(np.append(cticks,np.arange(0,clrmax,ctspc))))
      else:
        cticks = np.arange(clrmin,clrmax+ctspc,ctspc)

    cbax = plt.axes((0.15,0.07,0.7,0.015))
    cb = fig.colorbar(imgp, orientation="horizontal",
                      cax=cbax,
                      ticks=cticks)
    
    if args.label:
      skymapName = args.label
    else:
      if re.match("significance", skymapName):
          skymapName=r"significance [$\sigma$]"
    if args.gumbel:
      skymapName="log10(P-Value)"
    cb.set_label(skymapName)

    alpha = 1.0
    if textcolor != "#000000":
      alpha = 0.7
    
    # Either output the image to a file and quit or plot it in a window
    if args.output:
      fig.savefig(args.output, dpi=300)
    else:
      plt.show()

if __name__ == "__main__":
    main()
