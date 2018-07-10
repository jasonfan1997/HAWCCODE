#!/usr/bin/env python
################################################################################
# Plot HEALPix map in FITS column format using a Ortographic projection.  The
# healpy, pyfits and numpy packages need to be installed.
#
# Meant to be used for local maps
################################################################################

try:
    import pyfits as pf
    import argparse
    import numpy as np
    import healpy as hp
    from math import log

    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import re

    import MapPalette
except ImportError as e:
    print(e)
    raise SystemExit

degree = np.pi/180.0
hawcLat=19.030
hawcLon=97.269

def smoothTophat(skymap,radius):
    nside=hp.get_nside(skymap)
    
    for i in np.where(skymap!=hp.UNSEEN)[0]:
        pixInRadius=hp.query_disc(nside,hp.pix2vec(nside,i),radius,True)
        pixInRadius = pixInRadius[np.array([skymap[j]!=hp.UNSEEN for j in pixInRadius])]
        skymap[i] = sum([skymap[j] for j in pixInRadius])/len(pixInRadius) 

    return skymap

def maskHAWC(skymap,maxZenith,isHADec=False):
    """Set pixels outside the HAWC mask region to UNSEEN"""
    npix = skymap.size
    nside = hp.npix2nside(npix)

    pxls = np.arange(skymap.size)
    pxVec = hp.pix2vec(nside,pxls)
    
    if not isHADec:
        zenith = hp.ang2vec(0,0)
    else:
        zenith = hp.ang2vec((90-hawcLat)*degree,hawcLon*degree)

    pxTh = np.arccos(np.dot(np.transpose(pxVec),zenith))
        
    # Mask outside FOV
    cut = (pxTh > maxZenith*degree) | (skymap == 1e20)
    skymap[cut] = hp.UNSEEN

    return skymap

def setupColorbar(fig, title=None, ticks=None, coord="C"):
    """Create the color bar for a HEALPix map image."""
    for ax in fig.get_axes():
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
    p = argparse.ArgumentParser(description="Plots maps using the ortographic projection. Meant to plot local maps binned in azimuth-zenith angle")
    p.add_argument("fitsfile", nargs="?", help="Input HEALPix FITS file")
    p.add_argument("-C", "--column", dest="fitscol", type=int, default=0,
                   help="FITS file column number of HEALPix map")
    p.add_argument("--operation", nargs=3 , default="",
                   help="Apply an opperation to each pixel based on a second map. The options are: +, -, x, / (set to unseen if division by 0). Syntax: operation fits_file  colum_number")
    p.add_argument("--normalize", action="store_true", dest="normalize",
                   default=False,
                   help="Normalize (before operation)")
    p.add_argument("--smooth", dest="sradius",
                   default=0, type=float,
                   help="Smooth map taking the average over a circle of this radius. Performed after operation and mask.")
    p.add_argument("-m", "--min", dest="min", type=float, default=None,
                   help="Plot minimum value")
    p.add_argument("-M", "--max", dest="max", type=float, default=None,
                   help="Plot maximum value")
    p.add_argument("-n", "--ncolors", dest="ncolors", type=int, default=256,
                   help="Number of colors to use in the color map")
    p.add_argument("-o", "--output", dest="output", default=None,
                   help="Output image file (.png preferred)")
    p.add_argument("-S", "--size", dest="xsize", type=int, default=1500,
                   help="Map x-axis pixel number; increase for high resolution")
    p.add_argument("-t", "--ticks", dest="ticks", nargs="+", type=float,
                   help="List of ticks for color bar, space-separated")
    p.add_argument("-T", "--title", dest="title", default="",
                   help="Plot title")
    p.add_argument( "--label", dest="label", default=None,
                   help="Colobar label")
    p.add_argument( "--gridDiv",type=float , dest="griddiv", nargs=4 ,default=[5, 45, 9, 12],help="min_zenith max_zenith number_zenith_divisions   number_azimuth_divisions (default: 5 45 9 12 )  ")
    p.add_argument("--nogrid", action="store_false", dest="grid",
                   default=True,
                   help="Draw the altazimuthal grid")
    p.add_argument("--nomask", action="store_false", dest="mask",
                   default=True,
                   help="Don't mask unseen sky")
    p.add_argument("--alpha", dest="alpha",
                   default=0.7,
                   help="Alpha used for drawn lines.")
    p.add_argument("--logz", action="store_true", dest="logz",
                   help="Set color scale in log.")

    args = p.parse_args()

    # Start normal processing
    minZenith = args.griddiv[0]
    maxZenith = args.griddiv[1]
    nZenith = args.griddiv[2]
    nAzimuth = args.griddiv[3]

    if args.fitsfile == None:
        print("Please specify an FITS file to plot.")
        raise SystemExit

    znorm = None
    if (args.logz):
        znorm = 'log'

    skymap, skymapHeader = hp.read_map(args.fitsfile, args.fitscol, h=True)
    nside1 = hp.get_nside(skymap)

    try:
        hdulist = pf.open(args.fitsfile)
        header = hdulist[0].header
        if(header['COORDS']=="HA / Dec"):
            inHADec=True
        else:
            inHADec=False
    except:
        inHADec=False

    #Normalize
    if args.normalize:
        skymap = skymap/sum(skymap)
        
    #Apply operation
    if args.operation!="":
        operSkymap = hp.read_map(args.operation[1], int(args.operation[2]));
        operName =  args.operation[0];
        
        if operName == '+':
            skymap += operSkymap
        elif operName == '-':
            skymap -= operSkymap
        elif operName == 'x':
            skymap *= operSkymap
        elif operName == '/':
            skymap /= operSkymap
        else:
            print "The value of argument operation must be either +, -, x or /"
            raise SystemExit
        
    # remove naughty values
    skymap[np.isnan(skymap)] = hp.UNSEEN

    # Mask outside FOV
    if args.mask:
        skymap = maskHAWC(skymap,maxZenith,inHADec)

    #Smooth
    if args.sradius>0:
        skymap = smoothTophat(skymap,args.sradius*degree)
        
    # Set up the map minimum and maximum values
    notMask = skymap[skymap > hp.UNSEEN]
    dMin, dMax = (notMask[np.argmin(notMask)], notMask[np.argmax(notMask)])
    if args.min is not None:
        dMin = args.min
        skymap[(skymap<dMin) & (skymap > hp.UNSEEN)] = dMin
    if args.max is not None:
        dMax = args.max
        skymap[(skymap>dMax) & (skymap > hp.UNSEEN)] = dMax

    mpl.rc("font", family="serif", size=14)

    # Set up the color map for plotting
    textcolor, colmap = MapPalette.setupDefaultColormap(args.ncolors)

    # Set up the figure frame and coordinate rotation angle
    fig = plt.figure(1, figsize=(10,10), dpi=100)
        
    if(inHADec):
        print "Map in HA/Dec coordinates, will rotate"
        rotation = (hawcLon,hawcLat,180)
    else:
        rotation = (0,90,-90)
    
    hp.orthview(skymap, fig=1,
                        xsize=args.xsize,
                        coord='C',
                        min=dMin,
                        max=dMax,
                        rot=rotation,
                        title=args.title.replace("\\n","\n"),
                        cmap=colmap,
                        cbar=False,
                        notext=True,
                        norm=znorm,
                        half_sky=True,
                        flip='geo'
                        )


    # Hide not used space of image
    ax = fig.gca()
    imgp = ax.get_images()[0]
    imgp.get_cmap().set_under('w',alpha=0.)
    imgp.get_cmap().set_over('w',alpha=0.)
    imgp.get_cmap().set_bad('w',alpha=0.)

    ax.set_ylim([-1,1])


    # Draw grid
    npoints=1000;#if this numbes is too low healpy gets buggy
    
    if(inHADec):
        rotation = (hawcLon,hawcLat,-90)
    else:
        rotation = (0,90,0)
    
    if args.grid :
        #Constant Latitude
        for lat in np.linspace(90-minZenith,90-maxZenith,nZenith):
            hp.projplot([ np.linspace(0,360,npoints),np.linspace(lat,lat,npoints) ], lonlat=True, linestyle='-',color = 'black',coord='C', alpha = float(args.alpha), rot=rotation )
            hp.projtext(0,lat,"$%d^\circ$" % (90-lat)  , lonlat=True , coord='C',size='x-small',alpha = float(args.alpha), rot=rotation)
            
        #Constant Longitude
        for lon in np.linspace(0,360,nAzimuth,endpoint=False):
            hp.projplot([np.linspace(lon,lon,npoints),np.linspace(90-maxZenith,90-minZenith,npoints)], lonlat=True, linestyle = '-',color = 'black',coord='C',  alpha = float(args.alpha), rot=rotation )
            hp.projtext(lon,90-maxZenith,"$%d^\circ$" % lon  , lonlat=True , coord='C', rot=rotation)  
        
        #Center
        hp.projscatter([0,90], lonlat=True, color='Black',
                       s=5, marker='o', alpha=float(args.alpha),
                       coord="C", rot=rotation)
          
    # Set up the color bar and tick axis
    if args.label:
        paletteLabel = args.label
    else:
        toFind = 'TTYPE' + str(args.fitscol+1)
        skymapName = dict(skymapHeader)[toFind]
        paletteLabel = skymapName

    setupColorbar(fig, title=paletteLabel, coord="C")

    if args.output:
        fig.savefig(args.output, dpi=400)
        print("File %s created" % args.output)
    else:
        plt.show()

if __name__ == "__main__":
    main()

