#!/usr/bin/env python
################################################################################
# Plot HEALPix Map in column of FITS file on a Mercator Projection.
# Requires healpy to be installed.
################################################################################

__version__ = "$Id: plotEquatorial.py 26571 2015-08-04 15:42:51Z bbaugh $"

try:
    import argparse
    import numpy as np
    import healpy as hp

    import matplotlib as mpl
    # mpl.use('Agg') # uncomment if error "couldn't connect to display"
    import matplotlib.pyplot as plt
    import re

    degree = np.pi / 180.
except ImportError, e:
    print e
    raise SystemExit

# Try importing TeVCat
try:
    import TeVCat

    haveTeVCat = True
except ImportError, e:
    haveTeVCat = False
    print e

# Try importing the Fermi catalogs
try:
    import FGLCatalog
    import FHLCatalog

    haveFermi = True
except ImportError, e:
    haveFermi = False
    print e


def Pix2Ang(p, nside):
    (x, y) = hp.pix2ang(nside, p)
    return (y / degree, 90 - x / degree)


def gal2equ(l, b):
    """Convert galactic coordinates to equatorial (J2000.0).  This is taken
    from the galeq.f program in SLALIB.
    """
    rmtx = np.matrix([[-0.054875539726, 0.494109453312, -0.867666135858],
                      [-0.873437108010, -0.444829589425, -0.198076386122],
                      [-0.483834985808, 0.746982251810, 0.455983795705]])
    cosl = np.cos(l * degree)
    sinl = np.sin(l * degree)
    cosb = np.cos(b * degree)
    sinb = np.sin(b * degree)
    gvec = np.matrix([[cosl * cosb], [sinl * cosb], [sinb]])
    cvec = rmtx * gvec

    x, y, z = (cvec.item(0), cvec.item(1), cvec.item(2))
    r = np.sqrt(x * x + y * y)
    ra = 0
    dec = 0
    if r != 0.:
        ra = np.arctan2(y, x) / degree
        if ra < 0:
            ra += 360.
    if z != 0:
        dec = np.arctan2(z, r) / degree
    return (ra, dec)


def equ2gal(ra, dec):
    """Convert equatorial (J2000.0) coordinates to galactic.  This is taken
    from the eqgal.f program in SLALIB.

    Inputs expected in degrees.
    """
    rmtx = np.matrix([[-0.054875539726, 0.494109453312, -0.867666135858],
                      [-0.873437108010, -0.444829589425, -0.198076386122],
                      [-0.483834985808, 0.746982251810, 0.455983795705]])
    cosr = np.cos(ra * degree)
    sinr = np.sin(ra * degree)
    cosd = np.cos(dec * degree)
    sind = np.sin(dec * degree)
    evec = np.matrix([[cosr * cosd], [sinr * cosd], [sind]])
    gvec = rmtx.transpose() * evec

    x, y, z = (gvec.item(0), gvec.item(1), gvec.item(2))
    r = np.sqrt(x * x + y * y)
    l = 0.
    b = 0.
    if r != 0.:
        l = np.arctan2(y, x) / degree
        if l < 0:
            l += 360.
    if z != 0:
        b = np.arctan2(z, r) / degree
    return (l, b)


def gal2equArray(l, b):
    row, col = np.shape(l)
    ra = np.zeros(row * col).reshape(row, col)
    dec = np.zeros(row * col).reshape(row, col)
    for i in range(row):
        for j in range(col):
            ra[i, j], dec[i, j] = gal2equ(l[i, j] / degree, b[i, j] / degree)
    return (ra * degree, dec * degree)


def main():
    p = argparse.ArgumentParser(description="Plot map with Mercator projection")
    p.add_argument("fitsfile", nargs="?", help="Input HEALPix FITS file")
    p.add_argument("-a", "--abthresh", dest="athresh", default=None, type=float,
                   help="Apply an absolute (lower/upper) threshold to the map")
    p.add_argument("-c", "--coords", dest="coords", default="C",
                   help="C=equatorial, G=galactic")
    p.add_argument("-D", "--coldiff", dest="coldif", default=-1, type=int,
                   help="FITS column in which healpix map to be subtracted "
                        "(if any) resides")
    p.add_argument("-l", "--threshold", dest="thresh", type=float, default=None,
                   help="Apply a lower threshold to the plot.")
    p.add_argument("--norm-sqdeg", dest="sqdeg", action="store_true",
                   default=False,
                   help="Normalize values to square degree, according to nSides "
                        "(makes sense only certain dinds of maps)")
    p.add_argument("--interpolate", action="store_true", dest="interp",
                   default=False,
                   help="Use bilinear interpolation to smear HEALPix pixels")
    p.add_argument("-s", "--scale", dest="scale", type=float, default=1.,
                   help="scale up or down values in map")

    # Mutually exclusive option to specify plot xy range OR center + WxH
    argRange = p.add_mutually_exclusive_group(required=False)
    argRange.add_argument("--xyrange", dest="xyrange", type=float, nargs=4,
                          help="Xmin Xmax Ymin Ymax plot range [degree]")
    argRange.add_argument("--origin", dest="origin", type=float, nargs=4,
                          help="Central (x,y) coords + width + height [degree]")

    args = p.parse_args()

    # Start normal processing
    if args.fitsfile == None:
        print "Please specify an FITS file to plot."
        raise SystemExit

    # Fill 2D array
    xmin = 0.
    xmax = 360.
    ymin = 64
    ymax = -26

    if (args.xyrange):
        xmin, xmax, ymin, ymax = args.xyrange
    elif (args.origin):
        xC, yC, w, h = args.origin
        xmin = xC - 0.5 * w
        xmax = xC + 0.5 * w
        ymin = yC - 0.5 * h
        ymax = yC + 0.5 * h
    else:
        print("Using default zoom window. "
              "To customize the zoom window you must specify either "
              "(xyrange) or (origin and width and height).")

    # Set up a 2D grid of phi,theta points to evaluate the sky map
    if args.interp:
        N = 300
    else:
        N = 1000
    phi = np.linspace(xmax * degree, xmin * degree, 2 * N)
    theta = np.linspace((ymax + 90) * degree, (ymin + 90) * degree, N)
    Phi, Theta = np.meshgrid(phi, theta)

    # Rotate the phi,theta grid if galactic coordinates are requested
    if args.coords == "G":
        Phi, Theta = gal2equArray(Phi, Theta - np.pi / 2.)
        Theta = Theta + np.pi / 2.

    # Read in the flux map and mask out empty pixels
    fluxmap, fluxmapHeader = hp.read_map(args.fitsfile, 1, h=True)
    # remove naughty values
    #    fluxmap[np.isnan(fluxmap)] = hp.UNSEEN
    fluxmap *= args.scale
    nside1 = hp.get_nside(fluxmap)
    npix1 = hp.nside2npix(nside1)

    # I did not find header handler, thats all I came up with...
    toFind = 'TTYPE' + str(2)
    fluxmapName = dict(fluxmapHeader)[toFind]

    frot = 0.

    # Set up the figure frame and coordinate rotation angle
    coords = ["C", "C"]
    gratcoord = "C"
    if args.coords == "G":
        coords = ["C", "G"]
        gratcoord = "G"

    # Get extrema
    angverts = [[xmin, ymin], [xmax, ymin], \
                [xmax, ymax], [xmin, ymax]]
    vertlist = []
    cRot = hp.Rotator(coord=coords[::-1], rot=np.deg2rad(frot))
    for x, y in angverts:
        ctht, cph = np.deg2rad((y, x))
        ctht = 0.5 * np.pi - ctht
        if cph < 0:
            cph = 2. * np.pi + cph
        vertlist.append(cRot(hp.ang2vec(ctht, cph)))

    # get pixels inside origin
    ROIpix = hp.query_polygon(nside1, vertlist, inclusive=True)

    # generate mask with only pixels inside ROI
    masked = np.ma.array(fluxmap);
    masked.mask = True
    for index in ROIpix:
        masked[index] = fluxmap[index]

    # Read in map of flux errors
    fluxermap, fluxermapHeader = hp.read_map(args.fitsfile, 2, h=True)
    # remove naughty values
    #    fluxermap[np.isnan(fluxermap)] = hp.UNSEEN
    fluxermap *= args.scale
    nside2 = hp.get_nside(fluxermap)
    npix2 = hp.nside2npix(nside2)

    # I did not find header handler, thats all I came up with...
    toFind = 'TTYPE' + str(3)
    fluxermapName = dict(fluxermapHeader)[toFind]

    if args.interp:
        values = hp.get_interp_val(fluxmap, 180 * degree - Theta, Phi)
    else:
        pixels = hp.ang2pix(nside1, 180 * degree - Theta, Phi)
        values = fluxmap[pixels]

    print "Flux units: TeV^-1 cm^-2 s^-1"
    print "Coord units: deg"

    # Minimum flux with uncertainty
    p = np.argmin(masked)  # pixel containing min value
    (x, y) = Pix2Ang(p, nside1)  # location of min
    print "Min:                 %11.2e +/- %11.2e (%6.2f, %5.2f)" % (masked[p], fluxermap[p], x, y);

    # Maximum flux with uncertainty
    p = np.argmax(masked)  # pixel containing max value
    (x, y) = Pix2Ang(p, nside1)  # location of max
    print "Max:                 %11.2e +/- %11.2e (%6.2f, %5.2f)" % (masked[p], fluxermap[p], x, y);
    #    print masked[ROIpix].min();
    #    print fluxmap[ROIpix].max();

    if args.origin:
        # Indicate the center of the plot with a thick black cross
        idx = hp.ang2pix(nside1, (90. - yC) * degree, xC * degree)
        print "Map value at origin: %11.2e +/- %11.2e" % (fluxmap[idx], fluxermap[idx])


if __name__ == "__main__":
    main()
