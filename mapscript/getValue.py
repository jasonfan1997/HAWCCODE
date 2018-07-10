#!/usr/bin/env python
################################################################################
# Reads a single value from HEALPix Map.
# Requires healpy to be installed.
################################################################################

__version__ = "$Id: getValue.py 20811 2014-07-03 16:15:46Z fiorino $"

try:
    import argparse
    import numpy as np
    import healpy as hp
    import precess
    degree = np.pi / 180.
    import pyfits as pf
except ImportError,e:
    print e
    raise SystemExit


def readvaluefrommap(fitsfile,col=0,coldif=-1,lon=None,lat=None,outputtime=False,interpole=True,mjd=None,coord='C'):

    if lon is None or lat is None:
        print "Please specify both coordinates"
        raise SystemExit

    # Read in the skymap
    if coldif > -1:
        skymap,skymapdif = hp.read_map(fitsfile, (col,coldif), verbose=False, memmap=True)
        skymap -= skymapdif
    else:
        skymap = hp.read_map(fitsfile, col, verbose=False, memmap=True)

    th, ph = np.deg2rad([90-lat, lon])

    if mjd is not None or coord!='C':
        if mjd is not None:
            prot = precess.mjd2J2000ang(mjd=mjd,coord=['C',coord])
        else:
            prot = 0
            
        Rot = hp.Rotator(coord=['C',coord],rot=prot)
        th, ph = Rot.I(th,ph) 
    
    if interpole:
        value = hp.get_interp_val(skymap, th, ph)
    else:
        value = skymap[hp.ang2pix(hp.get_nside(skymap), th, ph)]
    
    if outputtime:
        # Not elegant, but healpy only reads half of the header
        # so I have to use another module pyfits
        hdun = 0
        hdulist = pf.open(fitsfile)
        start = hdulist[hdun].header['STARTMJD']
        stop = hdulist[hdun].header['STOPMJD']
        duration = hdulist[hdun].header['TOTDUR']
        return value, start, stop, duration
    else:
        return value


if __name__ == "__main__":

    p = argparse.ArgumentParser(description="Read map value at a given (RA,Dec)")
    p.add_argument("fitsfile", nargs="?", help="Input HEALPix FITS file")
    p.add_argument("-C", "--column", dest="col", default=0, type=int,
                   help="FITS column in which healpix map resides")
    p.add_argument("-D", "--coldiff", dest="coldif", default=-1, type=int,
                   help="FITS column in which healpix map to be subtracted "
                   "(if any) resides")
    
    p.add_argument("-r","--ra",dest="ra",type=float,
                   help="RA (deg)")

    p.add_argument("-d","--dec",dest="dec",type=float,
                   help="Dec (deg)")

    p.add_argument("--origin", dest="origin", type=float, nargs=2,
                        help="Additional option to input coordinate as lon lat pair")
    p.add_argument("-c", "--coords", dest="coord", default="C",
                   help="C=equatorial (default), G=galactic, E=Ecliptic")
    
    p.add_argument("-t","--time", action='store_true',
                   help="Flag to output start/stop dates and duration")
    p.add_argument("--mjd", default=None,type=float,
                   help="MJD of the map. Not needed if header has EPOCH key")
    p.add_argument("--nointerp", dest="interpole", default=True ,action='store_false',
                   help="Return value without interpolation")
    
    args = p.parse_args()

    #Check coordinates
    if args.origin is not None:
        if args.ra is not None or args.dec is not None:
            print "argument --origin: not allowed with argument -r/--ra and -d/--dec"
            raise SystemExit
        else:
            lon=args.origin[0]
            lat=args.origin[1]
    else:
        lon=args.ra
        lat=args.dec

    #Check epoch
    mjd=args.mjd

    if mjd is None:
        hdulist = pf.open(args.fitsfile)
        header = hdulist[0].header
        if 'EPOCH' in header:
            epoch = header['EPOCH']
            if epoch=='current':
                if 'STARTMJD' in header and 'STOPMJD' in header:
                    startmjd=header['STARTMJD']
                    stopmjd = header['STOPMJD']
                    if (startmjd>0 and stopmjd>0) or startmjd>stopmjd:
                        mjd = (stopmjd+startmjd)/2
            
    # Start normal processing
    if args.fitsfile == None:
        print "Please specify an FITS file to read."
        raise SystemExit

    print readvaluefrommap(args.fitsfile,
                           col=args.col,
                           coldif=args.coldif,
                           lon=lon,
                           lat=lat,
                           outputtime=args.time,
                           mjd=mjd,
                           interpole=args.interpole,
                           coord=args.coord)
