#!/usr/bin/env python

import argparse

import os

import pyfits as pf

import healpy as hp
import matplotlib as mpl

if not os.getenv("DISPLAY"):
    mpl.use('AGG')
    
import matplotlib.pyplot as plt

from hawc import hawcnest
from hawc.hawcnest import HAWCUnits as U


def main():

    p = argparse.ArgumentParser(description="Plot histogram")
    
    p.add_argument("fitsfile",
                   help="Input HEALPix FITS file. ")
    
    p.add_argument("-o", "--output", default=None, type=str,
                   help="Save to disk")
    p.add_argument("-T", "--title", default="Exposure", type=str,
                   help="Plot title")
    p.add_argument("-n", "--nBins", default=24, type=int,
                   help="Binning. Must be multiple of 4, "\
                   "and not greater than 4*nside")
    p.add_argument("-r", "--raw", default=False, action="store_true",
                   help="Do not divide by bin size")
    p.add_argument("-u", "--timeUnit", default="hour", type=str,
                   help="Time unit (one of HAWCUnits)")
    p.add_argument("--edgecolor", default="black", type=str,
                   help="Color of bin edges (you can use 'none')")
    p.add_argument("--facecolor", default="blue", type=str,
                   help="Face color of bins")
    
    args = p.parse_args()

    #Same unit for both axes
    timeUnit = getattr(U, args.timeUnit)
    timeUnitStr = args.timeUnit

    #Find what's the column of the exposure map
    hdulist = pf.open(args.fitsfile)
    header = hdulist[1].header

    fitscol = 1
    exposureFound = False
    while True:

        fieldName=("TTYPE"+str(fitscol))
        
        if fieldName in header:
                        
            if header[fieldName] == 'exposure':
                exposureFound = True
                break

            fitscol += 1

        else:
            break

    if not exposureFound:
        raise NameError("Exposure map not found")

    #Load exposure maps
    exposureMap = hp.read_map(args.fitsfile, fitscol-1)

    nside = hp.get_nside(exposureMap)

    #Binnin info
    nBins = args.nBins
    
    if nBins%4 != 0 or nBins > 4*nside:
        raise NameError("Invalid nBins")

    binSize = U.siderealDay/nBins/timeUnit

    binEdges = [i*binSize for i in range(0,nBins+1)]
    binCenters = [i*binSize+binSize/2 for i in range(0,nBins)]

    #Get histogram info from appropiate ring    
    ring = nBins/4-1
    
    pixStart = 2 * ring * (1 + ring)
    pixStop = pixStart + nBins
    
    exposure = exposureMap[pixStart:pixStop]/timeUnit

    totEx = sum(exposure)
    meanEx = sum(exposure)/nBins
    minEx = min(exposure)
    maxEx = max(exposure)
    
    if not args.raw:
        exposure /= binSize

    #Print info
    print ("Total          = %f "+timeUnitStr) % totEx
    print ("Mean           = %f "+timeUnitStr) % meanEx
    print ("Min            = %f "+timeUnitStr) % minEx
    print ("Max            = %f "+timeUnitStr) % maxEx
    print ("bin_size       = %f "+timeUnitStr) % binSize
    print ("Mean/bin_size  = %f") % (meanEx/binSize)
    print ("Min/bin_size   = %f") % (minEx/binSize)
    print ("Max/bin_size   = %f") % (maxEx/binSize)
        
    #Create histogram
    plt.hist(binCenters, binEdges, weights=exposure,
             edgecolor=args.edgecolor,
             facecolor=args.facecolor)

    #Pimp it
    plt.title(args.title)
    
    plt.xlabel("Greenwich mean sidereal time (GMST) [" + timeUnitStr + "]")
    
    plt.ylabel("Exposure / bin size")
    if args.raw:
        plt.ylabel("Exposure [" + timeUnitStr + "]")
        
    plt.gca().set_xlim(0,U.siderealDay/timeUnit)
    
    #Save/show plot
    if args.output is not None:

        plt.savefig(args.output, dpi=400, transparent=True)
            
    else:
        
        plt.show()
    
if __name__ == "__main__":
    main()




