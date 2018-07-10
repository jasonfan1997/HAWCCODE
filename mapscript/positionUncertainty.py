#!/usr/bin/env python
################################################################################
# Finds the position uncertainty of a source:
#  - Assumes it is given the coordinates of a local maximum.
#  - Get the TS contours corresponding to the reference position TS minus the
#    ad-hoc DeltaTS.
#  - Find the closest contour WRT the reference point
#  - Returns the maximum distance between the closest contour and the reference
#    point.
################################################################################

__version__ = "$Id: positionUncertainty.py 43026 2018-05-04 23:06:09Z criviere $"

try:
    import argparse
    import numpy as np
    import healpy as hp
    import os
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    degree = np.pi / 180.
    import scipy.stats
    import pyfits
except ImportError as e:
    print(e)
    raise SystemExit


def GetPositionUncertainty(fitsfile, origin, confidence=1., padding=3.,
                           column=0, interpolation=True, xsize=1000, freeposition=False):
    
    # Fill 2D array
    xmin=-180.
    xmax=180.
    ymax=90
    ymin=-90

    xC, yC = origin
    xmin = xC - padding / np.cos(yC*degree)
    xmax = xC + padding / np.cos(yC*degree)
    ymin = yC - padding
    ymax = yC + padding
    

    # Move to range expected by healpy
    while xmin < -180:
        xmin += 360
    while xmin > 180:
        xmin -= 360
    while xmax < -180:
        xmax += 360
    while xmax > 180:
        xmax -= 360

    if xmax < xmin:
        tmax = xmax
        tmin = xmin
        xmin = tmax
        xmax = tmin

    cxmin = xmin
    cxmax = xmax
    frot =0.
    if xmax > 90. and xmin < -90.:
        frot = 180.
        cxmin = xmax - 180.
        cxmax = xmin + 180.
    
    while xC>180:
        xC-=360

    # Read in the skymap and mask out empty pixels
    skymap = hp.read_map(fitsfile, column)
    # remove naughty values
    skymap[np.isnan(skymap)] = hp.UNSEEN
    nside = hp.get_nside(skymap)
    pixelResolution = hp.pixelfunc.nside2resol(nside) / degree
    print 'Nside %d, pixel resolution: %.4f degree' %(nside, pixelResolution)
    
    # Find FOV
    pxls = np.arange(skymap.size)
    nZeroPix = pxls[(skymap != 0)]
    pxTh, pxPh = hp.pix2ang(nside,pxls)
    
    # Mask outside FOV
    values = np.ma.masked_where((pxTh > pxTh[nZeroPix].max())
                                | (pxTh < pxTh[nZeroPix].min())
                                | (skymap == hp.UNSEEN)
                                | (skymap == 1e20) , skymap)

    # Load the skymap as an image
    faspect = abs(cxmax - cxmin)/abs(ymax-ymin)
    # Set up the figure frame and coordinate rotation angle
    coords = ["C","C"]

    # Get extrema
    angverts = [[xmin,ymin],[xmax,ymin],\
                [xmax,ymax],[xmin,ymax]]
    vertlist = [ ]
    cRot = hp.Rotator(coord=coords[::-1],rot=np.deg2rad(frot))
    for x,y in angverts:
        ctht,cph = np.deg2rad((y,x))
        ctht = 0.5*np.pi  - ctht
        if cph < 0:
            cph = 2.*np.pi + cph
        vertlist.append(cRot(hp.ang2vec(ctht,cph)))

    # Get pixels in image
    imgpix = hp.query_polygon(nside, vertlist, inclusive=True)

    # Set or Find min/max value of map
    dMin, dMax = (values[imgpix].min(), values[imgpix].max())
    print("Min: %g Max: %g" % (dMin, dMax))

    if interpolation:
        cRot = hp.Rotator(coord=coords[::-1])
        phi   = np.linspace(np.deg2rad(xmax), np.deg2rad(xmin), xsize)
        if xmin < 0 and xmax > 0 and (xmax - xmin) > 180.:
            phi   = np.linspace(np.deg2rad(xmin)+2.*np.pi,np.deg2rad(xmax), xsize)
            phi[(phi>2.*np.pi)] -= 2.*np.pi
        theta = 0.5*np.pi  - np.linspace(np.deg2rad(ymin), np.deg2rad(ymax),
                                         xsize/faspect)
        Phi, Theta = np.meshgrid(phi, theta)
        rTheta,rPhi = cRot(Theta.reshape(phi.size*theta.size),
                           Phi.reshape(phi.size*theta.size))
        rotimg = hp.get_interp_val(values,
                                   rTheta.reshape(Phi.shape),
                                   rPhi.reshape(Theta.shape))
    else:
        rotimg = hp.cartview(values,
                             coord=coords,
                             lonra=[cxmin,cxmax],
                             latra=[ymin,ymax],rot=frot,
                             notext=True,
                             xsize=xsize,
                             return_projected_map=True)

    contourmapcoord = 'C'

    fvalues = np.ma.masked_where(skymap == 0, skymap)
    if interpolation:
        cRot = hp.Rotator(coord=[coords[-1],contourmapcoord])
        phi   = np.linspace(np.deg2rad(xmax), np.deg2rad(xmin), xsize)
        if xmin < 0 and xmax > 0 and (xmax - xmin) > 180.:
            phi   = np.linspace(np.deg2rad(xmin)+2.*np.pi,np.deg2rad(xmax), xsize)
            phi[(phi>2.*np.pi)] -= 2.*np.pi
        theta = 0.5*np.pi - np.linspace(np.deg2rad(ymin), np.deg2rad(ymax),
                                        xsize/faspect)
        Phi, Theta = np.meshgrid(phi, theta)
        rTheta,rPhi = cRot(Theta.reshape(phi.size*theta.size),\
                            Phi.reshape(phi.size*theta.size))
        frotimg = hp.get_interp_val(fvalues,rTheta.reshape(Phi.shape),\
                                    rPhi.reshape(Theta.shape))
    else:
        frotimg = hp.cartview(fvalues,
                              coord=[contourmapcoord, coords[-1]],
                              lonra=[xmin,xmax],
                              latra=[ymin,ymax],
                              rot=frot,
                              xsize=1000,
                              min=dMin, max=dMax,
                              return_projected_map=True)


    def sigmaToDeltaTS(sigma, df=1):
        cdf = scipy.stats.chi2.cdf(sigma*sigma,1)
        dts = scipy.stats.chi2.ppf(cdf,df)
        return dts

    deltaTS = sigmaToDeltaTS(confidence,2)
    
    if freeposition:
        # update xC, yC to the values of the maximum of the map (within imgpix pixels)
        argmax = np.argmax(values[imgpix])
        pixmax = imgpix[argmax]
        yC, xC = np.rad2deg(hp.pix2ang(nside,pixmax))
        yC = 90-yC

    th, ph = np.deg2rad([90-yC, xC])
    centralPixelValue = skymap[hp.ang2pix(nside, th, ph)]
    contourSigmaValue = np.sqrt(np.power(centralPixelValue,2)-deltaTS)
    print 'Central pixel value: %.3f' %centralPixelValue
    print 'Sigma value for %.2f sigma uncertainty contour: %.3f' \
          %(confidence,contourSigmaValue)

    contp = plt.contour(frotimg,
                        levels=[contourSigmaValue,],
                        origin='upper',
                        extent=[cxmax, cxmin, ymax, ymin])
    
    def distance(lon0,lat0,lon1,lat1):
        dLon = lon1 - lon0
        dLat = lat1 - lat0
        mLat = (lat1+lat1) / 2.
        return np.sqrt(np.power(dLon*np.cos(mLat*degree),2)+np.power(dLat,2))
    
    minmax = []
    if len(contp.collections[0].get_paths()) == 0:
        print 'No contour, hopefully the uncertainty is too small to measure'
        print 'Returning pixel resolution.'
        return pixelResolution
    for i,cont in enumerate(contp.collections[0].get_paths()):
        print 'Contour %i:' %i
        distances = [distance(xy[0],xy[1],xC,yC) for xy in cont.vertices]
        distMin = np.min(distances)
        distMax = np.max(distances)
        print 'Min %.3f max %.3f ' %(distMin,distMax)
        minmax.append([distMin,distMax])
    print 'min:', min(minmax)
    
    uncertainty = float(min(minmax)[1])
    totUncertainty = float(np.sqrt(np.power(uncertainty,2)+np.power(pixelResolution,2)))
    print 'resolution: %.4f, plus pixel size: %.4f' %(uncertainty,totUncertainty)
    if freeposition:
        return xC, yC, totUncertainty, dMax, pyfits.open(fitsfile)[0].header['STARTMJD']
    else:
        return totUncertainty

if __name__ == "__main__":

    p = argparse.ArgumentParser(description="Compute position statistical "
                                            "uncertainty")
    p.add_argument("fitsfile", nargs='?', help="Input HEALPix FITS file")
    p.add_argument("--column", default=0, type=int,
                   help="FITS column in which healpix map resides")
    p.add_argument("--confidence", type=float, default = 1.,
                   help="Confidence level contour value, expressed in sigma.")
    p.add_argument("--nointerpolation", dest="interpolation",
                   action='store_false',
                   help="Do not uses bilinear interpolation in data map.")
    p.add_argument("--xsize", type=int, default=1000,
                   help="Number of X pixels, Y will be scaled acordingly.")
    p.add_argument("--origin", type=float, nargs=2, default=None,
                   help="Central (x,y) coords [degree]")
    p.add_argument("--padding", type=float, nargs=1, default=3.,
                   help="Map to load around the origin [degree]")
    p.add_argument("--max", action='store_true',
                   help="Use the maximum")

    args = p.parse_args()

    # Start normal processing
    if args.fitsfile == None:
        print("Please provide input FITS file.")
        raise SystemExit
    if args.origin is None:
        print("Please provide origin.")
        raise SystemExit

    GetPositionUncertainty(fitsfile=args.fitsfile,
                           column=args.column,
                           confidence=args.confidence,
                           interpolation=args.interpolation,
                           xsize=args.xsize,
                           origin=args.origin,
                           padding=args.padding,
                           freeposition=args.max)
