#!/usr/bin/env python
################################################################################
#  getHotSpots.py
#
#
#  Created by Brian Baughman on 11/14/13.
################################################################################
__version__ = "$Id: getLocalHotSpots.py 24903 2015-04-22 19:48:35Z bbaugh $"

try:
    import argparse
    import healpy as hp
    import numpy  as np

except ImportError,e:
    print e
    raise SystemExit
noephem = False
try:
  import ephem as epm
except:
  noephem = True
def equ2gal(rra, rdec, indegrees=True):
    """Convert equatorial (J2000.0) coordinates to galactic.  This is taken
    from the eqgal.f program in SLALIB.

    Inputs expected in degrees.
    """
    rmtx = np.matrix([[-0.054875539726,  0.494109453312, -0.867666135858],
                     [-0.873437108010, -0.444829589425, -0.198076386122],
                     [-0.483834985808,  0.746982251810,  0.455983795705]])
    if indegrees == True:
      rra = np.deg2rad(rra)
      rdec = np.deg2rad(rdec)
    cosr = np.cos(rra)
    sinr = np.sin(rra)
    cosd = np.cos(rdec)
    sind = np.sin(rdec)
    evec = np.matrix([[cosr*cosd], [sinr*cosd], [sind]])
    gvec = rmtx.transpose() * evec

    x, y, z =  (gvec.item(0), gvec.item(1), gvec.item(2))
    r = np.sqrt(x*x + y*y)
    l = 0.
    b = 0.
    if r != 0.:
        l = np.arctan2(y, x)
        if l < 0:
            l += 2.*np.pi
    if z != 0:
        b = np.arctan2(z, r)
    if indegrees == True:
      l = np.deg2rad(l)
      b = np.deg2rad(b)
    return (l, b)

if __name__ == "__main__":
    # Set up input
    p = argparse.ArgumentParser(description="Make significance map histogram")
    p.add_argument("-i","--input", dest="input", type=str,
                   help="Input HEALPix FITS file")
    p.add_argument("-o","--output", dest="output", type=str, default="",
                   help="output list of hotspots file")
    ohstr = "std - degrees, rad - radians"
    if noephem == False:
      ohstr += ", hours - hour angles"
    p.add_argument("--otype", dest="otype", type=str, default="std",
                     help="Output type: %s"%ohstr)
    p.add_argument("-d","--delimiter", dest="delimiter", type=str, default='\t',
                   help="Delimiter used in output file.")
    p.add_argument("--excludeRadius", dest="excludeRadius",
                   type=float, default=1.5,
                   help="Exclude a circular area around of size excludeRadius "
                   "around each peak. [degrees]")
    p.add_argument("-m", "--minSig", dest="minSig", type=float, default=5.,
                   help="Minimum sigma to print out")
    p.add_argument( "--ra", dest="ra", type=float, default=0,
                   help="RA of center of map in degrees.")
    p.add_argument( "--dec", dest="dec", type=float, default=0,
                   help="Declination of center of map in degrees.")
    p.add_argument( "--windowRadius", dest="windowRadius",
                   default=5, type=float,
                   help="Declination of center of map in degrees.")

    args = p.parse_args()

    if args.otype != 'std' and args.otype != 'rad' and args.otype != 'hours':
      print 'Invalid otype given: %s'%args.otype
      raise SystemExit
    if args.otype == 'hours' and noephem == True:
      print 'Cannot use otype hours without pyephem!'
      raise SystemExit

    try:
      sigmap, hdr = hp.read_map(args.input,h=True)
    except:
      print 'Cannot open input file: %s'%args.input
      raise SystemExit
    nside = hp.npix2nside(sigmap.size)
    centDec = np.pi/2. - np.deg2rad(args.dec)
    centRA = np.deg2rad(args.ra)
    winRad = np.deg2rad(args.windowRadius)
    centVec = hp.ang2vec(centDec,centRA)
    localidx = hp.query_disc(nside,centVec,winRad)
    srtlocalidx = sigmap[localidx].argsort()[::-1]
    rexcld = np.deg2rad(args.excludeRadius)
    srtidx = localidx[srtlocalidx]
    msked = np.array([])
    hotspots = []
    fhdr = ['#%s'%'Sigma'.rjust(9),\
        '%12s'%'Ra','%12s'%'Dec',\
        '%10s'%'gLong','%10s'%'gLat']
    hotspots.append(fhdr)
    print "\n%s\n%10s\t%10s\t%6s"%(''.center(26,'-'),'Ra','Dec','Sigma')
    for i in srtidx:
      if np.any(i == localidx) == False:
        continue
      if sigmap[i] < args.minSig:
        break
      if any(i == msked):
          continue
      nnghrs = hp.query_disc(nside,hp.pix2vec(nside,i),rexcld)
      msked = np.append(msked,nnghrs)
      dec, ra = hp.pix2ang(nside, i)
      ra  = ra
      dec = np.pi/2. - dec
      l, b = equ2gal(ra, dec, indegrees=False)
      if args.otype == 'std':
        carr=['%10.2f'%sigmap[i],\
              '%12s'%np.rad2deg(ra),'%12s'%np.rad2deg(dec),\
              '%10.2f'%np.rad2deg(l),'%10.2f'%np.rad2deg(b)]
      elif args.otype == 'rad':
        carr=['%10.2f'%sigmap[i],\
              '%12s'%(ra),'%12s'%(dec),\
              '%10.2f'%(l),'%10.2f'%(b)]
      elif args.otype == 'hours' and noephem==False:
        carr=['%10.2f'%sigmap[i],\
              '%12s'%epm.hours(ra),'%12s'%epm.degrees(dec),\
              '%10.2f'%np.rad2deg(l),'%10.2f'%np.rad2deg(b)]

      hotspots.append(carr)
      print "%s\t%s\t%s"%(carr[1],carr[2],carr[0])

if args.output != "":
  np.savetxt(args.output,hotspots,delimiter=args.delimiter,fmt='%s')
