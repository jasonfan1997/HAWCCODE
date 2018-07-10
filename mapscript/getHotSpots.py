#!/usr/bin/env python
################################################################################
#  getHotSpots.py
#
#
#  Created by Brian Baughman on 11/14/13.
################################################################################
__version__ = "$Id: getHotSpots.py 33834 2016-08-18 16:24:54Z jbecerra $"

try:
    import argparse, os, datetime
    import healpy as hp
    import numpy  as np

except ImportError,e:
    print e
    raise SystemExit

# Try importing TeVCat
try:
  import TeVCat
  haveTeVCat = True
except ImportError as e:
  haveTeVCat = False
  print(e)

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
    p.add_argument("-d","--delimiter", dest="delimiter", type=str, default=' ',
                   help="Delimiter used in output file.")
    p.add_argument("--excludeRadius", dest="excludeRadius",
                   type=float, default=1.5,
                   help="Exclude a circular area around of size excludeRadius "
                   "around each peak. [degrees]")
    p.add_argument("-m", "--minSig", dest="minSig", type=float, default=5.,
                   help="Minimum sigma to print out")
    p.add_argument("--tevcat", action="store_true", dest="tevcat", default=False,
                   help="Draw TeVCat sources.")
    p.add_argument("--liff", action="store_true", dest="liff", default=False,
                   help="For significance maps generated with liff, which in addition to significance values it also includes the flux and its error")

    args = p.parse_args()

    if args.otype != 'std' and args.otype != 'rad' and args.otype != 'hours':
      print 'Invalid otype given: %s'%args.otype
      raise SystemExit
    if args.otype == 'hours' and noephem == True:
      print 'Cannot use otype hours without pyephem!'
      raise SystemExit

    tevcat = None
    if haveTeVCat and args.tevcat:
      haveTeVCat = False
      tday = datetime.date.today()
      tevcatnmb = 'tevcat_data_%04i-%02i-%02i.txt'%(tday.year,tday.month,tday.day)
      if not os.path.isfile(tevcatnmb):
        tevcat = TeVCat.TeVCat()
      else:
        tevcat = TeVCat.TeVCat(tevcatnmb)
      if os.path.isfile(tevcatnmb):
        haveTeVCat = True


    print args.liff

    try:
        if args.liff:
            sigmap, fluxmap, fluxerrmap, hdr = hp.read_map(args.input,field=None,h=True)
        else:
            sigmap, hdr = hp.read_map(args.input,h=True)
    except:
      print 'Cannot open input file: %s'%args.input
      raise SystemExit
    nside = hp.npix2nside(sigmap.size)
    rexcld = np.deg2rad(args.excludeRadius)
    srtidx = sigmap.argsort()[::-1]
    msked = np.array([])
    hotspots = []

    if args.liff:
        fhdr = ['#%s'%'Sigma'.rjust(7),\
                '%8s'%'Flux','%8s'%'Flux_err',\
                '%8s'%'Ra','%8s'%'Dec',\
                '%8s'%'gLong','%8s'%'gLat','%7s'%'Name']
    else:
        fhdr = ['#%s'%'Sigma'.rjust(7),\
                '%8s'%'Ra','%8s'%'Dec',\
                '%8s'%'gLong','%8s'%'gLat','%7s'%'Name']

    hotspots.append(fhdr)
    for i in srtidx:
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
          if args.liff:
              carr=['%8.2f'%sigmap[i],\
                        '%2.5E'%fluxmap[i],\
                        '%2.5E'%fluxerrmap[i],\
                        '%8.2f'%np.rad2deg(ra),'%8.2f'%np.rad2deg(dec),\
                        '%8.2f'%np.rad2deg(l),'%8.2f'%np.rad2deg(b),'HHS%04i'%(i)]
          else:
              carr=['%8.2f'%sigmap[i],\
                        '%8.2f'%np.rad2deg(ra),'%8.2f'%np.rad2deg(dec),\
                        '%8.2f'%np.rad2deg(l),'%8.2f'%np.rad2deg(b),'HHS%04i'%(i)]

      elif args.otype == 'rad':
          if args.liff:
              carr=['%8.2f'%sigmap[i],\
                        '%2.5E'%fluxmap[i],\
                        '%2.5E'%fluxerrmap[i],\
                        '%8.2f'%(ra),'%8.2f'%(dec),\
                        '%8.2f'%(l),'%8.2f'%(b),'HHS%04i'%(i)]
          else:
              carr=['%8.2f'%sigmap[i],\
                        '%8.2f'%(ra),'%8.2f'%(dec),\
                        '%8.2f'%(l),'%8.2f'%(b),'HHS%04i'%(i)]

      elif args.otype == 'hours' and noephem==False:
          if args.liff:
              carr=['%8.2f'%sigmap[i],\
                        '%2.5E'%fluxmap[i],\
                        '%2.5E'%fluxerrmap[i],\
                        '%8.2f'%epm.hours(ra),'%8.2f'%epm.degrees(dec),\
                        '%8.2f'%np.rad2deg(l),'%8.2f'%np.rad2deg(b),'HHS%04i'%(i)]
          else:
              carr=['%8.2f'%sigmap[i],\
                        '%8.2f'%epm.hours(ra),'%8.2f'%epm.degrees(dec),\
                        '%8.2f'%np.rad2deg(l),'%8.2f'%np.rad2deg(b),'HHS%04i'%(i)]

      hotspots.append(carr)
      print carr


if haveTeVCat and args.tevcat:
  for cId in (1,2):
    if cId != 1:
      continue
    catalog = tevcat.GetCatalog(cId)
    ra = catalog.GetRA()
    dec = catalog.GetDec()
    assoc = catalog.GetCanonicalName()
    sThr = 3.
    if cId == 1:
      print 'TeVCat sources above', sThr, 'sigma:'
      for r, d, s in zip(ra, dec, assoc):
        s = s.replace(' ','_')
        valueTeVCat = hp.get_interp_val(sigmap, np.pi/2. - d, r)
        if valueTeVCat > sThr:
          l, b = equ2gal(r, d, indegrees=False)
          if args.otype == 'std':
              if args.liff:
                  carr=['%8.2f'%sigmap[i],\
                        '%2.5E'%fluxmap[i],\
                        '%2.5E'%fluxerrmap[i],\
                        '%8.2f'%np.rad2deg(ra),'%8.2f'%np.rad2deg(dec),\
                        '%8.2f'%np.rad2deg(l),'%8.2f'%np.rad2deg(b),'HHS%04i'%(i)]
              else:
                  carr=['%8.2f'%sigmap[i],\
                        '%8.2f'%np.rad2deg(ra),'%8.2f'%np.rad2deg(dec),\
                        '%8.2f'%np.rad2deg(l),'%8.2f'%np.rad2deg(b),'HHS%04i'%(i)]

          elif args.otype == 'rad':
              if args.liff:
                  carr=['%8.2f'%sigmap[i],\
                        '%2.5E'%fluxmap[i],\
                        '%2.5E'%fluxerrmap[i],\
                        '%8.2f'%(ra),'%8.2f'%(dec),\
                        '%8.2f'%(l),'%8.2f'%(b),'HHS%04i'%(i)]
              else:
                  carr=['%8.2f'%sigmap[i],\
                        '%8.2f'%(ra),'%8.2f'%(dec),\
                        '%8.2f'%(l),'%8.2f'%(b),'HHS%04i'%(i)]

          elif args.otype == 'hours' and noephem==False:
              if args.liff:
                  carr=['%8.2f'%sigmap[i],\
                        '%2.5E'%fluxmap[i],\
                        '%2.5E'%fluxerrmap[i],\
                        '%8.2f'%epm.hours(ra),'%8.2f'%epm.degrees(dec),\
                        '%8.2f'%np.rad2deg(l),'%8.2f'%np.rad2deg(b),'HHS%04i'%(i)]
              else: 
                  carr=['%8.2f'%sigmap[i],\
                        '%8.2f'%epm.hours(ra),'%8.2f'%epm.degrees(dec),\
                        '%8.2f'%np.rad2deg(l),'%8.2f'%np.rad2deg(b),'HHS%04i'%(i)]

          hotspots.append(carr)

if args.output != "":
  np.savetxt(args.output,hotspots,delimiter=args.delimiter,fmt='%s')
