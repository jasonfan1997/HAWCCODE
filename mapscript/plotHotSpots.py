#!/usr/bin/env python
#  hotspot-analysis.py
#  pyanalysis
#
#  Created by Brian Baughman on 3/18/15.
#  Copyright (c) 2015 University of Maryland. All rights reserved.

try:
  import argparse
  import numpy as np
  import os
  import subprocess as subp
  import datetime
except ImportError,e:
  print e
  raise SystemExit

# Try importing TeVCat
try:
  import TeVCat
  haveTeVCat = True
  import healpy as hp
except ImportError,e:
  haveTeVCat = False
  print e


# program description
dcrpt = ['Given an input skymap it will find hotspots and plot zooms of each.']
dcrpt.append("Output will be written to current directory with names from input.")
dcrpt.append(" Expects HAWC_INSTALL to point to an AERIE installation location")


parser = argparse.ArgumentParser(description='\n'.join(dcrpt))

parser.add_argument("-i","--input",\
                    action="store", dest="input",\
                    type=str,\
                    help="Input signficance map file.")
parser.add_argument("-r","--excludeRadius",\
                    action="store", dest="excludeRadius",\
                    default=3, type=float,\
                    help="Exclusion radius around which to exclude peaks in "
                          "signficance map.")
parser.add_argument("-s","--minSig",\
                    action="store", dest="minSig",\
                    default=5, type=float,\
                    help="Minimum sigma to call a hotspot.")
parser.add_argument("-e","--extraArgs",\
                    action="store", dest="extraArgs",\
                    default="--interpolation --gamma --min 0", type=str,\
                    help="Extra arguments to be passed to plotter")
parser.add_argument("--windowWidth",\
                    action="store", dest="wwidth",\
                    default=10, type=float,\
                    help="Window width in degrees.")
parser.add_argument("--windowHeight",\
                    action="store", dest="wheight",\
                    default=10, type=float,\
                    help="Window height in degrees.")
parser.add_argument("-g","--galactic", dest="galactic", \
               default=False, action='store_true',
               help="Plot in galactic coordinates.")
parser.add_argument("--tevcat", action="store_true", dest="tevcat", default=False,
               help="Draw TeVCat sources.")
args = parser.parse_args()


try:
  exdir = os.environ['HAWC_INSTALL']
except:
  print 'Cannot find HAWC_INSTALL in environment!'
  raise SystemExit


if not os.path.exists(exdir):
  print 'HAWC_INSTALL does not exist at: %s'%exdir
  raise SystemExit


hsexe = '%s/bin/getHotSpots.py'%(exdir)

if not os.path.exists(hsexe):
  print 'HAWC_INSTALL does not contain expected scripts please update AERIE.'
  raise SystemExit


pltexe = '%s/bin/plotMercator.py'%(exdir)

if not os.path.exists(pltexe):
  print 'HAWC_INSTALL does not contain expected scripts please update AERIE.'
  raise SystemExit

if not os.path.exists(args.input):
  print 'INPUT does not exist at %s.'%args.input
  raise SystemExit

input = os.path.abspath(args.input)

tevcatfn = ''
if haveTeVCat and args.tevcat:
  tevcatfn = subp.Popen('ls -t ${PWD}/tevcat_data_*.txt | head -1',\
                        shell=True,stdout=subp.PIPE).communicate()

  tevcatfn = tevcatfn[0].strip()
  tvdl = True
  if os.path.exists(tevcatfn):
    tvdstr = os.path.basename(tevcatfn)
    tvdstr = tvdstr.replace('.txt','')
    tvdstr = tvdstr.replace('tevcat_data_','')
    tvdarr = tvdstr.split('-')
    if len(tvdarr) == 3:
      tvdate = datetime.date(int(tvdarr[0]),int(tvdarr[1]),int(tvdarr[2]))
      today = datetime.date.today()
      if tvdate == today:
        tvdl = False
  if tvdl == True:
    print "Fetching data from tevcat.uchicago.edu."
    tevcat = TeVCat.TeVCat()
    tevcatfn = subp.Popen('ls -t ${PWD}/tevcat_data_*.txt | head -1',\
                          shell=True,stdout=subp.PIPE,stderr=None).communicate()
    tevcatfn = tevcatfn[0].strip()
  else:
    tevcat = TeVCat.TeVCat(tevcatfn)
  if not os.path.exists(tevcatfn):
    print 'Could not determine TeVCat file name path. Will not use TeVCat'
    haveTeVCat = False
else:
  haveTeVCat = False

ofnbase = os.path.basename(input)
ofnbase = ofnbase[0:ofnbase.find('.fits')]

hsfname = '%s-hotspots-sig.gt.%g.dat'%(ofnbase,args.minSig)
# Get hotspots.
hscmd = [ hsexe,\
         '--input','%s'%input ,\
         '--output','%s'%hsfname ,\
         '--minSig','%f'%args.minSig ,\
         '--excludeRadius','%f'%args.excludeRadius ]

hsout = subp.Popen(' '.join(hscmd),shell=True,stdout=subp.PIPE).communicate()
hsf = open(hsfname,'r')

tvposs = []
tvlocs = []
tvassocs = []
if haveTeVCat:
  for cId in (1,2):
    catalog = tevcat.GetCatalog(cId)
    ra = catalog.GetRA()
    dec = catalog.GetDec()
    assoc = catalog.GetCanonicalName()
    for r, d, s in zip(ra, dec, assoc):
      tvposs.append([np.rad2deg(d),np.rad2deg(r)])
      tvlocs.append(hp.ang2vec(np.pi*0.5 - d, r))
      tvassocs.append(s.replace(' ','_'))

tvposs = np.asarray(tvposs)
tvlocs = np.asarray(tvlocs)
tvassocs = np.asarray(tvassocs)
# Loop over hotspots
hsidx = 0
for line in hsf:
  line = line.strip()
  # Skip comments
  if line[0] == '#':
    continue
  if hsidx == 0:
    print '%-8s %-8s %-8s %-8s %-8s %-s'%('Index','Sigma','AngSep',\
                                'd(dec)','d(ra)','Association')
  carr = line.split()
  ra  = carr[1]
  dec = carr[2]
  cofname = '%s-hs%i.png'%(ofnbase,hsidx)
  xcent = ra
  ycent = dec
  if args.galactic:
    xcent = carr[3]
    ycent = carr[4]
  pltcmd = [ pltexe,\
            args.extraArgs,\
            '--origin %s %s %f %f'%(xcent,ycent, args.wwidth, args.wheight),\
            '--hotspots','%s'%hsfname ]
  if args.galactic:
    pltcmd.append('-c')
    pltcmd.append('G')

  cpos = [float(dec),float(ra)]
  cassoc = 'None'
  casep = 180.
  ctitle = 'Hotspot %i Sig=%s'%(hsidx,carr[0])
  if haveTeVCat:
    pltcmd.append('--tevcat %s'%tevcatfn)
    # check if nearby TeVCat source
    cvec = hp.ang2vec(np.deg2rad(90 - float(dec)),np.deg2rad(float(ra)))
    for i in xrange(tvassocs.size):
      asep = np.rad2deg(np.arccos(np.dot(cvec,tvlocs[i])))
      if asep < args.excludeRadius and asep < casep:
        casep = asep
        cpos = tvposs[i]
        cassoc = tvassocs[i]
        cofname = '%s-%s.png'%(ofnbase,cassoc.replace(' ','_'))
        ctitle = '%s Sig=%s'%(cassoc,carr[0])


  print '%-8i %-8s %-8.2g %-8.2g %-8.2g %-s'%(hsidx, carr[0], casep,\
                                cpos[0]-float(dec),cpos[1] - float(ra),\
                                cassoc)
  pltcmd.append('--title "%s"'%ctitle)
  pltcmd.append('--output %s'%cofname)
  pltcmd.append(input)
  pltout = subp.Popen(' '.join(pltcmd),shell=True,stdout=subp.PIPE).communicate()
  hsidx += 1





