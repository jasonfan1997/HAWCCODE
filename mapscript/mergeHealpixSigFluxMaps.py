#!/usr/bin/env python

try:
    import numpy as np
    import healpy as hp
    import argparse
    from sys import exit
    import os.path
except ImportError, e:
    print e
    raise SystemExit

p = argparse.ArgumentParser(description="Subtract average map from data map")
p.add_argument(dest="input", type=str, nargs='*',
               help="Input map parts in fits files")
p.add_argument("-o", "--output", dest="output", type=str,
               default="combinedMapParts.fits",
               help="Output file name")
p.add_argument("--ncols", dest="ncols", type=int,
               default=3,
               help="Number of columns in the fits file "
                    "(usually 3: sigma, flux, fluxError)")

args = p.parse_args()

output = args.output
if os.path.exists(output):
    print 'Error: output file %s already present, aborting' %output
    exit(1)

input = args.input
nfiles = len(input)
ncols = args.ncols

if ncols == 3:
    colnames = ["significance", "flux", "flux_error"]
elif ncols == 4:
    colnames = ["significance", "flux", "flux_lower_bound", "flux_upper_bound"]
elif ncols == 5:
    colnames = ["significance", "flux", "flux_error", "index", "index_error"]
else:
    print "I only know how to deal with 3, 4 or 5 columns (got %d)" % ncols
    exit(1)

print "Combining %d files, %d columns [%s]..." % (nfiles, ncols, ', '.join(colnames))

print '============ Reading map 1 of %d ============' % (nfiles)

maps = []

print 'Adding map %s' % input[0]
for i in range(ncols):
    maps.append(hp.read_map(input[0], i))

count = 2
for m in input[1:]:
    print '============ Reading map %d of %d ============' % (count, nfiles)
    print 'Adding map %s' %m
    for i in range(ncols):
        mapToAdd = hp.read_map(m, i)
        seenPix = mapToAdd != hp.UNSEEN
        if all(seenPix): #Maps used 0 for unseen
            maps[i] += mapToAdd
        else:
            maps[i][seenPix] = mapToAdd[seenPix]
    count += 1

hp.write_map(output,
             np.array(maps),
             column_names=colnames)
print('Wrote file %s' %output)
