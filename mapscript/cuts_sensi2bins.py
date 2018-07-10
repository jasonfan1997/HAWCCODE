#! /usr/bin/python

try:
    import sys
    import argparse
except ImportError, e:
    print e
    raise SystemExit

p = argparse.ArgumentParser(description="Convert sensi cuts file to liff bins-file")
p.add_argument("sensifile", help="Sensi input file")
p.add_argument("-o", "--output", dest="output", type=str, help="Output .bins file name.")
p.add_argument("-t", "--totalCh", default=-1, dest="totalCh", type=int,
               help="Total number of channels. If not given, last max-nCh will be used.")
args = p.parse_args()

binid = []
nhl = []
nhu = []
comcut = []
gw1 = []
gw2 = []
frac1 = []
corec2n = []
anglc2n = []

file = open(args.sensifile, "r")
for line in file:
    try:
        f = line.split()
        if (f[0][0] == "#"): continue
        binid.append(int(f[0]))
        nhl.append(int(f[1]))
        nhu.append(int(f[2]))
        comcut.append(float(f[8]))
        gw1.append(float(f[13]))
        gw2.append(float(f[14]))
        frac1.append(float(f[15]))
        corec2n.append(float(f[19]))
        anglc2n.append(float(f[20]))
    except:
        continue
file.close()

totalCh = args.totalCh
if (totalCh == -1):  totalCh = nhu[-1]

out = open(args.output, "w")
out.write("# ID  cuts \n")
upperFrac = -1
for b in xrange(0, len(binid)):
    out.write('%3d   ' % (binid[b]))
    if (b == 0):
        lowerFrac = float(nhl[b]) / float(totalCh)
    else:
        lowerFrac = upperFrac
    upperFrac = float(nhu[b]) / float(totalCh)
    out.write('"(1.*rec.nHit/rec.nChAvail>=%.3f)' % (lowerFrac))
    out.write(' && (1.*rec.nHit/rec.nChAvail<%.3f)' % (upperFrac))
    out.write(' && (1.*rec.angleFitChi2/rec.angleFitNdof<%.3f)' % (anglc2n[b]))
    out.write(' && (1.*rec.coreFitChi2/rec.nHit<%.3f)' % (corec2n[b]))
    out.write(' && (1.*rec.nHit/rec.CxPE40>=%.3f)"\n' % (comcut[b]))
out.close()

print "Double Gaussian PSF values:"
gw1s = ""
gw2s = ""
frac1s = ""
for b in xrange(0, len(binid)):
    gw1s += '%4.2f,' % (gw1[b])
    gw2s += '%4.2f,' % (gw2[b])
    frac1s += '%4.2f,' % (frac1[b])
print "  GW1:   ", gw1s[:-1]
print "  GW2:   ", gw2s[:-1]
print "  Frac1: ", frac1s[:-1]
