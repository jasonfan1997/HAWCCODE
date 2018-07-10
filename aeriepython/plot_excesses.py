#!/usr/bin/env python

try:
    from numpy import *
    from pylab import *
    import matplotlib.pyplot as plt
    from matplotlib import gridspec
    import argparse
except ImportError, e:
    print e
    raise SystemExit

##############################################################
# Plot the results produced with executable 'liff-BinExcess
##############################################################

p = argparse.ArgumentParser(description="Plot light curve")
p.add_argument("inputfile", help="input file with excess table")
p.add_argument("-o", "--output", dest="output", type=str,
               help="output picture file name (e.g. excesses.png")
p.add_argument("-p", "--output2", dest="output2", type=str,
               help="output picture file name (e.g. excesses2.png")
p.add_argument("-t", "--transits", dest="transits", type=float, default=0,
               help="give number of full transits to plot excess/transit, "
                    "otherwise total excess is shown")
p.add_argument("--show-expected", dest="showExpected", action="store_true", default=False,
               help="Show expected counts for input normalization before fit (default: False)")
p.add_argument("--show-per-bin", dest="showPerBin", action="store_true", default=False,
               help="Show individual single-bin best-fit expectations (default: False)")
args = p.parse_args()

try:
    file = open(args.inputfile)
except:
    raise

bins = []
aperture = []
ts = []
expSig = []
expBkg = []
bestFitSig = []
bestFitBinSig = []
bestFitBinErr = []
obsSig = []
obsErr = []
obsBkg = []

for line in file:
    try:
        if (line[0] == "#"): continue
        c = line.split()
        bins.append(int(c[0]))
        aperture.append(float(c[1]))
        ts.append(float(c[2]))
        expSig.append(float(c[3]))
        expBkg.append(float(c[4]))
        bestFitSig.append(float(c[5]))
        bestFitBinSig.append(float(c[6]))
        bestFitBinErr.append(float(c[7]))
        obsSig.append(float(c[8]))
        obsErr.append(float(c[9]))
        obsBkg.append(float(c[10]))
    except:
        raise

file.close()

bins = array(bins)
aperture = array(aperture)
ts = array(ts)
expSig = array(expSig)
expBkg = array(expBkg)
bestFitSig = array(bestFitSig)
bestFitBinSig = array(bestFitBinSig)
bestFitBinErr = array(bestFitBinErr)
obsSig = array(obsSig)
obsErr = array(obsErr)
obsBkg = array(obsBkg)

# f, axarr = plt.subplots(2, sharex=True, figsize=(16,8))
f = figure(figsize=(12, 8), tight_layout=True)
gs = gridspec.GridSpec(3, 1)
axarr = [0, 1]
axarr[0] = plt.subplot(gs[:2, :])

# fig = plt.figure(figsize=(14,8))
# ax = fig.add_subplot(2,1,1,figsize=(14,4))
axarr[0].set_yscale('log')

transits = args.transits
if (transits == 0):
    plt.ylabel(r'excess', fontsize=16)
    transits = 1.
    axarr[0].axis([float(bins[0]) - 0.4, float(bins[-1]) + 0.4,
                   10 ** (floor(log10(min(obsSig.min(), expSig.min())))),
                   10 ** (ceil(log10(max(obsSig.max(), expSig.max()))))])
else:
    plt.ylabel(r'excess per transit', fontsize=16)
    ymin = max(min(obsSig.min(), expSig.min()), 0.01)
    axarr[0].axis([float(bins[0]) - 0.4, float(bins[-1]) + 0.4,
                   10 ** (floor(log10(ymin / transits))),
                   10 ** (ceil(log10(max(obsSig.max(), expSig.max()) / transits)))])

if (args.showExpected):
    axarr[0].plot(bins, expSig / transits, 'x', color="green", label="expected excess", markersize=14., mew=2.0, ls=":",
                  lw=1.5)
axarr[0].plot(bins, bestFitSig / transits, '*', color="blue", markeredgecolor="blue", label="best fit excess",
              markersize=14.)
axarr[0].errorbar(bins, obsSig / transits, xerr=zeros(len(bins)), yerr=obsErr / transits,
                  label="data", markeredgecolor="red", markerfacecolor="none",
                  marker="o", ms=7, mew=1.5,
                  color="red", ls="", lw=1.5, capsize=6, elinewidth=2)
if (args.showPerBin):
    (_, _, elines) = axarr[0].errorbar(bins, bestFitBinSig / transits, xerr=zeros(len(bins)), yerr=bestFitBinErr / transits,
                                       label="per-bin best fit excess", markeredgecolor="orange",
                                       markerfacecolor="none", marker="s", ms=7, mew=1.5,
                                       color="orange", ls="", lw=1.5, capsize=5, elinewidth=2)
    for e in elines:
        e.set_linestyle('--')

axarr[0].legend(loc='upper right', fontsize=16, numpoints=1)
xticks(bins, fontsize=15)
yticks(fontsize=15)

plt.grid()
axarr[1] = plt.subplot(gs[2, :])
if (args.showExpected):
    axarr[1].errorbar(bins, obsSig / expSig, xerr=zeros(len(bins)), yerr=obsErr / expSig,
                      label="data / expected", markeredgecolor="green", markerfacecolor="green", marker="o", ms=7., mew=1.5,
                      color="green", lw=1.5, capsize=5)
axarr[1].errorbar(bins, obsSig / bestFitSig, xerr=zeros(len(bins)), yerr=obsErr / bestFitSig,
              label="data / best-fit", markeredgecolor="blue", markerfacecolor="blue", marker="o", ms=7., mew=1.5,
              color="blue", lw=1.5, capsize=4)
if (args.showPerBin & args.showExpected):
    axarr[1].errorbar(bins, bestFitBinSig / expSig, xerr=zeros(len(bins)), yerr=bestFitBinErr / expSig,
                      label="per-bin-best-fit / expected", markeredgecolor="green", markerfacecolor="none", marker="s",
                      ms=7., mew=1.1,
                      color="green", lw=1.5, ls="--", capsize=5)
if (args.showPerBin):
    axarr[1].errorbar(bins, bestFitBinSig / bestFitSig, xerr=zeros(len(bins)), yerr=bestFitBinErr / bestFitSig,
                      label="per-bin-best-fit / best-fit", markeredgecolor="blue", markerfacecolor="none", marker="s",
                      ms=7., mew=1.1,
                      color="blue", lw=1.5, ls="--", capsize=4)

plt.ylabel('event ratio data/sim', fontsize=16)
plt.xlabel('analysis bin #', fontsize=16)

plt.grid()
axarr[1].plot([-0.5, bins[-1]+1], [1., 1.], color="grey", ls="-", lw=1)
axarr[1].legend(loc='upper left', fontsize=12, numpoints=1, markerscale=0.4, handlelength=3)

axarr[1].axis([float(bins[0]) - 0.4, float(bins[-1]) + 0.4, 0., 3.])
xticks(bins, fontsize=15)
yticks(array(arange(0., 3., 0.5)), fontsize=15)
f.subplots_adjust(hspace=0)
plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)

if args.output:
    # f.savefig(args.output, dpi=300, tight_layout=True)
    f.savefig(args.output)
else:
    plt.show()

# Then background plots
f = figure(figsize=(12, 8), tight_layout=True)
gs = gridspec.GridSpec(3, 1)
axarr = [0, 1]
axarr[0] = plt.subplot(gs[:2, :])

# fig = plt.figure(figsize=(14,8))
# ax = fig.add_subplot(2,1,1,figsize=(14,4))
axarr[0].set_yscale('log')

transits = args.transits
if (transits == 0):
    plt.ylabel(r'background', fontsize=16)
    transits = 1.
    axarr[0].axis([float(bins[0]) - 0.4, float(bins[-1]) + 0.4,
                   10 ** (floor(log10(min(obsBkg.min(), expBkg.min())))),
                   10 ** (ceil(log10(max(obsBkg.max(), expBkg.max()))))])
else:
    plt.ylabel(r'background per transit', fontsize=16)

    ymin = min(obsBkg.min(), expBkg.min())
    axarr[0].axis([float(bins[0]) - 0.4, float(bins[-1]) + 0.4,
                   10 ** (floor(log10(ymin / transits))),
                   10 ** (ceil(log10(max(obsBkg.max(), expBkg.max()) / transits)))])

axarr[0].plot(bins, expBkg / transits, 'x', color="green", label="expected background", markersize=14., mew=2.0, ls=":",
              lw=1.5)
axarr[0].errorbar(bins, obsBkg / transits, xerr=zeros(len(bins)), yerr=sqrt(obsBkg) / transits,
                  label="observed background", markeredgecolor="red", markerfacecolor="none",
                  marker="o", ms=7, mew=1.5,
                  color="red", ls="", lw=1.5, capsize=6, elinewidth=2)


axarr[0].legend(loc='upper right', fontsize=16, numpoints=1)
xticks(bins, fontsize=15)
yticks(fontsize=15)

plt.grid()
axarr[1] = plt.subplot(gs[2, :])
axarr[1].errorbar(bins, obsBkg / expBkg, xerr=zeros(len(bins)), yerr=sqrt(obsBkg) / expBkg,
                  label="data / expected", markeredgecolor="green", markerfacecolor="green", marker="o", ms=7., mew=1.5,
                  color="green", lw=1.5, capsize=5)

plt.ylabel('background ratio data/sim', fontsize=16)
plt.xlabel('analysis bin #', fontsize=16)

plt.grid()
axarr[1].plot([float(bins[0]) - 0.4, float(bins[-1]) + 0.4], [1., 1.], color="grey", ls="-", lw=1)
axarr[1].legend(loc='upper left', fontsize=12, numpoints=1, markerscale=0.4, handlelength=3)

axarr[1].axis([float(bins[0]) - 0.4, float(bins[-1]) + 0.4, 0., 3.])
xticks(bins, fontsize=15)
yticks(array(arange(0., 3., 0.5)), fontsize=15)
f.subplots_adjust(hspace=0)
plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)

if args.output2:
    # f.savefig(args.output, dpi=300, tight_layout=True)
    f.savefig(args.output2)
else:
    plt.show()
