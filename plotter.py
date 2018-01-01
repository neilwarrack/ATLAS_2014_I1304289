#!/usr/local/opt/python/libexec/bin/python
import yoda, random
hs=yoda.read("Rivet.yoda", asdict=False)
#h = yoda.Histo1D(20, 0.0, 1.0, "/foo")
h=hs[0]
for b in h.bins:
    print b
ys=[b.numEntries for b in h.bins]
xs=[b.xMid for b in h.bins]
print hs
#for _ in range(1000):
#    h.fill(random.uniform(-0.5, 1.5))
#print h

import pylab as pl
pl.xlabel('p_T of hadronically decaing top quark')
pl.ylabel('No. of Entries')
pl.title(r'Inspire:1304289 ; pp at 7 TeV')
#pl.axis([40, 160, 0, 0.03])
pl.grid(True)

pl.plot(xs,ys)
pl.savefig("out.pdf")
