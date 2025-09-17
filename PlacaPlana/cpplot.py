# from pylab import *
from scipy.interpolate import interp1d

nc = 10
ns = 68

cp = loadtxt('placaplana.dat').reshape((ns,nc))

xteo, cpteo = loadtxt('teoria2.csv', unpack=True)
xt = linspace(min(xteo), max(xteo), 100)
xn = linspace(0.5/nc, 1.0 - 0.5/nc, nc)
cpt = interp1d(xteo, cpteo)

figure(1)
clf()
pcolor(cp)
colorbar()
show()

figure(2)
clf()
plot(xt, cpt(xt), '--', xn, cp[ns/2,:], 'o')
show()

