# from pylab import *

AR = 1
disc = 4
n = 2*disc
chord = 1.0 * n
span = AR * chord
area = span * chord
da = (span / n) * (chord / n)
cp = loadtxt('placaplana.dat')
cn = sum(cp,0) * da / area
cn_teo = loadtxt('valores_ref.csv', delimiter=',')

alfa = array(range(5))
alfa *= 5

figure(1)
clf()
plot( alfa, cn, 'ro')
plot( cn_teo[:,0], cn_teo[:,1], 'b--')
plot( cn_teo[:,0], cn_teo[:,3], 'ms')
plot( cn_teo[:,0], cn_teo[:,2], 'g^')
grid(linestyle='--')
xlabel('Angle of Attack (deg)')
ylabel(r'$C_{n}$')
legend(['simulacion', 'sergio', 'ceballos','gerbhart'], loc='upper left')
show()
#
#print 'blablabla'
