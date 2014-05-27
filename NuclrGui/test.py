__author__ = 'const.maslov'

from skimage import data
from numpy import *
from skimage import transform as tr
from pylab import *
from scipy.misc import factorial as fact
w = 30
l = 50
t = 50
im = data.coins()[l:l + w, t:t + w]
small = tr.rescale(im, 0.5)

imshow(small)
show()

x = linspace(0.0, 1.0, small.shape[0], endpoint=False)
y = linspace(0.0, 1.0, small.shape[1],endpoint=False)
t =[]
for i in x:
    for j in y:
        t.append([i,j])

t = array(t)

X = linspace(0.0, 1.0, 2*small.shape[0],endpoint=False)
Y = linspace(0.0, 1.0, 2*small.shape[1],endpoint=False)
T =[]
for i in X:
    for j in Y:
        T.append([i,j])

T = array(T)

m = 3
d = 2

def len(x1, x2):
    return sqrt((x1[0]-x2[0])**2 + (x1[1] - x2[1])**2)

Econst = (-1)**(d/2 + m + 1)/(2**(2*m-1)*pi**(d/2)*fact(m-1)*fact(m - d/2))

def E(s, t):
    l = len(s,t)
    if l == 0:
        l = 1
    res = Econst*l**(2*m-d) * log(l)
    return res

K = np.zeros((T.shape[0], t.shape[0]))

print K.shape

for i in xrange(K.shape[0]):
    for j in xrange(K.shape[1]):
        print i,j
        K[i,j] = E(T[i],t[j])

print 'done'

def phi(X, n):
    if n == 0:
        return 1
    elif n == 1:
        return X[0]
    elif n == 2:
        return X[1]
    elif n == 3:
        return X[0]*X[1]
    elif n==4:
        return X[0]**2
    elif n==5:
        return X[1]**2

U = zeros((T.shape[0], 6))

for i, elem in enumerate(T):
    for j in xrange(6):
        U[i, j] = phi(elem, j)

print 'done'

Ksmall = np.zeros((t.shape[0], t.shape[0]))

print Ksmall.shape

for i in xrange(Ksmall.shape[0]):
    for j in xrange(Ksmall.shape[1]):
        print i,j,E(t[i],t[j])

        Ksmall[i,j] = E(t[i],t[j])

print 'done'

Usmall = zeros((t.shape[0], 6))

for i, elem in enumerate(t):
    for j in xrange(6):
        Usmall[i, j] = phi(elem, j)

print 'done'
print small
g = []
for item in t:
    print item
    g.append(small[small.shape[0]*item[0], small.shape[1]*item[1]])

g = array(g)
print g
print g.shape
print 'DONE'
l = 1e-15
print Ksmall
W = Ksmall + Ksmall.shape[0]*l*np.eye(Ksmall.shape[0],Ksmall.shape[1])
print W

Winv = linalg.inv(W)
print 'INVERTED'
G = Usmall.transpose().dot(Winv).dot(Usmall)
Ginv = linalg.inv(G)
M1 = U.dot(Ginv).dot(Usmall.transpose()).dot(Winv)
M2 = Usmall.dot(Ginv).dot(Usmall.transpose()).dot(Winv)
print 'Part 1 done'
M = M1 + K.dot(Winv).dot(eye(M2.shape[0], M2.shape[1]) - M2)

res = M.dot(g.transpose())
print res.shape
print res
fig, ax = plt.subplots(1, 6)
out1 = tr.rescale(small, 2.0, mode = 'reflect',order = 1)
out2 = tr.rescale(small, 2.0, mode = 'reflect',order = 2)
out3 = tr.rescale(small, 2.0, mode = 'reflect', order = 3)

ax[0].imshow(im, cmap=plt.cm.gray).set_interpolation('nearest')
ax[1].imshow(small, cmap=plt.cm.gray).set_interpolation('nearest')
ax[2].imshow(out1, cmap=plt.cm.gray).set_interpolation('nearest')
ax[3].imshow(out2, cmap=plt.cm.gray).set_interpolation('nearest')
ax[4].imshow(out3, cmap=plt.cm.gray).set_interpolation('nearest')
ax[5].imshow(res.reshape(w,w), cmap=plt.cm.gray).set_interpolation('nearest')
plt.show()