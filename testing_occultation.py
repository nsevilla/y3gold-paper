import numpy as np
import os,sys
import sklearn
from sklearn.neighbors import NearestNeighbors as NN
import matplotlib
from matplotlib import pyplot as plt
matplotlib.use("Agg")
from astropy.io import fits

def ring(i,intr,extr):
    area = np.pi*(extr*extr-intr*intr)
    return area

distance = 20
nbins = int(distance)
bins,step = np.linspace(0.5, distance+0.5, nbins+1, retstep = True)
midbins =  np.linspace(1, distance, nbins)

#fake data to test method
#xstar = (np.random.rand(100000)*100+50)/2
#ystar = (np.random.rand(100000)*100+50)/2
#star = np.vstack((xstar,ystar)).T
#xgala = np.random.rand(100000)*100
#ygala = np.random.rand(100000)*100
#gala = np.vstack((xgala,ygala)).T

#real data
filename = '/Volumes/NO NAME/field_vvds_y3gold.fits'
hdulist = fits.open(filename,memmap=True)
tdata = hdulist[1].data

# Additional cut to be applied to the stars
cutname = 'sof_cm_mag_i'
cut     = [[18, 19], [19, 20], [20, 21], [21, 22.5]]

# STAR
cutS = (tdata['extended_class_mash_sof'] < 2)
xS = tdata['ra'][cutS]*3600. #degrees to arcseconds
yS = tdata['dec'][cutS]*3600. #degrees to arcseconds
star = np.vstack((xS,yS)).T
var = tdata['sof_cm_mag_i'][cutS]
nsta = len(tdata[cutS])

# GALAXY
cutG = (tdata['extended_class_mash_sof'] == 3)
xG = tdata['ra'][cutG]*3600. #degrees to arcseconds
yG = tdata['dec'][cutG]*3600. #degrees to arcseconds
gala = np.vstack((xG,yG)).T
ngal = len(tdata[cutG])

print('Calculating NN')
neigh = NN(radius = distance+1, metric = 'euclidean')
neigh.fit(gala)
distS = [] #distance from the star
distS, ind = neigh.radius_neighbors(star) #around star sample, what are the NN distances and indices in gala

print('Calculating histos')
d, derr, cnt = np.zeros((nbins)), np.zeros((nbins)), np.zeros((nbins))
for j in range(len(star)):
    tmphist, _ = np.histogram(distS[j], bins)
    #print tmphist,len(tmphist)
    d += tmphist
    derr += tmphist
    cnt += 1

hist = []
area = np.zeros((nbins))
for i in range(nbins):
    area[i] = ring(i+1,_[i],_[i+1]) # i=0 will give you funny results as area() is defined

d = d/cnt/area
derr = np.sqrt(derr)/cnt/area
hist.append(d/d[nbins-1])
print(midbins)
print(hist)
print(len(midbins),len(hist))
plt.figure()
plt.errorbar(midbins, hist[0],yerr=derr)
plt.savefig('test.png')

