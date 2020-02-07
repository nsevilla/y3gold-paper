#!/usr/bin/env python
'''
This script will plot the stellar occultation (or obscuration) in the survey
Author: Nacho Sevilla, based on work from Jelena Aleksic
Usage: python testing_occultation.py
'''
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
filename_galaxies = '/Users/nsevilla/y3gold-paper/data/field_300_redmagic.fits'  
hdulist = fits.open(filename_galaxies,memmap=True)
tdata_galaxies = hdulist[1].data
filename_stars = '/Users/nsevilla/y3gold-paper/data/field_300.fits'  #field_300.fits #field_vvds_y3gold.fits
hdulist = fits.open(filename_stars,memmap=True)
tdata_stars = hdulist[1].data

# Additional cut to be applied to the stars
cutname = 'sof_cm_mag_i'
#cut     = [[16, 18], [18, 19], [19, 20], [20, 21], [21, 22]]
cut     = [[16, 20], [20, 22]]
ncut = len(cut)
cutnamez = 'ZREDMAGIC'
cutz     = [[0.2, 0.4], [0.4, 0.6], [0.6, 0.8], [0.8,1.0]]
ncutz = len(cutz)

# STAR
cutS = (tdata_stars['extended_class_mash_sof'] < 1)
xS = tdata_stars['ra'][cutS]*3600. #degrees to arcseconds
yS = tdata_stars['dec'][cutS]*3600. #degrees to arcseconds
#star = np.vstack((xS,yS)).T
var = tdata_stars[cutname][cutS]
nsta = len(tdata_stars[cutS])

# GALAXY
#cutG = (tdata_galaxies['extended_class_mash_sof'] > 2)
cutG = (tdata_galaxies['ZREDMAGIC'] < 0.4)
xG = tdata_galaxies['ra'][cutG]*3600. #degrees to arcseconds
yG = tdata_galaxies['dec'][cutG]*3600. #degrees to arcseconds
#gala = np.vstack((xG,yG)).T
ngal = len(tdata_galaxies[cutG])

print(' * Applying additional cut to stars in', cutname)
xStar, yStar, nstar = [], [], []
for i in range(ncut):
    xStar.append(xS[(cut[i][0] <= var) & (var < cut[i][1])])
    yStar.append(yS[(cut[i][0] <= var) & (var < cut[i][1])])
    nstar.append(len(xStar[i]))
    print('     ** Cut', i+1, ': (', cut[i][0], '<=', cutname, '<', cut[i][1], ') --> Stars:', nstar[i])

star, gala = [], []
for i in range(ncut):
    tmpStar = []
    for j in range(nstar[i]): tmpStar.append([xStar[i][j], yStar[i][j]])
    star.append(np.asarray(tmpStar))
    del tmpStar
for i in range(ngal):
    gala.append([xG[i], yG[i]])
gala = np.asarray(gala)

print('Calculating NN')
neigh = NN(radius = distance+1, metric = 'euclidean')
neigh.fit(gala)
distS = [] #distance from the star
for i in range(ncut):
    dist, ind = neigh.radius_neighbors(star[i])
    distS.append(dist)
#distS, ind = neigh.radius_neighbors(star) #around star sample, what are the NN distances and indices in gala

print('Calculating histos')
hist,histerr = [],[]
area = np.zeros((nbins))

plt.figure()
for i in range(ncut):
    d, derr, cnt = np.zeros((nbins)), np.zeros((nbins)), np.zeros((nbins))
    #for j in range(len(star)):
    print(nstar[i])
    for j in range(nstar[i]):
        tmphist, _ = np.histogram(distS[i][j], bins)
        d += tmphist
        derr += tmphist
        cnt += 1
    for k in range(nbins):
        area[k] = ring(k+1,_[k],_[k+1]) # i=0 will give you funny results as area() is defined
    d = d/cnt/area
    derr = np.sqrt(derr)/cnt/area
    hist.append(d/d[nbins-1])
    histerr.append(derr/d[nbins-1])
    print(midbins)
    print(hist)
    print(derr)
    print(len(midbins),len(hist))
    plt.errorbar(midbins, hist[i], yerr=histerr[i], fmt='o', label=cutname+'='+str(cut[i]))
plt.ylabel('Relative abundance of galaxies with respect to 20 asec')
plt.xlabel('Distance (arcsec)') 
plt.legend()
plt.savefig('figs/y3gold_stellar_occultation_test.png')

