### photoz_validation.py
### Author: I.Sevilla-Noarbe
### This code is based on E.Gaztanaga's VIPERS validation
from astropy.table import Table
import pandas as pd
import numpy as np
import matplotlib #needed in MacOS
matplotlib.use("TkAgg") #needed in MacOS
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.special import erf
from scipy.spatial import KDTree
from scipy.special import factorial
from scipy.special import gammaln
from scipy import interpolate
from scipy import stats
import scipy.stats as st
from scipy.optimize import curve_fit as cu
from astropy.convolution import convolve, Gaussian1DKernel, Box1DKernel
from scipy.optimize import curve_fit as cu
from astroML.resample import bootstrap
from astroML.stats import sigmaG


def readVIPERS(directoryName):
    vipersFilename = 'VIPERS_ALL_SPECTROPHOT_DR2.txt'
    hdu = fits.open(directoryName+vipersSpectraFilename)
    return hdu[1].data

def delta_z(z_phot, z_spec):
    return z_phot - z_spec

def delta_z_1pz(z_phot, z_spec):
    return delta_z(z_phot, z_spec) / (1 + z_spec)

def sigma_68(arr, axis=None):
    """Input: an (multi-dimensional) array
    Optional input: the axis along which to calculate the metric
    Outputs: the 68% spread of data about the median value of the array
    """
    upper, lower = np.percentile(arr, [84.075, 15.825], axis=axis)
    return (upper - lower) / 2.0

def sigma_68_1pz(arr, z_spec, axis=None):
    """Input: an (multi-dimensional) array
    Optional input: the axis along which to calculate the metric
    Outputs: the 68% spread of data about the median value of the array
    """
    upper, lower = np.percentile(arr, [84.075, 15.825], axis=axis)
    return (upper - lower) / 2.0 / (1 + np.mean(z_spec))

def sigma_95(arr, axis=None):
    upper, lower = np.percentile(arr, [97.7, 2.3], axis=axis)
    return (upper - lower) / 2.0

def outlier_rate(arr, outR=None):
    """assumes frac outliers >0.15
    """
    if outR is None:
        outR = 0.15
    return np.sum(np.abs(arr) > outR)*1.0/len(arr)

def plotFoM(validationData,zbins,pzs,metric):
    refZ = 'Z' #'ZSPEC' #'Z_1'
    plt.figure()
    for refPz in pzs:
        metric_list = []
        errmetric_list = []
        metric_vip_list = []
        errmetric_vip_list = []
        zmid_list = []
        for zmin,zmax in zbins:
            print('Bin ',zmin,zmax)
            basicsel = (validationData[refPz] > zmin) & (validationData[refPz] < zmax) & (validationData['FLAGS_GOLD'] < 1) & (validationData['FLAGS_FOOTPRINT'] > 0) & (validationData['FLAGS_FOREGROUND'] < 2) & (validationData['EXTENDED_CLASS_MASH_SOF'] > 2)
            zsel = (validationData[refPz] > zmin) & (validationData[refPz] < zmax) & (validationData['FLAGS_GOLD'] < 1) & (validationData['FLAGS_FOOTPRINT'] > 0) & (validationData['FLAGS_FOREGROUND'] < 2) & (validationData['EXTENDED_CLASS_MASH_SOF'] > 2) & (validationData['SOF_CM_MAG_CORRECTED_I'] > 17.5) & (validationData['SOF_CM_MAG_CORRECTED_I'] < 18 + 4*validationData['DNF_ZMEAN_SOF'])
            zsel_vip = (validationData[refPz] > zmin) & (validationData[refPz] < zmax) & (validationData['FLAGS_GOLD'] < 1) & (validationData['FLAGS_FOOTPRINT'] > 0) & (validationData['FLAGS_FOREGROUND'] < 2) & (validationData['EXTENDED_CLASS_MASH_SOF'] > 2) & (validationData['source'] == "VIPERS" ) & (validationData['SOF_CM_MAG_CORRECTED_I'] > 17.5) & (validationData['SOF_CM_MAG_CORRECTED_I'] < 18 + 4*validationData['DNF_ZMEAN_SOF'])#& (validationData[refZ] > 0.01) & (validationData['zflg'] > 2.9) & (validationData['zflg'] < 9) & (validationData['classFlag'] > 0) 
        #zsel = (validationData[refPz] > zmin) & (validationData[refPz] < zmax) & (validationData[refZ] > 0.01) & (validationData['ZFLG'] > 2.9) & (validationData['ZFLG'] < 9) & (validationData['CLASSFLAG'] > 0)
            selection = validationData[zsel]
            selection_vip = validationData[zsel_vip]
            print(len(validationData[basicsel]),len(selection),len(selection_vip))
            zmid = (zmin + zmax)/2.0
            zmid_list.append(zmid)
            if metric == 'bias':
                val = np.mean(delta_z(selection[refPz],selection[refZ]))
                errval = np.std(delta_z(selection[refPz],selection[refZ]))/np.sqrt(len(selection[refPz]))
                val_vip = np.mean(delta_z(selection_vip[refPz],selection_vip[refZ]))
                errval_vip = np.std(delta_z(selection_vip[refPz],selection_vip[refZ]))/np.sqrt(len(selection_vip[refPz]))
                ylab = 'Photo-z - z_spec bias'
            if metric == 's68':
                val = sigma_68(selection[refPz]-selection[refZ])
                val_vip = sigma_68(selection_vip[refPz]-selection_vip[refZ])
                n = 50
                errval = np.std(bootstrap(selection[refPz]-selection[refZ], n, sigma_68, kwargs=dict(axis=1)))
                errval_vip = np.std(bootstrap(selection_vip[refPz]-selection_vip[refZ], n, sigma_68, kwargs=dict(axis=1)))
                #print(selection[refPz]-selection[refZ],len(selection[refPz]-selection[refZ]),n,errval)
                ylab = 'Sigma_68'
            if metric == 's681pz':
                val = sigma_68_1pz(selection[refPz]-selection[refZ],selection[refZ])
                val_vip = sigma_68_1pz(selection_vip[refPz]-selection_vip[refZ],selection[refZ])
                n = 50
                errval = np.std(bootstrap(selection[refPz]-selection[refZ], n, sigma_68_1pz, kwargs=dict(axis=1, z_spec = selection[refZ])))
                errval_vip = np.std(bootstrap(selection_vip[refPz]-selection_vip[refZ], n, sigma_68_1pz, kwargs=dict(axis=1, z_spec = selection[refZ])))
                ylab = 'Sigma_68/1+z'
            #print(len(selection),zmin,zmax,'bias',bias,'+-',errbias,'s68',s68,'s68/1+z',s681pz)
        #if metric == 'bias':
            print(metric,'val',val,'errval',errval)
            errmetric_list.append(errval)
            metric_list.append(val)
            errmetric_vip_list.append(errval_vip)
            metric_vip_list.append(val_vip)
        plt.errorbar(zmid_list,metric_list,yerr=errmetric_list,label=refPz,marker='o',ls='')
        plt.errorbar(zmid_list,metric_vip_list,yerr=errmetric_vip_list,label=refPz+' VIPERS',marker='o',ls='')
        #plt.scatter(zmid_list,metric_list,label=refPz)
        plt.xlabel('Photo-z bin')
        plt.ylabel(ylab)
        plt.savefig(metric+'_pz.png')
        #else:
        #    plt.scatter(zmid_list,metric_list,label=refPz)
        #    plt.xlabel('Photo-z bin')
        #    plt.ylabel('Sigma_68')
    plt.legend()
    plt.show()

def plotNz_histos(validationData,zbins,pzs):

    refZ = 'Z' #'ZSPEC' #'Z_1'
    listpzs = list(pzs)
    colors = ['red','blue']
    fig = plt.figure()
    # Turn off axis lines and ticks of the big subplot
    axmain = fig.add_subplot(111)
    axmain.spines['top'].set_color('none')
    axmain.spines['bottom'].set_color('none')
    axmain.spines['left'].set_color('none')
    axmain.spines['right'].set_color('none')
    axmain.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)
    # Create subplots
    ax = []
    for f in range(4):
        ax.append(fig.add_subplot(221+f))
    # Fill subplots with N(z)
    for b,(zmin,zmax) in enumerate(zbins):
        for i,(refPz,pzNz) in enumerate(listpzs):
            zsel = (validationData[refPz] > zmin) & (validationData[refPz] < zmax) & (validationData['FLAGS_GOLD'] < 1) & (validationData['FLAGS_FOOTPRINT'] > 0) & (validationData['FLAGS_FOREGROUND'] < 2) & (validationData['EXTENDED_CLASS_MASH_SOF'] > 2) & (validationData['SOF_CM_MAG_CORRECTED_I'] > 17.5) & (validationData['SOF_CM_MAG_CORRECTED_I'] < 18 + 4*validationData['DNF_ZMEAN_SOF'])
            selection = validationData[zsel]
            #print(len[selection])
            ax[b].hist(selection[pzNz],bins=50,histtype='step',range=(0,1.5),color=colors[i],label=refPz[0:3])
            ax[b].hist(selection[refZ],bins=50,histtype='step',range=(0,1.5),color=colors[i],label='Spec. with \n '+refPz[0:3]+' binning',ls='dashed')
    axmain.set_xlabel('Photometric redshift')
    fig.subplots_adjust(wspace=0.3,left = 0.1, right = 0.9, bottom = 0.3, top = 0.95)
    #ax[2].legend(bbox_to_anchor=(1.04,1), loc="upper center")
    ax[2].legend(loc='upper left', bbox_to_anchor=(0.5, -0.28), ncol=2)
    fig.savefig('nz_comp.png')

def main():
    
    plotbias = True
    plots68 = False
    plots681pz = True
    plotNz = True
    
    directoryName = ''
    validationFilename = 'data/validsample_may2018_2_2_v2.fits'
    hdu = fits.open(directoryName+validationFilename,memmap=True)
    validationData = hdu[1].data
    #zbins = [(0.6,0.65),(0.65,0.7),(0.7,0.75),(0.75,0.8),(0.8,0.85),(0.85,0.9),(0.9,0.95),(0.95,1.0)]
    #zbins = [(0.2,0.3),(0.3,0.4),(0.4,0.5),(0.5,0.6),(0.6,0.7),(0.7,0.8),(0.8,0.9),(0.9,1.0),(1.0,1.1),(1.1,1.2),(1.2,1.3)]
    zbins = [(0.2,0.43),(0.43,0.63),(0.63,0.9),(0.9,1.3)]

    pzs = ['DNF_ZMEAN_SOF'] #'Z_MEAN'
    print('Analyzing',pzs,'from',validationFilename)
    if plotbias:
        plotFoM(validationData,zbins,pzs,'bias')
    if plots68:
        plotFoM(validationData,zbins,pzs,'s68')
    if plots681pz:
        plotFoM(validationData,zbins,pzs,'s681pz')
       
    pzsnz = ['BPZ_ZMEAN_SOF','DNF_ZMEAN_SOF'] 
    print('Analyzing',pzsnz,'from',validationFilename)
    zbins = [(0.2,0.43),(0.43,0.63),(0.63,0.9),(0.9,1.3)]
    if plotNz:
        plotNz_histos(validationData,zbins,zip(pzs,pzsnz))
            
if __name__ == "__main__":
    main()
