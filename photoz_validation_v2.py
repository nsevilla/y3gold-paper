### photoz_validation.py
### Author: I.Sevilla-Noarbe
### This code is based on E.Gaztanaga's VIPERS validation

from astropy.table import Table
import pandas as pd
import numpy as np
import matplotlib 
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
from optparse import OptionParser
matplotlib.use("Agg") #needed in MacOS
matplotlib.style.use('des_dr1')

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
        print(refPz)
        metric_list = []
        errmetric_list = []
        metric_nocorr_list = []
        #errmetric_nocorr_list = []
        zmid_list = []
        for zmin,zmax in zbins:
            print('Bin ',zmin,zmax)
            basicsel = (validationData[refPz] > zmin) & (validationData[refPz] < zmax) #& (validationData['FLAGS_GOLD'] < 1) & (validationData['FLAGS_FOOTPRINT'] > 0) & (validationData['FLAGS_FOREGROUND'] < 2) & (validationData['EXTENDED_CLASS_MASH_SOF'] > 2)
            #zsel = (validationData[refPz] > zmin) & (validationData[refPz] < zmax) & (validationData['FLAGS_GOLD'] < 1) & (validationData['FLAGS_FOOTPRINT'] > 0) & (validationData['FLAGS_FOREGROUND'] < 2) & (validationData['EXTENDED_CLASS_MASH_SOF'] > 2) & (validationData['SOF_CM_MAG_CORRECTED_I'] > 17.5) & (validationData['SOF_CM_MAG_CORRECTED_I'] < 18 + 4*validationData['DNF_ZMEAN_SOF_v2_2'])
            #zsel_vip = (validationData[refPz] > zmin) & (validationData[refPz] < zmax) & (validationData['FLAGS_GOLD'] < 1) & (validationData['FLAGS_FOOTPRINT'] > 0) & (validationData['FLAGS_FOREGROUND'] < 2) & (validationData['EXTENDED_CLASS_MASH_SOF'] > 2) & (validationData['source'] == "VIPERS" ) & (validationData['SOF_CM_MAG_CORRECTED_I'] > 17.5) & (validationData['SOF_CM_MAG_CORRECTED_I'] < 18 + 4*validationData['DNF_ZMEAN_SOF_v2_2'])#& (validationData[refZ] > 0.01) & (validationData['zflg'] > 2.9) & (validationData['zflg'] < 9) & (validationData['classFlag'] > 0) 
            zsel = (validationData[refPz] > zmin) & (validationData[refPz] < zmax) #& (validationData['FLAGS_GOLD'] < 1) & (validationData['FLAGS_FOOTPRINT'] > 0) & (validationData['FLAGS_FOREGROUND'] < 2) & (validationData['EXTENDED_CLASS_MASH_SOF'] > 2) & (validationData['SOF_CM_MAG_CORRECTED_I'] > 17.5) & (validationData['SOF_CM_MAG_CORRECTED_I'] < 18 + 4*validationData[refPz])
        #zsel = (validationData[refPz] > zmin) & (validationData[refPz] < zmax) & (validationData[refZ] > 0.01) & (validationData['ZFLG'] > 2.9) & (validationData['ZFLG'] < 9) & (validationData['CLASSFLAG'] > 0)
            selection = validationData[zsel]
            #selection_vip = validationData[zsel_vip]
            #selection_nocorr = validationData[zsel_nocorr]
            #print(len(validationData[basicsel]),len(selection),len(selection_nocorr))
            print(len(validationData[basicsel]),len(selection))
            zmid = (zmin + zmax)/2.0
            zmid_list.append(zmid)
            if metric == 'bias':
                val = np.mean(delta_z(selection[refPz],selection[refZ]))
                errval = np.std(delta_z(selection[refPz],selection[refZ]))/np.sqrt(len(selection[refPz]))
                #val_nocorr = np.mean(delta_z(selection_nocorr[refPz],selection_nocorr[refZ]))
                #errval_nocorr = np.std(delta_z(selection_nocorr[refPz],selection_nocorr[refZ]))/np.sqrt(len(selection_nocorr[refPz]))
                ylab = '$z_{photo}-z_{spec}$'
            if metric == 's68':
                val = sigma_68(selection[refPz]-selection[refZ])
                #val_nocorr = sigma_68(selection_nocorr[refPz]-selection_nocorr[refZ])
                n = 50
                errval = np.std(bootstrap(selection[refPz]-selection[refZ], n, sigma_68, kwargs=dict(axis=1)))
                #errval_nocorr = np.std(bootstrap(selection_nocorr[refPz]-selection_nocorr[refZ], n, sigma_68, kwargs=dict(axis=1)))
                #print(selection[refPz]-selection[refZ],len(selection[refPz]-selection[refZ]),n,errval)
                ylab = 'Sigma_68'
            if metric == 's681pz':
                val = sigma_68_1pz(selection[refPz]-selection[refZ],selection[refZ])
                #val_nocorr = sigma_68_1pz(selection_nocorr[refPz]-selection_nocorr[refZ],selection[refZ])
                n = 50
                errval = np.std(bootstrap(selection[refPz]-selection[refZ], n, sigma_68_1pz, kwargs=dict(axis=1, z_spec = selection[refZ])))
                #errval_nocorr = np.std(bootstrap(selection_nocorr[refPz]-selection_nocorr[refZ], n, sigma_68_1pz, kwargs=dict(axis=1, z_spec = selection[refZ])))
                ylab = '$\sigma_{68}/(1+z_{spec})$'
            #print(len(selection),zmin,zmax,'bias',bias,'+-',errbias,'s68',s68,'s68/1+z',s681pz)
        #if metric == 'bias':
            print(metric,'val',val,'errval',errval)
            errmetric_list.append(errval)
            metric_list.append(val)
            #errmetric_nocorr_list.append(errval_nocorr)
            #metric_nocorr_list.append(val_nocorr)
        if refPz=='Z_MEAN':
	        zmid_list = list(np.array(zmid_list)+0.02)
        plt.errorbar(zmid_list,metric_list,yerr=errmetric_list,label=label_dict[refPz],marker='o',ls='')
        #plt.errorbar(zmid_list,metric_nocorr_list,yerr=errmetric_nocorr_list,label=refPz+' no chorm. corr.',marker='o',ls='')
        #plt.scatter(zmid_list,metric_list,label=refPz)
    plt.xlabel('Photo-z bin')
    plt.ylabel(ylab)
        #else:
        #    plt.scatter(zmid_list,metric_list,label=refPz)
        #    plt.xlabel('Photo-z bin')
        #    plt.ylabel('Sigma_68')
    
    plt.grid()
    plt.legend(loc="upper left")
    plt.savefig(metric+'_pz_test.png')
    plt.show()

def plotNz_histos(validationData,refZ,zbins,pzs):

    #refZ = 'Z' #'ZSPEC' #'Z_1'
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
    xtextpos = [0.5,0.5,0.55,0.1]
    for f in range(4):
        ax.append(fig.add_subplot(221+f))
    # Fill subplots with N(z)
    for b,(zmin,zmax) in enumerate(zbins):
        for i,(refPz,pzNz) in enumerate(listpzs):
            print(refPz)
            print(validationData[refPz])
            #zsel = (validationData[refPz] > zmin) & (validationData[refPz] < zmax) & (validationData['FLAGS_GOLD'] < 1) & (validationData['FLAGS_FOOTPRINT'] > 0) & (validationData['FLAGS_FOREGROUND'] < 2) & (validationData['EXTENDED_CLASS_MASH_SOF'] > 2) & (validationData['SOF_CM_MAG_CORRECTED_I'] > 17.5) & (validationData['SOF_CM_MAG_CORRECTED_I'] < 18 + 4*validationData['DNF_ZMEAN_SOF_v2_2'])
            zsel = (validationData[refPz] > zmin) & (validationData[refPz] < zmax) #& (validationData['FLAGS_GOLD'] < 1) & (validationData['FLAGS_FOOTPRINT'] > 0) & (validationData['FLAGS_FOREGROUND'] < 2) & (validationData['EXTENDED_CLASS_MASH_SOF'] > 2) & (validationData['SOF_CM_MAG_CORRECTED_I'] > 17.5) & (validationData['SOF_CM_MAG_CORRECTED_I'] < 18 + 4*validationData[refPz])
            selection = validationData[zsel]
            #print(len[selection])
            #print(refPz[0:3])
            if 'MEAN' in pzNz:
                labelpz = 'MEAN'
            elif 'MC' in pzNz:
                labelpz = '1NN'
            else:
                labelpz = 'ZPHOT' 
            ax[b].hist(selection[pzNz],bins=50,histtype='step',range=(0,1.5),color=colors[i],label=labelpz) 
            ax[b].hist(selection[refZ],bins=50,histtype='step',range=(0,1.5),color=colors[i],label=refZ,ls='dashed')
            ks,pval = st.ks_2samp(selection[pzNz],selection[refZ])
            ax[b].text(xtextpos[b],0.7-(i+1)*0.1,str(pval),transform=ax[b].transAxes, fontsize=14, color=colors[i])
        ax[b].text(xtextpos[b],0.7,str(zmin)+' $<$ z $<$ '+str(zmax),transform=ax[b].transAxes, fontsize=14)
    axmain.set_xlabel('Photometric redshift')
    fig.subplots_adjust(wspace=0.3,left = 0.1, right = 0.9, bottom = 0.3, top = 0.95)
    ax[2].legend(loc='upper left', bbox_to_anchor=(0.25, -0.28), ncol=2)
    fig.savefig('nz_comp_test.png')

def plotNz_single_histos(validationData,refZ,zbins,hrange,pzsnz):
    #refZ = 'Z'#'ZSPEC'
    colors = ['red','blue']
    fig = plt.figure()
    zmin = zbins[0]
    zmax = zbins[1]
    hdu = fits.open('/Users/nsevilla/des/data/y3_gold_2_2_sample01.fits',memmap=True)
    sampleData = hdu[1].data
    vselmask = (validationData[refZ] > zmin) & (validationData[refZ] < zmax) & (validationData['SOF_CM_MAG_CORRECTED_I'] > 17.5) & (validationData['SOF_CM_MAG_CORRECTED_I'] < 18 + 4*validationData[pzsnz[0]]) #& (validationData['SOF_CM_MAG_CORRECTED_I'] < 23.5) & (validationData['FLAGS_GOLD'] < 1) & (validationData['FLAGS_FOOTPRINT'] > 0) & (validationData['FLAGS_FOREGROUND'] < 2) & (validationData['EXTENDED_CLASS_MASH_SOF'] > 2) 
    vselection = validationData[vselmask]
    sselmask = (sampleData[pzsnz[0]] > zmin) & (sampleData[pzsnz[0]] < zmax) & (sampleData['FLAGS_GOLD'] < 1) & (sampleData['FLAGS_FOOTPRINT'] > 0) & (sampleData['FLAGS_FOREGROUND'] < 2) & (sampleData['EXTENDED_CLASS_MASH_SOF'] > 2) & (sampleData['SOF_CM_MAG_CORRECTED_I'] > 17.5)  & (sampleData['SOF_CM_MAG_CORRECTED_I'] < 18 + 4*sampleData[pzsnz[0]]) #& (sampleData['SOF_CM_MAG_CORRECTED_I'] < 23.5)
    sselection = sampleData[sselmask]
    for i,pzNz in enumerate(pzsnz):
        if 'MEAN' in pzNz:
            labelpz = 'MEAN'
        elif 'MC' in pzNz:
            labelpz = '1NN'
            plt.hist(vselection[pzNz],bins=25,histtype='step',range=hrange,color=colors[i],density=True,label=labelpz+' estimate, validation data')
            #plt.hist(sselection[pzNz],bins=25,histtype='step',range=(zmin,zmax),color=colors[i],density=True,label=labelpz+' estimate, all data',ls='dashed')
        else:
            labelpz = 'ZPHOT' 
    plt.hist(vselection[refZ],bins=25,histtype='step',range=hrange,color='black',density=True,label='Spectroscopic redshift, validation data')
    plt.xlabel('N(z) distribution')
    plt.ylabel('PDF')
    plt.ylim(0,13)
    plt.legend(loc='upper right')
    fig.savefig('nz_comp_single_test.pdf')

label_dict = {'DNF_ZMEAN_SOF_v2_2':'Corrected','Z_MEAN':'No SED-extinction / chrom. corrections','ZREDMAGIC':'redMaGiC','DNF_ZMEAN_SOF':'Mean Z DNF','DNF_ZMC_SOF':'1NN Z DNF'}

def main():
    '''
    Run code with options
    '''
    #myfile = '/pool/pcae75/data1/des/monroy/matched_VALIDSAMPLE_MAY2018_2_2_fluxfiducial_dnf.fits'
    usage = "%prog [options]"
    parser = OptionParser(usage=usage)
    parser.add_option("--plot_bias",action="store_true",dest="plot_bias",help="Toggle bias test",default=False)
    parser.add_option("--plot_s68",action="store_true",dest="plot_s68",help="Toggle sigma68 test",default=False)
    parser.add_option("--plot_s681pz",action="store_true",dest="plot_s681pz",help="Toggle sigma68/(1+z) test",default=False)
    parser.add_option("--plot_Nz",action="store_true",dest="plot_Nz",help="Toggle N(z) test",default=False)
    parser.add_option("--datapath",type="string",dest="datapath",help="Directory containing validation data",
                          default="/Users/nsevilla/y3gold-paper/data/validsample_may2018_2_2_v2.fits")#matched_VALIDSAMPLE_MAY2018_2_2_fluxfiducial_dnf.fits")
    #validsample_may2018_2_2_v2.fits
    #parser.add_option("--datapath",type="string",dest="datapath",help="Directory containing validation data",
                          #default=myfile)
    (options, args) = parser.parse_args()
         
    hdu = fits.open(options.datapath,memmap=True)
    validationData = hdu[1].data

    if options.plot_bias or options.plot_s68 or options.plot_s681pz:
        zbins = [(0.2,0.3),(0.3,0.4),(0.4,0.5),(0.5,0.6),(0.6,0.7),(0.7,0.8),(0.8,0.9)]
        pzs = ['DNF_ZMC_SOF','DNF_ZMEAN_SOF']#_v2_2'] 
        print('Analyzing',pzs,'from',options.datapath)
         
        if options.plot_bias:
            plotFoM(validationData,zbins,pzs,'bias')
        if options.plot_s68:
            plotFoM(validationData,zbins,pzs,'s68')
        if options.plot_s681pz:
            plotFoM(validationData,zbins,pzs,'s681pz')
       
    if options.plot_Nz:
        refZ = 'Z'
        pzs = ['Z','Z']#_v2_2'] #for binning purposes
        pzsnz = ['DNF_ZMEAN_SOF','DNF_ZMC_SOF']#_v2_2']
        print('Analyzing',pzsnz,'from',options.datapath,'for N(z)')
        zbins = [(0.2,0.43),(0.43,0.63),(0.63,0.9),(0.9,1.3)]
        plotNz_histos(validationData,refZ,zbins,zip(pzs,pzsnz))
        #pzsnz = 'DNF_ZMC_SOF'       
        zbins = [0.55,0.70]#[0.,1.5]
        hrange = [0.2,0.8]
        plotNz_single_histos(validationData,refZ,zbins,hrange,pzsnz)
            
if __name__ == "__main__":
    main()

    
