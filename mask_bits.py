from collections import OrderedDict as odict
from astropy import coordinates as coo

class mask_bits:
    def __init__(self):
        self.FOREGROUND_BITS = odict([
            (128, dict(name='Total area')),
            #(128, dict(name="Aggressive MW cut",filename='mwmasks/MC_mask.fits',mag=None,m1=None,m2=None,r1=None,r2=None,maxrad=None,minrad=None,cushion=None)),
            (64, dict(name='Famous stars',filename='famous_stars_footprint.fits',mag=None,m1=None,m2=None,r1=None,r2=None,maxrad=None,minrad=None,cushion=30./3600.)),
            (32, dict(name='Globular clusters',filename='globular_clusters_footprint.fits',mag=None,m1=None,m2=None,r1=None,r2=None,maxrad=None,minrad=None,cushion=30./3600.)),
            (16,  dict(name="Yale bright stars",filename='yale5th.fits',mag='Vmag',m1=4.45,m2=5.47,r1=0.4,r2=0.15,maxrad=0.4,minrad=0.15,cushion=30./3600.)),
            (8,  dict(name="Near the LMC",filename=None,mag=None,m1=None,m2=None,r1=None,r2=None,maxrad=None,minrad=None,cushion=None)),
            (4,   dict(name="2MASS bright star region (5<J<8)",filename='2mass_stars_jlt8.fits',mag='jmag',m1=5.0,m2=7.75,r1=0.075,r2=0.03,maxrad=0.075,minrad=0.03,cushion=10./3600.)),
            (2,   dict(name="Large nearby galaxy (HyperLEDA catalog)",filename='leda_galaxies.fits',mag='bt',m1=13.17,m2=13.74,r1=0.1,r2=0.03,maxrad=0.1,minrad=0.03,cushion=20./3600.)),
            (1,   dict(name="2MASS fainter star region (8<J<12)",filename='2mass_stars_jgt8.fits',mag='jmag',m1=8.0,m2=10.89,r1=0.03,r2=0.008,maxrad=0.03,minrad=0.008,cushion=10./3600.)),
])
        self.BAD_BITS = odict([
            (8,   dict(name='Total area')),
            (4,  dict(name="Crazy_colors",filename='y3a2_artifactsmask_v2.2.fits',tolerance=None,cut=None,dat=None,threshold=None)),
            (2,   dict(name="MOF_failure_regions",filename='y3a2_moffail_2_0.fits',tolerance=None,cut=None,dat=None,threshold=0)),#this is substituted by the temp fix in the main code
            (1,   dict(name="Coadd_PSF_failure_regions",filename='y3a2_psferr_2_0.fits',tolerance=0.2,cut=None,dat=None,threshold=70)),
            #(1,   dict(name="Astrometric_colors",filename='y3a2_badastrometry.fits',tolerance=coo.Angle('0d0m1s'),cut=0.2,dat=None,threshold=2)),
            #(1,   dict(name="Shallow_regions",filename='y3a2_oneexp.fits',tolerance=None,cut=None,dat=None,threshold=0)),
])

