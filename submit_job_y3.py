import numpy as np
import os
import sys

filters=['g','r','i','z']

def submit_job(fname_sql,do_preview=False) :
    command="easyaccess -s dessci -l"
    command+=" "+fname_sql
    print(command)
    os.system(command)

def add_filters(name,behind=True) :
    sthere=""
    for fi in filters :
        if behind :
            sthere+=", "+name+fi
        else :
            sthere+=", "+fi+name
    return sthere

def add_magcut(name,lo,hi) :
    sthere=""
    sthere+=" and %s between %.3lf and %.3lf"  %(name,lo,hi)
    return sthere

def writensub(stout,data_out,fname_sql,submit):
    stout+=";> %s\n" %(data_out)

    f=open(fname_sql,"w")
    f.write(stout)
    f.close()

    if submit :
        submit_job(fname_sql)

def write_sample(sampling=1,fname_sql="sample_y3gold.sql",data_out="sample_y3gold.fits",submit=False):
    stout="select coadd_object_id"
    stout+=", ra"
    stout+=", dec"
    stout+=", extended_class_coadd"
    stout+=", extended_class_sof"
    stout+=", extended_class_mash_sof"
    stout+=", sof_cm_t"
    stout+=", sof_cm_t_err"
    stout+=", spread_model_i"
    stout+=", spreaderr_model_i"
    stout+=add_filters("mag_auto_")
    stout+=add_filters("magerr_auto_")
    stout+=add_filters("sof_cm_mag_")
    stout+=add_filters("sof_cm_mag_err_")
    stout+=add_filters("sof_psf_mag_")
    stout+=add_filters("sof_psf_mag_err_")
    stout+=" FROM y3_gold_2_2 sample("+str(sampling)+") WHERE"
    stout+=" flags_footprint = 1"
    stout+=" and bitand(flags_gold,124) = 0"
    writensub(stout,data_out,fname_sql,submit)

def write_skydistribution(classifier='extended_class_mash_sof',extval=3,fname_sql="skydistribution_y3gold.sql",data_out="skydistribution_y3gold.fits",submit=False):
    stout="select hpix_4096 as PIXEL, count(*) as SIGNAL from Y3_GOLD_2_2 where"
    stout+=" flags_footprint = 1"
    stout+=" and "+classifier+" = "+str(extval)
    stout+=" and sof_psf_mag_r BETWEEN 16 AND 23 group by hpix_4096"
    stout+=";> %s\n" %(data_out)

    writensub(stout,data_out,fname_sql,submit)
    
def write_boxsearch(ra_range,dec_range,fname_sql="out.sql",data_out="out.fits",submit=False) :
    stout="select coadd_object_id"
    stout+=", ra"
    stout+=", dec"
    stout+=", extended_class_coadd"
    stout+=", extended_class_sof"
    stout+=", extended_class_mash_sof"
    stout+=", sof_cm_t"
    stout+=", sof_cm_t_err"
    stout+=", spread_model_i"
    stout+=", spreaderr_model_i"
    stout+=add_filters("mag_auto_")
    stout+=add_filters("magerr_auto_")
    stout+=add_filters("sof_cm_mag_")
    stout+=add_filters("sof_cm_mag_err_")
    stout+=add_filters("sof_psf_mag_")
    stout+=add_filters("sof_psf_mag_err_")
    stout+=" FROM y3_gold_2_2 WHERE"
    stout+=" flags_footprint = 1"
    stout+=" and bitand(flags_gold,124) = 0"
    stout+=add_magcut("mag_auto_i",16,25)
    stout+=" and ra between %.3lf and %.3lf" %(ra_range[0],ra_range[1])
    stout+=" and dec between %.3lf and %.3lf" %(dec_range[0],dec_range[1])
    writensub(stout,data_out,fname_sql,submit)

submit = False
boxsearch = False
skydistribution = True
sample = False
sampling = 0.1
classifier = 'extended_class_mash_sof'
extval=3

if boxsearch == True:
    field_names = ["sxds","deep2_3","vvds"]
    ra = [(33.6,35.5),(350.1,354.0),(336.9,341.4)]
    dec = [(-5.7,-3.9),(-1.7,1.1),(-0.7,2.2)]
    for field,ra_range,dec_range in zip(field_names,ra,dec):
        write_boxsearch(ra_range,dec_range,fname_sql="field_"+field+"_y3gold.sql",data_out="field_"+field+"_y3gold.fits",submit=submit)

if skydistribution == True:
    write_skydistribution(classifier=classifier,extval=extval,fname_sql="skydistribution_"+classifier+str(extval)+"_y3gold.sql",data_out="skydistribution_"+classifier+str(extval)+"_y3gold.fits",submit=submit)
        
if sample == True:
    write_sample(sampling,fname_sql="sample_"+str(sampling)+"_y3gold.sql",data_out="sample_"+str(sampling)+"_y3gold.fits",submit=submit)

exit(1)

