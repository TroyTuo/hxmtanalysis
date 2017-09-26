#!/usr/bin/env python

from function import *
import argparse

# find and create data dir list

parser = argparse.ArgumentParser()
parser.add_argument("-i","--input",help="data archived path")
parser.add_argument("-o","--output",help="products archived path")
parser.add_argument("-f","--fig",action="store_true",help="plot space craft data")
parser.add_argument("--rangefile",help="rangefile path of hegtigen")
parser.add_argument("--clobber",action="store_true",help="excute hegtigen and hescreen and overwrite existing gti or screens files")
parser.add_argument("--hxbary",action="store_true",help="carry out hxbary and copy Evt file to local directory")
args = parser.parse_args()
data_dir = args.input
product_dir = args.output
aux_dir = product_dir + "/AUX/" # AUX path
acs_dir = product_dir + "/ACS/" # ACS path
he_dir = product_dir + "/HE/"   # HE  path
ehkfilename = aux_dir + "EHK.fits"
print data_dir,args.input
print product_dir

if not os.path.isdir(product_dir):os.system('mkdir -p '+product_dir)
if not os.path.isdir(aux_dir):os.system('mkdir -p ' +aux_dir)
if not os.path.isdir(acs_dir):os.system('mkdir -p ' +acs_dir)
if not os.path.isdir(he_dir):os.system('mkdir -p '+he_dir)
print aux_dir
print acs_dir
print he_dir
print ehkfilename

#read filenames
v2_flag = commands.getoutput('if [ -f '+data_dir+'/HE/*HE-Evt_FFFFFF_V2* ];then ehco 1;else echo 0;fi')
if v2_flag == 1:
    filename = commands.getoutput('ls ' + data_dir + '/HE/*HE-Evt_FFFFFF_V2*')
else:
    filename = commands.getoutput('ls ' + data_dir + '/HE/*HE-Evt_FFFFFF_V1*')
print filename
hvfilename = commands.getoutput('ls ' + data_dir + '/HE/HXMT*HV_FFFFFF*')
pmfilename = commands.getoutput('ls ' + data_dir + '/HE/HXMT*PM_FFFFFF*')
deadfilename = commands.getoutput('ls ' + data_dir + '/HE/HXMT*DTime*')
tempfilename = commands.getoutput('ls ' + data_dir + '/HE/HXMT*TH*')

print hvfilename
print pmfilename
print deadfilename
print tempfilename

# generate ehk file using HXMT software
if not os.path.isfile(ehkfilename):
    ehkgen(data_dir,product_dir)

# generate data products utilizing HXMT software
if args.hxbary:
    if not os.path.isfile(he_dir+'/*HE-Evt*'):
        cp_text = 'cp -n '+filename+' '+he_dir
        os.system(cp_text)
        print cp_text #?????????????????hxbary is undo

if not os.path.isfile(he_dir+'he_pi.fits') or True:#???????????????????always True right now
    hepical_text = 'hepical evtfile='+filename+' outfile='+he_dir+'he_pi.fits clobber=yes'
    os.system(hepical_text)
    print hepical_text
if args.clobber:
    hegtigen_text = 'hegtigen hvfile='+hvfilename+\
    ' tempfile='+tempfilename+' outfile='+he_dir+'he_gti.fits ehkfile='+ehkfilename+\
    ' defaultexpr=no ELV=6 COR=10 SAA=yes T_SAA=10 TN_SAA=10 SUN_ANGLE=10'\
    ' MOON_ANGLE=10 ANG_DIST=100 clobber=yes'
    hegtigen_text = 'hegtigen hvfile='+hvfilename+\
    ' tempfile='+tempfilename+' outfile='+he_dir+'he_gti.fits ehkfile='+ehkfilename+\
    ' defaultexpr=yes rangefile='+args.rangefile
    os.system(hegtigen_text) 
    print hegtigen_text
    
    hescreen_text = 'hescreen evtfile="'+he_dir+'he_pi.fits" gtifile="'+he_dir+'he_gti.fits" outfile="'+he_dir+'he_screen.fits"'\
    ' baddetfile="" detid="det_id>=0" eventtype=1 anticoincidence=""'
    os.system(hescreen_text)
    print hescreen_text

# read Event data and space craft data 
raw_data,det_id,channel,pulse_width,acd,event_type = readdata(filename)
raw_data_screen,pi_screen,det_id_screen =readscreen(he_dir+'he_screen.fits')
time_ehk,sat_alt,sat_lon,sat_lat,elv,dye_elv,ang_dist,cor,saa_flag,sun_ang,moon_ang = readehk(ehkfilename)
time_pm,pm_cnt = readpm(pmfilename)
time_hv,pho_hv = readhv(hvfilename,phonum=0)

# plot light curve and all space craft data
if args.fig:
    plt.clf()
    plt.figure()
    plt.subplots_adjust(bottom=0)
    
    plt.subplot(9,1,1)
    lc_x,lc_y = genlc(raw_data,binsize=1,fig=False,rate=True,pannel=False)
    plt.step(lc_x,lc_y)
    
    plt.subplot(9,1,2)
    plt.plot(time_ehk,elv)
    plt.ylabel('elv')
    
    plt.subplot(9,1,3)
    plt.plot(time_ehk,dye_elv)
    plt.ylabel('dye_elv')
    
    plt.subplot(9,1,4)
    plt.plot(time_ehk,cor)
    plt.ylabel('cor')
    
    plt.subplot(9,1,5)
    plt.plot(time_pm,pm_cnt)
    plt.ylabel('pm')
    
    plt.subplot(9,1,6)
    plt.plot(time_ehk,saa_flag)
    plt.ylabel('saa_flag')
    
    plt.subplot(9,1,7)
    plt.plot(time_hv,pho_hv)
    plt.ylabel('hv')
    
    plt.subplot(9,1,8)
    plt.plot(time_ehk,sat_lat)
    plt.ylabel('lat')
    plt.savefig(he_dir+'craft.png')

    plt.figure()
    lc_x,lc_y = genlc(raw_data_screen,binsize=0.1,fig=False,rate=True,pannel=True)
    plt.step(lc_x,lc_y)
    plt.savefig(he_dir+'lc.png')
    
    





