from function import *

# find and create data dir list

if len(sys.argv) >= 3:
    data_dir = sys.argv[1]
    product_dir = sys.argv[2]
    aux_dir = product_dir + "/AUX/" # AUX path
    acs_dir = product_dir + "/ACS/" # ACS path
    ehkfilename = aux_dir + "EHK.fits"
else:
    data_dir = str(raw_input("Archive data path:"))
    product_dir = str(raw_input("Results out path:"))
    aux_dir = product_dir + "/AUX/" # AUX path
    acs_dir = product_dir + "/ACS/" # ACS path
    ehkfilename = aux_dir + "EHK.fits"

print aux_dir
print acs_dir
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
print hvfilename
print pmfilename

# generate ehk file using HXMT software
if not os.path.isfile(ehkfilename):
    ehkgen(product_dir)

# read Event data and space craft data 
raw_data,det_id,channel,pulse_width,acd,event_type = readdata(filename)
time_ehk,sat_alt,sat_lon,sat_lat,elv,dye_elv,ang_dist,cor,saa_flag,sun_ang,moon_ang = readehk(ehkfilename)
time_pm,pm_cnt = readpm(pmfilename)
time_hv,pho_hv = readhv(hvfilename,phonum=0)

# plot light curve and all space craft data
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

if sys.argv[3]=='fig=True':plt.show()






