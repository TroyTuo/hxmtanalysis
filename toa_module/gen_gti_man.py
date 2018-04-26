import numpy as np 
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i","--input",help="select instrument for HXMT(HE/ME/LE)")
parser.add_argument("-p","--epoch",help="select epoch before or after glitch(pre/post)")
args = parser.parse_args()
instru = args.input
pepoch = args.epoch

print instru
if instru == 'HE' or instru == 'he':
    gtilist = '/home/tuoyl/Data_analysis/Crab/gti_list_he.lst' # the time list of time intervals where counts rate is bad
    out_gti = '/home/tuoyl/Data_analysis/TimingTool/Glitch/gti_he.tim' # file name of output gti list for fselect 
    filename = 'he_screen_NaI_ang.fits' # input filename for fselect
    seloutname = filename # output filename for fselect
    selpioutname = filename # output filename for select channel
    seldetoutname = selpioutname # input filename for detector selection
    pathlistname = '/home/tuoyl/Data_analysis/TimingTool/Glitch/'+\
            'he_pre_glitch_path.lst_test' # file path for those carry fselect
    if pepoch == 'post':
        pathlistname = '/home/tuoyl/Data_analysis/TimingTool/Glitch/'+\
                'he_post_glitch_path.lst' # file path for those carry fselect

print instru
if instru == 'ME' or instru == 'me':
    print instru
    gtilist = '/home/tuoyl/Data_analysis/Crab/gti_list_me.lst' # the time list of time intervals where counts rate is bad
    out_gti = '/home/tuoyl/Data_analysis/TimingTool/Glitch/gti_me.tim' # file name of output gti list for fselect 
    detlist = '/home/tuoyl/Database/script/me_det_sel.expr' # me bad detector list for fselect
    detlist = '/home/tuoyl/Database/script/me_blind.expr' # me blind detector list for fselect
    filename = 'me_screen.fits' # input filename for fselect
    seloutname = 'me_screen_blind.fits' # output filename for fselect
    seldetoutname = 'me_screen_blind.fits' # input filename for detector selection
    pathlistname = '/home/tuoyl/Data_analysis/TimingTool/Glitch/'+\
            'me_pre_glitch_path.lst' # file path for those carry fselect
    if pepoch == 'post':
        pathlistname = '/home/tuoyl/Data_analysis/TimingTool/Glitch/'+\
                'me_post_glitch_path.lst' # file path for those carry fselect

if instru == 'LE' or instru == 'le':
    gtilist = '/home/tuoyl/Data_analysis/Crab/gti_list_le.lst_post' # the time list of time intervals where counts rate is bad
    out_gti = '/home/tuoyl/Data_analysis/TimingTool/Glitch/gti_le.tim' # file name of output gti list for fselect 
    filename = 'le_screen_ang.fits' # input filename for fselect
    seloutname = filename # output filename for fselect
    seldetoutname = seloutname # input filename for detector selection
    pathlistname = '/home/tuoyl/Data_analysis/TimingTool/Glitch/'+\
            'le_pre_glitch_path.lst' # file path for those carry fselect
    if pepoch == 'post':
        pathlistname = '/home/tuoyl/Data_analysis/TimingTool/Glitch/'+\
                'le_post_glitch_path.lst' # file path for those carry fselect



data = np.loadtxt(gtilist)
print data
with open(out_gti,'w')as f:
    for i in xrange(len(data)):
        if i == len(data)-1:
            string = "(Time <= "+str(data[i][0])+" || Time >= "+str(data[i][1])+ ")"
        else:
            string = "(Time <= "+str(data[i][0])+" || Time >= "+str(data[i][1])+ ") && "
        f.write(string)

with open(pathlistname,'r')as f:
    lines = f.read()
    pathlist = lines.split('\n')
    pathlist = pathlist[0:-1]
    
for i in xrange(len(pathlist)):
    fselect_text = 'fselect ' + pathlist[i] + filename + ' ' + pathlist[i] + seloutname + ' @'+out_gti+ ' clobber=yes > $glitch/ang_sel.log'
    print fselect_text
    os.system(fselect_text)
    if instru == 'HE' or instru == 'he':
    #    ### select PI 
    #    fselect_pi_text = 'fselect ' + pathlist[i] + filename + ' ' + pathlist[i] + selpioutname + ' "PI >= 20 && PI <=200" clobber=yes > $glitch/ang_sel.log'
    #    print fselect_pi_text
    #    os.system(fselect_pi_text)
        ### select detector 
        fselect_det_text = 'fselect ' + pathlist[i] + filename + ' ' + pathlist[i] + seldetoutname + ' "detid <=1 || (detid >= 3 && detid <= 8) || (detid >= 10 && detid <= 15) || detid == 17" clobber=yes > $glitch/ang_sel.log'
        print fselect_det_text
        os.system(fselect_det_text)
    if instru == 'ME' or instru == 'me':
        ### select detector
        fselect_det_text = 'fselect ' + pathlist[i] + seloutname + ' ' + pathlist[i] + seldetoutname + ' @'+detlist+ ' clobber=yes > $glitch/ang_sel.log'
        print fselect_det_text
        os.system(fselect_det_text)












