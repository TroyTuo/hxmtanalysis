 #! /bin/sh

'''
 Ver 0.1: basic script for HXMT data process
'''
'''
 Ver 0.2: Add configure file for bat processing HXMT data
'''
'''
 Ver 0.3: revise for 1R level data
'''
'''
Define the input and output parameters  for HXMT/HE data process
'''
'''
Basci script fromat for Dictionary
('script_comm_name', 'hepical')
('script_para_dict', {'outfile': 'outfile', 'evtfile': 'evtfile'})
('script_expr_dict', {'outfile': ' outfile=', 'evtfile': ' evtfile = '})
('script_valu_dict', {'outfile': './he_out.fits', 'evtfile': './hv.fits'})
('script_para_orde', ['evtfile', 'outfile'])
'''

'''
Only Generate the shell script for this Version
'''

hpipe_line_version='0.3'

import os
import sys
import commands

from xml.etree import ElementTree as ET

def Read_Configure(config_name, command_name):
    
    comm_str='./Command'
    com_namestr='com_name'
    expr_str='expr'
    valu_str='default'
    '''Load configure file'''
    root=ET.parse(config_name)
    comm=root.findall(comm_str)
    para_dict={}
    para_expr={}
    para_valu={}
    para_orde=[]
    
    '''Search the given command parameters'''
    for onecomm in comm:
        comm_attr=onecomm.attrib
        
        if command_name == comm_attr[com_namestr]:
            '''Format the parameters dictionary'''
            order=0
            for parameter in onecomm:
                attr= parameter.attrib
                para_dict[parameter.text]=parameter.text
                para_expr[parameter.text]=' '+ parameter.text + '='
                para_valu[parameter.text]=attr[valu_str]
                para_orde.append( parameter.text )
                
    comm_para=[command_name,para_dict,para_expr,para_valu,para_orde]
    
    return comm_para


def Format_Input():
    print 'Hello'

def Read_ParaInput():
    print 'Hello'


#para_dict=Read_Configure('./configure_he.xml','hegtigen')
#print len( para_dict )

'''Definition of basic parameters for the data process

Read the basic paramenters from the configure file
Format the script from the input parameters

'''

class script_tele:
    
    script_comm_name=''
    script_para_dict={}
    script_expr_dict={}
    script_valu_dict={}
    script_para_orde=[]
    hescript=''
    isprint=0
    '''Read configure content for special inst'''
    def __init__(self, configname, comm_name, para_list,isprint=0):
        
        result=Read_Configure(configname,comm_name)
        self.script_comm_name=result[0]
        self.script_para_dict=result[1]
        self.script_expr_dict=result[2]
        self.script_valu_dict=result[3]
        self.script_para_orde=result[4]
        self.isprint=isprint
        if (self.isprint == 1):
            print("Generating Script Name: ",result[0])
        tmppara=self.script_para_dict.keys()
        para_index=0
        for para_index in xrange(0,len(tmppara)):
            '''
            Check the input values:
            If the value is 'default', the script_valu_dict should be the content
            from the configure file of the given the command
            Or the values will be changed to the new ones
            '''
            if (self.isprint == 1):
                print('Para name= ',self.script_para_orde[para_index],' Value= ',para_list[para_index])
            if 'default' != para_list[para_index].lower():
                self.script_valu_dict[self.script_para_orde[para_index]]=para_list[para_index]
                
    
    '''Print informations'''
    def Print_Para(self):
        print self.script_para_dict
        print self.script_expr_dict
        print self.script_valu_dict
    '''Create the shell script according the congfigure file and input parameters'''
    def Script_Format(self):
        self.hescript=self.script_comm_name
        para_keys=self.script_expr_dict
        para_valu=self.script_valu_dict
        para_orde=self.script_para_orde
        
        for para_index in xrange(0,len(para_valu)):
            self.hescript += para_keys[para_orde[para_index]]+para_valu[para_orde[para_index]]
            
        return self.hescript


'''
Obtain the file information
for the given the file path
walk_dir: Basic function for obtaining the file information
get_fileinfo: Get all file names with given subdirections
get_filename: Get specific file name
'''
def walk_dir(dir,fileinfo,topdown=True):
    for root, dirs, files in os.walk(dir, topdown):
        for name in files:
            fileinfo.append(os.path.join(root,name))

def get_fileinfo(dirs,subdir,fileinfo,topdown=True):
    subdirs=dirs+subdir
    walk_dir(subdirs,fileinfo,topdown=True)
    del fileinfo[0]
def get_filename(fileinfo,suffix,outname,npos0):
    npos0=-1
    for fname in fileinfo:
        tmpstr=fname.split('/')
        tmpstr2=tmpstr[len(tmpstr)-1]
        npos=tmpstr2.find(suffix)
        if(npos >= 0):
            outname.append(fname)
            npos0=npos
    del outname[0]

'''
Search the data for the given path
'''

def hpipe_interface(basedir,outdir):
    para_print=0
    
    HE_fileinfo=['']
    ME_fileinfo=['']
    LE_fileinfo=['']
   
    he_evtname=['']
    he_temp=['']
    he_dead=['']
    he_hv=['']
    
    AUX_fileinfo=['']
    ehk_name=['']
    ACS_fileinfo=['']
    orb_name=['']
    att_name=['']
    
    npos=-1
    '''Get HE files and HouseKeeping files'''
    get_fileinfo(basedir,'/HE/',HE_fileinfo,topdown=True)
    get_filename(HE_fileinfo,'HE-Evt_FFFFFF_V',he_evtname,npos)
    print he_evtname
    for evtindex in range(len(he_evtname)):
        get_filename(HE_fileinfo,'TH',he_temp,npos)
        print he_temp
        get_filename(HE_fileinfo,'DTime',he_dead,npos)
        print he_dead 
        get_filename(HE_fileinfo,'HV',he_hv,npos)
        print he_hv
        
        '''Get auxiliary files and ACS and Orbit files'''
        get_fileinfo(basedir,'/AUX',AUX_fileinfo,topdown=True)
        get_filename(AUX_fileinfo,'EHK_man.fits',ehk_name,npos)
        
        get_fileinfo(basedir,'/ACS',ACS_fileinfo,topdown=True)
        get_filename(ACS_fileinfo,'Orbit',orb_name,npos)
        get_fileinfo(basedir,'/ACS',ACS_fileinfo,topdown=True)
        get_filename(ACS_fileinfo,'Att',att_name,npos)

        '''
        Format the paralist
        '''

        para_list=['default','default','default','default','default','default','default','default','default','default','default','default','default','default','default']
        print "Input dir: ",basedir
        print "Output dir: ",outdir
        
        fileid=he_evtname[evtindex][he_evtname[evtindex].find('HXMT_P0101'):he_evtname[evtindex].find('.FITS')]
        print 'fileid',fileid
        print heconfigname
        para_list_pical=[he_evtname[evtindex],outdir+'/'+fileid+'_pi.fits']
        
        para_list_hegti=[he_hv[0],he_temp[0],ehk_name[0],outdir+'/'+fileid+'_gti.fits','default']
        para_list_hescr=[para_list_pical[1],para_list_hegti[3],outdir+'/'+fileid+'_screen.fits','default','default','default','default','default','default','default','default']

        para_list_spec=[para_list_hescr[2],outdir+'/'+fileid+'_spec',he_dead[0],'default','default','default']
        para_list_lc=[para_list_hescr[2],outdir+'/'+fileid+'_lc',he_dead[0],'default','default','default','default','default']
#        para_list_rsp=[para_list_spec[1],ACS_fileinfo[0],outdir+'./hersp.fits']

        print "*********************************************************"
        print "*******************Mini script for HE.*******************"
        print "*********************************************************"
#        heconfigname='./configure_hpipe.xml'
        para=script_tele(heconfigname,'hepical',para_list_pical,para_print)
        hepical= para.Script_Format()
        para=script_tele(heconfigname,'hegtigen',para_list_hegti,para_print)
        hegtigen=para.Script_Format()
        para=script_tele(heconfigname,'hescreen',para_list_hescr,para_print)
        hescreen=para.Script_Format()
        para=script_tele(heconfigname,'hegenspec',para_list_spec,para_print)
        hegenspec=para.Script_Format()
        para=script_tele(heconfigname,'hegenlc',para_list_lc,para_print)
        hegenlc=para.Script_Format()
        #para=script_tele(heconfigname,'herspgen',para_list_rsp,para_print)
        #herspgen=para.Script_Format()
        
        [hepical_status,output]=commands.getstatusoutput(hepical)
        print hepical_status
        print output
        [hegtigen_status,output]=commands.getstatusoutput(hegtigen)
        print hegtigen
        print hepical_status
        print output
        [hescreen_status,output]=commands.getstatusoutput(hescreen)
        print hescreen
        print hepical_status
        print output
        #[hegenspec_status,output]=commands.getstatusoutput(hegenspec)
        print hegenspec
        print hepical_status
        print output
        [hegenlc_status,output]=commands.getstatusoutput(hegenlc)
        print hegenlc
        print hepical_status
        print output
#        [herspgen_status,output]=commands.getstatusoutput(herspgen)
#        print hepical_status
#        print output
        
#        fo=open('HE.sh','w')
#        fo.write(hepical)
#        fo.write('\n')
#        fo.write(hegtigen)
#        fo.write('\n')
#        fo.write(hescreen)
#        fo.write('\n')
#        fo.write(hegenspec)
#        fo.write('\n')
#        fo.write(hegenlc)
#        fo.write('\n')
#        fo.write(herspgen)
#        fo.close()
##        
##        print "*********************************************************"
##        print "*******************Mini script for ME.*******************"
##        print "*********************************************************"
##        me_evtname=['']
##        meconfigname='./configure_hpipe.xml'
##        para=script_tele(meconfigname,'mepical',para_list,para_print)
##        mepical=para.Script_Format()
##    
##        print "*********************************************************"
##        print "*******************Mini script for LE.*******************"
##        print "*********************************************************"
##        leconfigname='./configure_hpipe.xml'
##        le_evtname=['']
##        le_temp=['']
##        
##        get_fileinfo(basedir,'/LE',LE_fileinfo,topdown=True)
##        get_filename(LE_fileinfo,'SciEvt',le_evtname,npos)
##        get_filename(LE_fileinfo,'TH',le_temp,npos)
##        
##        para_list_lepi=[le_evtname[0],le_temp[0],outdir+'./le_pi.fits']
##        para_list_legti=[ehk_name[0],outdir+'./le_gti.fits','default']
##        para_list_ledet=[para_list_lepi[2],outdir+'./le_det.fits','default','default','default','default','default']
##        para_list_lerecon=[para_list_ledet[1],outdir+'./le_det_recon.fits']
##        para_list_lescr=[para_list_lerecon[1],para_list_legti[1],outdir+'./le_det_recon_scr.fits','default','default']
##        para_list_lespec=[para_list_lescr[2],'default','default','default','default','default',outdir+'./le_det_spec.fits']
##        para_list_lelc=[para_list_lescr[2],outdir+'./le_det_lc.fits','default','default','default','default','default','default']
##        para_list_lersp=[para_list_lespec[6],ACS_fileinfo[0],le_temp[0],outdir+'./le_rsp.fits']
##        
##        para=script_tele(leconfigname,'lepical',para_list_lepi,para_print)
##        lepical=para.Script_Format()
##        para=script_tele(leconfigname,'legtigen',para_list_legti,para_print)
##        legtigen=para.Script_Format()
##        para=script_tele(leconfigname,'ledetectorsel',para_list_ledet,para_print)
##        ledet=para.Script_Format()
##        para=script_tele(leconfigname,'lerecon',para_list_lerecon,para_print)
##        lerecon=para.Script_Format()
##        para=script_tele(leconfigname,'lescreen',para_list_lescr,para_print)
##        lescr=para.Script_Format()
##        para=script_tele(leconfigname,'legenspec',para_list_lespec,para_print)
##        lespec=para.Script_Format()
##        para=script_tele(leconfigname,'legenlc',para_list_lelc,para_print)
##        lelc=para.Script_Format()
##        para=script_tele(leconfigname,'lerspgen',para_list_lersp,para_print)
##        lersp=para.Script_Format()
##        
##        [lepical_status,output]=commands.getstatusoutput(lepical)
##        print lepical_status
##        print output
##        [legtigen_status,output]=commands.getstatusoutput(legtigen)
##        print legtigen_status
##        print output
##        [ledet_status,output]=commands.getstatusoutput(ledet)
##        print ledet_status
##        print output
##        [lerecon_status,output]=commands.getstatusoutput(lerecon)
##        print lerecon_status
##        print output
##        [lescr_status,output]=commands.getstatusoutput(lescr)
##        print lescr_status
##        print output
##        [lespec_status,output]=commands.getstatusoutput(lespec)
##        print lespec_status
##        print output
##        [lelc_status,output]=commands.getstatusoutput(lelc)
##        print lelc_status
##        print output
##        [lersp_status,output]=commands.getstatusoutput(lersp)
##        print lersp_status
##        print output
##        fo=open('LE.sh','w')
##        fo.write(lepical)
##        fo.write('\n')
##        fo.write(legtigen)
##        fo.write('\n')
##        fo.write(ledet)
##        fo.write('\n')
##        fo.write(lerecon)
##        fo.write('\n')
##        fo.write(lescr)
##        fo.write('\n')
##        fo.write(lespec)
##        fo.write('\n')
##        fo.write(lelc)
##        fo.write('\n')
##        fo.write(lersp)
##        fo.close()
##    
heconfigname='./configure_hpipe.xml'
basedir='../ana/P010100101001/P01010010100101-20170306-01-01/'
outdir='../ana/P010100101001/P01010010100101-20170306-01-01/stdpro/'

print "*********************************************************"
print "******************Running HXMT Pipeline******************"
print ("****************    Version %s       *******************" % hpipe_line_version)
print "************'python hpipeline.py -h' for help************"
print "*********************************************************"

if len(sys.argv)==2:
    if sys.argv[1]=='-h':
        print "Example 1: python hpipeline.py indir outdir"
        print "Example 2: Using interactive in prompt."
elif len(sys.argv)>=2:
    basedir=sys.argv[1]
    outdir=sys.argv[2]
    heconfigname=sys.argv[3]
    hpipe_interface(basedir,outdir)
else:
    basedir=str(raw_input("Archive data path:"))
    outdir =str(raw_input("Results out path:"))
    heconfigname=str(raw_input("Configure file path"))
    if len( basedir) == 0:
        print "Please input archive data path!"
        print "Type -h for help"
    else:
        hpipe_interface(basedir,outdir)

'''
TODO: Generate the values of the input parameters 
TODO: Search the basic the input files from Level 1 data
TODO: Excute the script
TODO: Add the message block
TODO: Log file design
TODO: Data file structure
'''














