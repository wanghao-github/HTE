# -*- coding: utf-8 -*-
import cPickle as pickle
import sys
import os
import commands
import shutil
#import gzip
import re
import glob
import warnings
import numpy as np
from itertools import permutations
from copy import deepcopy
from math import *
from ase import Atoms, Atom
from ase.calculators.vasp import Vasp
from ase.data import chemical_symbols
from ase.io.cif import *
from ase.io import write
from ase.io.wien2k import write_struct
from ase.io.wien2k import read_struct
from ase.visualize import view
from hte_dbentry import *
from transport import *
from ase.structure import bulk
from ase.data import *
from fplo import *
from ase.lattice.spacegroup import *
import ase.units as units
from magtools import *
from tools import Identify_Interstitial_Position
from tools.chull import *

class HTE(object):
 
    ############################################
    # routines to manage HTE                   #
    ############################################
    def __init__(self, dbfile=None, hostname='lichtenberg',silent=False):
        self.version="1.2"
        self.subversion="200514"
        self.main_directory=os.getcwd()
        self.cif_repository=".HTE_cif_repository"
        self.calc_schemes={}
        self.structureDB={}
        self.dbfile=dbfile
        self.job_commands={}
        self.searchpaths=['./',os.getcwd()]
        self.backupfile=None
        self.max_jobs_in_queue={'total':2000,'natoms_max':{20:600}}
        #dictionary to store tempory data which is read in only once during run time and removed afterwards
        self.tmpdata={'reference_energies':{},'prop_dict':{},'chulls':{}}
        self.dir_of_hte=os.path.split(os.path.abspath(__file__))[0]
        self.nsub_max=3
        self.log_messages=[]
        self.storage_options={'vasp':{'bzip2':['OUTCAR','out'],'cp':['POSCAR.INI','POSCAR','CONTCAR','KPOINTS','INCAR','hte_propdict.txt']},
                              'fplo':{'bzip2':['=.dens','out'],'cp':['=.in','hte_propdict.txt']}
                              }
        self.use_prop_dict=True
        self.use_job_array=False
        self.magnetic_atoms={'Mn':3.5,'Fe':2.5,'Co':1.5,'Ni':0.8}
        #directory to (permanently) store converged calculations in compressed form
        self.storage_directories=[]
        if silent==False:
            print "*** HTE version %s"%self.get_version()
            print "*** Author: Ingo Opahle"
        # local settings
        if hostname=='neptune': #ICAMS local host
            fplo_commands={'fedit':'fedit9.01-35-i386','job_command':'fplo9.01-35-i386 > out','job_file':None,'job_environment':None}
            self.job_commands['fplo']=fplo_commands
        elif  hostname=='lichtenberg': #TUD cluster with SLURM
            vasp_commands={'job_command':'bash prep-job-vasp.sh ; sbatch job-vasp.sh','job_file':'job-vasp/*','job_environment':'SLURM'}
            self.job_commands['vasp']=vasp_commands
            print "check_point49, vasp_commands should be:",vasp_commands      
            fplo_commands={'fedit':'fedit14.00-49-x86_64','job_command':'bash prep-job-fplo.sh ; sbatch job-fplo.sh','job_file':'job-fplo/*','job_environment':'SLURM'}
            self.job_commands['fplo']=fplo_commands
            self.job_commands['array_job']={'vasp': {'job_command':'bash prep-job-vasp.sh','array_job':True,'job_file':'jobarray-vasp/*','job_environment':'SLURM'}}        
            print "check_point51, self.job_commands['array_job']= ",self.job_commands['array_job']
        else: #ICAMS vulcan cluster
            fplo_commands={'fedit':'fedit9.01-35-x86_64','job_command':'qsub job-fplo.sh','job_file':'job-fplo.sh','job_environment':'SGE'}
            self.job_commands['fplo']=fplo_commands
            vasp_commands={'job_command':'qsub job.sh','job_file':'job.sh','job_environment':'SGE'}
            self.job_commands['vasp']=vasp_commands
        #
        if (dbfile==None) or (not  os.path.isfile(dbfile)):
            print 'Setting up default HTE'
            # VASP
            #ultrasoft pseudopotentials (by convention xc = 'PW91') for scan
            setups_uspp={'C': '_s', 'B': '_s', 'Ba': '_pv', 'N': '_s'}
            self.add_calc_scheme('vaspopt','vasp', ispin='auto', xc = 'PBE', encut = 500,  nsw = 60, lwave= False, lcharg = False, kspace_density=40, gamma=True, prec='High', ismear=0,ibrion=2, isif=2)
            self.add_calc_scheme('vasp-scf-hao','vasp', ispin='auto', xc = 'PBE', encut = 500,  nsw = 0, lwave= False, lcharg = False, kspace_density=50, gamma=True, prec='High', ismear=-5,init_structure='vaspopt')
            self.add_calc_scheme('us-k30-e250-opt','vasp', xc = 'PW91', encut = 250, enaug = 250, nsw = 40, lreal=True, lwave= False, lcharg = False, kspace_density=30, gamma=True, ibrion=2, isif=3,setups=setups_uspp)
            #self.add_calc_scheme('us-k30-e150-noopt','vasp', xc = 'PW91', encut = 150, enaug = 150, nsw = 0, lreal=True, lwave= False, lcharg = False, kspace_density=30, gamma=True)
            #self.add_calc_scheme('us-k30-e150-vscale','vasp', xc = 'PW91', encut = 150, enaug = 150, nsw = 0, lreal=True, lwave= False, lcharg = False, kspace_density=30, gamma=True, scale_volume='atomic_volume')
            #self.add_calc_scheme('us-k30-e150-vscale-opt','vasp', xc = 'PW91', encut = 150, enaug = 150, nsw = 40, lreal=True, lwave= False, lcharg = False, kspace_density=30, gamma=True, ibrion=2, isif=3, scale_volume='atomic_volume')
            #self.add_calc_scheme('us-k30-e150-opt','vasp', xc = 'PW91', encut = 150, enaug = 150, nsw = 40, lreal=True, lwave= False, lcharg = False, kspace_density=30, gamma=True, ibrion=2, isif=3)
            #PAW  pseudopotentials (by convention xc = 'PBE') for higher accuracy
            setups_paw={'C': '_s', 'Ba': '_sv', 'F': '_s', 'Sr': '_sv', 'Ca': '_sv', 'Fr': '_sv', 'O': '_s', 'N': '_s', 'Cs': '_sv', 'Ra': '_sv', 'Rb': '_sv', 'Y': '_sv', 'Nb': '_sv', 'K': '_sv', 'Zr': '_sv'}
            setups_hongbin={'Na': '_pv', 'Nb': '_pv', 'Mg': '', 'Pb': '_d', 'Sc': '_sv', 'Tl': '_d', 'Ra': '_sv', 'Rb': '_sv', 'Ti': '_sv', 'Fr': '_sv', 'Ba': '_sv', 'Bi': '_d', 'C': '_s', 'B': '_s', 'F': '_s', 'Sr': '_sv', 'K': '_sv', 'O': '_s', 'N': '_s', 'Sn': '_d', 'V': '_sv', 'Y': '_sv', 'Ca': '_sv', 'Ge': '_d', 'Ga': '_d', 'In': '_d', 'Cs': '_sv', 'Zr': '_sv'}
            self.add_calc_scheme('paw-k40-ecut350-noopt','vasp', xc = 'PBE', encut = 350, enaug = 350, nsw = 0, lwave= False, lcharg = False, kspace_density=40, gamma=True, prec='High')
            self.add_calc_scheme('paw-k40-ecut350-opt','vasp', xc = 'PBE', encut = 350, enaug = 350, nsw = 40, lwave= False, lcharg = False, kspace_density=40, gamma=True, prec='High', ibrion=2, isif=3,setups={ 'Y': '_sv'})
            self.add_calc_scheme('vasp-default','vasp', ispin='auto', xc = 'PBE', encut = 350, enaug = 350, nsw = 40, lwave= False, lcharg = False, kspace_density=40, gamma=True, prec='High', ibrion=2, isif=3,init_structure=['us-k30-e250-opt','paw-k40-ecut350-opt'], scale_volume_mag=1.1,setups=setups_paw)
            #high pressure schemes
            self.add_calc_scheme('vasp-default-P5GPa','vasp', pstress=50, ispin='auto', xc = 'PBE', encut = 350, enaug = 350, nsw = 40, lwave= False, lcharg = False, kspace_density=40, gamma=True, prec='High', ibrion=2, isif=3, setups=setups_paw, init_structure='vasp-default')
            self.add_calc_scheme('vasp-default-P10GPa','vasp', pstress=100, ispin='auto', xc = 'PBE', encut = 350, enaug = 350, nsw = 40, lwave= False, lcharg = False, kspace_density=40, gamma=True, prec='High', ibrion=2, isif=3, setups=setups_paw, init_structure='vasp-default-P5GPa')
            self.add_calc_scheme('vasp-default-P15GPa','vasp', pstress=150, ispin='auto', xc = 'PBE', encut = 350, enaug = 350, nsw = 40, lwave= False, lcharg = False, kspace_density=40, gamma=True, prec='High', ibrion=2, isif=3, setups=setups_paw, init_structure='vasp-default-P10GPa')
            #PAW  pseudopotentials with semi core (by convention xc = 'PBE') for higher accuracy
            #self.add_calc_scheme('paw-k50tet-ecut500-sc-noopt','vasp', xc = 'PBE', ismear=-5, encut = 500, enaug = 500, nsw = 0, lwave= False, lcharg = False, setups={'Mg':'_pv', 'Si':'_h'}, kspace_density=50, gamma=True, prec='High')
            self.add_calc_scheme('paw-k50tet-ecut500-sc-opt','vasp', xc = 'PBE', ismear=-5, encut = 500, enaug = 500, nsw = 40, lwave= False, lcharg = False, setups={'Mg':'_pv', 'Si':'_h'}, kspace_density=50, gamma=True, prec='High', ibrion=2, isif=3)
            setups_paw_hongbin={'Na': '_pv', 'Nb': '_pv', 'Mg': '', 'Pb': '_d', 'Sc': '_sv', 'Tl': '_d', 'Ra': '_sv', 'Rb': '_sv', 'Ti': '_sv', 'Fr': '_sv', 'Ba': '_sv', 'Bi': '_d', 'C': '_s', 'B': '_s', 'F': '_s', 'Sr': '_sv', 'K': '_sv', 'O': '_s', 'N': '_s', 'Sn': '_d', 'V': '_sv', 'Y': '_sv', 'Ca': '_sv', 'Ge': '_d', 'Ga': '_d', 'In': '_d', 'Cs': '_sv', 'Zr': '_sv'}
            self.add_calc_scheme('vasp-accurate-opt','vasp', ispin='auto', xc = 'PBE', encut = 500,  nsw = 40, lwave= False, lcharg = False, kspace_density=50, gamma=True, prec='High', ibrion=2, isif=3, setups=setups_paw_hongbin, init_structure='vasp-default')
            self.add_calc_scheme('vasp-accurate','vasp', ispin='auto', xc = 'PBE', encut = 500,  nsw = 0, lwave= False, lcharg = False, kspace_density=50, gamma=True, prec='High', ismear=-5, setups=setups_paw_hongbin, init_structure='vasp-accurate-opt')
            # FPLO
            self.add_calc_scheme('fplo-gga-default','fplo')
            self.add_calc_scheme('fplo-cont-vasp-default','fplo',init_structure='vasp-default',ispin='auto',init_nm=True,kspace_density=40)
        else:
            if not silent:
                print 'Loading database file <%s>'%dbfile
            fin=open(dbfile,"r")
            hte_loaded=pickle.load(fin)
            fin.close()
            #check if version of db file is up to date
            if (hasattr(hte_loaded,'version')==False):
                if silent==False:
                    print "Upgrading database from HTE version 0.x to %s"%self.version
                self.backupfile=self.dbfile+'_HTE0.x'
                if silent==False:
                    print "Backup file for old version: <%s>"%self.backupfile
                #upgrade dbentries
                for uid in hte_loaded.structureDB:
                    dbentry_old=hte_loaded.structureDB[uid]
                    atoms_initial=dbentry_old.atoms_initial
                    source_info=self.get_source_info_old_version(uid,hte_loaded.structureDB)
                    dbentry=HTEdbentry(atoms_initial,uid,_source_info=source_info)
                    dbentry.get_data_from_cif_source()
                    if hasattr(dbentry_old,'submitted_jobs'):
                        dbentry.submitted_jobs=dbentry_old.submitted_jobs
                    if hasattr(dbentry_old,'stored_calc_results'):
                        dbentry.stored_calc_results=dbentry_old.stored_calc_results
                    if hasattr(dbentry_old,'do_not_calc'):
                        dbentry.do_not_calc=dbentry_old.do_not_calc
                    self.add_structDB_entry(uid, structureDB_entry=dbentry,_source_info=source_info, silent=True)
            else:
                self.structureDB=hte_loaded.structureDB
            if hasattr(hte_loaded,'calc_schemes'):
                self.calc_schemes=hte_loaded.calc_schemes
            if hasattr(hte_loaded,'searchpaths'):
                self.searchpaths=hte_loaded.searchpaths
            if hasattr(hte_loaded,'nsub_max'):
                self.nsub_max=hte_loaded.nsub_max
            if hasattr(hte_loaded,'storage_directories'):
                self.storage_directories=hte_loaded.storage_directories
        if silent==False:
            print "Entries in database: %d"%len(self.structureDB)
            print "Main storage:",self.get_main_storage_directory()
            print "Storage directories:",self.get_storage_directories()
            print "Scratch directory:",self.get_scratch_directory()
            print "Searchpaths:",self.get_searchpaths()
                
    def set_storage_directories(self,dirlist):
        if not isinstance(dirlist,list):
            print "WARNING(set_storage_directories): argument must be list of directories, nothing done!"
            return self.storage_directories
        self.storage_directories=dirlist
        #for stdir in dirlist:
        #    if not os.path.isdir(stdir):
        #        exitcode, out = commands.getstatusoutput("mkdir -p %s"%stdir)
        #        #os.mkdir(stdir)
        #    self.storage_directories.append(stdir)
        return self.storage_directories

    def set_main_storage_directory(self,stdir):
        
        if not os.path.isdir(stdir):
            exitcode, out = commands.getstatusoutput("mkdir -p %s"%stdir)
        if not os.path.isdir(stdir):
            self.add_logmessage("WARNING(set_main_storage_directory): Failed to create <%s>, nothing done!"%stdir)
        else:
            stordirs=[stdir]
            for stdir in self.storage_directories:
                if not stdir in stordirs:
                    stordirs.append(stdir)
            self.storage_directories=stordirs
        return self.get_main_storage_directory

    def get_storage_directories(self):
        """Returns list of directories where calculations are (supposed to be) safely stored.
        Existence of directories is checked (access permissions still to be done).
        Relative paths are changed to absolute paths!
        """
        storage_directories=[]
        storage_ok=True
        for sdir in self.storage_directories:
            if (sdir!=None) and (not (os.path.abspath(sdir) in storage_directories)) and (os.path.isdir(sdir)): # TODO: test access permissions?
                storage_directories.append(os.path.abspath(sdir))
            elif (sdir!=None) and (not (os.path.isdir(sdir))):
                storage_ok=False
                self.add_logmessage("WARNING(get_storage_directories): Check storage directory <%s>!"%sdir)
        if (storage_ok==False):
            self.add_logmessage("1) Make a backup of your .hte file")
            self.add_logmessage("2) Copy the following lines to a python script and edit the storage directories in the first line, run the script")
            self.add_logmessage("#### start python script")
            self.add_logmessage("storage_directories=%s"%str(self.storage_directories))
            self.add_logmessage("from hte import *")
            self.add_logmessage("hte=HTE('%s', silent=True)"%self.dbfile)
            self.add_logmessage("hte.set_storage_directories(storage_directories)")
            self.add_logmessage("print 'updated storage directories:', hte.get_storage_directories()")
            self.add_logmessage("hte.write()")
            self.add_logmessage("#### end python script")
            self.write_logmessages()
            sys.exit("Failed to access storage_directories, check log file!")
        return storage_directories

    
    def get_main_storage_directory(self):
        """Returns directory where the present HTE database is storing (finished) calculations in a compressed format.
        Existence of directory is checked, but no access permissions.
        Returns None if no main storage is defined (or accessible).
        """
        try:
            if os.path.isdir(self.storage_directories[0]): # TODO: test access permissions?
                main_storage_directory=self.storage_directories[0]
            else:
             main_storage_directory=None
        except:
             main_storage_directory=None
        return main_storage_directory

    def get_storage_options(self,calc_scheme):
        storage_options={}
        if calc_scheme in self.storage_options:
            storage_options=self.storage_options[calc_scheme]
        elif calc_scheme in self.calc_schemes:
            calc_name,settings=self.calc_schemes[calc_scheme]
            if calc_name.lower() in self.storage_options:
                storage_options=self.storage_options[calc_name.lower()]
        return storage_options
    
    def get_scratch_directory(self):
        """Returns directory where calculations are submitted
        Returns None if no scratchdir is defined (or accessible).
        """
        try:
            if not (os.path.isdir(self.searchpaths[0])):
                exitcode, out = commands.getstatusoutput("mkdir -p %s"%self.searchpaths[0])
            scratchdir=os.path.abspath(self.searchpaths[0])
            if not (os.path.isdir(scratchdir)):
                scratchdir=None
                self.add_logmessage("WARNING(get_scratch_directory): no job submission because scratchdir not accessible")
        except:
            scratchdir=None
            self.add_logmessage("WARNING(get_scratch_directory): no job submission because scratchdir not accessible")
        return scratchdir


    def set_scratch_directory(self,scratchdir):
        """Set scratch directory where calculations are submitted
        Returns None if no scratchdir is defined (or accessible).
        """
        searchpaths=[scratchdir]
        for sp in self.searchpaths:
            if not (sp in searchpaths):
                searchpaths.append(sp)
        self.searchpaths=searchpaths
        return self.get_scratch_directory()


    def add_logmessage(self, message):
        if not (message in self.log_messages):
            self.log_messages.append(message)
            print message
        return len(self.log_messages)

    def write_logmessages(self, filename="hte-logmessages.txt"):
        outf=open(filename,"w")
        for mes in self.log_messages:
            outf.write("%s\n"%mes)
        outf.close()
    
        

    def write(self, name=None):
        """Save HTE object, should be done after every call to avoid data loss.
        """
        if self.backupfile!=None:
            shutil.copy(self.dbfile,self.backupfile)
        if (name==None) and (self.dbfile==None):
            print 'No db filename given, nothing done'
        elif name!=None:
            self.dbfile=name
        if self.dbfile!=None:
            print 'Removing temporary data...'
            self.tmpdata={'reference_energies':{},'prop_dict':{},'chulls':{}}
            for uid in self.structureDB:
                self.structureDB[uid].atoms=None
            print 'Writing file ',self.dbfile
            fout = open(self.dbfile, "w")
            pickle.dump(self, fout, protocol=0)
            fout.close()
            main_storage=self.get_main_storage_directory()
            if (main_storage!=None):
                shutil.copy(self.dbfile,main_storage)
            self.write_logmessages()
        ### write 这个函数, name参数里写的定义为名字 如果不是空的就把这个文件复制到主存储目录 把db的参考能量 性质词典和chull三个属性都清空了 然后atoms属性设置为空 然后字节化
    
    def merge(self, dbfile, **kwargs):
        """Merge present HTE object with HTE object from file dbfile
        (e.g. in a different user's directory).
        Structure data and calc_schemes are imported and results of
        other HTE object's calculations are avaiable (make sure that
        file permissions are given).
        kwargs are the same as in select(), e.g.:
        merge('xxx.hte', 'nelements'=2) will upload all binaries.
        """
        parentdir=os.getcwd()
        cif_repository=self.get_cif_repository()
        if os.path.isfile(dbfile):
            print "* Uploading HTE database <%s> with kwargs %s"%(dbfile,str(kwargs))
            dbpath,dbname=os.path.split(dbfile)
            if os.path.isdir(dbpath):
                os.chdir(dbpath)
        else:
            print "* Failed to upload database <%s> (kwargs=%s)"%(dbfile,str(kwargs))
            return None
        mhte=HTE(dbname, silent=True)
        mhte_storage=mhte.get_storage_directories()
        print mhte_storage,mhte.storage_directories
        if mhte_storage!=[]:
            storage=[]
            if self.storage_directories==[]:
                storage=[None]
            for sdir in self.storage_directories:
                storage.append(sdir)
            for sdir in mhte_storage:
                if not (os.path.abspath(sdir) in storage):
                    storage.append(os.path.abspath(sdir))
            if storage!=self.storage_directories:
                print "updating storage directories:"
                print "old:",self.get_storage_directories()
                self.storage_directories=storage
                print "new:",self.get_storage_directories()
                print "main storage:",self.get_main_storage_directory()
        elif hasattr(mhte,'searchpaths'):
            for sp in mhte.get_searchpaths():
                if not (os.path.abspath(sp) in self.searchpaths):
                    self.searchpaths.append(os.path.abspath(sp))
        #else:
        #    self.searchpaths.append(os.getcwd())
        for cs in mhte.calc_schemes:
            if not (cs in self.calc_schemes):
                self.calc_schemes[cs]=mhte.calc_schemes[cs]
            #todo check if cs is defined the same
        selection=mhte.select(**kwargs)
        for uid in selection:
            dbentry=mhte.structureDB[uid]
            cif_file,cif_id,cif_index=dbentry.get_cif_source(methods=['internal'])
            if (cif_file!=None) and (cif_repository!=None):
                shutil.copy(cif_file,os.path.join(parentdir,cif_repository))
            dbentry.comment="** uploaded from <%s> **\n%s"%(dbfile,dbentry.comment)
            if uid in self.structureDB:
                #update submitted job data, todo: check if structures are really identical
                if (hasattr(self.structureDB[uid],'submitted_jobs')) and (hasattr(dbentry,'submitted_jobs')) and (self.structureDB[uid].atoms_initial==dbentry.atoms_initial):
                    for sjob in dbentry.submitted_jobs:
                        if not (sjob in self.structureDB[uid].submitted_jobs):
                            self.structureDB[uid].submitted_jobs[sjob]=deepcopy(dbentry.submitted_jobs[sjob])
                else:
                    self.add_logmessage("WARNING(merge): check %s in <%s>"%(uid,dbfile))
            self.add_structDB_entry(uid, structureDB_entry=dbentry,silent=True, exclude_similar=False)
        os.chdir(parentdir)
        
        
        ### 比较老的HTE文件和新的HTE文件然后合并包括CS 以及每个UID条目等等

    def collect_calculations(self,uid_list='all',calc_schemes='all', include_transport=True, reduce_searchpath='auto'):
        """collect (finished) calculations from other HTE repositories.
        May be used to create a central HTE archive where all calculations
        are stored (default setting: uid_list='all',calc_schemes='all').
        Alternatively this may be used to continue calculations from a
        remote HTE directory, e.g. to run transport calculations (set
        uid_list and calc_schemes for the calculations to be copied).
        reduce_searchpath=True: reduce HTE searchpath to present HTE directory
        reduce_searchpath=False: keep HTE searchpath
        reduce_searchpath='auto': reduce HTE searchpath if all calculations are archieved"""
        print "START collect_calculations"
        sp=self.get_searchpaths()[0]
        if uid_list=='all':
            uid_list=self.select()
        if calc_schemes=='all':
            calc_schemes=self.calc_schemes
        i=0
        n=len(uid_list)
        for uid in uid_list:
            dbentry=self.structureDB[uid]
            i=i+1
            for cs in calc_schemes:
                print "* %d/%d %s %s (%s)"%(i,n,uid,self.structureDB[uid].get_composition(reduce=True),cs)
                convpaths=self.get_converged_calculation_paths(uid, cs)
                if convpaths!=[]:
                    targetpath=os.path.join(sp,dbentry.calcdir,cs)
                    if not (targetpath in convpaths):
                        print "Copying %s to %s"%(convpaths[0],targetpath)
                        exitcode, out = commands.getstatusoutput("mkdir -p %s"%targetpath)
                        rsynccmd="rsync -Caz --exclude '*/' %s/ %s/"%(convpaths[0],targetpath)
                        exitcode, out = commands.getstatusoutput(rsynccmd)
                        print rsynccmd,exitcode, out
                    transppaths=[]
                    if include_transport==True:
                        for cpath in convpaths:
                            transppaths.extend(glob.glob(os.path.join(cpath,'transp-*/')))
                    for tpath in transppaths:
                        kdir=tpath.split('/')[len(tpath.split('/'))-2]
                        targetpath=os.path.join(dbentry.calcdir,cs,kdir)
                        if not (os.path.isdir(targetpath)):
                            if (kdir=='transp-ini'):
                                calc=self.setup_calculator(uid, cs, bandstructure_init=True)
                            elif (kdir.startswith('transp-k')):
                                kdens=int(kdir.split('-k')[1])
                                calc=self.setup_calculator(uid, cs, transport=True, transport_kspace_density=kdens)
                            else:
                                break
                            if (dbentry.check_convergency(tpath,calc)):
                                print "Copying %s to %s"%(tpath,targetpath)
                                exitcode, out = commands.getstatusoutput("mkdir -p %s"%targetpath)
                                if os.path.isfile(os.path.join(tpath,'boltzgen','hte.trace')):
                                    rsynccmd="rsync -Caz %s/ %s/"%(tpath,targetpath)
                                else: #boltzrap still needs to be run
                                    rsynccmd="rsync -Caz --exclude '*/' %s/ %s/"%(tpath,targetpath)
                                exitcode, out = commands.getstatusoutput(rsynccmd)
                                print rsynccmd,exitcode, out
                                    

        self.searchpaths=["./",os.getcwd()]
        print "END collect_calculations"
        
    def get_searchpaths(self,check=True):
        """return paths where HTE objects looks for results of calculations,
        first entry is used for new calculations.
        """
        if check:
            sps=[]
            for sp in self.searchpaths:
                if (sp!=None) and (os.path.isdir(sp)) and (not (os.path.abspath(sp) in sps)):
                    sps.append(os.path.abspath(sp))
            return sps
        return self.searchpaths
    
    ### 这个函数是从searchpaths属性中获得不重复的路径
    
    def get_version(self, main_only=False):
        """return present HTE version
        """
        if main_only:
            return self.version
        return self.version+'-'+self.subversion
    
    def get_cif_repository(self):
        """return dirctory for internal storage of cif sources
        """
        if self.cif_repository==None:
            return None
        if os.path.isdir(self.cif_repository):
            return self.cif_repository
        else:
            try:
                os.mkdir(self.cif_repository)
                return self.cif_repository
            except:
                return None
            
    def self_check(self, uid_list='all', symprec=1e-3):
        """Self check for entries in HTE structure database:
        Checks with different methods if symmetry information and composition are consistent
        returns a tuple (warnings,errors) with uids and diagnostics.
        uid_list: list with uids to be checked, default: all database entries
        """
        if uid_list=='all':
            uidlist_in=self.select()
        else:
            uidlist_in=uid_list
        warnings={}
        errors={}
        print '*** Starting self_check for %d data base entries...'%len(uidlist_in)
        nchecks=0
        for uid in uidlist_in:
            nchecks=nchecks+1
            sys.stdout.write('..%d'%nchecks)
            if nchecks % 20==0:
                sys.stdout.flush()
            warn,err=self.structureDB[uid].self_check()
            if err!=None:
                errors[uid]=err
            if warn!=None:
                warnings[uid]=warn
        print '\n** self_check finished:  found ',len(errors),' errors and ',len(warnings),' warnings in ',len(uidlist_in),' data base entries'
        return warnings,errors

    def reset_jobstatus(self, uid_list, calc_schemes=[],nsubmit=0):
        """Reset the number of submitted jobs for unconverged calculations.
        Useful if jobs did not converge because of cluster failures.
        uid_list: list with uids to be checked, default: empty list
        calc_schemes: list of calculation schemes, default: empty list
        todo: include subdirectories
        """
        resets=[]
        for uid in uid_list:
            for cs in calc_schemes:
                calcdir=os.path.join(self.structureDB[uid].calcdir,cs,"")
                if (calcdir in self.structureDB[uid].submitted_jobs):
                    if self.get_energy(uid,cs)==None:
                        nsub=self.structureDB[uid].submitted_jobs[calcdir]['nsubmit']
                        self.structureDB[uid].submitted_jobs[calcdir]['nsubmit']=nsubmit
                        print uid,cs," reset from ",nsub," to ",self.structureDB[uid].submitted_jobs[calcdir]['nsubmit']
                        resets.append((uid,cs))
        return resets
    ## 如果获得能量为None 就把uid cs添加到reset列表里 然后把那个任务的submitted_jobs[calcdir]['nsubmit']设置成0
    def get_cif_source(self,uid):
        """returns the cif source as tuple (filename,cif_uid,cif_index)
        or a tuple (None,None,None) if no source is accessible.
        """
        if uid in self.structureDB:
            return self.structureDB[uid].get_cif_source()
        
    ############################################
    # routines to manage calculational schemes #
    ############################################
    def add_calc_scheme(self, name, calculator_name, **settings):
        if name in self.calc_schemes:
            print 'calc_scheme ',name,' already exists, nothing done!'
        else:
            if (calculator_name.lower()=='vasp') or (calculator_name.lower()=='fplo'):
                self.calc_schemes[name]= calculator_name.lower(), settings
                print 'Added calc_scheme ',name 
            else:
                print 'Calculator',calculator_name,'not yet implemented in HTE, nothing done!'
        return name

    def check_vasp_pp_files(self, element, settings, return_dict=False):
        ppdict={}
        if 'VASP_PP_PATH' in os.environ:
            pppaths = os.environ['VASP_PP_PATH'].split(':')
        else:
            print 'check_vasp_pp_files(): VASP_PP_PATH not given, abort!'
            return False
        if 'xc' in settings:
            if settings['xc'] == 'PW91':
                xc = '_gga/'
            elif settings['xc'] == 'PBE':
                xc = '_pbe/'
            else:
                print 'check_vasp_pp_files(): xc not given, abort!'
                return False
        if 'setups' in settings:
            setups=settings['setups']
        else:
            setups={}
        pp=element
        if element in setups:
            pp=element+setups[element]
        found=False
        for path in pppaths:
            ppfile=os.path.join(path,'potpaw'+xc.upper()+str(pp),'POTCAR')
            if isfile(ppfile) or islink(ppfile) or isfile(ppfile+'.Z') or islink(ppfile+'.Z'):
                found=True
                if return_dict==True:
                    if isfile(ppfile+'.Z') or islink(ppfile+'.Z'):
                        ppdict['filename']=ppfile+'.Z'
                        os.system('gunzip -c %s > hte_tmp.txt'%ppdict['filename'])
                        infile=open('hte_tmp.txt','r')
                    else:
                        ppdict['filename']=ppfile
                        infile=open(ppdict['filename'],'r')
                    line=infile.readline()
                    while line:
                        if (line.strip().startswith('ENMAX')):
                            ppdict['ENMAX']=float(line.split("=")[1].split(";")[0])
                        if (line.strip().startswith('EAUG')):
                            ppdict['EAUG']=float(line.split("=")[1])
                        if (line.strip().startswith('RWIGS')):
                            ppdict['RWIGS']=round(float(line.split(";")[0].split("=")[1].strip().split(";")[0].strip())*0.5291772,3)
                        line=infile.readline()
                    infile.close()
                break
        if return_dict==True:
            return found,ppdict
        return found

    
    def update_calc_scheme(self, name, setups=None, init_structure=None, lorbit=None):
        if name in self.calc_schemes:
            calc_name,settings=self.calc_schemes[name]
            if ('init_structure' in settings) and (init_structure!=None):
                print "old:",settings['init_structure'],"new:",init_structure
                settings['init_structure']=init_structure
            if (calc_name=='vasp') and (lorbit!=None):
                #lorbit tag should have no influence results
                settings['lorbit']=lorbit
            if (calc_name=='vasp') and (setups!=None):
                #allow change of setups if PP doesn't exist for the element
                if 'VASP_PP_PATH' in os.environ:
                    pppaths = os.environ['VASP_PP_PATH'].split(':')
                else:
                    print 'update_calc_scheme(): VASP_PP_PATH not given, abort!'
                    return False
                if 'xc' in settings:
                    if settings['xc'] == 'PW91':
                        xc = '_gga/'
                    elif settings['xc'] == 'PBE':
                        xc = '_pbe/'
                    else:
                        print 'update_calc_scheme(): xc not given, abort!'
                        return False
                if 'setups' in settings:
                    setups_old=settings['setups']
                else:
                    setups_old={}
                setups_new=setups_old
                print 'old setups:',setups_old
                for el in setups:
                    if el in setups_old:
                        oldpp=setups_old[el]
                    else:
                        oldpp=el
                    found=False
                    for path in pppaths:
                        ppfile=os.path.join(path,'potpaw'+xc.upper()+str(oldpp),'POTCAR')
                        if isfile(ppfile) or islink(ppfile) or isfile(ppfile+'.Z') or islink(ppfile+'.Z'):
                            found=True
                            break
                    if found:
                        print 'PP file for old setup ',oldpp,' exists, keeping old setup!'
                    else:
                        found=False
                        for path in pppaths:
                            ppfile=os.path.join(path,'potpaw'+xc.upper()+str(el)+str(setups[el]),'POTCAR')
                            if isfile(ppfile) or islink(ppfile) or isfile(ppfile+'.Z') or islink(ppfile+'.Z'):
                                found=True
                                break
                        if found:
                            setups_new[el]=setups[el]
                settings['setups']=setups_new
                print 'new setups:',settings['setups']
                return True
            return False
        
        
        #### 这个函数看上去只升级了下赝势 同时也加入了赝势还有

    def get_calc_schemes(self,substring=None, calc_names=[]):
        """returns dictionary with calc_schemes defined for present HTE object
        """
        cs_list=[]
        for cs in self.calc_schemes:
            add=True
            calculator, settings=self.calc_schemes[cs]
            if (substring!=None) and (not(substring in cs)):
                add=False
            elif (calc_names!=[]) and (not(calculator in calc_names)):
                add=False
            if add:
                cs_list.append(cs)
        return cs_list
    
    def print_calc_schemes(self,substring=None):
        for i in self.calc_schemes:
            if (substring==None) or (substring in i):
                calculator, settings=self.calc_schemes[i]
                print i,': Calculator=', calculator
                print '   Settings=', settings

    def get_job_commands(self, calc_scheme):
        """Returns dictionary with commands for job submission for the given calc_scheme (or calculator name)
        'job_file': file(s) to be copied before job submission (e.g. PBS/SLURM script); relative path is changed to absulute path
        'job_command': shell command to submit job
        'job_environment': SLURM/PBS
        """
        
        job_commands={}
        if self.use_job_array==True:
            print "check_point50, use_job_array is True"
            try:
                if calc_scheme in self.job_commands['array_job']:
                    job_commands=self.job_commands['array_job'][calc_scheme]
                elif calc_scheme in self.calc_schemes:
                    calc_name, settings=self.calc_schemes[calc_scheme]
                    if calc_name.lower() in self.job_commands['array_job']:
                        job_commands=self.job_commands['array_job'][calc_name.lower()]
                if ('job_file' in job_commands) and (not (job_commands['job_file'].startswith("/"))):
                    jobfile=os.path.join(self.main_directory,job_commands['job_file'])
                    job_commands['job_file']=jobfile
            except:
                self.add_logmessage("WARNING(get_job_commands): Failed to get job commands for %s"%str(calc_scheme))
        if job_commands!={}:
            return job_commands
        try:
            if calc_scheme in self.job_commands:
                job_commands=self.job_commands[calc_scheme]
            elif calc_scheme in self.calc_schemes:
                calc_name, settings=self.calc_schemes[calc_scheme]
                if calc_name.lower() in self.job_commands:
                    job_commands=self.job_commands[calc_name.lower()]
            if ('job_file' in job_commands) and (not (job_commands['job_file'].startswith("/"))):
                jobfile=os.path.join(self.main_directory,job_commands['job_file'])
                job_commands['job_file']=jobfile
        except:
            self.add_logmessage("WARNING(get_job_commands): Failed to get job commands for %s"%str(calc_scheme))
        print "check_point52, final job_commands is :", job_commands
        return job_commands
    
    ### use_job_array==True的话就是嵌套词典了

    def set_job_commands(self, calc_scheme, **kwargs):
        """Sets the job commands for the given calc_scheme (or calculator name).
        (calculator names in lower case)
        """
        job_commands={}
        if calc_scheme in self.job_commands:
            job_commands=self.job_commands[calc_scheme]
        for args in kwargs:
            job_commands[args]=kwargs[args]
        self.job_commands[calc_scheme]=job_commands
        return job_commands
    
    ############################################
    # routines to manage structure db          #
    ############################################
    def add_structDB_entry(self, uid, structureDB_entry=None, atoms_obj=None, comment=None, cif_file=None, cif_index=None, cif_uid=None, exclude_similar=True, silent=True, _source_info=None):
        """ add_structDB_entry(self, uid, structureDB_entry=None, atoms_obj=None, comment=None, cif_file=None, cif_index=None, cif_uid=None, exclude_similar=True, silent=False)
        add a new entry to HTE structureDB (unless uid or similar structure is already in DB)
        """
        do_not_add=False
        similar=None
        if uid in self.structureDB:
            do_not_add=True
        elif (structureDB_entry==None) and (atoms_obj==None):
            do_not_add=True
        elif exclude_similar==True:
            if structureDB_entry!=None:
                similar=self.check_similar_entries(structureDB_entry) #.atoms_initial) 
            elif atoms_obj!=None:
                similar=self.check_similar_entries(atoms_obj)
            if similar!=None:
                do_not_add=True
        if do_not_add==True:
            if (silent==False):
                print ' %s not added (because similar structure already in DB)' % uid
            return False
        elif structureDB_entry!=None:
            self.structureDB[uid]=structureDB_entry
        else:
            if (_source_info==None):
                source_info=[('non_HTE',{'comment':comment})]
            else:
                source_info=_source_info
            self.structureDB[uid]=HTEdbentry(calcdir=uid, atoms_obj=atoms_obj, _source_info=source_info)
        if (silent==False):
            print ' %s added to DB' % uid
        return True
####使用 HTEdbentry这个函数添加新的计算条目 主要内容是计算路径uid 原子对象
    def print_source_info(self,uidlist):
        """print source information of structure database entry (e.g. cif source)
        """
        if isinstance(uidlist,list):
            uid_list=uidlist
        else:
            uid_list=[uidlist]
        for uid in uid_list:
            print "*** ",uid," ***"
            if uid in self.structureDB:
                for i in range(len(self.structureDB[uid].source_info)):
                    method,args=self.structureDB[uid].source_info[i]
                    print " %d) %-25s: %s"%(i, method,args)
            else:
                print "Warning(print_source_info):",uid," not found!"
            
    def output(self, uid_list, *args, **kwargs):
        """ output(self, uid_list, *args, **kwargs)
        args: list of output order
          'energy','energy_per_atom','formation_energy','formation_energy_difference',
          'jobstatus','uid','cputime',comment,'cell_composition','elnum_max',
          'spacegroup','spacegroup_no','spacegroup_symbol','structure_type',
          'atfraction_','cell_volume', 'wyckoff_information', 'sg_international'
         kwargs: additional options, e.g. to set reference structure for formation energy
           calc_scheme='paw-k40-ecut350-cont-us-k30-e250-opt',file="pip.dat"
          'sep1'],'sep1':'&'
        """
        # set some defaults and if they are redefined in kwargs:
        sort_output=True
        if 'sort_output' in kwargs:
            sort_output=kwargs['sort_output']
        debug=False
        if 'debug' in kwargs:
            debug=kwargs['debug']
        format_strings={'lattice_parameter_a':"%.2f",'lattice_parameter_b':"%.2f",'lattice_parameter_c':"%.2f",'energy':"%.3f",'energy_per_atom':"%.3f",'formation_energy':"%.3f",'dE_FM_per_atom':"%.3f",'MAE_BAND':"%.3f",'distance_convex_hull':"%.3f"
                        ,'magnetic_moment':"%.1f",'magnetic_moments':"%.1f",'magnetic_moment_per_atom':"%.1f",'cell_volume':"%.1f"}
        if 'format_strings' in kwargs:
            for argument in kwargs['format_strings']:
                format_strings[argument]=kwargs['format_strings'][argument]
        sort_chem_formula='alpha'
        if 'sort_chem_formula' in kwargs:
            sort_chem_formula=kwargs['sort_chem_formula']
        calc_scheme=None
        if 'calc_scheme' in kwargs:
            calc_scheme=kwargs['calc_scheme']
        magsettings={'submitted':True}
        if 'magsettings' in kwargs:
            magsettings=kwargs['magsettings']
        if uid_list=='all':
            uid_list=self.select(all=True)
        update=False
        if 'update' in kwargs:
            update=kwargs['update']
        comment='!'
        if 'comment' in kwargs:
            comment=kwargs['comment']
        nsub_max=self.nsub_max
        if 'nsub_max' in kwargs:
            nsub_max=kwargs['nsub_max']
        header=comment
        for argument in args:
            header=header+argument+' '
        if 'file' in kwargs:
            outfile=open(kwargs['file'],'w')
            outfile.write(header+'\n')
        else:
            outfile=None
            print header
        latex=False
        if 'latex' in kwargs:
            latex=kwargs['latex']
        #
        if ('energy' in args) or ('energy_per_atom' in args) or ('formation_energy' in args) or ('jobstatus' in args):
            if 'calc_scheme' in kwargs:
                if (not (kwargs['calc_scheme'] in self.calc_schemes)) and (not (isinstance(calc_scheme,tuple))): #check this (and if still needed)
                    print 'hte.output: calc_scheme ',kwargs['calc_scheme'],' undefined, abort'
                    return
            else:
                print 'hte.output: No calc_scheme defined, abort'
                return
        if 'formation_energy_difference' in args:
            statistics_E_form_diff={}
            for cs in kwargs['calc_schemes']:
                statistics_E_form_diff[cs]=(0.0,0,0.0)
        ######################################
        #formation enthalpies require reference structure for elemental compound:
        if 'formation_enthalpy' in args:
            #get energies for reference structures
            H_ref={}
            for uid in uid_list:
                comp={}
                comp=self.structureDB[uid].get_composition(reduce=False)
                for element in comp:
                    if not (element in H_ref):
                        H_ref[element]=None
            #determine reference energies automatically if no reference is given:
            for element in H_ref:
                #probe if reference structure was given (better way)
                #keyw='reference_'+element
                #if keyw in kwargs:
                #    uid=kwargs[keyw]
                #    H_ref[element]=self.get_energy_per_atom(uid,calc_scheme,update=update)
                #else:
                H_ref[element],uid_ref=self.get_reference_enthalpy_per_atom(element,calc_scheme=kwargs['calc_scheme'],update=update,return_uid=True)
                if uid_ref==None:
                    continue
                if outfile!=None:
                    outfile.write("! %s: reference_uid=%s H_ref=%.3f\n"%(element,uid_ref,H_ref[element]))
                else:
                    print "! %s: reference_uid=%s H_ref=%.3f"%(element,uid_ref,H_ref[element])
        #formation energies require reference structure for elemental compound:
        E_ref={}
        if False: #'formation_energy' in args:
            #get energies for reference structures
            E_ref={}
            for uid in uid_list:
                comp={}
                comp=self.structureDB[uid].get_composition(reduce=False)
                for element in comp:
                    if not (element in E_ref):
                        E_ref[element]=None
            #determine reference energies automatically if no reference is given:
            for element in E_ref:
                #probe if reference structure was given (better way)
                keyw='reference_'+element
                if keyw in kwargs:
                    uid=kwargs[keyw]
                    E_ref[element]=self.get_energy_per_atom(uid,calc_scheme,update=update,nsub_max=nsub_max)
                else:
                    E_ref[element]=self.get_reference_energy_per_atom(element,calc_scheme=kwargs['calc_scheme'],update=update,nsub_max=nsub_max)
        ############################################
        outstr=[] #header+"\n"
        subdirectories=[]
        for uid in uid_list:
            if ('all' in magsettings) and (magsettings['all']==True):
                magconfigs=self.setup_magnetic_structures(uid,calc_scheme, magsettings=magsettings)
                for magconf in magconfigs:
                    subdir=os.path.join(calc_scheme,magconf)
                    subdirectories.append((uid,{subdir:magconfigs[magconf]}))
            else:
                subdirectories.append((uid,{}))

        print "check_point106 , check subdirectories in output", subdirectories
        for uid,subdir in subdirectories:
            status=None
            if debug==True:
                print "DEBUG(output): uid=",uid,subdir
            print "check_point160, uid and subdir are:",uid,subdir
            line=''
            print "check_point111, before converged_ao"
            converged_ao=self.get_atoms_object(uid,calc_scheme=calc_scheme,magsettings=magsettings, sub_directories=subdir)
            print "cooonvergeed ao :", converged_ao
            commentline=False
            for argument in args:
                if argument=='uid':
                    line=line+uid+' '
                elif (argument.startswith('sep')) and (argument in  kwargs):
                    line=line+kwargs[argument]+' '
                elif argument=='cell_composition':
                    line=line+self.structureDB[uid].get_composition(reduce=True,latex=latex, sort=sort_chem_formula)+' '
                elif argument=='setups':
                    comp=self.structureDB[uid].get_composition(reduce=False)
                    calc_name,settings=self.calc_schemes[kwargs['calc_scheme']]
                    for el in comp.keys():
                        line=line+el
                        if ('setups' in settings) and (el in settings['setups']):
                            if latex==True:
                                line=line+settings['setups'][el]
                            else:
                                line=line+settings['setups'][el]
                elif argument=='cell_composition_parent':
                    uidp=self.get_parent_uid(uid)
                    line=line+self.structureDB[uidp].get_composition(reduce=True,latex=latex, sort=sort_chem_formula)+' '
                elif argument=='interstitial_atom':
                    line=line+self.get_interstitial_atom(uid)+' '
                elif argument=='interstitial_position':
                    line=line+self.get_interstitial_position(uid)+' '
                elif argument=='interstitial_concentration':
                    iat=self.get_interstitial_atom(uid)
                    if iat!="--":
                        conc=self.structureDB[uid].get_atfraction(iat)
                    else:
                        conc=0.
                    line=line+"%.1f"%(conc*100)+' '
                elif argument=='cell_composition_cif_org':
                    line=line+self.structureDB[uid].get_composition(reduce=True,methods=['cif_org'])+' '
                elif argument=='sg_international':
                    line=line+self.get_sg_international(uid)+' '
                elif argument=='data_path':
                    pd=self.get_properties(uid,kwargs['calc_scheme'], magsettings=magsettings, sub_directories=subdir)
                    if 'data_path' in pd:
                        line=line+pd['data_path']+' '
                    elif subdir!={}:
                        line=line+' %s '%str(subdir.keys())
                    else:
                        line=line+'None '
                elif argument=='wyckoff_information':
                    try:
                        line=line+str(self.structureDB[uid].get_wyckoff(sort=True))+' '
                    except:
                        line=line+'None'+' '
                elif argument=='elnum_max':
                    elnum=0
                    for Z in self.structureDB[uid].atoms_initial.arrays['numbers']:
                        if Z>elnum:
                            elnum=Z
                    line=line+str(elnum)+' '
                elif argument=='spacegroup_ini':
                    line=line+str(self.structureDB[uid].atoms_initial.info['spacegroup'].no)+'['+self.structureDB[uid].atoms_initial.info['spacegroup'].symbol+'] '
                elif argument=='spacegroup_no_ini':
                    line=line+str(self.get_spacegroup_number(uid))+' '
                elif argument=='spacegroup_no':
                    if 'calc_scheme' in kwargs:
                        line=line+str(self.get_spacegroup_number(uid,calc_scheme=kwargs['calc_scheme'], magsettings=magsettings))+' '
                    else:
                        line=line+str(self.get_spacegroup_number(uid))+' '
                elif argument=='spacegroup_symbol':
                    line=line+str(self.structureDB[uid].atoms_initial.info['spacegroup'].symbol)+' '
                elif argument=='structure_type':
                    strtype=self.get_structure_type(uid)
                    if (latex==True) and (strtype==None):
                        line=line+"-- "
                    else:
                        line=line+str(self.get_structure_type(uid))+' '
                elif argument=='gle_label':
                    line=line+"\""+self.structureDB[uid].get_composition(reduce=True, gle=True, sort=sort_chem_formula)+"("
                    strtype=self.get_structure_type(uid)
                    if (strtype==None):
                        line=line+"-,-,"
                        if 'calc_scheme' in kwargs:
                            line=line+str(self.get_spacegroup_number(uid,calc_scheme=kwargs['calc_scheme'], magsettings=magsettings))
                        else:
                            line=line+str(self.get_spacegroup_number(uid))
                    else:
                        line=line+str(self.get_structure_type(uid))
                    line=line+")\" "
                elif argument=='lattice_parameter_a':
                    if converged_ao!=None:
                        cell=converged_ao.get_cell()
                        line=line+"%5.2f "%norm(cell[0])
                elif argument=='lattice_parameter_a_ini':
                    line=line+"%5.2f "%self.get_lattice_parameters(uid)[0]
                elif argument=='lattice_parameter_b':
                    if converged_ao!=None:
                        cell=converged_ao.get_cell()
                        line=line+"%5.2f "%norm(cell[1])
                elif argument=='lattice_parameter_c':
                    if converged_ao!=None:
                        cell=converged_ao.get_cell()
                        line=line+"%5.2f "%norm(cell[2])
                #
                elif argument.startswith('atfraction_'):
                    keyw, element=argument.split('_')
                    line=line+"%5.3f " % self.structureDB[uid].get_atfraction(element)
                elif argument=='energy':
                    energy=self.get_energy(uid,kwargs['calc_scheme'], update=update, nsub_max=nsub_max, magsettings=magsettings, sub_directories=subdir)
                    line=line+value2string(energy,latex=latex, format_string=format_strings['energy'])
                    if energy==None:
                        commentline=True
                elif argument=='energy_per_atom':
                    energy=self.get_energy_per_atom(uid,kwargs['calc_scheme'], update=update, nsub_max=nsub_max, magsettings=magsettings, sub_directories=subdir)
                    line=line+value2string(energy,latex=latex, format_string=format_strings[argument])
                    if energy==None:
                        commentline=True
                elif argument=='dE_FM_per_atom':
                    E_FM=self.get_energy_per_atom(uid,kwargs['calc_scheme'], update=update, nsub_max=nsub_max, magsettings={'default_aoini*':True})
                    energy=self.get_energy_per_atom(uid,kwargs['calc_scheme'], update=update, nsub_max=nsub_max, magsettings=magsettings, sub_directories=subdir)
                    if (E_FM!=None) and (energy!=None):
                        energy=energy-E_FM
                    else:
                        energy=None
                    line=line+value2string(energy,latex=latex, format_string=format_strings[argument])
                    if energy==None:
                        commentline=True
                elif argument.startswith('MAE_BAND'):
                    argx=argument
                    unitmae='meV/cell'
                    qaxis=[[0,0,1],[1,0,0]]
                    kmesh=[]
                    if '_per_atom' in argx:
                        unitmae='meV/atom'
                        argx=argx.replace('_per_atom','')
                    if '_kp' in argx:
                        kspl=argx.split('_kp')[1].split('_')[0].split(',')
                        kmesh=[int(kspl[0]),int(kspl[1]),int(kspl[2])]
                        rplc=argx.split('_kp')[1].split('_')[0]
                        argx=argx.replace('_kp'+rplc,'')
                    qspl=argx.split('_BAND')[1].split('_')
                    if len(qspl)>2:
                        qaxis=[]
                        for i in range(1,len(qspl)):
                            qaxs=qspl[i].split(',')
                            qax=[]
                            for j in range(len(qaxs)):
                                qax.append(float(qaxs[j]))
                            qaxis.append(qax)
                    energy=self.get_band_MAE(uid,kwargs['calc_scheme'], update=update, qaxis=qaxis, k_mesh=kmesh, nsub_max=nsub_max,magsettings=magsettings, unit=unitmae)
                    line=line+value2string(energy,latex=latex, format_string=format_strings['MAE_BAND'])
                    if energy==None:
                        commentline=True
                #obsolete, remove after test of upper part
                elif argument=='MAE_BAND':
                    energy=self.get_band_MAE(uid,kwargs['calc_scheme'], update=update, nsub_max=nsub_max,magsettings=magsettings)
                    line=line+value2string(energy,latex=latex, format_string=format_strings[argument])
                    if energy==None:
                        commentline=True
                elif argument.startswith('MAE_BAND_'):
                    kmesh=[]
                    if argument.startswith('MAE_BAND_kp'):
                        kspl=argument.split('_kp')[1].split('_')[0].split(',')
                        kmesh=[int(kspl[0]),int(kspl[1]),int(kspl[2])]
                        qspl=argument.split('_kp')[1].split('_')
                    else:
                        qspl=argument.split('_BAND')[1].split('_')
                    if len(qspl)<3:
                        qaxis=[[0,0,1],[1,0,0]]
                    else:
                        qaxis=[]
                        for i in range(1,len(qspl)):
                            qaxs=qspl[i].split(',')
                            qax=[]
                            for j in range(len(qaxs)):
                                qax.append(float(qaxs[j]))
                            qaxis.append(qax)
                    if kmesh==[]:
                        energy=self.get_band_MAE(uid,kwargs['calc_scheme'], update=update, qaxis=qaxis, nsub_max=nsub_max,magsettings=magsettings)
                    else:
                        energy=self.get_band_MAE(uid,kwargs['calc_scheme'], update=update, qaxis=qaxis, k_mesh=kmesh, nsub_max=nsub_max,magsettings=magsettings)
                    line=line+value2string(energy,latex=latex, format_string=format_strings['MAE_BAND'])
                    if energy==None:
                        commentline=True
                elif argument=='distance_convex_hull':
                    energy=self.get_distance_from_convex_hull(uid,kwargs['calc_scheme'], magsettings=magsettings, sub_directories=subdir)
                    line=line+value2string(energy,latex=latex, format_string=format_strings[argument])
                    if energy==None:
                        commentline=True
                elif argument=='magnetic_moment':
                    val=self.get_magnetic_moment(uid,kwargs['calc_scheme'],magsettings=magsettings, sub_directories=subdir)
                    line=line+value2string(val,latex=latex, format_string=format_strings[argument])
                elif argument=='magnetic_moment_per_atom':
                    val=self.get_magnetic_moment_per_atom(uid,kwargs['calc_scheme'],magsettings=magsettings, sub_directories=subdir)
                    line=line+value2string(val,latex=latex, format_string=format_strings[argument])
                elif argument=='magnetic_moments':
                    pd=self.get_properties(uid,kwargs['calc_scheme'], magsettings=magsettings, sub_directories=subdir)
                    print "check_point163, pd is", pd
                    moms = []
                    if 'magnetic_moments' in pd:
                        line=line+pd['magnetic_moments']+' '
                    elif subdir!={}:
                        unique_key = next(iter(pd))
                        nested_pd = pd[unique_key] 
                        print "check_point162, subdir!={}! and is", subdir
                        print "check_point164, nested_pd is:  ", nested_pd
                        # for each_subdir in subdir.keys():
                        if ('magnetic_moments' in nested_pd) and ('chemical_symbols' in nested_pd):
                            for i in range(len(nested_pd['chemical_symbols'])):
                                moms.append((nested_pd['chemical_symbols'][i],nested_pd['magnetic_moments'][i]))
                        # line=line+' %s '%str(subdir.keys())
                        if 'energy' in nested_pd:
                            energy = nested_pd['energy']
                        print "check_point161, moms is ", moms
                        print "unique_key is: ",unique_key
                        line = line + str(energy) + ' ' + str(moms)
                        if subdir != None:
                            # for key, value in subdir.items():
                            #     print  "check_point180, ", key + ": " + str(value)
                            # # for sub in subdir:
                            # #     if subdir
                            if unique_key in subdir and subdir[unique_key] is not None and bool(subdir[unique_key]):
                                print "check_point179, atom_object is, ", subdir[unique_key]
                        
                            # if subdir[unique_key] != None:
                                
                        # #### write mcif files by hao
                        file_tmp = str(unique_key+'.mcif')
                        print "file_tmp is ",file_tmp
                        # cell_tmp=nested_pd['cell']
                        # a_tmp = norm(cell_tmp[0])
                        # b_tmp = norm(cell_tmp[1])
                        # c_tmp = norm(cell_tmp[2])
                        # alpha_tmp = arccos(dot(cell_tmp[1], cell_tmp[2])/(b_tmp*c_tmp))*180./pi
                        # beta_tmp  = arccos(dot(cell_tmp[0], cell_tmp[2])/(a_tmp*c_tmp))*180./pi
                        # gamma_tmp = arccos(dot(cell_tmp[0], cell_tmp[1])/(a_tmp*b_tmp))*180./pi
                        # lines_tmp=['data_%s'%uid,'_cell_angle_alpha              %5.2f'%alpha_tmp,'_cell_angle_beta               %5.2f'%beta_tmp,'_cell_angle_gamma              %5.2f'%gamma_tmp,'loop_','_space_group_symop_magn_operation.id','_space_group_symop_magn_operation.xyz','1 x,y,z,+1','']
                        # abc_tmp=['a','b','c']
                        # for i in range(3):
                        #     lines_tmp.append("_cell_length_%s\t %.3f"%(abc_tmp[i],norm(nested_pd['cell'][i])))
                        # lines_tmp=lines_tmp+['loop_','_atom_site_label','_atom_site_type_symbol','_atom_site_fract_x','_atom_site_fract_y','_atom_site_fract_z']
                        # for i in range(len(nested_pd['chemical_symbols'])):
                        #     lines_tmp.append("%s%d %s %.8f %.8f %.8f"%(nested_pd['chemical_symbols'][i],i+1,nested_pd['chemical_symbols'][i],nested_pd['scaled_positions'][i][0],nested_pd['scaled_positions'][i][1],nested_pd['scaled_positions'][i][2]))
                        # lines_tmp=lines_tmp+['','loop_','_atom_site_moment.label','_atom_site_moment.crystalaxis_x','_atom_site_moment.crystalaxis_y','_atom_site_moment.crystalaxis_z']
                        
                        # magnetic_moments_tmp=nested_pd['magnetic_moments']
                        
                        # for i in range(len(nested_pd['chemical_symbols'])):
                        #     if isinstance(magnetic_moments_tmp[i],float):
                        #         moms_tmp=np.dot(magnetic_moments_tmp[i],qaxis/norm(qaxis))
                        #     else:
                        #         moms_tmp=magnetic_moments_tmp[i]

                        #     lines_tmp.append("%s%d %.4f %.4f %.4f"%(nested_pd['chemical_symbols'][i],i+1,moms_tmp[0],moms_tmp[1],moms_tmp[2]))
                        # outfile_tmp=open(file_tmp,"w")
                        # for line in lines_tmp:
                        #     outfile_tmp.write("%s\n"%line)
                        # outfile_tmp.close()        
                    else:
                        line=line+'None '
                    
                    
                    # print "check_point151, entering output magnetic_moments"
                    # print "check_point159, print subdir", subdir
                    # vallist=self.get_magnetic_moments(uid,kwargs['calc_scheme'],magsettings=magsettings, sub_directories=subdir)
                    # print "check_point152,vallist is: ",vallist
                    # for (valel,valmom) in vallist:
                    #     line=line+valel+" ("
                    #     if (isinstance(valmom,np.ndarray)) or (isinstance(valmom,list)):
                    #         for vmom in valmom:
                    #             line=line+value2string(vmom,latex=latex, format_string=format_strings[argument])
                    #     else:
                    #         line=line+value2string(valmom,latex=latex, format_string=format_strings[argument])
                    #     line=line+") "
                elif argument=='magnetic_type':
                    if latex==True:
                        line=line.strip()
                    line=line+str(self.get_magnetic_structure(uid,kwargs['calc_scheme'],latex=latex,magsettings=magsettings, sub_directories=subdir))+' '
                elif argument=='jobstatus_magnetic_structures':
                    magsettings_jobst=deepcopy(magsettings)
                    magsettings_jobst['return_all']=True
                    prop_dict=self.get_properties(uid, calc_scheme, magsettings=magsettings_jobst) #{'submitted':True,'return_all':True})
                    nmag=len(prop_dict)
                    nconverged=0
                    for pd in prop_dict:
                        if ('energy' in prop_dict[pd]) and (prop_dict[pd]['energy']!=None):
                            nconverged=nconverged+1
                    line=line+" %d/%d "%(nmag,nconverged)
                    if latex==True:
                        line=line.strip()
                    line=line+str(self.get_magnetic_structure(uid,kwargs['calc_scheme'],latex=latex,magsettings=magsettings, sub_directories=subdir))+' '
                elif argument.startswith('bandgap'):
                    if argument=='bandgap':
                        energy=self.get_bandgap(uid,kwargs['calc_scheme'], update=update)
                    else:
                        cs=argument.split('_')[1]
                        energy=self.get_bandgap(uid,cs, update=update)
                    if energy!=None:
                        fstr="%4.1f "
                        if 'fstr_bandgap' in kwargs:
                            fstr=kwargs['fstr_bandgap']
                        line=line+fstr%energy
                    else:
                        line=line+str(energy)+' '
                    #if energy==None:
                    #    commentline=True
                elif argument.startswith('DOS_at_E_F'):
                    if argument=='DOS_at_E_F':
                        energy=self.get_DOS_at_E_F(uid,kwargs['calc_scheme'], update=update)
                    else:
                        cs=argument.split('_')[1]
                        energy=self.get_DOS_at_E_F(uid,cs, update=update)
                    if energy!=None:
                        fstr="%4.1f "
                        if 'fstr_bandgap' in kwargs:
                            fstr=kwargs['fstr_bandgap']
                        line=line+fstr%energy
                    else:
                        line=line+str(energy)+' '
                #
                elif (argument.startswith('formation_energy')) and (argument!='formation_energy_difference'):
                    if argument=='formation_energy':
                        cs=kwargs['calc_scheme']
                        reference_energy=E_ref
                    else:
                        cs=argument.split('_')[2]
                        reference_energy={}
                    val=self.get_formation_energy(uid,cs, reference_energy=reference_energy, update=update, nsub_max=nsub_max, magsettings=magsettings, sub_directories=subdir)
                    line=line+value2string(val,latex=latex, format_string=format_strings[argument])
                    if val==None:
                        commentline=True
                #
                elif argument.startswith('formation_enthalpy'):
                    if argument=='formation_enthalpy':
                        cs=kwargs['calc_scheme']
                    else:
                        cs=argument.split('_')[2]
                    E_form=self.get_formation_enthalpy(uid,cs, update=update, reference_enthalpy=H_ref)
                    if (E_form!=None) and ('fstr_formation_energy' in kwargs):
                        line=line+kwargs['fstr_formation_energy']%E_form
                    else:
                        line=line+str(E_form)+' '
                    if E_form==None:
                        commentline=True
                #
                elif (argument=="elastic constants") or (argument=="elastic_constants"):
                    C=self.get_elastic_constants(uid, kwargs['calc_scheme'], update=update, magsettings=magsettings, sub_directories=subdir)
                    line=line+str(C)
                #
                elif argument=='formation_energy_difference':
                    E_form=self.get_formation_energy(uid, kwargs['calc_scheme'], update=update)
                    for calc_scheme in kwargs['calc_schemes']:
                        E_form_diff=self.get_formation_energy(uid, calc_scheme, update=update)
                        if (E_form!=None) and (E_form_diff!=None):
                            E_form_diff=E_form_diff-E_form
                            line=line+"%.3f "%abs(E_form_diff)
                        else:
                            E_form_diff=None
                            line=line+str(E_form_diff)+' '
                        if E_form_diff==None:
                            commentline=True
                        else:
                            sumefd,nstat,maxefd=statistics_E_form_diff[calc_scheme]
                            absefd=sqrt(E_form_diff*E_form_diff)
                            nstat=nstat+1
                            sumefd=sumefd+absefd
                            if absefd>maxefd:
                                maxefd=absefd
                            statistics_E_form_diff[calc_scheme]=sumefd,nstat,maxefd
                #
                elif argument=='jobstatus':
                    print "check_point47, the output have jobstatus"
                    if status==None:
                        print "check_point48, jobstatus is none"
                        job_commands=self.get_job_commands(kwargs['calc_scheme'])
                        if 'job_environment' in job_commands:                            
                            job_environment=job_commands['job_environment']
                            print "check_point53, job_environment is: ",job_environment
                        else:
                            job_environment=None
                        status=self.structureDB[uid].get_jobstatus( kwargs['calc_scheme'], None,job_environment=job_environment)
                        print "check_point63, kwargs['calc_scheme'] is: ", kwargs['calc_scheme']
                        print "check_point54, status is: ", status
                    line=line+str(status['qstat'])+'('+str(status['nsubmit'])
                    if status['converged']:
                        print "check_point104, converged",status['converged']
                        line=line+', converged111) '
                    else:
                        line=line+', not converged1) '
                        print "check_point105, not converged",status['converged']
                elif argument=='cputime':
                    if status==None:
                        status=self.structureDB[uid].get_jobstatus( kwargs['calc_scheme'], calc)
                    line=line+str(status['cputime'])+' '
                elif argument=='comment':
                    line=line+comment
                elif argument=='cell_volume':
                    if converged_ao!=None:
                        line=line+"%8.3f " %converged_ao.get_volume()
                else:
                    line=line+" %s "%argument
            if commentline==True:
                line=comment+line
            if (outfile==None):
                print line
            outstr.append(line)
        if (sort_output==True):
            outstr=sorted(outstr)
        if 'formation_energy_difference' in args:
            for cs in kwargs['calc_schemes']:
                sumefd,nstat,maxefd=statistics_E_form_diff[cs]
                if nstat>2:
                    line="%s statistics for formation_energy_difference to %s: abssum=%8.3f N=%d max=%5.3f mean=%5.3f"%(comment,cs,sumefd,nstat,maxefd,sumefd/nstat)
                    if (outfile==None):
                        print line
                    else:
                        outstr.append(line)
        if (outfile!=None):
            for line in outstr:
                outfile.write(line+'\n')
            outfile.close()
        return outstr


    def check_similar_entries(self, dbstructure):
        similar=[]
        for dbstruct in self.structureDB:
            if self.struct_diff(self.structureDB[dbstruct],dbstructure)==False:
                similar.append(dbstruct)
        if len(similar)==0:
            similar=None
        else:
            print similar
        return similar
    
    ###看着像是使用struct_diff把相似结构的uid存到similar列表里面    dbstructure是要比较的结构

    def import_mcif(self,cifnames, primitive_cell=True, exclude_similar=True, exclude_fractional=True, silent=True, check=True, filename_as_uid=False, methods=['ase','cif2struct'],modify_cif=False,dry_run=False, symprec=1e-3):
        """ import_cif(cifnames): import structures from cif files in HTE structure database
        cifnames: shell expression like '*.cif'
        check: True=self_check directly after importing, 'all'=self_check at the end for all imported cif-structures         at once, False=no self_check.
        """
        # import structures from cif file(s)
        print '* Importing structure data from ',cifnames
        cif_list=glob.glob(cifnames)
        print "check_point110, print cif_list",cif_list
        import_statistics={'imported':[],'already_in':[],'excluded':[],'failed':[]}
        nentry_total=0
        index=0 #for the moment: only single mcifs
        for ciffile in cif_list:
            nentry_total=nentry_total+1
            try:
                print "check_point109, entering read_mcif"
                uid,pd,ciftags=MSG().read_full_mcif(ciffile)
                mcif_short=ciffile.split('/')[-1].split('.mcif')[0]
                print "check_point108, read_mcif", mcif_short
                if (filename_as_uid==True):
                    uid=mcif_short
                ao=Atoms(symbols=pd['chemical_symbols'],scaled_positions=pd['scaled_positions'],cell=pd['cell'],pbc=True)
                print ao,spglib.get_symmetry_dataset(ao, symprec=symprec)['number']
                lattice, scaled_positions, numbers=spglib.find_primitive(ao, symprec=symprec)
                aop=Atoms(numbers=numbers,cell=lattice,scaled_positions=scaled_positions,pbc=True)
                print "primitive:",aop,spglib.get_symmetry_dataset(aop, symprec=symprec)['number']
                if  (uid in self.structureDB) and (dry_run==False):
                    import_statistics['already_in'].append((ciffile,uid))
                    if silent==False:
                        print "already in database:",uid
                else:
                    cif_sources={'cif_file_external':os.path.abspath(ciffile),
                                 'cif_index_external':index,
                                 'cif_id':uid}
                    source_info=[('import_mcif',cif_sources)]
                    if (primitive_cell==True):
                        dbentry=HTEdbentry(aop,uid,_source_info=source_info)
                    else:
                        dbentry=HTEdbentry(ao,uid,_source_info=source_info)
                    if not ('mcif_structure' in dbentry.cif_info):
                        dbentry.cif_info['mcif_structure']={}
                    dbentry.cif_info['mcif_structure']['mcif_%s'%mcif_short]={'atoms_object':ao,'magmom':pd['initial_magnetic_moments'],'lnoncollinear':True} #TODO: reduce to coll
                    #print "xx",dbentry.cif_info
                    if self.add_structDB_entry(uid, structureDB_entry=dbentry,_source_info=source_info, exclude_similar=exclude_similar, silent=silent)==True:
                        import_statistics['imported'].append((ciffile,uid))
                        if self.get_cif_repository()!=None:
                            shutil.copy(ciffile,os.path.join(self.get_cif_repository(),uid+'.mcif'))
            except:
                import_statistics['failed'].append((ciffile,uid))
        print " Imported %d structures from %d [%d already in db/%d failed/%d excluded]"%(
            len(import_statistics['imported']),nentry_total,len(import_statistics['already_in']),
            len(import_statistics['failed']),len(import_statistics['excluded']))
        if (check=='all'):
            import pprint
            pprint.pprint(self.self_check())
        if dry_run==True:
            return None,None
        return import_statistics


    def import_cif(self,cifnames, primitive_cell=True, exclude_similar=True, exclude_fractional=True, silent=True, check=True, filename_as_uid=False, methods=['ase','cif2struct'],modify_cif=False,dry_run=False):
        """ import_cif(cifnames): import structures from cif files in HTE structure database
        cifnames: shell expression like '*.cif'
        check: True=self_check directly after importing, 'all'=self_check at the end for all imported cif-structures         at once, False=no self_check.
        """
        # import structures from cif file(s)
        print '* Importing structure data from ',cifnames
        cif_list=glob.glob(cifnames)
        import_statistics={'imported':[],'already_in':[],'excluded':[],'failed':[]}
        nentry_total=0
        for ciffile in cif_list:
            if os.path.isfile(ciffile):
                if silent==False:
                    print 'Evaluating cif file ',ciffile
                infile=open(ciffile,'r')
                cif_blocks=[]
                cif_block=''
                uid=None
                symloop=False
                line=infile.readline()
                while line:
                    #if silent==False:
                    #    print line
                    #check again for multiple cifs
                    #if (((line.startswith('# End of data set')) or (line.startswith('_cod_database_code'))) and (uid!=None)): 
                    if ((line.startswith('# End of data set')) or (line.startswith('#End of data'))) and (uid!=None):
                        cif_block=cif_block+line
                        cif_blocks.append((uid,cif_block))
                        cif_block=''
                        uid=None
                    elif (line.startswith('data_')):
                        nentry_total=nentry_total+1
                        if uid!=None:
                            cif_blocks.append((uid,cif_block))
                            cif_block=line
                        else:
                            cif_block=cif_block+line
                        uid=line.strip()
                        if (filename_as_uid==True):
                            uid=ciffile.split('/')[-1].split('.cif')[0]
                    elif '_symmetry_equiv_pos_as_xyz' in line:
                        symloop=True
                        cif_block=cif_block+line
                    elif 'loop_' in line:
                        symloop=False
                        cif_block=cif_block+line
                    else:
                        if (symloop==True) and (modify_cif==True):
                            #replace 0.333 with 1/3, TODO: check
                            num=""
                            replaces=[]
                            for a in line:
                                if (a.isdigit()) or (a=='.'):
                                    num=num+a
                                else:
                                    if num!="":
                                        for ystr,y in [("1/3",1./3.),("2/3",2./3.),("1/2",0.5)]:
                                            if abs(float(num)-y)<0.001:
                                                replaces.append((num,ystr))
                                    num=""
                            for r,y in replaces:
                                line=line.replace(r,y)
                        if (modify_cif==True) and ("'" in line):
                            #filter line with ' in author names etc.
                            nx=0
                            linef=""
                            for x in line:
                                if x=="'":
                                    if nx==0:
                                        nx=1
                                    else:
                                        nx=0
                                else:
                                    linef=linef+x
                            if nx==1:
                                line=linef
                        cif_block=cif_block+line
                    line=infile.readline()
                if (cif_blocks==[])  and (uid!=None):#and (cif_block!='')
                    cif_block=cif_block+line
                    cif_blocks.append((uid,cif_block))
                    cif_block=''
                    uid=None
                #print cif_blocks
                for index in range(len(cif_blocks)):
                    (uid,block)=cif_blocks[index]
                    if  (uid in self.structureDB) and (dry_run==False):
                        import_statistics['already_in'].append((ciffile,uid))
                        if silent==False:
                            print "already in database:",uid
                    else:
                        cif_sources={'cif_file_external':os.path.abspath(ciffile),
                                     'cif_index_external':index,
                                     'cif_id':uid}
                        if self.get_cif_repository()!=None:
                            cif_sources['cif_file']=os.path.join(self.get_cif_repository(),
                                                                 uid+'.cif')
                        fname='tmp_hte_'+uid+'.cif'
                        outfile=open(fname,"w")
                        outfile.write(block)
                        outfile.close()
                        is_disordered=False
                        for method in methods:
                            if method=='ase':
                                try:
                                    cifstruct=read_cif(fname, primitive_cell=primitive_cell)
                                except:
                                    cifstruct=None
                                    if silent==False:
                                        print "ase_cif_reader: failed to read ",ciffile,uid
                            elif method=='cif2struct':
                                try:
                                    exitcode, mes=commands.getstatusoutput("cif2struct %s"%fname)
                                    structfile=fname.split(".cif")[0]+".struct"
                                    cifstruct=read_struct(structfile)
                                    if silent==False:
                                        print "succes: %s %s (%s)"%(ciffile,uid,fname)
                                except:
                                    cifstruct=None
                                    if silent==False:
                                        print "cif2struct: failed to read ",ciffile,uid
                            dbentry=None
                            if cifstruct!=None:
                                print "check_point178, cifstruct!=None",cifstruct
                                source_info=[('import_cif',cif_sources)]
                                dbentry=HTEdbentry(cifstruct,uid,_source_info=source_info)
                                print "check_point179, dbentry is ",dbentry
                                dbentry.get_data_from_cif_source(fname)
                                if (dbentry.is_disordered()):
                                    is_disordered=True
                                if (check==True):
                                    warn,err=dbentry.self_check()
                                    print "check_point180, self_check err is ", err
                                    if err!=None:
                                        dbentry=None
                                        if (silent==False):
                                            print uid,method,err
                            if (dbentry!=None) or ((is_disordered==True) and (exclude_fractional==True)):
                                break
                        if ((is_disordered==True) and (exclude_fractional==True)):
                            import_statistics['excluded'].append((ciffile,uid))
                            if silent==False:
                                print "import_cif:",ciffile,uid," excluded because of fractional occupations"
                        elif dbentry==None:
                            import_statistics['failed'].append((ciffile,uid))
                            if silent==False:
                                print "WARNING(import_cif): failed to import ",ciffile,uid
                        elif dry_run==True:
                            return dbentry,source_info
                        elif self.add_structDB_entry(uid, structureDB_entry=dbentry,_source_info=source_info, exclude_similar=exclude_similar, silent=silent)==True:
                            import_statistics['imported'].append((ciffile,uid))
                            if self.get_cif_repository()!=None:
                                shutil.copy(fname,os.path.join(self.get_cif_repository(),uid+'.cif'))
                        else:
                            import_statistics['excluded'].append((ciffile,uid))
                        try:
                            for fname in glob.glob('tmp_hte_'+uid+'.*'):
                                #print fname
                                os.remove(fname)
                        except:
                            if silent==False:
                                print "WARNING(import_cif): failed to remove tempory files"
        print " Imported %d structures from %d [%d already in db/%d failed/%d excluded]"%(
            len(import_statistics['imported']),nentry_total,len(import_statistics['already_in']),
            len(import_statistics['failed']),len(import_statistics['excluded']))
        if (check=='all'):
            import pprint
            pprint.pprint(self.self_check())
        if dry_run==True:
            return None,None
        return import_statistics
        
    def select(self, **kwargs):
        """return list with structure uids matching conditions in kwargs, e.g.:
        "'nelements'=2,'elements'=['Fe','Co','Ni']" -> all binary FeCo, FeNi and CoNi alloys
        "'nelements'=2,'elements'=['Fe']" -> all binary FeX alloys, X=any element
        possible keywords (combinations with logical AND):
        nelements: number of different elements in alloy (integer)
        nelements_max: maximum number of different elements in alloy
        nelements_min: minimum number of different elements in alloy
        uid_list: uid has to be in given list
        uid_matches: substring must be in uid
        uid_matches_not: substring must not be in uid
        atfraction_max: alloys with at most at% of elements in dict, e.g:
                        atfraction_max={'Si':0.5} -> alloys with at most 50% Si content
        atfraction_min: alloys with at least at% of elements in dict
        structure_type: structure type of uid matches full('CsCl,cP2,221') or short('CsCl')
                        structure type
                        
        ### unknown1:structure_type?
        structure_types: same as structure_type for a list of types
        spacegroup_no: spacegroup number of uid matches given spacegroup number
        calculated    : Is there a non empty stored_calc_results-dictionary? If 
                        yes: select
        """
        # return list with structures uids according to keyword list
        magsettings={'submitted':True}
        if 'magsettings' in kwargs:
            magsettings=kwargs['magsettings']
        if 'uid_list' in kwargs:
            uidlist_in=kwargs['uid_list']
        else:
            uidlist_in=self.structureDB.keys()
            ### uid就是每个条目的键值
        uidlist=[]
        if 'update' in kwargs:
            update=kwargs['update']
        else:
            update=False
            
            ###默认FALSE 除非给定了
        if 'formation_energy_max' in kwargs:
            if ('calc_scheme' in kwargs) and (kwargs['calc_scheme'] in self.calc_schemes):
                calc_scheme=kwargs['calc_scheme']
                E_ref={}
                
                ##如果参数里有formation_energy_max 并且给定了存在的calc_scheme
            else:
                print 'select(): formation_energy_max needs calc_scheme, nothing done!'
                return None
                ### 否则直接返回None
        for uid in uidlist_in:
            if ('all' in kwargs) and (kwargs['all']==True):
                uidlist.append(uid)
            else:
                add=True
                #if 'calc_scheme' in kwargs:
                #    calc=self.setup_calculator(uid, kwargs['calc_scheme'])
                if ('do_not_calc' in kwargs) and (self.structureDB[uid].do_not_calc==False):
                    add=False
                elif ('calculated' in kwargs and kwargs['calculated'] == True 
                      and not self.structureDB[uid].stored_calc_results):
                    add=False
                elif (not ('do_not_calc' in kwargs)) and (self.structureDB[uid].do_not_calc==True):
                    add=False
                elif ('nelements' in kwargs) and (self.structureDB[uid].get_number_of_elements()!=kwargs['nelements']):
                    
                    ###通过get_number_of_elements这个函数判断每个条目包含的元素个数
                    add=False
                elif ('nelements_max' in kwargs) and (self.structureDB[uid].get_number_of_elements()>kwargs['nelements_max']):
                    add=False
                elif ('nelements_min' in kwargs) and (self.structureDB[uid].get_number_of_elements()<kwargs['nelements_min']):
                    add=False
                elif ('natoms_max' in kwargs) and (len(self.structureDB[uid].atoms_initial)>kwargs['natoms_max']):
                    add=False
                if (add==True) and ('elements' in kwargs):
                    elements=self.structureDB[uid].get_composition()
                    if self.structureDB[uid].get_number_of_elements()>=len(kwargs['elements']):
                        #require that all elements in the list are in compound
                        for el in kwargs['elements']:
                            if not (el in elements):
                                add=False
                            ###这个是通过成分判断加不加到选择列表里
                    else:
                        #require that all elements in the compound are in the list
                        for el in elements:
                            if not (el in kwargs['elements']):
                                add=False
                if (add==True) and ('exclude_elements' in kwargs):
                    elements=self.structureDB[uid].get_composition()
                    for el in elements:
                        if (el in kwargs['exclude_elements']):
                            add=False
                if (add==True) and ('uid_matches' in kwargs):
                    if isinstance(kwargs['uid_matches'],list):
                        matchlist=kwargs['uid_matches']
                    else:
                        matchlist=[kwargs['uid_matches']]
                    for ml in matchlist:
                        if not (ml in uid):
                            add=False
                if (add==True) and ('uid_matches_not' in kwargs):
                    if isinstance(kwargs['uid_matches_not'],list):
                        matchlist=kwargs['uid_matches_not']
                    else:
                        matchlist=[kwargs['uid_matches_not']]
                    for ml in matchlist:
                        if ml in uid:
                            add=False
                            break
                if (add==True) and ('atfraction_max' in kwargs):
                    atfrac_list=kwargs['atfraction_max']
                    for el in atfrac_list:
                        if self.structureDB[uid].get_atfraction(el)>atfrac_list[el]:
                            add=False
                if (add==True) and ('atfraction_min' in kwargs):
                    atfrac_list=kwargs['atfraction_min']
                    for el in atfrac_list:
                        if self.structureDB[uid].get_atfraction(el)<atfrac_list[el]:
                            add=False
                if (add==True) and ('structure_type' in kwargs):
                    structure_type=self.get_structure_type(uid)
                    short_type=structure_type
                    if structure_type!=None:
                        short_type=structure_type.split(',')[0]
                    if (kwargs['structure_type']!=structure_type) and (kwargs['structure_type']!=short_type):
                        add=False
                #
                if (add==True) and ('structure_types' in kwargs):
                    structure_type=self.get_structure_type(uid)
                    short_type=structure_type
                    if structure_type!=None:
                        short_type=structure_type.split(',')[0]
                    if (not(structure_type in kwargs['structure_types'])) and (not(short_type in kwargs['structure_types'])):
                        add=False
                        
                        ###上面我们可以看出对列表 和字符串两种数据的不同处理 方式 !=和not in
                #
                if (add==True) and ('spacegroup_no' in kwargs):
                    if (self.get_spacegroup_number(uid)!=kwargs['spacegroup_no']):
                        add=False
                if (add==True) and ('formation_energy_max' in kwargs):
                    comp=self.structureDB[uid].get_composition(reduce=False)
                    for element in comp:
                        if not (element in E_ref):
                            #probe if reference structure was given (better way)
                            keyw='reference_'+element
                            if keyw in kwargs:
                                uid=kwargs[keyw]
                                E_ref[element]=self.get_energy_per_atom(uid,calc_scheme,update=update, silent=True)
                            else:
                                E_ref[element]=self.get_reference_energy_per_atom(element,calc_scheme=calc_scheme,update=update)
                    E_form=self.get_formation_energy(uid, calc_scheme, reference_energy=E_ref)
                    if (E_form==None) or (E_form>kwargs['formation_energy_max']):
                        add=False
                if (add==True) and ('dehull_max' in kwargs):
                    deh=self.get_distance_from_convex_hull(uid,kwargs['calc_scheme'],magsettings=magsettings)
                    if (deh==None) or (deh>kwargs['dehull_max']):
                        add=False
                if add==True:
                    uidlist.append(uid)
        if ('ground_state' in kwargs) and (kwargs['ground_state']==True):
            #select structure with lowest energy for each composition
            uidlist_gs=[]
            for uid in uidlist:
                add=True
                comp=self.structureDB[uid].get_composition(reduce=False)
                for uid2 in uidlist:
                    if uid==uid2:
                        continue
                    samecomp=True
                    for el in comp:
                        if (self.structureDB[uid].get_atfraction(el)-self.structureDB[uid2].get_atfraction(el))>0.0001:
                            samecomp=False
                    if (samecomp==False):
                        continue
                    if (self.get_energy_per_atom(uid,kwargs['calc_scheme'],magsettings=magsettings)<self.get_energy_per_atom(uid2,kwargs['calc_scheme'],magsettings=magsettings)):
                        continue
                    add=False
                    break
                if (add==True):
                    uidlist_gs.append(uid)
            return uidlist_gs
        return uidlist

    def mark_do_not_calc(self, uid, comment=None):
        if uid in self.structureDB:
            self.structureDB[uid].do_not_calc=True
            if self.structureDB[uid].comment!=None:
                if comment!=None:
                    self.structureDB[uid].comment=self.structureDB[uid].comment+'\n  '+comment
            else:
                self.structureDB[uid].comment=comment
            print '* do_not_calc flag set for uid ',uid
            
    def exclude_uids(self, uid_list={}, query='all'):
        #uid_list: dictionary with uid and comments
        nexcluded=0
        print 'About to set do_not_calc flag for %d data base entries'%len(uid_list)
        if query=='once':
            yesno=raw_input('Shall I do this [y(es)/a(bort)] ? ')
            if yesno.startswith('y'):
                print 'O.k., starting...'
            else:
                print 'Nothing done, abort..'
                return nexcluded
        abort=False
        delentry=True
        for uid in uid_list:
            print '**********',uid,':',self.structureDB[uid].get_composition(reduce=True)
            print '- original source:',self.get_cif_source(uid)
            print '- atoms.info:',self.structureDB[uid].atoms_initial.info
            print '- joblist:',self.structureDB[uid].submitted_jobs
            print '-> comment:',uid_list[uid]
            if query=='all':
                yesno=raw_input('Shall I exclude this entry from DB [y(es)/n(o)/a(bort)/!(exclude all)] ?')
                if yesno.startswith('a'):
                    abort=True
                elif yesno.startswith('y'):
                    delentry=True
                elif yesno=='!':
                    delentry=True
                    query='once'
                else:
                    delentry=False
            if abort:
                break
            if delentry:
                self.mark_do_not_calc(uid,comment=uid_list[uid])
                print 'do_not_calc flag set for ',uid
                nexcluded=nexcluded+1
                
                
                ### 询问是不是把这个uid添加到排除列表里面 也就是直接添加到了 do_not_calc条目属性里面 也会在comment里面加上uid信息
            else:
                print 'nothing done for ',uid
        print '** do_not_calc flag set for ',nexcluded,' data base entries'
        return nexcluded
    
    def create_interstitial_defect(self, uid, interstitial_atom, calc_scheme=None, rmin='auto', scale_radii=0.7, radii='auto', Nmax=10, nscan=24, symprec=1e-3,distortion=[],exclude_similar=False
                                   ,std_lattice=False, supercell=[], interstitial_blocks=[[0,0,0]],prefer_high_sym=False, append_to_uid="",exclude_positions=[],silent=False):
        """ Identify potential interstitial positions and create a defect cells
        with one interstiatial atom at this position
        uid: unique identifier of initial structure in HTE database
        interstitial_atom: atom to be put at interstitial positions ('N', 'C', etc.)
        rmin: minimal radius for interstitial position in the structure
              (default is 'auto', value is taken from covalent radii in ASE and scaled with scale_radii)
        Nmax: maximum number of interstitial positions
        returns a list with newly created uids
    """
        new_uids=[]
        ao_ini=self.get_atoms_object(uid, calc_scheme=calc_scheme)
        if ao_ini==None:
            print 'create_interstitial_defect: No atoms object for uid=%s and calc_scheme=%s, abort'%(uid,str(calc_scheme))
            return new_uids
        if std_lattice==True:
            spglib_info=spglib.get_symmetry_dataset(ao_ini, symprec=symprec)
            ao_ini=Atoms(numbers=spglib_info['std_types'],cell=spglib_info['std_lattice'],scaled_positions=spglib_info['std_positions'],pbc=True)
            print 'std_lattice',ao_ini.get_cell(),spglib_info['number'],spglib_info['international']
            for el,pos in zip(ao_ini.get_chemical_symbols(),ao_ini.get_scaled_positions()):
                print el,pos
        # make sure that scaled positions are in [0,1):
        cell=ao_ini.get_cell()
        new_symbols=ao_ini.get_chemical_symbols()
        new_scaled_positions=[]
        for x in ao_ini.get_scaled_positions():
            new_scaled_positions.append(x%1.)
        newatoms=Atoms(symbols=new_symbols,scaled_positions=new_scaled_positions,cell=cell,pbc=True)
        spglib_info=spglib.get_symmetry_dataset(newatoms, symprec=symprec)
        # find potential sites for interstitial atoms
        interst_pos={}
        if (rmin=='auto'):
            Z=atomic_numbers[interstitial_atom]
            rmin=covalent_radii[Z]*scale_radii
        Nist_sites=0
        N_radii=0 #len(newatoms)
        while (Nist_sites<=Nmax):
            newpos=Identify_Interstitial_Position(newatoms,rmin, nscan=nscan, symprec=symprec,interstitial_atom=interstitial_atom,scale_radii=scale_radii, radii=radii,N_radii=N_radii,prefer_high_sym=prefer_high_sym)
            if newpos==[]:
                break
            Nist_sites=Nist_sites+1
            # create new atoms object with interstitials at all symmetry eqivalent positions (for symmetry label)
            new_symbols=newatoms.get_chemical_symbols()
            new_scaled_positions=[]
            for x in newatoms.get_scaled_positions():
                new_scaled_positions.append(x)
            multiplicity=0
            for rot,trans in zip(spglib_info['rotations'],spglib_info['translations']):
                gpos=(np.dot(rot,newpos)+trans)%1.
                #tmp
                for j in range(3):
                    if np.fabs(gpos[j]-1.0)<0.0001:
                        gpos[j]=gpos[j]-1.0
                is_new=True
                for pos in new_scaled_positions:
                    if np.sum(np.fabs(gpos-pos))<0.0001:
                        is_new=False
                        break
                if (is_new==True):
                    new_symbols.append(interstitial_atom)
                    new_scaled_positions.append(gpos)
                    multiplicity=multiplicity+1
            newatoms=Atoms(symbols=new_symbols,scaled_positions=new_scaled_positions,cell=cell,pbc=True)
            spglib_new=spglib.get_symmetry_dataset(newatoms, symprec=symprec)
            if "%d%s"%(multiplicity,spglib_new['wyckoffs'][-1]) in exclude_positions:
                continue
            label="%s_%d%s"%(interstitial_atom,multiplicity,spglib_new['wyckoffs'][-1])
            i=1
            while label in interst_pos:
                label="%s_%d%s_%d"%(interstitial_atom,multiplicity,spglib_new['wyckoffs'][-1],i)
                i=i+1
            interst_pos[label]=newpos
        if (Nist_sites>Nmax):
            print 'create_interstitial_defect: More than Nmax=%d interstitial positions found for uid=%s, abort'%(Nmax,uid)
            return new_uids
        #create structures with one interstitial at the positions identified and put them in database
        if (distortion==[]):
            new_cell=cell
        else:
            new_cell=np.dot(distortion,cell)
        for label in interst_pos:
            new_symbols=ao_ini.get_chemical_symbols()
            new_scaled_positions=[]
            for x in ao_ini.get_scaled_positions():
                new_scaled_positions.append(x)
            new_symbols.append(interstitial_atom)
            new_scaled_positions.append(interst_pos[label])
            new_atoms_object=Atoms(symbols=new_symbols,scaled_positions=new_scaled_positions,cell=new_cell,pbc=True)
            new_uid=uid+'_intst_%s'%label+append_to_uid
            if (supercell!=[]):
                # create a supercell
                aosc=Create_Interstitial_Supercell(new_atoms_object, interstial_atoms=[len(new_atoms_object)-1], interstitial_blocks=interstitial_blocks, N=supercell, silent=silent)
                ib=0
                iblabel=""
                for N_0 in range(supercell[0]):
                    for N_1 in range(supercell[1]):
                        for N_2 in range(supercell[2]):
                            if [N_0,N_1,N_2] in interstitial_blocks:
                                iblabel=iblabel+"_%d"%ib
                            ib=ib+1
                new_uid=uid+'_intst_%s_sc_%dx%dx%d_ib%s'%(label+append_to_uid,supercell[0],supercell[1],supercell[2],iblabel)
                new_atoms_object=aosc
            source_info=[] #deepcopy(self.structureDB[uid].source_info)
            source_info.append(('create_interstitial_defect',{'interstitial':interstitial_atom
                                                              ,'scaled position':interst_pos[label]
                                                              ,'wyckoff':label.split("_")[-1]
                                                              ,'distortion':distortion
                                                              ,'supercell':supercell
                                                              ,'interstitial_blocks':interstitial_blocks
                                                              ,'initial_uid':uid}))
            dbentry=HTEdbentry(new_atoms_object,new_uid,_source_info=source_info)
            #dbentry.cif_info=self.structureDB[uid].cif_info.copy()
            #dbentry.user_info=self.structureDB[uid].user_info.copy()
            if self.add_structDB_entry(new_uid, structureDB_entry=dbentry, silent=silent,exclude_similar=exclude_similar):
                new_uids.append(new_uid)
        self.self_check(uid_list=new_uids)
        return new_uids

    
    def get_parent_uid(self,uid):
        """returns uid of parent compound for interstitial alloy
        """
        uidp=uid
        try:
            meth,info=self.structureDB[uid].source_info[-1]
            if (meth=='create_interstitial_defect') and ('initial_uid' in info) and (info['initial_uid'] in self.structureDB):
                uidp=info['initial_uid']
        except:
            uidp=uid
        return uidp
                
    def get_interstitial_atom(self,uid):
        """returns interstitial atom for interstitial alloy
        """
        iatom="--"
        try:
            meth,info=self.structureDB[uid].source_info[-1]
            if (meth=='create_interstitial_defect') and ('interstitial' in info):
                iatom=info['interstitial']
        except:
            iatom="--"
        return iatom
                
    def get_interstitial_position(self,uid):
        """returns interstitial position for interstitial alloy
        """
        iatom="--"
        try:
            meth,info=self.structureDB[uid].source_info[-1]
            if (meth=='create_interstitial_defect') and ('wyckoff' in info):
                iatom=info['wyckoff']
        except:
            iatom="--"
        return iatom
                


    def create_test_set(self, elements=[], structure_type_file='', cif_dir='', primitive_cell=True, silent=True, exclude_similar=True, Nmax=-1, atfractions={},max_no_of_atoms=-1):
        if structure_type_file=='':
            cif_dir=os.path.join(self.dir_of_hte,'structure_types_%d'%len(elements))
            structure_type_file=os.path.join(cif_dir,'structure_types_%d.txt'%len(elements))
        new_uids=[]
        if not os.path.isfile(structure_type_file):
            print "WARNING(create_test_set): Could not find ",structure_type_file
            return new_uids
        infile=open(structure_type_file,'r')
        line=infile.readline()
        N_loaded=0
        while (line):
            if (Nmax>0) and (N_loaded>=Nmax):
                break
            n,strtype,uid=line.split()
            line=infile.readline()
            if silent==False:
                print  n,strtype,uid, elements
            cif_file=os.path.join(cif_dir,uid+".cif")
            dbentry,source_info=self.import_cif(cif_file,dry_run=True,silent=silent)
            if dbentry==None:
                print "Check:",cif_file
                continue
            cif_struct=dbentry.atoms_initial
            if (max_no_of_atoms>0) and (len(cif_struct)>max_no_of_atoms):
                continue
            cif_elements=cif_struct.get_chemical_symbols()
            element_types={}
            nel=0
            dummy_symbols=[]
            cif_symbols=[]
            for el in cif_elements:
                if el in element_types:
                    dummy_symbols.append(element_types[el])
                else:
                    element_types[el]=nel
                    nel=nel+1
                    dummy_symbols.append(element_types[el])
                    cif_symbols.append(el)
            if len(elements)!=nel:
                print "WARNING(create_test_set): structure type %s not suitable for %s"%(strtype,str(elements))
            new_symbols=[]
            for iel in dummy_symbols:
                new_symbols.append(elements[iel])
            cif_struct.set_chemical_symbols(new_symbols)
            #check if alloy has requested composition
            skip=False
            for el in atfractions:
                atfr=dbentry.get_atfraction(el)
                if (atfr<atfractions[el][0]) or (atfr>atfractions[el][1]):
                    skip=True
                    break
            if skip==True:
                continue
            N_loaded=N_loaded+1
            new_uid='htets'
            for x in strtype.split(','):
                new_uid=new_uid+'_'+x
            for el in elements:
                new_uid=new_uid+"_%s"%el
            dbentry.calcdir=new_uid
            source_info.append(('create_test_set',{'initial_uid':uid, 'elements':elements, 'cif_elements':cif_symbols}))
            if self.add_structDB_entry(new_uid, structureDB_entry=dbentry, silent=silent,exclude_similar=exclude_similar):
                new_uids.append(new_uid)
        return new_uids

    
    def create_test_set_new(self, elements=[], fix_elements=[], dbfile="", primitive_cell=True, silent=True, exclude_similar=True, Nmax=-1, atfractions={},max_no_of_atoms=-1):
        new_uids=[]
        subels=[]
        for el in elements:
            if not (el in fix_elements):
                subels.append(el)
        if len(subels)>2:
            print "create_test_set: too many elements"
            return new_uids
        if dbfile=='':
            dbfile=os.path.join(self.dir_of_hte,'structure_types/allcifs.hte')
        if os.path.isfile(dbfile):
            hte=HTE(dbfile, silent=True)
        sel=hte.select(nelements=len(elements),elements=fix_elements)
        typelist={}
        for uid in sel:
            strt=hte.get_structure_type(uid)
            if strt!=None:
                if strt in typelist:
                    typelist[strt].append(uid)
                else:
                    typelist[strt]=[uid]
        for strt in sorted(typelist,reverse=True,  key=lambda k: len(typelist[k])):
            uid=typelist[strt][0]
            comp=hte.structureDB[uid].get_composition(reduce=False)
            atom_lists=[]
            i=0
            for el in comp:
                if el in fix_elements:
                    continue
                if atom_lists==[]:
                    atom_lists=[{el:subels[0]},{el:subels[-1]}]
                else:
                    atom_lists[0][el]=subels[-1]
                    atom_lists[1][el]=subels[0]
            for atom_list in atom_lists:
                print comp,atom_list
                new_uids=new_uids+hte.substitute_atoms(uid_list=[uid],atom_list=atom_list,exclude_similar=False)
            if (Nmax>0) and (len(new_uids)>Nmax):
                break
        for uid in new_uids:
            dbentry=hte.structureDB[uid]
            self.add_structDB_entry(uid, structureDB_entry=dbentry, silent=silent,exclude_similar=False)
        hte.output(new_uids,"uid","cell_composition","structure_type")
        hte=None
        return new_uids
        new_uids=[]
        if not os.path.isfile(structure_type_file):
            print "WARNING(create_test_set): Could not find ",structure_type_file
            return new_uids
        infile=open(structure_type_file,'r')
        line=infile.readline()
        N_loaded=0
        while (line):
            if (Nmax>0) and (N_loaded>=Nmax):
                break
            n,strtype,uid=line.split()
            line=infile.readline()
            if silent==False:
                print  n,strtype,uid, elements
            cif_file=os.path.join(cif_dir,uid+".cif")
            dbentry,source_info=self.import_cif(cif_file,dry_run=True,silent=silent)
            if dbentry==None:
                print "Check:",cif_file
                continue
            cif_struct=dbentry.atoms_initial
            if (max_no_of_atoms>0) and (len(cif_struct)>max_no_of_atoms):
                continue
            cif_elements=cif_struct.get_chemical_symbols()
            element_types={}
            nel=0
            dummy_symbols=[]
            cif_symbols=[]
            for el in cif_elements:
                if el in element_types:
                    dummy_symbols.append(element_types[el])
                else:
                    element_types[el]=nel
                    nel=nel+1
                    dummy_symbols.append(element_types[el])
                    cif_symbols.append(el)
            if len(elements)!=nel:
                print "WARNING(create_test_set): structure type %s not suitable for %s"%(strtype,str(elements))
            new_symbols=[]
            for iel in dummy_symbols:
                new_symbols.append(elements[iel])
            cif_struct.set_chemical_symbols(new_symbols)
            #check if alloy has requested composition
            skip=False
            for el in atfractions:
                atfr=dbentry.get_atfraction(el)
                if (atfr<atfractions[el][0]) or (atfr>atfractions[el][1]):
                    skip=True
                    break
            if skip==True:
                continue
            N_loaded=N_loaded+1
            new_uid='htets'
            for x in strtype.split(','):
                new_uid=new_uid+'_'+x
            for el in elements:
                new_uid=new_uid+"_%s"%el
            dbentry.calcdir=new_uid
            source_info.append(('create_test_set',{'initial_uid':uid, 'elements':elements, 'cif_elements':cif_symbols}))
            if self.add_structDB_entry(new_uid, structureDB_entry=dbentry, silent=silent,exclude_similar=exclude_similar):
                new_uids.append(new_uid)
        return new_uids

    def get_substitutional_elements(self, element, type='column'):
        element_blocks={
            '3d_transition_metals':['Sc', 'Ti', 'V', 'Cr','Mn','Fe','Co','Ni','Cu','Zn'],
            '4d_transition_metals':['Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd'],
            '5d_transition_metals':['La','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg'],
            'sp_row3':['Na','Mg','Al','Si','P','S'],
            'sp_row4':['K','Ca','Ga','Ge','As','Se'],
            'sp_row5':['Rb','Sr','In','Sn','Sb','Te'],
            'sp_row6':['Cs','Ba','Tl','Pb','Bi','Po']
            }
        #
        substitutions={}
        if isinstance(element, list):
            elements=element
        else:
            elements=[element]
        for el in elements:
            substitutions[el]=[el]
            for block in element_blocks:
                if el in element_blocks[block]:
                    index=element_blocks[block].index(el)
                    if 'transition_metals' in block:
                        substitutions[el]=[element_blocks['4d_transition_metals'][index],element_blocks['3d_transition_metals'][index],element_blocks['5d_transition_metals'][index]]
                    elif block=='sp_row3':
                        substitutions[el]=[element_blocks['sp_row3'][index],element_blocks['sp_row4'][index],element_blocks['sp_row5'][index]]
                    elif block=='sp_row4':
                        substitutions[el]=[element_blocks['sp_row4'][index],element_blocks['sp_row5'][index],element_blocks['sp_row3'][index]]
                    elif block=='sp_row5':
                        substitutions[el]=[element_blocks['sp_row5'][index],element_blocks['sp_row6'][index],element_blocks['sp_row4'][index]]
                    elif block=='sp_row6':
                        substitutions[el]=[element_blocks['sp_row6'][index],element_blocks['sp_row5'][index],element_blocks['sp_row4'][index]]
                    break
        return substitutions

        
    def substitute_atoms(self, uid_list=[], atom_list={}, append_to_uid='', silent=False, exclude_similar=True):
        """Substitute atoms for structures in uid_list and add the new structure to HTE structureDB:
        uid_list: list with uids of HTE database 
        atom_list: dictionary with substitutions,
                   e.g. {'Fe':'Co','Pd':'Pt'} will create CoPt from FePd
                   new since HTE version 1.1:
                       {'Fe':'Co','Co':'Ni'} will replace all Fe atoms by Co and all Co atoms by Ni as expected
                       {'Fe_g':'Co'} will replace Fe atoms on Wyckoff position g by Co
        append_to_uid: name appended to uid of new structure, default chosen by HTE
        returns a list with uid of newly created structures
        """
        new_uids=[]
        if not isinstance(atom_list,dict):
            print "WARNING(substitute_atoms): from HTE version 1.1 atom_list must be a dictionary (see help(HTE.substitute_atoms)), nothing done!"
        for els in sorted(atom_list):
            if '_' in els:
                el=els.split('_')[0]
            else:
                el=els
            if atom_list[els]==el:
                del atom_list[els]
        for uid in uid_list:
            if (uid in self.structureDB) and (self.structureDB[uid].is_runnable()):
                source_info=deepcopy(self.structureDB[uid].source_info)
                source_info.append(('substitute_atoms',{'atom_list':[atom_list], #for compatibility with old version
                                                        'initial_uid':uid,
                                                        'hte_version':self.get_version()}))
                new_atoms_object=self.structureDB[uid].atoms_initial.copy()
                new_chem_symbols=new_atoms_object.get_chemical_symbols()
                try:
                    wyckoffs=spglib.get_symmetry_dataset(new_atoms_object)['wyckoffs']
                except:
                    wyckoffs=len(new_chem_symbols)*['']
                uidxts={}
                for i in range(len(new_chem_symbols)):
                    el_w="%s_%s"%(new_chem_symbols[i],wyckoffs[i])
                    if new_chem_symbols[i] in atom_list:
                        uidxts[new_chem_symbols[i]]=atom_list[new_chem_symbols[i]]
                        new_chem_symbols[i]=atom_list[new_chem_symbols[i]]
                    elif el_w in atom_list:
                        uidxts[el_w]=atom_list[el_w]
                        new_chem_symbols[i]=atom_list[el_w]
                        
                if uidxts=={}:
                    continue
                if append_to_uid=='':
                    append_to_uid='_subst'
                    for el in sorted(uidxts):
                        append_to_uid=append_to_uid+'_%s-%s'%(el,uidxts[el])
                new_uid=uid+append_to_uid
                new_atoms_object.set_chemical_symbols(new_chem_symbols)
                dbentry=HTEdbentry(new_atoms_object,new_uid,_source_info=source_info)
                dbentry.cif_info=self.structureDB[uid].cif_info.copy()
                if 'wyckoff_elem' in dbentry.cif_info:
                    for i in range(len(dbentry.cif_info['wyckoff_elem'])):
                        if dbentry.cif_info['wyckoff_elem'][i] in atom_list:
                            dbentry.cif_info['wyckoff_elem'][i]=atom_list[dbentry.cif_info['wyckoff_elem'][i]]
                dbentry.user_info=self.structureDB[uid].user_info.copy()
                if self.add_structDB_entry(new_uid, structureDB_entry=dbentry, silent=silent,exclude_similar=exclude_similar):
                    new_uids.append(new_uid)
        return new_uids

    
    def psubstitute_atoms(self, uid, atom_list={}, append_to_uid='', silent=False, std_lattice=False, exclude_similar=False):
        """Partial substitution of one atoms for structures in uid_list and add the new structure to HTE structureDB:
        uid_list: list with uids of HTE database 
        atom_list: dictionary with substitutions,
                   e.g. {'Fe':'Co','Pd':'Pt'} will create CoPt from FePd
                   new since HTE version 1.1:
                       {'Fe':'Co','Co':'Ni'} will replace all Fe atoms by Co and all Co atoms by Ni as expected
        append_to_uid: name appended to uid of new structure, default chosen by HTE
        returns a list with uid of newly created structures
        """
        new_uids=[]
        if not isinstance(atom_list,dict):
            print "WARNING(substitute_atoms): from HTE version 1.1 atom_list must be a dictionary (see help(HTE.substitute_atoms)), nothing done!"
        for el in sorted(atom_list):
            if atom_list[el]==el:
                del atom_list[el]
        if (uid in self.structureDB) and (self.structureDB[uid].is_runnable()):
            ao_ini=self.get_atoms_object(uid)
            spglib_info=spglib.get_symmetry_dataset(ao_ini)
            if std_lattice==True:
                ao_ini=Atoms(numbers=spglib_info['std_types'],cell=spglib_info['std_lattice'],scaled_positions=spglib_info['std_positions'],pbc=True)
                spglib_info=spglib.get_symmetry_dataset(ao_ini)
            chem_symbols=ao_ini.get_chemical_symbols()
            uidxts={}
            for i in range(len(chem_symbols)):
                if chem_symbols[i] in atom_list:
                    uidxt="_psubst_%s_%s-1%s"%(chem_symbols[i],spglib_info['wyckoffs'][i],atom_list[chem_symbols[i]])
                    if (uidxt in uidxts) and (spglib_info['equivalent_atoms'][i]==i):
                        uidxt="_psubst_%s%d_%s-1%s"%(chem_symbols[i],i,spglib_info['wyckoffs'][i],atom_list[chem_symbols[i]])
                    if not (uidxt in uidxts):
                        new_chem_symbols=deepcopy(chem_symbols)
                        new_chem_symbols[i]=atom_list[chem_symbols[i]]
                        uidxts[uidxt]=new_chem_symbols
            for uidxt in uidxts:
                new_atoms_object=ao_ini.copy()
                for i in range(len(uidxts[uidxt])):
                    if uidxts[uidxt][i]=='E':
                        del new_atoms_object[i]
                        del uidxts[uidxt][i] # do not remove more than one atom!
                        break
                new_atoms_object.set_chemical_symbols(uidxts[uidxt])
                source_info=deepcopy(self.structureDB[uid].source_info)
                source_info.append(('psubstitute_atoms',{'atom_list':[atom_list], #for compatibility with old version
                                                         'initial_uid':uid,
                                                         'std_lattice':std_lattice,
                                                         'hte_version':self.get_version()}))
                if std_lattice==True:
                    new_uid=uid+'_std_latt'+uidxt
                else:
                    new_uid=uid+uidxt
                dbentry=HTEdbentry(new_atoms_object,new_uid,_source_info=source_info)
                #dbentry.cif_info=self.structureDB[uid].cif_info.copy()
                if self.add_structDB_entry(new_uid, structureDB_entry=dbentry, silent=silent,exclude_similar=exclude_similar):
                    new_uids.append(new_uid)
        return new_uids

    
    def substitute_atoms_old(self, uid_list=None, atom_list=None, append_to_uid=None, silent=False, exclude_similar=True):
        """Substitute atoms for structures in uid_list and add the new
        structure to HTE structureDB:
        uid_list: list with uids of HTE database 
        atom_list: dictionary or list of dictionaries with substitutions,
                   e.g. {'Fe':'Co','Pd':'Pt'}
                   will create CoPt from FePd
                   not allowed: {'Fe':'Co','Co':'Ni'} because it will produce
                   unpredictable results, use instead [{'Fe':'Co'},{'Co':'Ni'}]
                   to first replace all Fe with Co and than all Co with Ni
                   atoms.
        append_to_uid: name appended to uid of new structure, typically chosen
                       by HTE
        returns a list with uid of newly created structures
        """
        new_uids=[]
        if isinstance(atom_list,list):
            atlist=atom_list
        else:
            atlist=[atom_list]
        #check if atom_list is allowed expression
        uidxts=''
        for substitution in atlist:
            ellist=[]
            uidxts=uidxts+'subst_'
            for el in substitution:
                if (el in ellist) or (substitution[el] in ellist):
                    print 'check substitution: %s'%str(substitution)
                    return None
                else:
                    ellist.append(el)
                    ellist.append(substitution[el])
                uidxts=uidxts+'%s-%s'%(el,substitution[el])
        if append_to_uid==None:
            append_to_uid=uidxts
        for uid in uid_list:
            if (uid in self.structureDB) and (self.structureDB[uid].is_runnable()):
                source_info=deepcopy(self.structureDB[uid].source_info)
                source_info.append(('substitute_atoms',{'atom_list':atlist,
                                                        'initial_uid':uid}))
                new_uid=uid+append_to_uid
                new_chem_composition=self.structureDB[uid].get_composition(reduce=False,
                                                                      methods=['cif'])
                new_atoms_object=self.structureDB[uid].atoms_initial.copy()
                new_chem_symbols=new_atoms_object.get_chemical_symbols()
                new_struct_cpy=deepcopy(self.structureDB[uid])#.cif_info)#.copy()
                try:
                    new_wyckoff_elem=new_struct_cpy.get_wyckoff(sort=None)['wyckoff_elem']
                except:
                    new_wyckoff_elem=None
                for substitution in atlist:
                    for i in range(len(new_chem_symbols)):
                        if new_chem_symbols[i] in substitution:
                            new_chem_symbols[i]=substitution[new_chem_symbols[i]]
                    if new_wyckoff_elem!=None:
                        for i in range(len(new_wyckoff_elem)):
                            if new_wyckoff_elem[i] in substitution:
                                new_wyckoff_elem[i]=substitution[new_wyckoff_elem[i]]
                new_atoms_object.set_chemical_symbols(new_chem_symbols)
                dbentry=HTEdbentry(new_atoms_object,new_uid,_source_info=source_info)
                dbentry.cif_info=self.structureDB[uid].cif_info.copy()
                if new_wyckoff_elem!=None:
                    dbentry.cif_info['wyckoff_elem']=new_wyckoff_elem     
                dbentry.user_info=self.structureDB[uid].user_info.copy()
                if self.add_structDB_entry(new_uid, structureDB_entry=dbentry, silent=silent,exclude_similar=exclude_similar):
                    new_uids.append(new_uid)
        return new_uids

    
    def get_composition(self, ao, reduce=True):
        #return string with chemical composition of atoms object ao
        #atoms.get_chemical_symbols() behaves sometimes strange..
        comp={}
        for Z in ao.arrays['numbers']:
            sym=chemical_symbols[Z]
            if sym in comp:
                comp[sym]=comp[sym]+1
            else:
                comp[sym]=1
        formula=''
        for sym in sorted(comp.keys()):
            formula=formula+sym+str(comp[sym])
        return formula
        ###获取化学表达式
        
        
    def struc_diff_wyckoff(self, structA, structB):
        #experimental!
        #wyckoffA=zip(*structA.info['wyckoff_info'])
        #wyckoffB=zip(*structB.info['wyckoff_info'])
        liste=[]
        #'''(???IO?
        for pos in structA.get_wyckoff(sorted=True):
            print 'pos',pos
            print structB.get_wyckoff(sorted=True)
            if pos in structB.get_wyckoff(sorted=True):
                liste.append(True)
                print liste
        print 'liste', liste
        #'''
        liste=[True for pos in structA.get_wyckoff(sort=True) if pos in structB.get_wyckoff(sort=True)]
        #return any(True for pos in structA.get_wyckoff(sorted=True)
        #           if pos in structB.get_wyckoff(sorted=True))
        if (any(liste) and (len(liste)==len(structB.get_wyckoff(sort=True)))): 
            return True
        else:
            return False
        #return any(True for pos in wyckoffA if pos not in wyckoffB)

    def struct_diff(self, structA, structB, Tol=0.5, comp_sg=True):
        #compare two atoms objects and return True if objects have same SG and composition and structure type
        #(or differ by more than Tol if comp_sg=False)
        #(AND have same elements at same wyckoff positions, multiplicities.)
        #status: experimental
        structA_dbentry=None
        structB_dbentry=None
        if isinstance(structA, Atoms)==False:
            structA_dbentry=structA
            structA=structA.atoms_initial
        if isinstance(structB, Atoms)==False:
            structB_dbentry=structB
            structB=structB.atoms_initial
        #print structA.get_number_of_atoms(), structB
        if structA.get_number_of_atoms()!=structB.get_number_of_atoms():
            return True
        if ('spacegroup' in structA.info) and ('spacegroup' in structB.info):
            if structA.info['spacegroup'].no!=structB.info['spacegroup'].no:
                return True
            
            ### 很重要的逻辑 就是先判断有没有 再判断一不一样
            
        if self.get_composition(structA)!=self.get_composition(structB):
            return True
        #'''
        if ((structA_dbentry) and (structB_dbentry) and ('_atom_site_wyckoff_symbol' in structA) and 
            ('_atom_site_wyckoff_symbol' in structB) and
            (self.struc_diff_wyckoff(structA_dbentry, structB_dbentry)==False)):
            return True
        #'''
        if (structA_dbentry!=None) and (structB_dbentry!=None) and (structA_dbentry.get_structure_type()!=structB_dbentry.get_structure_type()):
            return True
        if ('structure_type' in structA.info) and ('structure_type' in structB.info):
            if (structA.info['structure_type']!=structB.info['structure_type']):
                return True
        if comp_sg==True:
            return False
        atnumsA=structA.arrays['numbers']
        atnumsB=structB.arrays['numbers']
        posA=structA.get_positions()
        posB=structB.get_positions()
        for i in range(len(atnumsA)):
            for j in range(len(atnumsB)):
                found=False
                if atnumsA[i]==atnumsB[j]:
                    #print 'compare ',posA[i],posB[j]
                    found=True
                    for k in range(3):
                        if np.abs(posA[i][k]-posB[j][k])>Tol:
                            found=False
                            
                    ####判断原子位置是否一致
                    
                if found==True:
                    #print posA[i],posB[j],' similar'
                    break
            if found==False:
                #print posA[i],' not found'
                return True
        return False
    
    
    
    ############################################
    # routines to run & evaluate calculations  #
    ############################################
    def setup_calculator(self, uid, calc_scheme, basic=False, update_atoms=True, transport=False, bandstructure=False, bandstructure_init=False, transport_kspace_density=None,elastic=False,emat=None,afm={},update=False, options={},return_settings=False):
        print "check_point62",uid, calc_scheme,afm
        # set up calculator for calc_scheme and uid:
        return_failed=None
        if return_settings==True:
            return_failed=None,{}
        if not (calc_scheme in self.calc_schemes):
            print 'setup_calculator(): %s undefined, nothing done!'%calc_scheme
            return return_failed
        if not (uid in self.structureDB):
            print 'setup_calculator(): unknown uid %s, nothing done!'%uid
            return return_failed
        calculator_name, settings=self.calc_schemes[calc_scheme]
        print "check_point65, calculator_name, settings=:", calculator_name, settings
        if (basic==True): #convergence tests etc. do not need to know about details
            print "check_point66, here show the basic is True"
            if calculator_name.lower()=='vasp':
                if 'ibrion' in settings:
                    return Vasp(ibrion=settings['ibrion'])
                else:
                    return Vasp()
            elif (calculator_name.lower()=='fplo'):
                job_settings=self.get_job_commands(calc_scheme)
                return FPLO(fedit=job_settings['fedit'])
            else:
                return None
        magatoms={'Cr':2.5,'Mn':3.5,'Fe':2.5,'Co':1.5,'Ni':0.8}
        # always start with the initial atoms object to ensure reproducabilty
        if (update_atoms==True):
            self.structureDB[uid].atoms=self.structureDB[uid].atoms_initial.copy()
        pass2calc=settings.copy()
        print "check_point61, pass2calc is :",pass2calc
        if options!={}:
            print "check_point129,options!={}",options
            for keyw in options:
                pass2calc[keyw]=options[keyw]
        if (elastic==True) and (update_atoms==True):
            tmpatoms=self.get_atoms_object(uid,calc_scheme)
            print "check_point130,(elastic==True) and (update_atoms==True),tmpatoms:",tmpatoms
            if tmpatoms==None:
                return None
            self.structureDB[uid].atoms=tmpatoms
            self.structureDB[uid].atoms.set_cell(np.dot(emat,tmpatoms.get_cell()), scale_atoms=True)
            if 'init_structure' in pass2calc:
                del pass2calc['init_structure']
            if 'scale_volume' in pass2calc:
                del pass2calc['scale_volume']
        if ('init_structure' in pass2calc) and (update_atoms==True):
            if (isinstance(pass2calc['init_structure'],list)==True):
                calc_schemes_ini=pass2calc['init_structure']
                print "check_point94"
            else:
                calc_schemes_ini=[pass2calc['init_structure']]
                print "check_point95"
            for calc_scheme_ini in calc_schemes_ini:
                #print "setup_calc:",uid,calc_scheme_ini
                if not calc_scheme_ini in self.calc_schemes:
                    continue
                calcdir=os.path.join(self.structureDB[uid].calcdir,calc_scheme_ini)
                print "check_point96"
                if not (calcdir in self.structureDB[uid].submitted_jobs):
                    calcdir=os.path.join(self.structureDB[uid].calcdir,calc_scheme_ini,"")
                    print "check_point97"
                tmpatoms=self.get_atoms_object(uid,calc_scheme_ini, magsettings={}) #exclude AF structures here
                print "check_point98, tmpatoms is:", tmpatoms
                if (tmpatoms!=None):
                    print "check_point99"
                    break       
                elif (calcdir in self.structureDB[uid].submitted_jobs) and ((self.structureDB[uid].submitted_jobs[calcdir]['nsubmit']>=self.nsub_max) or (self.structureDB[uid].submitted_jobs[calcdir]['nsubmit']<0)):
                    #print "Failed to get init_structure for ",uid,calc_scheme_ini
                    print "check_point100"
                    continue
                else:
                    self.get_energy(uid,calc_scheme_ini,update=update,nsub_max=self.nsub_max, sloppy_mode=False)
                    print "check_point101"
                    return return_failed
            if (tmpatoms!=None):
                self.structureDB[uid].atoms=tmpatoms
                del pass2calc['init_structure']
                print "check_point101"
            else:
                return return_failed
            if 'symprec' in pass2calc:
                try:
                    spglib_info=spglib.get_symmetry_dataset(tmpatoms, symprec=pass2calc['symprec'])
                    self.structureDB[uid].atoms=Atoms(numbers=spglib_info['std_types'],cell=spglib_info['std_lattice'],scaled_positions=spglib_info['std_positions'],pbc=True)
                except:
                    print "WARNING(setup_calculator): could not determine symmetry for uid ",uid," calc scheme=",calc_scheme
                    return None
        if ('ispin' in pass2calc) and (pass2calc['ispin']=='auto'):
            print "check_point131,('ispin' in pass2calc) and (pass2calc['ispin']=='auto')"
            pass2calc['ispin']=1
            magmom_ini=[]
            for el in self.structureDB[uid].atoms.get_chemical_symbols():
                if el in magatoms:
                    pass2calc['ispin']=2
                    magmom_ini.append(magatoms[el])
                else:
                    magmom_ini.append(0)
            if pass2calc['ispin']==2:
                pass2calc['magmom']=magmom_ini
            print "check_point132, magmom_ini is :",magmom_ini
            
        # if ('I_CONSTRAINED_M' in pass2calc):
        #     pass2calc['I_CONSTRAINED_M']=pass2calc['I_CONSTRAINED_M']

        print "check_point165, add I_CONSTRAINED_M into pass2calc:",pass2calc
        if ('scale_volume_mag' in pass2calc):
            if pass2calc['ispin']==2:
                pass2calc['scale_volume']=pass2calc['scale_volume_mag']
            del pass2calc['scale_volume_mag']
        if ('scale_volume' in pass2calc) and (update_atoms==True): 
            if pass2calc['scale_volume']=='atomic_volume':
                comp=self.structureDB[uid].get_composition(reduce=False)
                scaled_volume=0.0
                for element in comp:
                    if (scaled_volume!=None) and (self.get_atomic_volume(element)!=None):
                        scaled_volume=scaled_volume+comp[element]*self.get_atomic_volume(element)
                    else:
                        scaled_volume=None
                if (scaled_volume!=None):
                    scale=pow(scaled_volume/self.structureDB[uid].atoms.get_volume(),1./3)
                    self.structureDB[uid].atoms.set_cell(scale*self.structureDB[uid].atoms.get_cell(), scale_atoms=True)
                del pass2calc['scale_volume']
            elif isinstance(pass2calc['scale_volume'],float):
                scale=pow(pass2calc['scale_volume'],1./3.)
                self.structureDB[uid].atoms.set_cell(scale*self.structureDB[uid].atoms.get_cell(), scale_atoms=True)
                del pass2calc['scale_volume']    
        if afm!={}:
            print "check_point133, afm!={}",afm
            if 'ispin' in afm:
                if ('ispin' in pass2calc) and (pass2calc['ispin']==afm['ispin']):
                    print "check_point134, ('ispin' in pass2calc) and (pass2calc['ispin']==afm['ispin']) retrurn failed"
                    return return_failed # don't make double calculation of default
                pass2calc['ispin']=afm['ispin']
                #if (afm['ispin']==1) and ('magmom' in pass2calc):
                if ('magmom' in pass2calc):
                    del pass2calc['magmom']
            if 'magmom' in afm:
                pass2calc['magmom']=afm['magmom']
                pass2calc['ispin']=2
            if 'lnoncollinear' in afm:
                pass2calc['lnoncollinear']=afm['lnoncollinear']
            if 'atoms_object' in afm:
                self.structureDB[uid].atoms=afm['atoms_object']
        #vasp specific
        if calculator_name.lower()=='vasp':
            print "check_point168,calculator_name.lower()=='vasp' "
            # check if PP files are present
            comp=self.structureDB[uid].get_composition(reduce=False)
            for element in comp:
                print "check_point173,get_composition",element
                if self.check_vasp_pp_files(element, pass2calc)==False:
                    calcdir=os.path.join(self.structureDB[uid].calcdir,calc_scheme)
                    self.structureDB[uid].submitted_jobs[calcdir]={'nsubmit':-1,'jobid':-1}
                    print "check_point135, check_vasp_pp_files False,self.structureDB[uid].submitted_jobs[calcdir] is ",self.structureDB[uid].submitted_jobs[calcdir]
                    self.add_logmessage("WARNING(setup_calculator): Could not find PP files for %s and calc_scheme %s, check setups!"%(element,calc_scheme))
                    return return_failed
            # preprocess some arguments before sending them to ase-vasp interface:
            if ('ispin' in pass2calc) and (pass2calc['ispin']==2):
                print "check_point138,('ispin' in pass2calc) and (pass2calc['ispin']==2):"
                if not ('lorbit' in pass2calc):
                    pass2calc['lorbit']=11 #seems to be only way to get local atomic moments
            if ('lambda1' in pass2calc):
                pass2calc['lambda']=pass2calc['lambda1']
                del pass2calc['lambda1']
            # add this because lambda is a keyword in python, could not directly use this.
            if ('i_constrained_m' in pass2calc) and (pass2calc['lnoncollinear']==True):
                pass2calc['lorbit']=1
                pass2calc['m_constr']=pass2calc['magmom']
                
                rwigs_list = []
                elpot=""
                chem_sym_list = self.structureDB[uid].atoms.get_chemical_symbols()
                print "check_point175, chem_sym_dict is", chem_sym_list
                for el in chem_sym_list:
                    if el!=elpot:
                        elpot=el
                        isok,pp_dict=self.check_vasp_pp_files(el,settings,return_dict=True)
                        print "check_point174, pp_dict is", pp_dict
                        rwigs_list.append(pp_dict['RWIGS'])
                print  "check_point176, rwigs_list  is ", rwigs_list
                pass2calc['rwigs'] = rwigs_list
            #     for el in self.structureDB[uid].atoms.get_chemical_symbols():
            #         print "check_point172, each element", el
            # print "check_point170, add pass2calc['lambda'] and show pass2calc:", pass2calc    
            # if ('I_CONSTRAINED_M' in pass2calc):
            #     pass2calc['I_CONSTRAINED_M']=pass2calc['I_CONSTRAINED_M']
            #     del pass2calc['I_CONSTRAINED_M']
            # if ('LAMBDA' in pass2calc):
            #     pass2calc['LAMBDA']=pass2calc['LAMBDA']
            #     del pass2calc['LAMBDA']
            if ('LSDA_U' in pass2calc):
                lsdau=False
                ldauL=[]
                ldauU=[]
                ldauJ=[]
                el0=""
                for el in self.structureDB[uid].atoms.get_chemical_symbols():
                    if el==el0:
                        continue
                    else:
                        el0=el
                    if el in pass2calc['LSDA_U']:
                        lsdau=True
                        uset=pass2calc['LSDA_U'][el]
                        ldauL.append(uset['L'])
                        ldauU.append(uset['U'])
                        ldauJ.append(uset['J'])
                    else:
                        ldauL.append(-1)
                        ldauU.append(0.0)
                        ldauJ.append(0.0)
                if lsdau==True:
                    pass2calc['ldau']=True
                    if not ('ldautype' in pass2calc):
                        pass2calc['ldautype']=1
                    pass2calc['ldaul']=ldauL
                    pass2calc['ldauu']=ldauU
                    pass2calc['ldauj']=ldauJ
                del pass2calc['LSDA_U']
                #print "LSDA+U",pass2calc
                #print "TODO: LSDA+U"
                #bla    
            if (transport==True) or (bandstructure==True) or (bandstructure_init==True):
                #change settings for electronic structure:
                if ('nsw' in pass2calc) and (pass2calc['nsw']>0):
                    self.structureDB[uid].atoms=self.get_atoms_object(uid,calc_scheme=calc_scheme)
                    if self.structureDB[uid].atoms==None:
                        return None
                pass2calc['nsw']=0 # do not re-optimize
                if 'ibrion' in  pass2calc:
                    del pass2calc['ibrion']
                if 'isif' in  pass2calc:
                    del pass2calc['isif']
                if (transport==True) or (bandstructure_init==True):
                    pass2calc['ismear']=-5 # use tetrahedron method
                elif (bandstructure==True) and ('ismear' in pass2calc):
                    del pass2calc['ismear']
                pass2calc['lcharg']=True # write CHARGE file
                if bandstructure or transport:
                    pass2calc['icharg']=11
                if transport_kspace_density!=None:
                    pass2calc['kspace_density']=transport_kspace_density
            if 'kspace_density' in pass2calc:
                # adopt k-mesh to structure (AUTO len not possible in ase)
                kdens=pass2calc['kspace_density']
                rc=self.structureDB[uid].atoms.get_reciprocal_cell()
                N=[1,1,1]
                for i in range(rc.shape[0]):
                    x2=0.0
                    for j in range(rc.shape[1]):
                        x2=x2+rc[i][j]*rc[i][j]
                    N[i]=int(ceil(kdens*sqrt(x2)))
                del pass2calc['kspace_density']
                pass2calc['kpts']=N
            if elastic==True:
                if 'isif' in  pass2calc: #no reoptimization for the moment
                    del pass2calc['isif']
                if 'ibrion' in  pass2calc:
                    del pass2calc['ibrion']
                pass2calc['nsw']=0
            print "check_point169, before Vasp(**pass2calc)",pass2calc
            calc=Vasp(**pass2calc)
            print "check_point167, calc is :", calc
            #print pass2calc
        elif (calculator_name.lower()=='fplo') and (self.get_scratch_directory()!=None):
            if 'subdir' in afm:
                subdir=afm['subdir']
            else:
                subdir=""
            #if second run is needed increase number of iterations and choose LCIterat
            if (os.path.isfile(os.path.join(self.get_scratch_directory(),uid,calc_scheme,subdir,"=.dens"))) and (update_atoms==True):
                pass2calc['niter']=100
                pass2calc['iterat_version']=4
            if ('ispin' in pass2calc) and (pass2calc['ispin']==2) and (update_atoms==True):
                #test if spinpolarized calculation was already set up
                if FPLO_densfile_is_spinpolarized(os.path.join(self.get_scratch_directory(),uid,calc_scheme,subdir)):
                    pass2calc['initial_spin_split']='f'
                else:
                    grit=FPLO().read_convergence_level(calcdir=os.path.join(self.get_scratch_directory(),uid,calc_scheme,subdir,'ini_nm'))
                    if (grit!=None) and (grit<1.e-3):
                        inipath=os.path.join(self.get_scratch_directory(),uid,calc_scheme,subdir,'ini_nm')
                    else:
                        inipath=self.get_converged_calculation_path(uid,calc_scheme,special=os.path.join(subdir,'ini_nm'))
                    if inipath!=None:
                        pass2calc['initial_spin_split']='t'
                        pass2calc['initdir']=inipath
                        syminfo=FPLO().read_input_file(calcdir=inipath)
                        magmom_ini=[]
                        for el in syminfo['wyckoff_elements']:
                            if el in magatoms:
                                magmom_ini.append(magatoms[el])
                            else:
                                magmom_ini.append(0)
                        if not ('magmom' in afm):
                            pass2calc['initial_spin_split_sorts']=magmom_ini
                    else:
                        pass2calc['ispin']=1
                        pass2calc['initial_spin_split']='f'
                        pass2calc['subdir']='ini_nm'
            job_settings=self.get_job_commands(calc_scheme)
            if (job_settings!=None) and ('fedit' in job_settings):
                calc=FPLO(fedit=job_settings['fedit'],**pass2calc)
            else:
                print 'setup_calculator(): need fedit for FPLO()'
                calc=None
            if ('lnoncollinear' in afm) and (afm['lnoncollinear']==True):
                calc=None
        else:
            calc=None
            print 'setup_calculator:', calculator_name,'not yet implemented in HTE, nothing done!'
        if return_settings==True:
            print "check_point136,returned calc,pass2calc :",calc,pass2calc
            return calc,pass2calc
        return calc
    
                
    def delete_calculations(self, structure_list=None, calc_schemes=None):
        if (structure_list==None) or (calc_schemes==None):
            print 'No structures or calc_scheme selected, nothing done!'
            return None
        else:
            del_folders=[]
            for uid in structure_list:
                self.structureDB[uid].stored_calc_results={} #todo: only deleted
                for cs in calc_schemes:
                    folder=os.path.join(uid,cs)
                    if folder in self.structureDB[uid].submitted_jobs:
                        del self.structureDB[uid].submitted_jobs[folder]
                        print 'removed history of ',folder
                    if os.path.isdir(folder):
                        del_folders.append(folder)
            print 'List of folders to be deleted:'
            for folder in del_folders:
                print folder
            if raw_input('Marked %d folders for deletion, shall I continue (yes/no) ? ' %len(del_folders))=='yes':
                for folder in del_folders:
                    shutil.rmtree(folder)
                print 'Done!'
            else:
                print 'Nothing done!'


    def get_atoms_object(self, uid, calc_scheme=None, magsettings={'submitted':True}, sub_directories={}):
        """returns the atoms object of database entry with uid
        default (no calc_scheme given): initial, unconverged structure
        calc_scheme: if given, return converged structure or None if
                     calculation is not converged
        """
        parentdir=os.getcwd()
        if uid in self.structureDB:
            if calc_scheme==None:
                print "check_point29, in get_atoms_object, calc_schmem is none,return init_ao", self.structureDB[uid].atoms_initial
                return self.structureDB[uid].atoms_initial
            if calc_scheme in self.calc_schemes:
                print "check_point30,calc_scheme is existed in calc_schemes"
                if (self.use_prop_dict==True):
                    print "check_point31, in get_atoms_object, use_prop_dict=True"
                    mags=deepcopy(magsettings)
                    print "check_point32, in get_atoms_object, mags_setting is:",mags
                    print "check_point33, in get_atoms_object, sub_directories are:",sub_directories
                    mags['get_atoms']=False  #for the moment (recursion problem)
                    prop_dict=self.get_properties(uid, calc_scheme, magsettings=mags, sub_directories=sub_directories)
                    print "check_point34, in get_atoms_object, prop_dict is :", prop_dict
                    if ('chemical_symbols' in prop_dict) and ('cell' in prop_dict) and ('scaled_positions' in prop_dict):
                        print "check_point35, in get_atoms_object, pd_dict have scaled_positions"
                        atoms=Atoms(prop_dict['chemical_symbols'],cell=prop_dict['cell'],scaled_positions=prop_dict['scaled_positions'],pbc=True)
                        return atoms
                    elif ('chemical_symbols' in prop_dict) and ('cell' in prop_dict) and ('positions' in prop_dict):
                        print "check_point36, in get_atoms_object, pd_dict have positions"
                        atoms=Atoms(prop_dict['chemical_symbols'],cell=prop_dict['cell'],positions=prop_dict['positions'],pbc=True)
                        return atoms
                    else:
                        print "check_point37, this prop_dict is not include position"
                    return None
                #check if lower part is still necessary
                for sp in self.get_searchpaths():
                    try:
                        print "check_point38, change dir to sp" 
                        os.chdir(sp)
                    except OSError, e:
                        warnings.warn("Searchpaths: OSError {0}".format(e))
                        continue
                    calc=self.setup_calculator(uid,calc_scheme,basic=True)
                    print "check_point67, here is get_atoms_object and set the basic==True"
                    atoms=self.structureDB[uid].get_converged_structure(calc_scheme,calc)
                    print "check_point39, get the ao in sp,and the atoms is:", atoms 
                    os.chdir(parentdir)
                    if atoms!=None:
                        return atoms
        return None

    def get_lattice_parameters(self, uid, calc_scheme=None, magsettings={'submitted':True}):
        """returns lattice parameters [a,b,c] of database entry with uid
        default (no calc_scheme given): initial, unconverged structure
        calc_scheme: if given, return converged structure or None if
                     calculation is not converged
        """
        cell=None
        if calc_scheme==None:
            cell=self.get_atoms_object(uid).get_cell()
        else:
            prop_dict=self.get_properties(uid, calc_scheme, magsettings=magsettings)
            if 'cell' in prop_dict:
                cell=prop_dict['cell']
        if isinstance(cell,type(None)):
            return None
        return [norm(cell[0]),norm(cell[1]),norm(cell[2])]
                
            
            
    def get_converged_calculation_path(self, uid, calc_scheme,special='', subdir=''):
        """returns a path where a converged calculation for the given database uid and
        calc_scheme can be found."""
        calcpath=None
        calc_scheme_sp=calc_scheme
        print "check_point70, here is get_converged_calculation_path: and calc_scheme_sp is:", calc_scheme_sp
        if subdir!='':
            calc_scheme_sp=subdir
            print "check_point71, here is get_converged_calculation_path: and subdir!='', calc_scheme_sp is:", calc_scheme_sp
        elif special!='':
            calc_scheme_sp=os.path.join(calc_scheme,special)
            print "check_point72, special!='', calc_scheme_sp is:",calc_scheme_sp
        if calc_scheme in self.calc_schemes:
            parentdir=os.getcwd()
            calc=self.setup_calculator(uid,calc_scheme,basic=True)
            print "check_point68, here is get_converged_calculation_path and set the basic=True.parentdir is:",parentdir
            for sp in self.get_searchpaths():
                print "check_point69, here is self.get_searchpaths(): and sp is:",sp
                os.chdir(sp)
                dir_test = os.getcwd()
                print "check_point74, here we change the dir now in:",dir_test
                print "check_point75, before check the convergence, first we print the calc_scheme_sp: ", calc_scheme_sp, "calc is:",calc
                if self.structureDB[uid].check_convergency(calc_scheme_sp,calc):
                    print "check_point76, here the structureDB[uid].check_convergency is Ture:"
                    calcpath=os.path.join(sp,uid,calc_scheme_sp)
                os.chdir(parentdir)
                if calcpath!=None:
                    break
        print "check_point73, here is the calcpath:", calcpath
        return calcpath

    def get_effective_calc_scheme(self, uid, calc_scheme):
        """returns name of calculation scheme to be used for compound with unique identifier uid
        can be used to combine calculations which require special treatment for some elements (e.g. rare earths)
        with standard calculations e.g. for convex hull calculations
        usage:
        """
        if not (calc_scheme in self.calc_schemes):
            self.add_logmessage("WARNING(get_effective_calc_scheme): unknown calc_scheme %s, nothing done!"%calc_scheme)
            return None,None,None
        calculator_name, settings=self.calc_schemes[calc_scheme]
        #print calculator_name, settings
        if not (uid in self.structureDB):
            self.add_logmessage("WARNING(get_effective_calc_scheme): unknown uid %s, nothing done!"%uid)
            return calc_scheme,calculator_name, settings
        comp=self.structureDB[uid].get_composition(reduce=False)
        if 'default' in settings:
            eff_csname=settings['default']
            calculator_name, settings_default=self.calc_schemes[eff_csname]
            eff_settings=deepcopy(settings_default)
            for el in sorted(comp):
                if el in settings:
                    eff_csname=eff_csname+"_%s"%el+settings[el]['csext']
                    for tag in settings[el]['settings']:
                        if tag in eff_settings:
                            eff_settings[tag][el]=settings[el]['settings'][tag]
                        else:
                            eff_settings[tag]={el:settings[el]['settings'][tag]}
            if (eff_csname!=settings['default']) and ("general" in settings):
                for tag in settings['general']:
                    if (settings['general'][tag]==None) and (tag in eff_settings):
                        del eff_settings[tag]
                    else:
                        eff_settings[tag]=settings['general'][tag]
            return eff_csname,calculator_name,eff_settings
        return calc_scheme,calculator_name, settings
                    
                            
                
        #if 'if' in settings:
        #    for case in settings['if']:
                
                
                
        

    def get_properties(self, uid, calc_scheme, magsettings={'submitted':True}, sub_directories={}, reference_prop='energy', update=False, storage_options='default', nsub_max=-1):
        """returns a dictionary with properties if a converged calculation exists"""
        # searchpaths and storage options
        print "check_point1: start get_properties"
        searchpaths=self.get_storage_directories()+self.get_searchpaths()
        if storage_options=='default':
            st_files=self.get_storage_options(calc_scheme)
        else:
            st_files=storage_options
        if nsub_max<0:
            nsub_max=self.nsub_max
        prop_dict={}
        print "check_point2: check pd, should be empty dictionary", prop_dict
        if not uid in self.structureDB:
            return prop_dict
        if not (uid in self.tmpdata['prop_dict']):
            self.tmpdata['prop_dict'][uid]={}
        # subdirectories with AF calculations etc.
        subdirs={}
        print "check_point93, self.tmpdata is:",self.tmpdata
        if sub_directories!={}:
            subdirs=sub_directories
            print "check_point3, sub_directories is not empty, so what is seems like", sub_directories
        else:
            magconfigs=self.setup_magnetic_structures(uid,calc_scheme, magsettings=magsettings)
            print "check_point4, sub_directories is empty!! the magconfigs from the setup_magnetic_structures function:", magconfigs
            for magconf in magconfigs:
                print "check_point5, here is a cycle, out put of each magconfigs subdirs,", magconf
                subdir=os.path.join(calc_scheme,magconf)
                subdirs[subdir]=magconfigs[magconf]
        for subdir in subdirs:
            if subdir in self.tmpdata['prop_dict'][uid]:
                print "check_point6, subdir in self.tmpdata, and the subdir is : ", subdir
                prop_dict[subdir]=self.tmpdata['prop_dict'][uid][subdir]
                print "check_point40, subdir in self.tmpdata and prop_dict[subdir]", prop_dict[subdir]
                continue
            prop_dict[subdir]={}
            print "check_point7, search path is : ", searchpaths
            for sp in searchpaths:
                prop_file=os.path.join(sp,self.structureDB[uid].calcdir,subdir,'hte_propdict.txt')
                print "check_point8, prop_file is : ", prop_file
                if os.path.isfile(prop_file):
                    try:
                        infile=open(prop_file)
                        line=infile.readline()
                        while line:
                            prop,val=string2prop_dict(line)
                            prop_dict[subdir][prop]=val
                            line=infile.readline()
                        infile.close()
                        prop_dict[subdir]['data_path']=os.path.join(sp,self.structureDB[uid].calcdir,subdir)
                        print "check_point86,prop_dict[subdir]['data_path'] is : ", prop_dict[subdir]['data_path']
                    except:
                        print "WARNING(get_properties()): Failed to read ",prop_file
                        prop_dict[subdir]={}
                    print "check_point87,the dict prop_dict[subdir] is : ",prop_dict[subdir]                
                if (not ('hte_version' in prop_dict[subdir])) or (('errors' in  prop_dict[subdir]) and (prop_dict[subdir]['errors']!='False')): # or (float(prop_dict[subdir]['hte_version'].split('-')[0])<1.2):
                    prop_dict[subdir]={}
                if 'initial_magnetization' in prop_dict[subdir]: #tmp version
                    print "check_point85, prop_dict[subdir] have initial_magnetization "
                    prop_dict[subdir]={}
                if (not (reference_prop in prop_dict[subdir])):
                    prop_dict[subdir]={}
                if (prop_dict[subdir]!={}):
                    print "check_point89, self.get_main_storage_directory() is :",self.get_main_storage_directory()
                    print "check_point90, get_storage_directories() is :", self.get_storage_directories()
                    print "check_point91, sp is :", sp
                    if (self.get_main_storage_directory()!=None) and (not (sp in self.get_storage_directories())):
                        print "check_point88,self.get_main_storage_directory()!=None) and (not (sp in self.get_storage_directories())"
                        store_files(os.path.join(sp,self.structureDB[uid].calcdir,subdir),os.path.join(self.get_main_storage_directory(),self.structureDB[uid].calcdir,subdir),files=st_files)
                    break
            if prop_dict[subdir]=={}:
                #check if path with converged calculation exists
                calcpath=self.get_converged_calculation_path(uid, calc_scheme, subdir=subdir)
                print "check_point81, here prop_dict={}, first check_convergence and get the calcpath", calcpath
                if calcpath!=None:
                    print "get_properties(): converged calculation in ",calcpath
                    calc_name,settings=self.calc_schemes[calc_scheme]
                    if calc_name.lower()=='vasp':
                        pdict=get_properties_vasp(calcdir=calcpath)
                        print "check_point82, get_properties_vasp, the pdict is:",pdict
                    elif calc_name.lower()=='fplo':
                        pdict=get_properties_fplo(calcdir=calcpath)
                    else:
                        pdict={}
                    if ('errors' in pdict) and (pdict['errors']=='True'):
                        self.add_logmessage("WARNING(get_properties): Check <%s>, %s!"%(calcpath,pdict['errors']))
                        pdict={}
                    if (not (reference_prop in pdict)):
                        pdict={}
                    if pdict!={}:
                        lines=["hte_version = %s"%self.get_version()]
                        for prop in sorted(pdict):
                            line="%s = "%prop
                            if isinstance(pdict[prop],list):
                                for i in range(len(pdict[prop])):
                                    val=pdict[prop][i]
                                    if isinstance(val,list):
                                        line=line+"%s"%val[0]
                                        for j in range(1,len(val)):
                                            line=line+";%s"%val[j]
                                    else:
                                        line=line+"%s"%val
                                    if (i<len(pdict[prop])-1):
                                        line=line+","
                            else:
                                line=line+"%s"%pdict[prop]
                            lines.append(line)
                        for line in lines:
                            prop,val=string2prop_dict(line)
                            prop_dict[subdir][prop]=val
                        #store in self.storage_directories?
                        if (prop_dict[subdir]!={}):
                            if (self.get_scratch_directory()!=None):
                                if os.path.abspath(calcpath).startswith(self.get_scratch_directory()):
                                    outfile=open(os.path.join(calcpath,'hte_propdict.txt'),"w")
                                    for line in lines:
                                        outfile.write("%s\n"%line)
                                    outfile.close()
                            if (self.get_main_storage_directory()!=None):
                                stdir=os.path.join(self.get_main_storage_directory(),self.structureDB[uid].calcdir,subdir)
                                store_files(calcpath,stdir,files=st_files)
                                outfile=open(os.path.join(stdir,'hte_propdict.txt'),"w")
                                for line in lines:
                                    outfile.write("%s\n"%line)
                                outfile.close()
                                prop_dict[subdir]['data_path']=stdir
                            else:
                                prop_dict[subdir]['data_path']=calcpath
            self.tmpdata['prop_dict'][uid][subdir]=prop_dict[subdir]
        print "check_point83, prop_dict[subdir] = self.tmpdata and they are:",self.tmpdata['prop_dict'][uid][subdir]
        pd_ref={}
        print "check_point9: print all subdirs:", subdirs
        print "check_point140: print magsettings:", magsettings
        for subdir in subdirs:
            print "check_point143, entering subdir:", subdir
            print "check_point145, prop_dict[subdir]", prop_dict[subdir]
            if (prop_dict[subdir]!={}) and (not ('updating' in prop_dict[subdir])):
                print "check_point137: (prop_dict[subdir]==!{}) and (update==True), prop_dict[subdir]:",prop_dict[subdir]
                #tmp solution if DB file was not stored
                print "check_point146, self.structureDB[uid].submitted_jobs",self.structureDB[uid].submitted_jobs
                if not os.path.join(uid,subdir) in self.structureDB[uid].submitted_jobs:
                    print "check_point144, entering submitted_jobs:",self.structureDB[uid].submitted_jobs
                    self.structureDB[uid].submitted_jobs[os.path.join(uid,subdir)]={'jobid':'unknown','nsubmit':0, 'settings':subdirs[subdir]}
                    if subdir.strip("/")!=calc_scheme:
                        self.add_logmessage("WARNING(get_properties): Submission history lost for %s %s, pd=%s"%(uid,subdir,str(prop_dict[subdir])))
                if (isinstance(magsettings,dict)) and ('non_collinear' in magsettings) and (magsettings['non_collinear']==False):
                    if ('magnetic_moments' in prop_dict[subdir]) and (isinstance(prop_dict[subdir]['magnetic_moments'][0],list)):
                        print "check_point141:('magnetic_moments' in prop_dict[subdir]), collinear, jump out of the cycle"
                        continue
                elif (isinstance(magsettings,dict)) and ('non_collinear' in magsettings) and (magsettings['non_collinear']==True):
                    if (not ('magnetic_moments' in prop_dict[subdir])) or (isinstance(prop_dict[subdir]['magnetic_moments'][0],float)):
                        print "check_point142:('magnetic_moments' in prop_dict[subdir]), noncollinear, jump out of the cycle"
                        continue
                if (reference_prop in prop_dict[subdir]) and ('chemical_symbols' in prop_dict[subdir]):
                    if (pd_ref=={}) or (prop_dict[subdir][reference_prop]/float(len(prop_dict[subdir]['chemical_symbols']))<pd_ref[reference_prop]/float(len(pd_ref['chemical_symbols']))):
                        pd_ref=prop_dict[subdir]
                        if 'subdir' in subdirs[subdir]:
                            pd_ref['subdir']=subdirs[subdir]['subdir']
            if (prop_dict[subdir]=={}) and (update==True):
                print "check_point109: (prop_dict[subdir]=={}) and (update==True)"
                print "check_point110 uid,calc_scheme,  ", uid, calc_scheme
                try:
                    self.run_calculation(uid,calc_scheme, sub_directories={subdir:subdirs[subdir]}, nsub_max=nsub_max)
                    print "check_point113 after run_calculation"
                except:
                    self.add_logmessage("WARNING(get_properties): Failed to run uid %s (%s) in %s!"%(uid,calc_scheme,subdir))
                self.tmpdata['prop_dict'][uid][subdir]={'updating':True}
                print "check_point112 self.tmpdata['prop_dict'][uid][subdir]:", self.tmpdata['prop_dict'][uid][subdir]
                ### if (prop_dict[subdir]=={}) and (update==True):
        if (not ('return_all' in magsettings)) or (magsettings['return_all']==False):
            print "check_point114 (not ('return_all' in magsettings)) or (magsettings['return_all']==False)"
            prop_dict=pd_ref
            
        print "check_point10: the final prop_dict for each uid:", prop_dict
        return prop_dict

    def restore_calculation(self, uid, calc_scheme, target_dir="", magsettings={'submitted':True}, uncompress=True):
        prop_dict=self.get_properties(uid, calc_scheme, magsettings=magsettings)
        if 'data_path' in prop_dict:
            try:
                print prop_dict['data_path'],target_dir
                if not (os.path.isdir(target_dir)):
                    exitcode, out = commands.getstatusoutput("mkdir -p %s"%target_dir)
                exitcode, out = commands.getstatusoutput("cp %s/* %s/"%(prop_dict['data_path'],target_dir))
                if uncompress==True:
                    for fname in glob.glob("%s/*.bz2"%target_dir):
                        exitcode, out = commands.getstatusoutput("bunzip2 %s"%fname)
                    for fname in glob.glob("%s/*.gz"%target_dir):
                        exitcode, out = commands.getstatusoutput("gunzip %s"%fname)
                calculator_name, settings=self.calc_schemes[calc_scheme]
                if (calculator_name.lower()=='vasp') and (not (os.path.isfile("%s/POTCAR"))):
                    elpot=""
                    commandline=""
                    pipe=" > "
                    for el in prop_dict['chemical_symbols']:
                        if el!=elpot:
                            elpot=el
                            isok,pp_dict=self.check_vasp_pp_files(el,settings,return_dict=True)
                            print "check_point171: pp_dict is:",pp_dict
                            if 'filename' in pp_dict:
                                if pp_dict['filename'].endswith(".Z"):
                                    command="gunzip -c %s"%pp_dict['filename']
                                else:
                                    command="cat %s"%pp_dict['filename']
                                if commandline=="":
                                    commandline="%s > %s/POTCAR "%(command,target_dir)
                                else:
                                    commandline=commandline+"; %s >> %s/POTCAR "%(command,target_dir)
                            else:
                                self.add_logmessage("WARNING(restore_calculation): Failed to restore %s %s in %s"%(uid, calc_scheme, target_dir))
                    print commandline
                    os.system(commandline)
                return True
            except:
                self.add_logmessage("WARNING(restore_calculation): Failed to restore %s %s in %s"%(uid, calc_scheme, target_dir))
        return False
        
    def get_converged_calculation_paths(self, uid, calc_scheme,property=None, **kwargs):
        """Returns a list with paths where converged calculations for the given uid and
        calc_scheme can be found. Optional arguments:
        property=None: default (energy) calculation
        property='transport': Boltztrap transport calculations, minimal k-space density
                              for transport calculation may be given as additional
                              keyword, e.g. 'ksdmin=80'; paths are sorted by ascending
                              k-space density.
        """
        calcpaths=[]
        tcps=[]
        ksdmin=0
        if 'ksdmin' in kwargs:
            ksdmin=kwargs['ksdmin']
        if calc_scheme in self.calc_schemes:
            parentdir=os.getcwd()
            calc=self.setup_calculator(uid,calc_scheme)
            for sp in self.get_searchpaths():
                os.chdir(sp)
                if self.structureDB[uid].check_convergency(calc_scheme,calc):
                    if property==None:
                        calcpaths.append(os.path.join(sp,uid,calc_scheme))
                    elif property=='transport':
                        tpaths=glob.glob(os.path.join(sp,uid,calc_scheme,'transp-k*'))
                        for tpath in tpaths:
                            ksd=int(tpath.split('transp-k')[1])
                            if (ksd>=ksdmin) and (os.path.isfile(os.path.join(tpath,'boltzgen','hte.trace'))):
                                tcps.append((ksd,tpath))
                os.chdir(parentdir)
            if property=='transport':
                #return sorted array of transport calculations (best first):
                for (ksd,tpath) in sorted(tcps,reverse=True):
                    calcpaths.append(tpath)
        return calcpaths

    def get_band_MAE(self, uid, calc_scheme, afm_dir="",magsettings={}, qaxis=[[0,0,1],[1,0,0]], k_mesh=[], update=False, nsub_max=2, unit='meV/cell'):
        """Get magnetic anisotropy energy with magnetic force theorem for given uid and calc_scheme (only FPLO supported!)
        Default is MAE for ferro magnetic spin alignment, options:
            afm_dir: specify subdirectory with (scalarrelativistic) AFM calculations
            magsettings: Dictionary to specify spin structure, e.g. {'submitted': True} to calculate MAE for spin structure
            with lowest energy (of the ones which have been submitted)
        """
        calcname=""
        if calc_scheme in self.calc_schemes:
            calcname, settings=self.calc_schemes[calc_scheme]
        if calcname.lower()!='fplo':
            self.add_logmessage("WARNING(get_band_MAE): Only FPLO supported, check calc_scheme!")
            return None
        edict={}
        subdir=calc_scheme
        if afm_dir!="":
            subdir=os.path.join(calc_scheme,afm_dir)
        elif magsettings!={}:
            prop_dict=self.get_properties(uid, calc_scheme,magsettings=magsettings,update=update)
            if not ('energy' in prop_dict):
                return None
            if 'subdir' in prop_dict:
                subdir=os.path.join(calc_scheme,prop_dict['subdir'])
        calcpathini=None
        for qax in qaxis:
            qaxdir=os.path.join(subdir,"MAEBAND_%d_%d_%d"%(qax[0],qax[1],qax[2]))
            if k_mesh!=[]:
                qaxdir=os.path.join(subdir,"MAEBAND_%d_%d_%d_kp_%d_%d_%d"%(qax[0],qax[1],qax[2],k_mesh[0],k_mesh[1],k_mesh[2]))
            prop_dict=self.get_properties(uid, calc_scheme, sub_directories={qaxdir:{}}, reference_prop='band_energy', storage_options={'bzip2':['out'],'cp':['=.in','hte_propdict.txt']})
            if 'band_energy' in prop_dict:
                edict[str(qax)]=prop_dict['band_energy']
                natom=len(prop_dict['chemical_symbols'])
            elif (update==True) and (self.get_scratch_directory()!=None):
                if calcpathini==None:
                    calcpathini=self.get_converged_calculation_path(uid, calc_scheme,subdir=subdir)
                if calcpathini==None:
                    if self.restore_calculation(uid, calc_scheme,magsettings=magsettings,target_dir=os.path.join(self.get_scratch_directory(),uid,subdir))==True:
                        calcpathini=os.path.join(self.get_scratch_directory(),uid,subdir)
                    else:
                        return None
                noccu=FPLO().estimate_occupied_bands(calcdir=calcpathini)
                subdirmae=os.path.join(self.get_scratch_directory(),uid,qaxdir)
                if k_mesh==[]:
                    calc=FPLO(fedit=self.get_job_commands(calc_scheme)['fedit'], initdir=calcpathini,update_symmetry=False,use_str=False,ispin=2,initial_spin_split='f',niter=1,relativistic="fullrelat",quantization_axis=qax,occupied_bands=noccu, xc_vers=None)
                else:
                    calc=FPLO(fedit=self.get_job_commands(calc_scheme)['fedit'], initdir=calcpathini,update_symmetry=False,use_str=False,ispin=2,initial_spin_split='f',niter=1,relativistic="fullrelat",quantization_axis=qax,occupied_bands=noccu, xc_vers=None,kmesh=k_mesh)
                self.structureDB[uid].atoms=self.structureDB[uid].atoms_initial
                self.structureDB[uid].run_calculation(subdirmae,calc, job_commands=self.get_job_commands(calc_scheme),nsub_max=nsub_max)
        emin=None
        emax=None
        for qax in qaxis:
            if str(qax) in edict:
                if (emin==None) or (edict[str(qax)]<emin):
                    emin=edict[str(qax)]
                if (emax==None) or (edict[str(qax)]>emax):
                    emax=edict[str(qax)]
            else:
                return None
        if len(qaxis)==2: #signed value with respect to first axis given
            emax=edict[str(qaxis[0])]
            emin=edict[str(qaxis[1])]
        if unit=='meV/atom':
            return (emax-emin)*1000/natom #meV/atom
        return (emax-emin)*1000 #meV/cell
    
    
    def get_bandgap(self, uid, calc_scheme, update=False):
        """Returns an estimate of the band gap for the database entry with given uid.
        calc_scheme may be a single calculation scheme or a list of calculation schemes
        (then first entry != None is returned).
        The band gap is estimated with the k-point set of the self consistent calculation,
        check band structure for accurate value.
        Works only for FPLO and nonmagnetic Vasp calculations."""
        parentdir=os.getcwd()
        gap=None
        #check if uid is valid
        if uid in self.structureDB:
            dbentry=self.structureDB[uid]
        else:
            return gap
        #allow for multiple calc_schemes
        if isinstance(calc_scheme,list):
            calc_schemes=calc_scheme
        else:
            calc_schemes=[calc_scheme]
        #look if we can find some value for the band gap
        for calc_scheme in calc_schemes:
            if (calc_scheme in self.calc_schemes):
                #check for converged calculations in searchpaths:
                calc=self.setup_calculator(uid,calc_scheme)
                for sp in self.get_searchpaths():
                    os.chdir(sp)
                    try:
                        if (calc!=None) and (dbentry.check_convergency(calc_scheme,calc)):
                            gap=self.structureDB[uid].get_bandgap(calc_scheme,calc, update=False, job_commands=self.get_job_commands(calc_scheme))
                    except:
                        print "get_bandgap(): check calculation in %s/%s/%s"%(sp,uid,calc_scheme)
                    os.chdir(parentdir)
                    if gap!=None:
                        break
            if gap!=None:
                break
        if ((gap==None) and (update==True)):
            sp=self.searchpaths[0]
            calc_scheme=calc_schemes[0]
            os.chdir(sp)
            calc=self.setup_calculator(uid,calc_scheme)
            gap=self.structureDB[uid].get_bandgap(calc_scheme,calc, update=update, job_commands=self.get_job_commands(calc_scheme))
            os.chdir(parentdir)
        return gap

    def get_initial_magnetic_moment(self,el):
        if el in self.magnetic_atoms:
            return self.magnetic_atoms[el]
        return 0

    def setup_magnetic_structures(self,uid,calc_scheme, magsettings={'submitted':True}, max_configs=1000,silent=True, eps_mom=0.1, get_atoms=True, symprec=5e-2, AFatoms=['Mn'], forceMAG=[]):
        """returns a dictionary with magnetic configuration according to the options given in the dictionary 'settings', possible options are:
        settings={'submitted':True} (default): returns magnetic configurations which have been submitted before
        settings={}: only default FM/NM configuration
        settings={'init_structure':{'submitted':True}}: groundstate magnetic structure obtained in initial calc_scheme (e.g. if a Vasp optimization is continued with FPLO)
        settings={'sublattices','max_subgroups':[[0,0,0],[0,0,0.5]]}: AF coupling of sublattices and magnetic structures compatible with maximal subgroup symmetry, see examples
        max_configs: maximal number of spin configurations (only for 'max_subgroups' and 'sublattices'), can be overwritten by 'max_configs' in settings
        """
        settings=deepcopy(magsettings)
        print "check_point11, the start of setup_magnetic_structures, show magsetting:",settings
        AF_atoms=AFatoms
        if 'AFatoms' in settings:
            AF_atoms=settings['AFatoms']
        debug=False
        if ('debug' in settings) and (settings['debug']==True):
            debug=True
        if 'forceMAG' in settings: 
            forceMAG=settings['forceMAG']
        report_magnetic_structures=False
        if ('get_atoms' in settings) and (settings['get_atoms']==False): #TODO: clean 'get_atoms' in settings
            report_magnetic_structures=False
        elif ('report_magnetic_structures' in settings) and (settings['report_magnetic_structures']==True):
            report_magnetic_structures=True
        comp_r=self.structureDB[uid].get_composition(reduce=True)
        print "check_point12, comp_r is ", comp_r
        symdirhte=os.path.join(self.dir_of_hte,'verified-msg/')
        print "check_point13, symdirhte is :",symdirhte
        # determine symmetry
        if debug==True:
            ao_ini=self.get_atoms_object(uid)
            ao=self.get_atoms_object(uid,calc_scheme,magsettings={})
            print "check_point14, debug=True, ao_ini is:",ao_ini
            print "check_point15, debug=True, ao is:",ao
            for sympr in [symprec, 1e-4, 1e-3, 5e-2]:
                symspg=spglib.get_symmetry_dataset(ao, symprec=sympr)
                lattice_p, scaled_positions_p, numbers_p=spglib.find_primitive(ao, symprec=sympr)
                symspg_ini=spglib.get_symmetry_dataset(ao_ini, symprec=sympr)
                lattice_ini_p, scaled_positions_ini_p, numbers_ini_p=spglib.find_primitive(ao_ini, symprec=sympr)
                self.add_logmessage("  - SG(ao_ini) %s for symprec=%.1e primitive=%s"%(symspg_ini['number'],sympr,len(ao_ini)==len(scaled_positions_ini_p)))
                self.add_logmessage("  - SG(ao) %s for symprec=%.1e primitive=%s"%(symspg['number'],sympr,len(ao)==len(scaled_positions_p)))
                #if (symspg_ini['number']==symspg['number']):
                #    symprec=sympr
                #    #break
        # rm when tested
        #symdirhte='msg/'
        # print "check_point16, debug=False, ao is:",ao
        if ('auto' in settings):
            aset={}
            ao_ini=self.get_atoms_object(uid)
            ao=self.get_atoms_object(uid,calc_scheme,magsettings={})
            if ao==None:
                self.add_logmessage("* autoAF: initial FM calculation not ready for %s (%s/%s)"%(uid,comp_r,calc_scheme))
                print "check_point17, in auto, ao is None!"
                return {'':{}} 	
            symspg=spglib.get_symmetry_dataset(ao, symprec=symprec)
            symspg_ini=spglib.get_symmetry_dataset(ao_ini, symprec=symprec)
            print "check_point18, in auto, symspg is:", symspg
            print "check_point19, in auto, symspg_ini is:", symspg_ini
            if symspg['number']<symspg_ini['number']:
                self.add_logmessage("* autoAF: check symmetry of %s (%s/%s:%d/%d)"%(uid,comp_r,calc_scheme,symspg['number'],symspg_ini['number']))
                return {}
            ao_std=Atoms(numbers=symspg['std_types'],cell=symspg['std_lattice'],scaled_positions=symspg['std_positions'],pbc=True)
            afmconfigs=get_magnetic_sublattices(ao, symprec=symprec)
            aset['sublattices']=afmconfigs
            if (len(ao)!=len(ao_std)) and (len(ao_std)<30) and (not (symspg['number'] in [227,15])): #exclude some cases for the moment
                aset['maxmags']=[[0,0,0],[0,0,1]]
            else:
                comp=self.structureDB[uid].get_composition(reduce=False)
            if ('Mn' in comp) and (comp['Mn']==1):
                aset['maxmags']=[[0.5,0.5,0],[0,0,0],[0,0,0.5]]
            else:
                aset['maxmags']=[[0,0,0]]
            self.add_logmessage("* autoAF: %s %s"%(uid,str(aset)))
            settings['max_subgroups']=aset['maxmags']
            settings['sublattices']=True
            if len(afmconfigs)>5:
                settings['max_subgroups']={}
        if ('include_default' in settings) and (settings['include_default']==False):
            magconfigs={}
        else:
            magconfigs={'':{}} #always include default FM/NM calculation
        if ('default_aoini' in settings) and (settings['default_aoini']==True):
            #probe default calculation without pre-optimization
            ao_ini=self.get_atoms_object(uid)
            magmoms=[]
            for el in ao_ini.get_chemical_symbols():
                magmoms.append(self.get_initial_magnetic_moment(el))
            magconfigs['q_0_0_0_ao_ini']={'atoms_object':ao_ini,'magmom':magmoms}
        if ('default_aoini*' in settings) and (settings['default_aoini*']==True): #check
            for jobdir in self.structureDB[uid].submitted_jobs:
                if (jobdir.startswith(os.path.join(uid,calc_scheme,'q_0_0_0_ao_ini'))):
                    magconfigs[jobdir.split('/')[2]]=self.structureDB[uid].submitted_jobs[jobdir]['settings']
        if ('mcif_structure' in settings) and (settings['mcif_structure']==True) and ('mcif_structure' in self.structureDB[uid].cif_info):
            for mcifname in self.structureDB[uid].cif_info['mcif_structure']:
                magconfigs[mcifname]=self.structureDB[uid].cif_info['mcif_structure'][mcifname]
        if ('submitted' in settings) and (settings['submitted']==True):
            print "submitted jobs: ",self.structureDB[uid].submitted_jobs
            for jobdir in self.structureDB[uid].submitted_jobs:
                if ((jobdir.startswith(os.path.join(uid,calc_scheme,'afm_q'))) or (jobdir.startswith(os.path.join(uid,calc_scheme,'q_'))) or (jobdir.startswith(os.path.join(uid,calc_scheme,'ispin')))) and (len(jobdir.split('/'))==3):
                    name=jobdir.split('/')[2]
                    if 'settings' in self.structureDB[uid].submitted_jobs[jobdir]:
                        magconfigs[jobdir.split('/')[2]]=self.structureDB[uid].submitted_jobs[jobdir]['settings']
                    else:
                        magconfigs[jobdir.split('/')[2]]={}
                        self.add_logmessage("WARNING(setup_magnetic_structures): No settings found for %s, check magnetic structure!"%jobdir)
                    #if (report_magnetic_structures==True): 
                    #    self.add_logmessage("* %s not added for %s (already in, E=%s)."%(name,uid,str(self.get_energy_per_atom(uid,calc_scheme, sub_directories={os.path.join(calc_scheme,name):magconfigs[name]}))))
        
        ao=self.get_atoms_object(uid,magsettings={})
        # ao=self.get_atoms_object(uid,calc_scheme=calc_scheme,magsettings={})
        if ('max_subgroups' in settings):
            #symprec=1e-2
            print "check_point 24, begin_max_groups."
            if 'symprec' in settings:
                symprec=settings['symprec']
                symspg=spglib.get_symmetry_dataset(ao, symprec=settings['symprec'])
                ao=Atoms(numbers=symspg['std_types'],cell=symspg['std_lattice'],scaled_positions=symspg['std_positions'],pbc=True)
                if silent==False:
                    print 'std_lattice',ao.get_cell(),symspg['number'],symspg['international']
                    for el,pos in zip(ao.get_chemical_symbols(),ao.get_scaled_positions()):
                        print el,pos
            symspg=spglib.get_symmetry_dataset(ao)
            print "check_point25, symspg is: ",symspg
            lattice, scaled_positions, numbers=spglib.find_primitive(ao, symprec=symprec)
            if (len(scaled_positions)!=len(ao)):
                print "not primitive",len(scaled_positions),len(ao)
                #bla #TODO
            aop=Atoms(numbers=numbers,cell=lattice,scaled_positions=scaled_positions,pbc=True)
            ao_std=Atoms(numbers=symspg['std_types'],cell=symspg['std_lattice'],scaled_positions=symspg['std_positions'],pbc=True)
            sym_parent=spglib.get_symmetry_dataset(ao,symprec=1e-2) #symprec)
            ao_ini=self.get_atoms_object(uid)
            symspg_parent_ini=spglib.get_symmetry_dataset(ao_ini)
            for qvec in settings['max_subgroups']:
                qstr="q"
                for i in range(len(qvec)):
                    if fabs(qvec[i])<0.01:
                        qstr=qstr+"_0"
                    else:
                        qstr=qstr+"_%5.3f"%qvec[i]
                self.add_logmessage("setup_magnetic_structures: uid=%s, SG=%s, q=%s"%(uid,str(sym_parent['number']),qstr))
                if qstr=="q_0_0_0":
                    aoq=setup_supercell(ao,q=qvec)
                else:
                    aoq=setup_supercell(ao_std,q=qvec)
                symspg=spglib.get_symmetry_dataset(aoq,symprec=1e-2)
                print "**********new symspg with q vec*********", symspg
                msg=MSG(symelem=zip(symspg['rotations'],symspg['translations']),eps=symprec)
                #msg.get_maximal_subgroups(symdir='msg-par',sym_parent=sym_parent,q=qvec,symdirhte="XXX",silent=False)
                #bla
                print "msg is _________" , msg
                greymsg=msg.grey_msg()
                #print "XXX",len(msg.elements),len(greymsg.elements)
                #print symspg
                maxsubs=greymsg.get_maximal_subgroups(symdir='msg',sym_parent=sym_parent,q=qvec,silent=silent,symdirhte=symdirhte)
                print "*******maxsubs is *********", maxsubs
                if maxsubs==[]:
                    self.add_logmessage("WARNING(setup_magnetic_structures): failed to determine maximal subgroups for uid=%s, SG=%s, q=%s"%(uid,str(sym_parent['number']),qstr))
                exclude_ferromagnetic=True
                if 'exclude_ferromagnetic' in settings:
                    exclude_ferromagnetic=settings['exclude_ferromagnetic']
                reduce_collinear=True
                if 'reduce_collinear' in settings:
                    reduce_collinear=settings['reduce_collinear']
                for maxsub in maxsubs:
                    if (debug==True):
                        self.add_logmessage("-- %s,%s -- "%(maxsub,str(maxsubs[maxsub].name)))
                    try:
                        afmconfigs=maxsubs[maxsub].get_magnetic_configurations(aoq,name=qstr,exclude_ferromagnetic=exclude_ferromagnetic,reduce_collinear=reduce_collinear, AFatoms=AF_atoms)
                    except:
                        self.add_logmessage("WARNING(setup_magnetic_structures): failed to set up spin configurations for %s(%s)"%(uid,maxsub))
                        continue
                    if len(afmconfigs)>max_configs:
                        self.add_logmessage("WARNING(setup_magnetic_structures): too many spin configurations (%d) for %s(%s), increase max_configs=%d"%(len(afmconfigs),uid,maxsub,max_configs))
                        continue
                    for name in afmconfigs:
                        if not name in magconfigs:
                            magconfigs[name]=afmconfigs[name]
                            if (debug==True):
                                self.add_logmessage("* %s (%s,%s) added for %s."%(name,maxsub,str(maxsubs[maxsub].name),uid))
                        elif (debug==True):
                            self.add_logmessage("* %s (%s,%s) not added for %s (already in, E=%s)."%(name,maxsub,str(maxsubs[maxsub].name),uid,str(self.get_energy_per_atom(uid,calc_scheme, sub_directories={os.path.join(str(calc_scheme),name):afmconfigs[name]}))))
 
        
        
        
        
        
        
        if ('get_atoms' in settings) and (settings['get_atoms']==False):
            print "check_point147, 'get_atoms' in settings",settings
            print "check_point149, magconfigs", magconfigs
            return magconfigs
        
        ### 
        magconfigs_submitted={}
        if (debug==True) or (report_magnetic_structures==True):
            self.add_logmessage("----- setup_magnetic_structures(): %s (%s/%s/%s) -----"%(comp_r,uid,calc_scheme,str(settings)))
            # determine submitted magnetic structures:
            Nsubmitted=0
            Nconverged=0
            #E_GS=self.get_energy_per_atom(uid,calc_scheme)
            E_GS=10000
            magconfigs_submitted={}
            for jobdir in self.structureDB[uid].submitted_jobs:
                if ((jobdir.startswith(os.path.join(uid,calc_scheme,'afm_q'))) or (jobdir.startswith(os.path.join(uid,calc_scheme,'q_'))) or (jobdir.startswith(os.path.join(uid,calc_scheme,'ispin')))) and (len(jobdir.split('/'))==3):                
                    name=jobdir.split('/')[2]
                    if 'settings' in self.structureDB[uid].submitted_jobs[jobdir]:
                        magconfigs_submitted[jobdir.split('/')[2]]=self.structureDB[uid].submitted_jobs[jobdir]['settings']
                    else:
                        magconfigs_submitted[jobdir.split('/')[2]]={}
                    Nsubmitted=Nsubmitted+1
                    Epat=self.get_energy_per_atom(uid,calc_scheme, sub_directories={os.path.join(calc_scheme,name):magconfigs_submitted[name]})
                    if (Epat)!=None:
                        Nconverged=Nconverged+1
                        if (Epat<E_GS):
                            E_GS=Epat
            self.add_logmessage(" %d submitted structures, %d converged, E_min=%.3f"%(Nsubmitted,Nconverged,E_GS))
        if ('user' in settings):
            for name in settings['user']:
                use_parent_symmetry=False
                symprec=1e-3
                matches='exact'
                if 'matches' in settings['user'][name]:
                    if settings['user'][name]['matches'] in ['any_AF','exact','prototype','compound','any']:
                        matches=settings['user'][name]['matches']
                if 'symprec' in settings['user'][name]:
                    symprec=settings['user'][name]['symprec']
                if 'use_parent_symmetry' in settings['user'][name]:
                    use_parent_symmetry=settings['user'][name]['use_parent_symmetry']
                ao_ini=self.get_atoms_object(uid)
                if use_parent_symmetry==True:
                    source_info=self.structureDB[uid].get_source_info()
                    if source_info[-1][0]=='substitute_atoms':
                        ao_ini=self.get_atoms_object(source_info[-1][1]['initial_uid'])
                        atom_list_ini=source_info[-1][1]['atom_list'][0]
                    else:
                        use_parent_symmetry=False
                spglib_info_ini=spglib.get_symmetry_dataset(ao_ini,symprec=symprec)
                ao_mag=None
                if ('ao' in settings['user'][name]):
                    ao_mag=settings['user'][name]['ao']
                    spglib_info=spglib.get_symmetry_dataset(ao_mag,symprec=symprec)
                elif ('prototype' in settings['user'][name]):
                    ao_org=settings['user'][name]['prototype'].copy()
                    spglib_info=spglib.get_symmetry_dataset(ao_org,symprec=symprec)
                    if (spglib_info['number']==spglib_info_ini['number']) and (len(spglib_info['std_positions'])==len(spglib_info_ini['std_positions'])):
                        atom_list={}
                        for i in range(len(spglib_info_ini['std_positions'])):
                            for j in range(len(spglib_info['std_positions'])):
                                found=True
                                for l in range(3):
                                    if (abs(spglib_info_ini['std_positions'][i][l]-spglib_info['std_positions'][j][l])>0.001):
                                        found=False
                                if (found==True):
                                    atom_list[chemical_symbols[spglib_info['std_types'][j]]]=chemical_symbols[spglib_info_ini['std_types'][i]]
                                    break
                            if (found==False):
                                break
                        if (found==True):
                            ao_mag,uidxts=substitute_atoms(ao_org,atom_list)
                            if use_parent_symmetry==True:
                                ao_mag,uidxts=substitute_atoms(ao_mag,atom_list_ini)
                        else:
                            ao_mag=None
                    else:
                        ao_mag=None
                if (ao_mag==None):
                    self.add_logmessage("WARNING(setup_magnetic_structures): Failed to set up prototype spin structure %s for uid %s"%(name,uid))
                    continue
                ao_ini=self.get_atoms_object(uid)
                spglib_info=spglib.get_symmetry_dataset(ao_mag,symprec=symprec)
                origin_shift=spglib_info['origin_shift']
                spglib_info_ini=spglib.get_symmetry_dataset(ao_ini,symprec=symprec)
                if (spglib_info['number']==spglib_info_ini['number']) and (sorted(spglib_info['std_types'])==sorted(spglib_info_ini['std_types'])):
                    for i in range(len(spglib_info['std_types'])):
                        found=False
                        for j in range(len(spglib_info['std_types'])):
                            if (spglib_info_ini['std_types'][i]==spglib_info['std_types'][j]):
                                found=True
                                for l in range(3):
                                    if (abs(spglib_info_ini['std_positions'][i][l]-spglib_info['std_positions'][j][l]+origin_shift[l]+0.0005)%1>0.001):
                                        found=False
                            if (found==True):
                                break
                        if (found==False):
                            break
                else:
                    found=False
                if (found==True) and ('initial_spin_structure' in settings['user'][name]):
                        chem_symb=ao_mag.get_chemical_symbols() #use prop_dict instead?
                        spinstruct=settings['user'][name]['initial_spin_structure']
                        magmoms=[]
                        Nup=0
                        Ndn=0
                        namex=name
                        match_uid=True
                        match_proto=True
                        set_ispin=True
                        for iat,momf,el in zip(range(len(chem_symb)),spinstruct,chem_symb):
                            mom_el=self.get_initial_magnetic_moment(el)
                            if abs(mom_el)>0.1:
                                set_ispin=False
                            if ('non_default_magnetic_atoms' in settings['user'][name]) and (el in settings['user'][name]['non_default_magnetic_atoms']):
                                mom_el=settings['user'][name]['non_default_magnetic_atoms'][el]
                                if (mom_el!=self.get_initial_magnetic_moment(el)) and (not (el in namex.split('_'))):
                                    namex=namex+'_%s_%.1f'%(el,mom_el)
                            mom=momf*mom_el
                            magmoms.append(mom)
                            if (mom>0.1):
                                Nup=Nup+1
                            elif (mom<-0.1):
                                Ndn=Ndn+1
                            else:
                                if (momf!=0):
                                    match_proto=False
                                elif (abs(mom_el)>0.1):
                                    match_uid=False
                        if (matches=='any_AF'):
                            doit=True
                        elif (matches=='exact') and (match_proto==True) and (match_uid==True):
                            doit=True
                        elif  (matches=='prototype') and (match_proto==True):
                            doit=True
                        elif  (matches=='compound') and (match_uid==True):
                            doit=True
                        else:
                            doit=False
                        if ((Nup>0) and (Ndn>0)) and (doit==True):
                            magconfigs['q_user_'+namex]={'magmom':magmoms,'atoms_object':ao_mag}
                            if set_ispin==True:
                                magconfigs['q_user_'+namex]['ispin']=2
                        else:
                            self.add_logmessage("WARNING(setup_magnetic_structures): spin structure %s for uid %s excluded because it does not match condition <%s>."%(name,uid,matches))
                else:
                        self.add_logmessage("WARNING(setup_magnetic_structures): Failed to set up spin structure %s for uid %s, ao not o.k."%(name,uid))
        if ('parent_cell' in settings):
            pd_sc={}
            pd={}
            try:
                method,si=self.structureDB[uid].source_info[-1]
                if (method=='create_interstitial_defect') and ('initial_uid' in si):
                    pd=self.get_properties(si['initial_uid'],calc_scheme,magsettings=settings['parent_cell'])
                    #print "***",si
                    #print pd
                    interstial_atoms=[]
                    for ib in si['interstitial_blocks']:
                        interstial_atoms.append((si['interstitial'], si['scaled position'], ib, None))
                    N=[1,1,1]
                    if si['supercell']!=[]:
                        N=si['supercell']
                    distortion=[]
                    if 'distortion' in si:
                        distortion=si['distortion']
                    pd_sc=Create_Magnetic_Supercell(pd, interstial_atoms=interstial_atoms,distortion=distortion,N=N)
            except:
                pd_sc={}
            if pd_sc=={}:
                self.add_logmessage("WARNING(setup_magnetic_structures): Failed to set up magnetic structure from parent cell for %s(%s,settings=%s)"%(uid,calc_scheme, settings['parent_cell']))
            else:
                ao=Atoms(pd_sc['chemical_symbols'],cell=pd_sc['cell'],scaled_positions=pd_sc['scaled_positions'],pbc=True)
            if ('subdir' in pd) and (pd['subdir']!=''):
                name='%s_pc'%pd['subdir']
                if not name in magconfigs:
                    magconfigs[name]={'magmom':pd_sc['magnetic_moments'],'atoms_object':ao,'subdir':name}
                    if isinstance(pd_sc['magnetic_moments'][0],list):
                        magconfigs[name]['lnoncollinear']=True
        if ('init_structure' in settings):
            print "check_point42, 'init_structure' in settings."
            calc_name, calc_settings=self.calc_schemes[calc_scheme]
            prop_dicts=[]
            print "check_point43, calc_name is :",calc_name,"calc_settings is:", calc_settings
            if 'init_structure' in calc_settings:
                print "check_point44, 'init_structure' in calc_settings."
                cs_ini=calc_settings['init_structure']
                print "check_point45, cs_ini shold be like 'vaspopt':",cs_ini
                if ('all' in settings['init_structure']) and (settings['init_structure']['all']==True):
                    configs=self.setup_magnetic_structures(uid,cs_ini, magsettings=settings['init_structure'])
                    for config in configs:
                        prop_dict=deepcopy(self.get_properties(uid, cs_ini, magsettings={},sub_directories={os.path.join(cs_ini,config):configs[config]}))
                        if prop_dict!={}:
                            prop_dicts.append(prop_dict)
                else:
                    prop_dicts=[deepcopy(self.get_properties(uid,cs_ini, magsettings=settings['init_structure']))]
            print "init setting is here"
            for prop_dict in prop_dicts:
                add_conf=True
                for x in ['chemical_symbols','cell','scaled_positions']:
                    if not (x in prop_dict):
                        add_conf=False
                if not ('magnetic_moments' in prop_dict) and ('initial_magnetic_moments' in prop_dict):
                    prop_dict['magnetic_moments']=prop_dict['initial_magnetic_moments']
                if not ('magnetic_moments' in prop_dict):
                    add_conf=False
                if add_conf==True:
                    if 'symprec' in calc_settings:
                        #create a dummy atoms object with symmetry of magnetic configuration
                        dummy_symb_pool=[]
                        for el in ['Sc', 'Ti', 'V', 'Cr','Mn','Fe','Co','Ni','Cu','Zn','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','La','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg']:
                            if not el in prop_dict['chemical_symbols']:
                                dummy_symb_pool.append(el)
                        dummy_symb=[]
                        dumno=0
                        dummy_tab={}
                        for el,mom in zip(prop_dict['chemical_symbols'],prop_dict['magnetic_moments']):
                            newdum=True
                            for dum_el in dummy_tab:
                                elx,momx=dummy_tab[dum_el]
                                if (el==elx) and (abs(mom-momx)<eps_mom):
                                    dummy_symb.append(dum_el)
                                    newdum=False
                                    break
                            if newdum==True:
                                dum_el=dummy_symb_pool[dumno]
                                dummy_tab[dum_el]=(el,mom)
                                dummy_symb.append(dum_el)
                                dumno=dumno+1
                        #print "XXX",prop_dict['chemical_symbols'],dummy_symb,dummy_tab
                        ao=Atoms(dummy_symb,cell=prop_dict['cell'],scaled_positions=prop_dict['scaled_positions'],pbc=True)
                        #symspg=spglib.get_symmetry_dataset(ao, symprec=calc_settings['symprec']) #returns sometimes a larger cell
                        prop_dict['cell'], prop_dict['scaled_positions'], numbers = spglib.standardize_cell(ao, symprec=calc_settings['symprec'], to_primitive=True)
                        chem_symb_new=[]
                        magmom_new=[]
                        for atnum in numbers: #symspg['std_types']:
                            for dum_el in dummy_tab:
                                if atomic_numbers[dum_el]==atnum:
                                    el,mom=dummy_tab[dum_el]
                                    chem_symb_new.append(el)
                                    magmom_new.append(mom)
                                    break
                        prop_dict['chemical_symbols']=chem_symb_new
                        #prop_dict['cell']=symspg['std_lattice']
                        #prop_dict['scaled_positions']=symspg['std_positions']
                        prop_dict['magnetic_moments']=magmom_new
                    ao=Atoms(prop_dict['chemical_symbols'],cell=prop_dict['cell'],scaled_positions=prop_dict['scaled_positions'],pbc=True)
                    if ('subdir' in prop_dict) and (prop_dict['subdir']!=''):
                        name=prop_dict['subdir']+'cs_ini'
                        if not name in magconfigs:
                            magconfigs[name]={'magmom':prop_dict['magnetic_moments'],'atoms_object':ao,'subdir':name}
                            if isinstance(prop_dict['magnetic_moments'][0],list):
                                magconfigs[name]['lnoncollinear']=True
        # print "magconfigs is:",magconfigs
        # print "before ao return"
        # print "calc_sheme is", calc_scheme
        # print "uid is", uid
        # print "only uid ao is:", self.get_atoms_object(uid)
        print "check_point26, before get_atom_object"
        ao=self.get_atoms_object(uid, calc_scheme=calc_scheme, magsettings={})
        print "check_point27,after ao return"
        if (ao==None):
            print "check_point28,ao is None, will return empty dictionary"
            return magconfigs
        # determine site magnetic moments
        if isinstance(AF_atoms,dict):
            # determine AFatoms by site magnetic moments
            AF_atoms_new=[]
            magmomsFM=self.get_magnetic_moments(uid, calc_scheme, magsettings={})
            for (el,mom) in magmomsFM:
                if (el in AF_atoms) and (abs(mom)>AF_atoms[el]):
                    if not (el in AF_atoms_new):
                        AF_atoms_new.append(el)
            AF_atoms=AF_atoms_new
            print "AF atoms:",AF_atoms
        if forceMAG!=[]:
            print "XXX",forceMAG
            magmomsFM=self.get_magnetic_moments(uid, calc_scheme, magsettings={})
            print "forceMAG",magmomsFM
            for (el,mom) in magmomsFM:
                print el,mom,norm(mom)>0.5
        # FM coupling of non-default magnetic atoms
        if ('initial_magnetic_moments' in settings):
            magmoms=[]
            els_nondef=[]
            for el in ao.get_chemical_symbols():
                M_el=self.get_initial_magnetic_moment(el)
                if (el in settings['initial_magnetic_moments']) and (abs(settings['initial_magnetic_moments'][el]-M_el)>0.1):
                    M_el=settings['initial_magnetic_moments'][el]
                    if not (el in els_nondef):
                        els_nondef.append(el)
                magmoms.append(M_el)
            if len(els_nondef)!=0:
                name="q_0_0_0"
                for el in sorted(els_nondef):
                    name=name+"_%s%.1f"%(el,settings['initial_magnetic_moments'][el])
                if not name in magconfigs:
                    magconfigs[name]={'magmom':magmoms,'atoms_object':ao}
        # AF coupling of non-equivalent atoms 
        if ('sublattices' in settings):
            maxconfigs=max_configs
            print "check_point20, in sublattice, maxconfigs is :",maxconfigs
            if 'max_configs' in settings:
                maxconfigs=settings['max_configs']
            afmconfigs=get_magnetic_sublattices(ao, symprec=symprec)
            if afmconfigs==False:
                print "check_point21, in sublattice, afmconfigs is False:", afmconfigs
                self.add_logmessage("WARNING(setup_magnetic_structures): Failed to get sublattices for %s (%s/%s), check spglib!"%(uid,comp_r,calc_scheme))
            elif (pow(2,len(afmconfigs)-1)-1)>maxconfigs:
                print "check_point22, pow(2,len(afmconfigs)-1)-1 is larger than max config", 
                self.add_logmessage("WARNING(setup_magnetic_structures): %d sublattice configurations for %s (%s), increase max_configs!"%(pow(2,len(afmconfigs)-1)-1,uid,calc_scheme))
            else:
                afmconfigs=get_magnetic_sublattices(ao,return_afm=True, symprec=symprec)
                print "check_point23, print_afmconfigs:", afmconfigs
                if len(afmconfigs)>maxconfigs:
                    self.add_logmessage("WARNING(setup_magnetic_structures): %d sublattice configurations for %s, increase max_configs!"%(len(afmconfigs),uid))
                    afmconfigs={}
                for name in afmconfigs:
                    if not name in magconfigs:
                        magconfigs[name]=afmconfigs[name]
                        if (debug==True):
                            self.add_logmessage("* %s added for %s."%(name,uid))
                    elif (debug==True):
                        self.add_logmessage("* %s not added for %s (already in, E=%s)."%(name,uid,str(self.get_energy_per_atom(uid,calc_scheme, sub_directories={os.path.join(calc_scheme,name):afmconfigs[name]}))))
        if ('initspin' in settings):
            #user defined initial spin split
            if len(ao)==len(settings['initspin']):
                name="afm_q0_isp_"
                for mom in settings['initspin']:
                    name=name+"_%.1f"%mom
                if not name in magconfigs:
                    magconfigs[name]={'magmom':settings['initspin'],'atoms_object':ao}
        if ('initial_spin_splits' in settings):
            #user defined initial spin split
            for spinstr in settings['initial_spin_splits']:
                try:
                    if (isinstance(spinstr,list)) and (len(ao)==len(spinstr)):
                        name="q_0_0_0_isp_"
                        aox=ao
                        magmoms=spinstr
                    elif (isinstance(spinstr,dict)) and ('q' in spinstr):
                        name="q_"
                        for i in range(3):
                            if isinstance(spinstr['q'][i],int):
                                name=name+"%d_"%spinstr['q'][i]
                            else:
                                name=name+"%.2f_"%spinstr['q'][i]
                        name=name+"isp_"
                        aox=setup_supercell(ao, q=spinstr['q'])
                        magmoms=spinstr['spins']
                    chem_symb=aox.get_chemical_symbols() #use prop_dict instead?
                    for iat,mom,el in zip(range(len(chem_symb)),magmoms,chem_symb):
                        if norm(mom)>0.05:
                            name=name+"%s%d"%(el,iat)
                            if isinstance(mom,np.ndarray):
                                for i in range(len(mom)):
                                    name=name+"_%.1f"%mom[i]
                            else:
                                name=name+"_%.1f"%mom
                    if not name in  magconfigs:
                        magconfigs[name]={'magmom':magmoms,'atoms_object':aox}
                        if isinstance(magmoms[0],np.ndarray):
                            magconfigs[name]['lnoncollinear']=True
                except:
                    self.add_logmessage("WARNING(setup_magnetic_structures): Failed to set up spin structure %s for uid=%s and calc scheme %s (%s)"%(str(spinstr),uid,calc_scheme, str(sys.exc_info())))
        if ('initial_spin_structures' in settings):
            #user defined initial spin structure
            for spinstr in settings['initial_spin_structures']:
                try:
                    if (isinstance(spinstr,list)) and (len(ao)==len(spinstr)):
                        name="q_0_0_0_isp_"
                        aox=ao
                        spinstruct=spinstr
                    elif (isinstance(spinstr,dict)) and ('q' in spinstr):
                        name="q_"
                        for i in range(3):
                            if isinstance(spinstr['q'][i],int):
                                name=name+"%d_"%spinstr['q'][i]
                            else:
                                name=name+"%.2f_"%spinstr['q'][i]
                        name=name+"isp_"
                        aox=setup_supercell(ao, q=spinstr['q'])
                        spinstruct=spinstr['spins']
                    chem_symb=aox.get_chemical_symbols() #use prop_dict instead?
                    magmoms=[]
                    for iat,momf,el in zip(range(len(chem_symb)),spinstruct,chem_symb):
                        mom=momf*self.get_initial_magnetic_moment(el)
                        magmoms.append(mom)
                        if norm(mom)>0.05:
                            name=name+"%s%d"%(el,iat)
                            if isinstance(mom,np.ndarray):
                                for i in range(len(mom)):
                                    name=name+"_%.1f"%mom[i]
                            else:
                                name=name+"_%.1f"%mom
                    if not name in  magconfigs:
                        magconfigs[name]={'magmom':magmoms,'atoms_object':aox}
                        if isinstance(magmoms[0],np.ndarray):
                            magconfigs[name]['lnoncollinear']=True
                except:
                    self.add_logmessage("WARNING(setup_magnetic_structures): Failed to set up spin structure %s (%s)"%(str(spinstr),str(sys.exc_info())))
        if ('mcifs' in settings):
            #to be done: full mcif reader, this works only for q=[0 0 0]
            mcifs=glob.glob(settings['mcifs'])
            for mcif in mcifs:
                msg=MSG(mcif)
                afmconfigs=msg.get_magnetic_configurations(ao,name="q_0_0_0") #mcif)
                for name in afmconfigs:
                    if not name in magconfigs:
                        magconfigs[name]=afmconfigs[name]
        if ('user_magconf' in settings):
            #to be done: make sure this is consistent
            for name in settings['user_magconf']:
                namex=name+'_user_magconf'
                if not namex in  magconfigs:
                    magconfigs[namex]={'magmom':settings['user_magconf'][name]['magmom'],'atoms_object':settings['user_magconf'][name]['atoms_object']}
                    if isinstance(settings['user_magconf'][name]['magmom'][0],np.ndarray):
                        magconfigs[namex]['lnoncollinear']=True
        if ('autoprio' in settings):
            status_Ntot=0
            status_Nnew=0
            status_Nuncoverged=0
            status_failed=False
            Nmin=10
            Nmax=30
            dmax=3.5
            use_sublattices=True
            use_primitive=False
            engine_new=False
            ordmaxmsg=-1
            if 'use_sublattices' in settings['autoprio']:
                use_sublattices=settings['autoprio']['use_sublattices']
            if 'Nmin' in settings['autoprio']:
                Nmin=settings['autoprio']['Nmin']
            if 'Nmax' in settings['autoprio']:
                Nmax=settings['autoprio']['Nmax']
            if 'dmax' in settings['autoprio']:
                dmax=settings['autoprio']['dmax']
            if 'use_primitive' in settings['autoprio']:
                use_primitive=settings['autoprio']['use_primitive']
            if 'engine_new' in settings['autoprio']:
                engine_new=settings['autoprio']['engine_new']
            if 'ordmaxmsg' in settings['autoprio']:
                ordmaxmsg=settings['autoprio']['ordmaxmsg'] # debug option for msg tests to avoid msgs which require large CPU time
            if 'symdirhte' in settings['autoprio']:
                symdirhte=settings['autoprio']['symdirhte']
            sym_parent=spglib.get_symmetry_dataset(ao,symprec=symprec)
            print "sym_parent is ", sym_parent
            symspg=spglib.get_symmetry_dataset(ao,symprec=symprec)
            ao_std=Atoms(numbers=symspg['std_types'],cell=symspg['std_lattice'],scaled_positions=symspg['std_positions'],pbc=True)
	    #tmp-CHECK
            if use_primitive==True:
                ao_std=ao
                symdir='msgp'
            else:
                symdir='msg'
	    #ao_std=ao
            # check sublattices
            if use_sublattices==True:
                afmconfigs=get_magnetic_sublattices(ao_std, symprec=symprec,return_afm=True, Nmax=Nmax)
                if afmconfigs==False:
                    status_failed=True
                    self.add_logmessage("WARNING(setup_magnetic_structures):autoprio, failed to get sublattices for %s (%s/%s), check spglib!"%(uid,comp_r,calc_scheme))
                else:
                    self.add_logmessage(" -- probing AF coupling of sublattices")
                    for name in afmconfigs:
                        if len (magconfigs)>=Nmax:
                            self.add_logmessage(" Nmax reached, stopping here (more configs possible")
                        if (name in magconfigs_submitted) or (name in magconfigs):
                            if (debug==True) or (report_magnetic_structures==True):
                                self.add_logmessage(" * %s already in, E=%s)."%(name,str(self.get_energy_per_atom(uid,calc_scheme, sub_directories={os.path.join(calc_scheme,name):afmconfigs[name]}))))
                            if not name in magconfigs:
                                magconfigs[name]=afmconfigs[name]
                                status_Ntot=status_Ntot+1
                        else:
                            magconfigs[name]=afmconfigs[name]
                            status_Ntot=status_Ntot+1
                            status_Nnew=status_Nnew+1
                            status_Nuncoverged=status_Nuncoverged+1
                            if (debug==True) or (report_magnetic_structures==True):
                                self.add_logmessage(" * %s added for %s."%(name,uid))
                        if len (magconfigs)>=Nmax:
                            break
                self.add_logmessage(" --Ntot %d for %s (after get_magnetic_sublattices())"%(len(magconfigs),uid))
            #determine q-vectors from lattice parameters
            self.add_logmessage(" -- Setting up magnetic structures with maximal magnetic subgroup symmetry")
            qvecs=[[0,0,0]]
            if (len(ao_std)!=len(ao)):
                qvecs=[[0,0,0], [0,0,1]]
            logmes=" - autoq: lattice parameters: "
            for i in range(3):
                ai=norm(ao_std.get_cell()[i])
                logmes=logmes+"%5.2f "%ai
            spgno=symspg['number']
            if (spgno>195) and (norm(ao_std.get_cell()[0])<dmax):
                if (use_primitive==True) or (spgno<221):
                    qvecs=qvecs+[[0, 0, 0.5],[0.5, 0.5, 0],[0.5, 0.5, 0.5]]
                else:
                    qvecs=qvecs+[[0, 0, 0.5],[0.5, 0.5, 0]]
            else:
                qix=[[0],[0],[0]]
                for i in range(3):
                    ai=norm(ao_std.get_cell()[i])
                    for j in range(2,5):
                        if (ai*(j-1)<dmax):
                            qix[i].append(float(1./j))
                print qix
                for qx in qix[0]:
                    for qy in qix[1]:
                        for qz in qix[2]:
                            qv=[qx,qy,qz]
                            if qv!=[0,0,0]:
                                qvecs.append(qv)
                            print qv,qvecs
            logmes=logmes+" --> qvecs=%s"%str(qvecs)
            self.add_logmessage(logmes)
            for qvec in qvecs: #[[0, 0, 0], [0,0,1],[0, 0, 0.5],[0.5, 0, 0],[0.5, 0.5, 0],[0.5, 0, 0.5],[0.5, 0.5, 0.5]]:
                # qvec for cubic, tet, etc.
                if len(magconfigs)>=Nmax:
                    break
                # ignore centerings for the moment
                qstr="q"
                for i in range(len(qvec)):
                    if fabs(qvec[i])<0.01:
                        qstr=qstr+"_0"
                    else:
                        qstr=qstr+"_%5.3f"%qvec[i]
                if qstr=="q_0_0_0":
                    aoq=setup_supercell(ao,q=qvec)
                else:
                    aoq=setup_supercell(ao_std,q=qvec)
                if (qvec==[0,0,1]) and (len(aoq)==len(ao)):
                    continue
                if check_ao_equivalence(ao,aoq)==False:
                    print "SCSC: supercell inequivalent to original ao",ao,aoq
                    bla
                    
                symspg=spglib.get_symmetry_dataset(aoq, symprec=symprec)
                if (debug==True):
                    print "XXX",symspg['number']
                    for i in zip(symspg['rotations'],symspg['translations']):
                        print i
                msg=MSG(symelem=zip(symspg['rotations'],symspg['translations']),eps=symprec)
                print "Parent SG done"
                greymsg=msg.grey_msg()
                print "Grey MSG done"
                if (ordmaxmsg<0) or (greymsg.get_order()<=ordmaxmsg):
                    maxsubs=greymsg.get_maximal_subgroups(symdir=symdir,sym_parent=sym_parent,q=qvec,silent=silent,symdirhte=symdirhte,engine_new=engine_new)
                print "Maximal subgroups done"
                if maxsubs==[]:
                    self.add_logmessage("WARNING(setup_magnetic_structures): failed to determine maximal subgroups for uid=%s, SG=%s, q=%s"%(uid,str(sym_parent['number']),qstr))
                    status_failed=True                    
                Ndim=100000
                for maxsub in sorted(maxsubs, key=lambda(k):len(maxsubs[k].get_elements()),reverse=True):
                    if len (magconfigs)>=Nmax:
                        self.add_logmessage("  - setup_magnetic_structures(): Nmax reached, more configurations possible")
                        break
                    if (len(maxsubs[maxsub].get_elements())<Ndim) and (len(magconfigs)>Nmax):
                        self.add_logmessage("--Ntot skipping %s (%d elements) for %s)"%(maxsub,len(maxsubs[maxsub].get_elements()),uid))
                        continue
                    Ndim=len(maxsubs[maxsub].get_elements())
                    try:
                        afmconfigs,status=maxsubs[maxsub].get_magnetic_configurations(aoq,name=qstr, return_status=True, AFatoms=AF_atoms)
                    except:
                        self.add_logmessage("WARNING(setup_magnetic_structures): failed to set up spin configurations for %s(%s)"%(uid,maxsub))
                        status_failed=True
                        continue
                    self.add_logmessage("-- %s: # %d (%s)"%(maxsub,len(afmconfigs),str(status)))
                    for name in afmconfigs:
                        if len (magconfigs)>=Nmax:
                            self.add_logmessage("  - setup_magnetic_structures(): Nmax reached, more configurations possible")
                            break
                        if (name in magconfigs_submitted) or (name in magconfigs):
                            if (debug==True) or (report_magnetic_structures==True):
                                self.add_logmessage(" * %s already in, E=%s)."%(name,str(self.get_energy_per_atom(uid,calc_scheme, sub_directories={os.path.join(calc_scheme,name):afmconfigs[name]}))))
                            if not name in magconfigs:
                                magconfigs[name]=afmconfigs[name]
                                status_Ntot=status_Ntot+1
                        else:
                            magconfigs[name]=afmconfigs[name]
                            status_Ntot=status_Ntot+1
                            status_Nnew=status_Nnew+1
                            status_Nuncoverged=status_Nuncoverged+1
                            if (debug==True) or (report_magnetic_structures==True):
                                self.add_logmessage("* %s added for %s."%(name,uid))
                    self.add_logmessage("--Ntot %d for %s (after %s (%s), %d elements, %s)"%(len(magconfigs),uid,maxsub,str(maxsubs[maxsub].name),len(maxsubs[maxsub].get_elements()),str(status)))
                logmes="   status_autoprio(%s (%s/%s): "%(comp_r,uid,calc_scheme)
                if (status_Nuncoverged==0)  and (status_failed==False):
                    logmes=logmes+" complete (%d magnetic configurations converged)"%status_Ntot
                else:
                    logmes=logmes+" incomplete (%d magnetic configurations"%status_Ntot
                    if (status_Nuncoverged>0):
                        logmes=logmes+", %d not converged"%status_Nuncoverged
                    if (status_Nnew>0):
                        logmes=logmes+", %d new configurations"%status_Nuncoverged
                    if (status_failed==True):
                        logmes=logmes+", setup failed"
                logmes=logmes+")"
                self.add_logmessage(logmes)
                #self.add_logmessage("--NTOT %d for %s (after %s)"%(len(magconfigs),uid,qstr))
        # if ('max_subgroups' in settings):
        #     #symprec=1e-2
        #     print "check_point 24, begin_max_groups."
        #     if 'symprec' in settings:
        #         symprec=settings['symprec']
        #         symspg=spglib.get_symmetry_dataset(ao, symprec=settings['symprec'])
        #         ao=Atoms(numbers=symspg['std_types'],cell=symspg['std_lattice'],scaled_positions=symspg['std_positions'],pbc=True)
        #         if silent==False:
        #             print 'std_lattice',ao.get_cell(),symspg['number'],symspg['international']
        #             for el,pos in zip(ao.get_chemical_symbols(),ao.get_scaled_positions()):
        #                 print el,pos
        #     symspg=spglib.get_symmetry_dataset(ao)
        #     print "check_point25, symspg is: ",symspg
        #     lattice, scaled_positions, numbers=spglib.find_primitive(ao, symprec=symprec)
        #     if (len(scaled_positions)!=len(ao)):
        #         print "not primitive",len(scaled_positions),len(ao)
        #         #bla #TODO
        #     aop=Atoms(numbers=numbers,cell=lattice,scaled_positions=scaled_positions,pbc=True)
        #     ao_std=Atoms(numbers=symspg['std_types'],cell=symspg['std_lattice'],scaled_positions=symspg['std_positions'],pbc=True)
        #     sym_parent=spglib.get_symmetry_dataset(ao,symprec=1e-2) #symprec)
        #     ao_ini=self.get_atoms_object(uid)
        #     symspg_parent_ini=spglib.get_symmetry_dataset(ao_ini)
        #     for qvec in settings['max_subgroups']:
        #         qstr="q"
        #         for i in range(len(qvec)):
        #             if fabs(qvec[i])<0.01:
        #                 qstr=qstr+"_0"
        #             else:
        #                 qstr=qstr+"_%5.3f"%qvec[i]
        #         self.add_logmessage("setup_magnetic_structures: uid=%s, SG=%s, q=%s"%(uid,str(sym_parent['number']),qstr))
        #         if qstr=="q_0_0_0":
        #             aoq=setup_supercell(ao,q=qvec)
        #         else:
        #             aoq=setup_supercell(ao_std,q=qvec)
        #         symspg=spglib.get_symmetry_dataset(aoq,symprec=1e-2)
        #         print "**********new symspg with q vec*********", symspg
        #         msg=MSG(symelem=zip(symspg['rotations'],symspg['translations']),eps=symprec)
        #         #msg.get_maximal_subgroups(symdir='msg-par',sym_parent=sym_parent,q=qvec,symdirhte="XXX",silent=False)
        #         #bla
        #         print "msg is _________" , msg
        #         greymsg=msg.grey_msg()
        #         #print "XXX",len(msg.elements),len(greymsg.elements)
        #         #print symspg
        #         maxsubs=greymsg.get_maximal_subgroups(symdir='msg',sym_parent=sym_parent,q=qvec,silent=silent,symdirhte=symdirhte)
        #         print "*******maxsubs is *********", maxsubs
        #         if maxsubs==[]:
        #             self.add_logmessage("WARNING(setup_magnetic_structures): failed to determine maximal subgroups for uid=%s, SG=%s, q=%s"%(uid,str(sym_parent['number']),qstr))
        #         exclude_ferromagnetic=True
        #         if 'exclude_ferromagnetic' in settings:
        #             exclude_ferromagnetic=settings['exclude_ferromagnetic']
        #         reduce_collinear=True
        #         if 'reduce_collinear' in settings:
        #             reduce_collinear=settings['reduce_collinear']
        #         for maxsub in maxsubs:
        #             if (debug==True):
        #                 self.add_logmessage("-- %s,%s -- "%(maxsub,str(maxsubs[maxsub].name)))
        #             try:
        #                 afmconfigs=maxsubs[maxsub].get_magnetic_configurations(aoq,name=qstr,exclude_ferromagnetic=exclude_ferromagnetic,reduce_collinear=reduce_collinear, AFatoms=AF_atoms)
        #             except:
        #                 self.add_logmessage("WARNING(setup_magnetic_structures): failed to set up spin configurations for %s(%s)"%(uid,maxsub))
        #                 continue
        #             if len(afmconfigs)>max_configs:
        #                 self.add_logmessage("WARNING(setup_magnetic_structures): too many spin configurations (%d) for %s(%s), increase max_configs=%d"%(len(afmconfigs),uid,maxsub,max_configs))
        #                 continue
        #             for name in afmconfigs:
        #                 if not name in magconfigs:
        #                     magconfigs[name]=afmconfigs[name]
        #                     if (debug==True):
        #                         self.add_logmessage("* %s (%s,%s) added for %s."%(name,maxsub,str(maxsubs[maxsub].name),uid))
        #                 elif (debug==True):
        #                     self.add_logmessage("* %s (%s,%s) not added for %s (already in, E=%s)."%(name,maxsub,str(maxsubs[maxsub].name),uid,str(self.get_energy_per_atom(uid,calc_scheme, sub_directories={os.path.join(str(calc_scheme),name):afmconfigs[name]}))))
        # # if ('get_atoms' in settings) and (settings['get_atoms']==False):
        #     print "check_point147, 'get_atoms' in settings",settings
        #     return magconfigs
        # #先把这个挪到这试试
        # 试了下是不行的
        if ('ispin' in settings):
            if isinstance(settings['ispin'],list):
                ispins=settings['ispin']
            else:
                ispins=[settings['ispin']]
            for isp in ispins:
                name="ispin%d"%isp
                magconfigs[name]={'ispin':isp}
        delconfs=[]
        for name in magconfigs:
            if ('non_collinear' in settings):
                if (settings['non_collinear']==False) and ('lnoncollinear' in magconfigs[name]) and (magconfigs[name]['lnoncollinear']==True):
                    delconfs.append(name)
                else:
                    self.add_logmessage(" added %s"%name) #%(%s), settings=%s"%(name,maxsub,str(magconfigs[name])))
            magconfigs[name]['subdir']=name
        for name in delconfs:
            del magconfigs[name]
        if (debug==True):
            self.add_logmessage("setup_magnetic_structures(): %d magnetic structures for uid=%s (%s, magsettings=%s)"%(len(magconfigs),uid,calc_scheme, str(settings)))
        if ('auto' in settings) and (settings['auto']==False):
            for x in magconfigs.keys():
                self.add_logmessage(x)
            self.write_logmessages()
            magconfigs={}
        return magconfigs
    
    def get_energy(self, uid, calc_scheme, update=False, sloppy_mode=False, return_dict=False, magsettings={'submitted':True}, sub_directories={}, nsub_max=2, special=''):
        energy=None
        use_prop_dict=self.use_prop_dict
        if (use_prop_dict==True) and (return_dict==False):
            prop_dict=self.get_properties(uid, calc_scheme, magsettings=magsettings, sub_directories=sub_directories, update=update)
            if 'energy' in prop_dict:
                return prop_dict['energy']
            else:
                return None
        edict={}
        Ts={'C11':np.array([[1.0, 0.0, 0.0],[0.0, 0.0, 0.0],[0.0, 0.0,0.0]])
            ,'C22':np.array([[0.0, 0.0, 0.0],[0.0, 1.0, 0.0],[0.0, 0.0,0.0]])
            ,'C33':np.array([[0.0, 0.0, 0.0],[0.0, 0.0, 0.0],[0.0, 0.0,1.0]])
            ,'C44':np.array([[0.0, 0.0, 0.0],[0.0, 0.0, 0.5],[0.0, 0.5,0.0]])
            ,'C55':np.array([[0.0, 0.0, 0.5],[0.0, 0.0, 0.0],[0.5, 0.0,0.0]])
            ,'C66':np.array([[0.0, 0.5, 0.0],[0.5, 0.0, 0.0],[0.0, 0.0,0.0]])
            ,'C12':np.array([[1.0, 0.0, 0.0],[0.0, 1.0, 0.0],[0.0, 0.0,0.0]])
            ,'C13':np.array([[1.0, 0.0, 0.0],[0.0, 0.0, 0.0],[0.0, 0.0,1.0]])
            ,'C23':np.array([[0.0, 0.0, 0.0],[0.0, 1.0, 0.0],[0.0, 0.0,1.0]])
            #,'Cp':np.array([[1.,0.,0.],[0.,-0.5,0.],[0.,0.,-0.5]])
           }
        if special.startswith("V_"): #bulkmodulus
            subdir=os.path.join(calc_scheme,special)
            elastic=True
            eps=pow(float(special.split("_")[1]),1./3.)
            emat=eps*np.identity(3)
        elif (special.startswith('C')) and (special.split("_")[0] in Ts):
            subdir=os.path.join(calc_scheme,special)
            elastic=True
            index=special.split("_")[0]
            T=Ts[index]
            eps=float(special.split("_")[1])
            emat=np.identity(3)+eps*T
        else:
            subdir=calc_scheme
            elastic=False
            emat=None
        subdirs={subdir:{}}
        if isinstance(magsettings,dict):
            if ('include_default' in magsettings) and (magsettings['include_default']==False):
                subdirs={}
            magconfigs=self.setup_magnetic_structures(uid,calc_scheme, magsettings=magsettings)
            for magconf in magconfigs:
                subdir=os.path.join(calc_scheme,magconf)
                subdirs[subdir]=magconfigs[magconf]
        if (calc_scheme in self.calc_schemes) and (uid in self.structureDB):
            for subdir in subdirs:
                prop_dict=self.get_properties(uid, calc_scheme, magsettings={},sub_directories={subdir:subdirs[subdir]})
                if 'energy' in prop_dict:
                    edict[subdir]=prop_dict['energy']
                    continue
                #check if energy is already stored in database
                if (sloppy_mode==True):
                    energy=self.structureDB[uid].get_energy(subdir, None, update=False)
                    if (energy!=None):
                        edict[subdir]=energy
                        continue
                parentdir=os.getcwd()
                #check for converged calculations in searchpaths:
                calc=self.setup_calculator(uid,calc_scheme,elastic=elastic,emat=emat,afm=subdirs[subdir]) #check afm
                for sp in self.get_searchpaths():
                    try:
                        os.chdir(sp)
                    except OSError, e:
                        warnings.warn("Searchpaths: OSError {0}".format(e))
                        continue
                    energy=self.structureDB[uid].get_energy(subdir,calc, sloppy_mode=sloppy_mode, update=False, job_commands=self.get_job_commands(calc_scheme))
                    os.chdir(parentdir)
                    if energy!=None:
                        #edict[subdir]=energy
                        break
                if ((energy==None) and (update==True)):
                    calc=self.setup_calculator(uid,calc_scheme,elastic=elastic,emat=emat,afm=subdirs[subdir],update=update)
                    sp=self.get_scratch_directory()
                    if sp!=None:
                        os.chdir(sp)
                        energy=self.structureDB[uid].get_energy(subdir,calc, sloppy_mode=sloppy_mode, update=update, job_commands=self.get_job_commands(calc_scheme),nsub_max=nsub_max) 
                        os.chdir(parentdir)
                edict[subdir]=energy
            energy=None
            reference_state=''
            for subdir in edict:
                if (edict[subdir]!=None) and ((energy==None) or (edict[subdir]<energy)):
                    energy=edict[subdir]
                    reference_state=subdir
            if reference_state!='':
                edict['reference_state']=reference_state
        if return_dict==False:
            return energy
        else:
            return energy,edict
    
    def get_energy_per_atom(self, uid, calc_scheme, sloppy_mode=False, update=False, nsub_max=2, magsettings={'submitted':True}, sub_directories={}):
        E=None
        prop_dict=self.get_properties(uid, calc_scheme, magsettings=magsettings, sub_directories=sub_directories, nsub_max=nsub_max, update=update)
        if ('energy' in prop_dict) and ('chemical_symbols' in prop_dict):
            E=prop_dict['energy']/float(len(prop_dict['chemical_symbols']))
        return E
    
    def get_formation_energy(self, uid, calc_scheme, sloppy_mode=False, update=False, nsub_max=2, reference_energy={},returns=[], magsettings={'submitted':True}, sub_directories={}):
        """get the formation energy per atom for compound with respect to decomposition
        into the elements for specified calc_scheme. 
        """
        used_calc_scheme=None
        E_form=None
        if isinstance(calc_scheme,tuple):
            for cs in calc_scheme:
                Ef=self.get_energy_per_atom(uid, cs, sloppy_mode=sloppy_mode, update=update, nsub_max=nsub_max, magsettings=magsettings, sub_directories=sub_directories)
            if (E_form==None) or ((Ef!=None) and (Ef<E_form)):
                E_form=Ef
                used_calc_scheme=cs
        else:
            E_form=self.get_energy_per_atom(uid, calc_scheme, sloppy_mode=sloppy_mode, update=update, nsub_max=nsub_max, magsettings=magsettings, sub_directories=sub_directories)
        used_calc_scheme=calc_scheme
        comp=self.structureDB[uid].get_composition(reduce=False)
        for element in comp:
            if element in reference_energy:
                E_ref=reference_energy[element]
            elif 'reference_energy' in magsettings:
                #allow to reference FM states to FM reference elements etc.
                E_ref=None
                for uidx in self.select(nelements=1,elements=[element]):
                    Er=self.get_energy_per_atom(uidx, calc_scheme, update=update, nsub_max=nsub_max, magsettings=magsettings['reference_energy'])
                    if (E_ref==None) or ((Er!=None) and (Er<E_ref)):
                        E_ref=Er
            else:
                E_ref=None
                if isinstance(calc_scheme,tuple):
                    for cs in calc_scheme:
                        Er=self.get_reference_energy_per_atom(element,calc_scheme=cs, sloppy_mode=sloppy_mode, update=update, nsub_max=nsub_max)
                        if (E_ref==None) or ((Er!=None) and (Er<E_ref)):
                            E_ref=Er
                else:
                    E_ref=self.get_reference_energy_per_atom(element,calc_scheme=calc_scheme, update=update, nsub_max=nsub_max)
            if (E_form!=None) and (E_ref!=None):
                E_form=E_form-E_ref*self.structureDB[uid].get_atfraction(element)
            else:
                E_form=None
        if returns!=[]:
            retval=[]
            for x in returns:
                if x=='calc_scheme':
                    retval.append(used_calc_scheme)
            return E_form,retval
        else: 
            return E_form
    
    
    def get_formation_enthalpy(self, uid, calc_scheme, sloppy_mode=True, update=False, reference_enthalpy={},returns=[]):
        """get the formation enthalpy per atom for compound with respect to decomposition
        into the elements for specified calc_scheme. 
        """
        used_calc_scheme=None
        E_form=None
        if isinstance(calc_scheme,tuple):
            css=calc_scheme
        else:
            css=[calc_scheme]
        for cs in css:
            prop_dict=self.get_properties(uid,calc_scheme, update=update)
            Ef=None
            if ('enthalpy' in prop_dict) and ('chemical_symbols' in prop_dict):
                Ef=prop_dict['enthalpy']/float(len(prop_dict['chemical_symbols']))
            if (E_form==None) or ((Ef!=None) and (Ef<E_form)):
                E_form=Ef
                used_calc_scheme=cs
        comp=self.structureDB[uid].get_composition(reduce=False)
        for element in comp:
            if element in reference_enthalpy:
                E_ref=reference_enthalpy[element]
            else:
                E_ref=None
                for cs in css:
                    Er=self.get_reference_enthalpy_per_atom(element,calc_scheme=cs, sloppy_mode=sloppy_mode, update=update)
                    if (E_ref==None) or ((Er!=None) and (Er<E_ref)):
                        E_ref=Er
            if (E_form!=None) and (E_ref!=None):
                E_form=E_form-E_ref*self.structureDB[uid].get_atfraction(element)
            else:
                E_form=None
        if returns!=[]:
            retval=[]
            for x in returns:
                if x=='calc_scheme':
                    retval.append(used_calc_scheme)
                return E_form,retval
        else: 
            return E_form
    
    def get_elastic(self, uid, calc_scheme, update=False, sloppy_mode=True, nsub_max=2, tau=["-0.03","-0.02","-0.01","0","0.01","0.02","0.03"],special='C11'):
        """returns elastic constants
        """
        B=None
        ao=self.get_atoms_object(uid, calc_scheme=calc_scheme)
        if ao==None:
            return None
        etab=[]
        ax=[]
        ex=[]
        for alpha in tau:
            a=float(alpha) #.split("_")[1])
            if alpha=="0":
                e=self.get_energy(uid, calc_scheme, update=update, sloppy_mode=sloppy_mode, nsub_max=nsub_max)
            else:
                e=self.get_energy(uid, calc_scheme, update=update, sloppy_mode=sloppy_mode, nsub_max=nsub_max, special='%s_%s'%(special,alpha))
            if e!=None:
                etab.append([a,e*get_conversion('eV','Hartree')])
                ax.append(a)
                ex.append(e)
        if len(etab)>4:
            parentdir=os.getcwd()
            os.chdir(os.path.join(uid,calc_scheme))
            fname="e_%s.dat"%special
            outfile=open(fname,"w")
            for (v,e) in etab:
                outfile.write("%8.5f %.8e\n"%((1.0+v),e/(3*ao.get_volume()))) #check
            outfile.close()
            exitcode, out = commands.getstatusoutput("Bulkmodul e_%s.dat -oe_%s.fit -V -A -d2"%(special,special))
            B=float(out.split("Bulkmodul=")[1].split("=")[1].split('[GPa')[0])
            os.chdir(parentdir)
            poly3 = np.poly1d(np.polyfit(ax, ex, 3))
            poly3d2 = np.polyder(poly3, 2)
            GPa=(units.kJ*10**-3)/(10**10)**3*10**9
            B=poly3d2(0.)/(ao.get_volume())/GPa
        if B!=None:
            if special=="C12":
                C11=self.get_elastic(uid, calc_scheme, update=update, sloppy_mode=sloppy_mode, nsub_max=nsub_max, special='C11')
                C22=self.get_elastic(uid, calc_scheme, update=update, sloppy_mode=sloppy_mode, nsub_max=nsub_max, special='C22')
                if (C11!=None) and (C22!=None):
                    B=0.5*(B-C11-C22)
            elif special=="C13":
                C11=self.get_elastic(uid, calc_scheme, update=update, sloppy_mode=sloppy_mode, nsub_max=nsub_max, special='C11')
                C33=self.get_elastic(uid, calc_scheme, update=update, sloppy_mode=sloppy_mode, nsub_max=nsub_max, special='C33')
                if (C11!=None) and (C33!=None):
                    B=0.5*(B-C11-C33)
            elif special=="C23":
                C22=self.get_elastic(uid, calc_scheme, update=update, sloppy_mode=sloppy_mode, nsub_max=nsub_max, special='C22')
                C33=self.get_elastic(uid, calc_scheme, update=update, sloppy_mode=sloppy_mode, nsub_max=nsub_max, special='C33')
                if (C22!=None) and (C33!=None):
                  B=0.5*(B-C22-C33)
        return B ##,poly3,poly3d2
        

    def get_elastic_constants(self, uid, calc_scheme, elements='tensor', magsettings={'submitted':True}, sub_directories={}, update=False):
        """Calculation of elastic constants adopted from D. Ohmer
        in /home/groups/da_tmm/HTE_Test/hte-v1.0/hte.py
        status: untested
        """
        C=None
        prop_dict=self.get_properties(uid, calc_scheme, magsettings=magsettings, sub_directories=sub_directories, update=update)
        if ('total_elastic_moduli' in prop_dict):
            C=prop_dict['total_elastic_moduli']
        return C
    
    def get_bulkmodulus(self, uid, calc_scheme, update=False, sloppy_mode=True, Vrel=["V_0.85","V_0.9","V_0.95","V_0.98","V_1.0","V_1.02","V_1.05","V_1.1","V_1.15"]):
        """returns bulkmodulus
        """
        B=None
        ao=self.get_atoms_object(uid, calc_scheme=calc_scheme)
        if ao==None:
            return None
        V_at=ao.get_volume()/len(ao)
        etab=[]
        for V in Vrel:
            e=self.get_energy(uid, calc_scheme, update=update, sloppy_mode=sloppy_mode, special=V)
            if e!=None:
                etab.append([V_at*float(V.split("_")[1]),e/len(ao)*get_conversion('eV','Hartree')])
        if len(etab)>6:
            parentdir=os.getcwd()
            os.chdir(os.path.join(uid,calc_scheme))
            outfile=open("e_vat.dat","w")
            for (v,e) in etab:
                outfile.write("%8.4f %12.8f\n"%(v,e))
            outfile.close()
            exitcode, out = commands.getstatusoutput("Bulkmodul e_vat.dat -oe_vat.fit -d5 -V -A")
            B=float(out.split("Bulkmodul=")[1].split("=")[1].split('[GPa')[0])
            os.chdir(parentdir)
        return B
    
    def get_magnetic_moment(self, uid, calc_scheme, magsettings={'submitted':True}, sub_directories={}, update=False):
        """return the total magnetic moment"""
        mom=None
        prop_dict=self.get_properties(uid, calc_scheme, magsettings=magsettings, sub_directories=sub_directories, update=update)
        if ('magnetic_moment' in prop_dict):
            mom=abs(prop_dict['magnetic_moment'])
        return mom

    def get_magnetic_moment_per_atom(self, uid, calc_scheme, magsettings={'submitted':True}, sub_directories={}, update=False):
        """return the total magnetic moment per atom"""
        mom=None
        prop_dict=self.get_properties(uid, calc_scheme, magsettings=magsettings, sub_directories=sub_directories, update=update)
        if ('magnetic_moment' in prop_dict) and ('chemical_symbols' in prop_dict):
            mom=abs(prop_dict['magnetic_moment'])/float(len(prop_dict['chemical_symbols']))
        return mom

    def get_magnetic_moments(self, uid, calc_scheme, magsettings={'submitted':True}, sub_directories={}):
        """return local magnetic moments (as list: [(atom,mom),...])"""
        print "check_point158, sub_directories are:", sub_directories
        prop_dict=self.get_properties(uid, calc_scheme, magsettings=magsettings, sub_directories=sub_directories)
        # print "check_point157, sub_directories are:", sub_directories
        # print "check_point153, get_magnetic_moments, prop_dict is :",prop_dict
        # print_dict = []
        # for key in prop_dict.keys():
            # print_dict.append(key)
        # print "check_point156, all prodict key is:", print_dict
        moms=[]
        
        # if prop_dict!={}:
        ### 这个地方明天再弄 标记一下 获得磁矩在有子目录的情况下有问题
 ##       for sub_dir in print_dict:
        # if ('magnetic_moments' in sub_dir) and ('chemical_symbols' in sub_dir):
        #     print "check_point154,entering condition"
        #     for i in range(len(prop_dict['chemical_symbols'])):
        #         moms.append((prop_dict['chemical_symbols'][i],prop_dict['magnetic_moments'][i]))
        # print "check_point155,moms is :",moms
        # return moms
         
    # def get_magnetic_moments(self, uid, calc_scheme, magsettings={'submitted':True}, sub_directories={}):
    #     """return local magnetic moments (as list: [(atom,mom),...])"""
    #     prop_dict=self.get_properties(uid, calc_scheme, magsettings=magsettings, sub_directories=sub_directories)
    #     print "check_point153,get_magnetic_moments,prop_dict is :",prop_dict
    #     moms=[]
    #     if ('magnetic_moments' in prop_dict[uid]) and ('chemical_symbols' in prop_dict[uid]):
    #         print "check_point154,entering condition"
    #         for i in range(len(prop_dict[uid]['chemical_symbols'])):
    #             moms.append((prop_dict[uid]['chemical_symbols'][i],prop_dict[uid]['magnetic_moments'][i]))
    #     print "check_point155,moms is :",moms
    #     return moms
            
    def get_magnetic_structure(self, uid, calc_scheme, latex=False, latex_symbols={"--":"","FM":"","AF":"$^{\dagger}$","NC":"$^{*}$"}, magsettings={'submitted':True}, sub_directories={}):
        """return magnetic structure type (NM, FM, AF, NC)"""
        prop_dict=self.get_properties(uid, calc_scheme, magsettings=magsettings, sub_directories=sub_directories)
        magtype="--"
        af=False
        if ('magnetic_moments' in prop_dict) and ('chemical_symbols' in prop_dict):
            for i in range(len(prop_dict['chemical_symbols'])):
                momi=np.array(prop_dict['magnetic_moments'][i])
                if (magtype=="--") and (np.dot(momi,momi)>0.01):
                    magtype="FM"
                eli=prop_dict['chemical_symbols'][i]
                for j in range(len(prop_dict['chemical_symbols'])):
                    momj=np.array(prop_dict['magnetic_moments'][j])
                    elj=prop_dict['chemical_symbols'][j]
                    if (isinstance(prop_dict['magnetic_moments'][i],list)) and (norm(np.cross(momi,momj))>0.1):
                        magtype="NC"
                    elif (ismagnetic(eli)) and (ismagnetic(elj)) and (np.dot(momi,momj)<-0.1):
                        magtype="AF"
                if magtype=="NC":
                    break
        if (latex==True) and (magtype in latex_symbols):
            return latex_symbols[magtype]
        return magtype
            
    def get_reference_energy_per_atom(self, element, sloppy_mode=False, calc_scheme=None, update=False, nsub_max=2):
        energy=None
        found=False
        if (element in self.tmpdata['reference_energies']):
            if calc_scheme in self.tmpdata['reference_energies'][element]:
                energy=self.tmpdata['reference_energies'][element][calc_scheme]
                if (update==False):
                    return energy
        else:
            self.tmpdata['reference_energies'][element]={}
        for uid in self.structureDB:
            comp=self.structureDB[uid].get_composition(reduce=False)
            if (len(comp)==1) and (element in comp):
                found=True
                e=self.get_energy_per_atom(uid, calc_scheme, sloppy_mode=sloppy_mode,update=update,nsub_max=nsub_max)
                if (e!=None) and ((energy==None) or (e<energy)):
                    energy=e
        if found==False:
            #try to get reference structure from ase data
            crystr=self.get_reference_state(element)
            if crystr!=None:
                source_info=[('get_reference_state',{'comment':'reference structure from ase data'})]
                uid='reference_ase_'+element
                self.add_structDB_entry(uid, atoms_obj=crystr, _source_info=source_info)
                calc=self.setup_calculator(uid,calc_scheme)
                self.structureDB[uid].get_energy_per_atom(calc_scheme,calc,update=update, job_commands=self.get_job_commands(calc_scheme))
        self.tmpdata['reference_energies'][element][calc_scheme]=energy
        return energy

    def get_reference_enthalpy_per_atom(self, element, sloppy_mode=True, calc_scheme=None, update=False,return_uid=False):
        enthalpy=None
        found=False
        uid_ref=None
        if not ('reference_enthalpies' in self.tmpdata):
            self.tmpdata['reference_enthalpies']={}
        if (element in self.tmpdata['reference_enthalpies']):
            if calc_scheme in self.tmpdata['reference_enthalpies'][element]:
                enthalpy,uid_ref=self.tmpdata['reference_enthalpies'][element][calc_scheme]
                if (update==False):
                    if return_uid==True:
                        return enthalpy,uid_ref
                    return enthalpy
        else:
            self.tmpdata['reference_enthalpies'][element]={}
        for uid in self.structureDB:
            comp=self.structureDB[uid].get_composition(reduce=False)
            if (len(comp)==1) and (element in comp):
                found=True
                H=None
                prop_dict=self.get_properties(uid,calc_scheme, update=update)
                if 'enthalpy' in prop_dict:
                    H=prop_dict['enthalpy']
                if (H!=None):
                    H=H/len(self.get_atoms_object(uid))
                    if (enthalpy==None) or (H<enthalpy):
                        enthalpy=H
                        uid_ref=uid
        if not ('reference_ase_'+element in self.structureDB): #found==False:
            #try to get reference structure from ase data
            crystr=self.get_reference_state(element)
            if crystr!=None:
                source_info=[('get_reference_state',{'comment':'reference structure from ase data'})]
                uid='reference_ase_'+element
                if (self.add_structDB_entry(uid, atoms_obj=crystr, _source_info=source_info)==True) and (update==True):
                    calc=self.setup_calculator(uid,calc_scheme)
                    self.structureDB[uid].get_energy_per_atom(calc_scheme,calc,update=update, job_commands=self.get_job_commands(calc_scheme))
        self.tmpdata['reference_enthalpies'][element][calc_scheme]=enthalpy,uid_ref
        if return_uid==True:
            return enthalpy,uid_ref
        return enthalpy

        

    def get_cif_composition(self,uid):
        if 'chemical_formula_sum' in self.structureDB[uid].atoms_initial.info:
            cfsumx=self.structureDB[uid].atoms_initial.info['chemical_formula_sum'] #.strip('\'').strip(' ').strip('[').strip(']')
            cfsum=''
            for i in range(len(cfsumx)):
                if (cfsumx[i].isalpha()) or (cfsumx[i].isdigit()):
                    cfsum=cfsum+cfsumx[i]
            #print cfsum
            i=0
            n=len(cfsum)
            comp={}
            while (i<n):
                #get element
                if (i+1<n) and (cfsum[i+1].isalpha()) and (cfsum[i+1]==cfsum[i+1].lower()):
                    el=cfsum[i:i+2]
                    i=i+2
                else:
                    el=cfsum[i:i+1]
                    i=i+1
                #get multiplicity
                j=i
                while (j<n) and (cfsum[i:j+1].isdigit()):
                    j=j+1
                if j==i:
                    mult=1
                else:
                    mult=int(cfsum[i:j])
                    i=j
                comp[el]=mult
            return comp
        else:
            return None
                
    def get_atomic_volume(self, element):
        # atomic volumes from experimental bulk structure
        at_volumes={'Mg':23.237, 'Si':20.021, 'Sn':34.189}
        if element in at_volumes:
            return at_volumes[element]
        # try to get atomic volume from ase data for reference structures
        Z=atomic_numbers[element]
        ref=reference_states[Z]
        if ref!=None:
            if 'symmetry' in ref:
                sym=reference_states[Z]['symmetry'].lower()
            else:
                sym=None
            if 'a' in ref:
                a=reference_states[Z]['a']
            else:
                a=None
            if 'c/a' in ref:
                coa=reference_states[Z]['c/a']
            else:
                coa=None
            if sym in ['sc', 'fcc', 'bcc', 'hcp', 'diamond']:
                crystruct=bulk(element, sym, a=a, covera=coa)
                atvol=crystruct.get_volume()/len(crystruct.arrays['numbers'])
                return atvol
        return None

    def get_reference_state(self, element):
        # try to get reference state form ase data and return atoms object
        Z=atomic_numbers[element]
        ref=reference_states[Z]
        if ref!=None:
            if 'symmetry' in ref:
                sym=reference_states[Z]['symmetry'].lower()
            else:
                sym=None
            if 'a' in ref:
                a=reference_states[Z]['a']
            else:
                a=None
            if 'c/a' in ref:
                coa=reference_states[Z]['c/a']
            else:
                coa=None
            if sym in ['sc', 'fcc', 'bcc', 'hcp', 'diamond']:
                sg={'sc':221, 'fcc':225,'bcc':229, 'hcp':194, 'diamond':227}
                crystruct=bulk(element, sym, a=a, covera=coa)
                crystruct.info['spacegroup']=Spacegroup(sg[sym])
                crystruct.info['_cell_length_a']=a
                crystruct.info['_cell_length_b']=a
                if coa!=None:
                    crystruct.info['_cell_length_c']=a*coa
                else:
                    crystruct.info['_cell_length_c']=a
                crystruct.info['_cell_angle_alpha']=90
                crystruct.info['_cell_angle_beta']=90
                if sym=='hcp':
                    crystruct.info['_cell_angle_gamma']=120
                else:
                    crystruct.info['_cell_angle_gamma']=90
                crystruct.info['_atom_site_type_symbol']=[element]
                if sym in ['sc', 'fcc', 'bcc']:
                    crystruct.info['_atom_site_fract_x']=[0.0]
                    crystruct.info['_atom_site_fract_y']=[0.0]
                    crystruct.info['_atom_site_fract_z']=[0.0]
                elif sym=='diamond':
                    crystruct.info['_atom_site_fract_x']=[1./8.]
                    crystruct.info['_atom_site_fract_y']=[1./8.]
                    crystruct.info['_atom_site_fract_z']=[1./8.]
                elif sym=='hcp':
                    crystruct.info['_atom_site_fract_x']=[1./3.]
                    crystruct.info['_atom_site_fract_y']=[2./3.]
                    crystruct.info['_atom_site_fract_z']=[1./4.]
                return crystruct
        return None

    def get_spacegroup_number(self, uid,calc_scheme=None, magsettings={'submitted':True}, symprec=1e-2):
        """returns the spacegroup number for database entry with given uid.
        """
        no=None
        if uid in self.structureDB:
            if calc_scheme==None:
                no=self.structureDB[uid].get_spacegroup_number()
            elif calc_scheme in self.calc_schemes:
                calcname,settings=self.calc_schemes[calc_scheme]
                if calcname=='fplo':
                    prop_dict=self.get_properties(uid,calc_scheme, magsettings=magsettings)
                    if 'data_path' in prop_dict:
                        no,sym=FPLO().get_spacegroup(calcdir=prop_dict['data_path'])
                else:
                    ao=self.get_atoms_object(uid,calc_scheme=calc_scheme)
                    if ao!=None:
                        try:
                            spglib_info=spglib.get_symmetry_dataset(ao, symprec=symprec)
                            no=spglib_info['number']
                        except:
                            no=None
        return no
        
    def get_structure_type(self, uid):
        """returns the structure type for database entry with given uid.
        NOTE: Changes of the structure type due to relaxations or substitutions
              are NOT considered, for instance if 'Cl' in 'NaCl' was substituted
              for 'Na' the output is still a 'NaCl' structure.
        """
        if uid in self.structureDB:
            return self.structureDB[uid].get_structure_type()
        return None

    def get_sg_international(self, uid):
        """
        Returns the space group in the international notation, (Hermann-Mauguin)
        e.g: hP8.
        """
        try:
            return self.structureDB[uid].spglib_info['international']
        except KeyError as ke:
            if not silent:
                warnings.warn("KeyError: {0} in 'get_sg_international' for uid {1}".format(ke, uid))
        except (AttributeError, TypeError):
            try:
                return spglib.get_symmetry_dataset(self.structureDB[uid].converged_ao)['international']
            except AttributeError as ae2:
                return spglib.get_symmetry_dataset(self.structureDB[uid].atoms_initial)['international']
                
    
    def get_structure_types(self, uid_list):
        """returns a list with structure types for database entries in uid_list
        """
        type_list=[]
        for uid in uid_list:
            strtype=self.get_structure_type(uid)
            if (strtype!=None) and (not(strtype in type_list)):
                type_list.append(strtype)
        return type_list
    
        
    def get_bandstructure(self, uid_list, calc_scheme, sympoints='FPLO', outtype='FPLO', update=False, results_dir=None, silent=False):
        if not (calc_scheme in self.calc_schemes):
            print 'get_bandstructure(): calc_scheme %s undefined, nothing done!'%calc_scheme
            return False
        parentdir=os.getcwd()
        if (results_dir!=None):
            if (not os.path.isdir(results_dir)):
                os.mkdir(results_dir)
            if outtype=='FPLO':
                fout=open(os.path.join(results_dir,'bands.tex'),'w')
                fout.write("\\documentclass{article}\n\\usepackage[dvips]{color}\n\\usepackage{graphicx}\n\\begin{document}\n")
                fout.close()
        trafo=None
        btrafo=None
        for uid in uid_list:
            statusout=uid+'('+calc_scheme+'): '
            if uid in self.structureDB:
                dbentry=self.structureDB[uid]
            else:
                print 'get_bandstructure(): unknown uid %s, abort'%uid
                return False
            if sympoints=='FPLO':
                bsdir='bs-symp-fplo'
            else:
                bsdir='bs'
                for ikp in range(len(sympoints)):
                    bsdir=bsdir+'-'+sympoints[ikp]['name'].strip()
            bspath=os.path.join(calc_scheme,bsdir)
            bs_ini_path=os.path.join(calc_scheme,'bs-ini')
            bs_fplo_file=os.path.join(dbentry.calcdir,bspath,'band.ps')
            bs_xny_file=os.path.join(dbentry.calcdir,bspath,'bands.dat')
            #
            bs_ready=False
            # check if we are already done:
            if ((outtype=='FPLO') and (os.path.isfile(bs_fplo_file))) or ((outtype=='xny') and (os.path.isfile(bs_xny_file))):
                bs_ready=True
                statusout=statusout+'band structure ready'
            else:
                # we still need to do something...
                calc_parent=self.setup_calculator(uid, calc_scheme)
                parent_converged=dbentry.check_convergency(calc_scheme,calc_parent)
                if (parent_converged==False) and (calc_parent!=None) and update:
                    dbentry.run_calculation(calc_scheme, calc_parent, job_commands=self.get_job_commands(calc_scheme))
                    statusout=statusout+'updating parent'
                if parent_converged:
                    statusout=statusout+'parent converged'
                    bs_ini_ready=False
                    fplo_files_ready=False
                    bs_calc_ready=False
                    #
                    if calc_parent.name.lower()=='vasp':
                        calc_bs_ini=self.setup_calculator(uid, calc_scheme, bandstructure_init=True)
                        if dbentry.check_convergency(bs_ini_path,calc_bs_ini):
                            bs_ini_ready=True
                            statusout=statusout+'/bsini converged'
                        elif (calc_bs_ini!=None) and update:
                            dbentry.run_calculation(bs_ini_path, calc_bs_ini, job_commands=self.get_job_commands(calc_scheme))
                            statusout=statusout+'/updating '+bs_ini_path
                    elif calc_parent.name.lower()=='fplo':
                        bs_ini_ready=True
                    #
                    calc_bs=self.setup_calculator(uid, calc_scheme, bandstructure=True)
                    if (bs_ini_ready) and (calc_bs!=None):
                        fplodir=os.path.join(dbentry.calcdir,calc_scheme,'tmp-fplo'+bsdir)
                        if (calc_bs.name.lower()=='fplo') or ((sympoints!='FPLO') and (outtype!='FPLO')):
                            fplo_files_ready=True
                        elif (os.path.isfile(os.path.join(fplodir,'+symminfo'))) and (os.path.isfile(os.path.join(fplodir,'=.in'))):
                            fplo_files_ready=True
                            statusout=statusout+'/fplo files ready'
                        else:
                            statusout=statusout+'/updating '+fplodir
                            if os.path.isdir(fplodir)==False:
                                os.mkdir(fplodir)
                            os.chdir(fplodir)
                            FPLO(job_settings=self.get_job_commands('fplo')).initialize(dbentry.atoms)
                            os.chdir(parentdir)
                            if (os.path.isfile(os.path.join(fplodir,'+symminfo'))) and (os.path.isfile(os.path.join(fplodir,'=.in'))):
                                fplo_files_ready=True
                        #
                        if (fplo_files_ready) and (dbentry.check_convergency(bspath,calc_bs)):
                            bs_calc_ready=True
                        else:
                            if (sympoints=='FPLO') and (calc_bs.name.lower()!='fplo'):
                                rec_latt=get_reciprocal_fplo(pathname=os.path.join(fplodir,'+tmp'), outfilename='+symminfo') #+symminfo only updated after FPLO run
                                if rec_latt!=None:
                                    symp=get_sympoints_fplo(pathname=fplodir, outfilename='=.in')
                                    trafo=inv(rec_latt)
                                else:
                                    symp=None
                            else:
                                symp=sympoints
                                trafo=None
                            if (symp!=None) and update:
                                dbentry.setup_calculation(bspath,calc_bs, bs_calc=True)
                                os.chdir(dbentry.calcdir)
                                if calc_bs.name.lower()=='vasp':
                                    shutil.copy(os.path.join(bs_ini_path,'CHGCAR'),bspath)
                                    shutil.copy(os.path.join(bs_ini_path,'CHG'),bspath)
                                    export_sympoints_to_vasp(symp,pathname=bspath,ktrafo=trafo)
                                    export_points_file(symp, pathname=bspath)
                                elif calc_bs.name.lower()=='fplo':
                                    shutil.copy(os.path.join(calc_scheme,'=.dens'),bspath)
                                os.chdir(parentdir)
                                if dbentry.run_calculation(bspath, calc_bs, job_commands=self.get_job_commands(calc_bs.name.lower())):
                                    statusout=statusout+'/updating '+bspath
                                else:
                                     statusout=statusout+'/update of '+bspath+' failed'
                                     print self.get_job_commands(calc_bs.name.lower())
                    if (bs_ini_ready) and (fplo_files_ready) and (bs_calc_ready):
                        # we have all necessary input to create a bandstructure
                        statusout=statusout+'/updating band structure'
                        if calc_bs.name.lower()=='vasp':
                            bs=get_vasp_bandstructure(pathname=os.path.join(dbentry.calcdir,bspath), pathname_sc=os.path.join(dbentry.calcdir,bs_ini_path))
                            if outtype=='FPLO':
                                rec_latt=get_reciprocal_fplo(pathname=os.path.join(fplodir,'+tmp'), outfilename='+symminfo')
                                btrafo=rec_latt
                                xnyfile=os.path.join(dbentry.calcdir,bspath,'+band')
                            else:
                                xnyfile=os.path.join(dbentry.calcdir,bspath,'bands.dat')
                                btrafo=None
                            bs_ready=write_bandstructure_xny(bs, xscale=1.0, unit='eV', E_Fermi_zero=True, shifte=0.0, ktrafo=btrafo, filename=xnyfile)
                        if outtype=='FPLO':
                            title=self.structureDB[uid].get_composition(reduce=True)
                            bs_ready=run_fplo_bandplot(calcdir=os.path.join(dbentry.calcdir,bspath),title=title,job_commands=self.get_job_commands('fplo'))
            if bs_ready and (results_dir!=None):
                if outtype=='FPLO':
                    shutil.copy(os.path.join(dbentry.calcdir,bspath,'band.ps'),os.path.join(results_dir,'bs-'+uid+'-'+calc_scheme+'.ps'))
                    fout=open(os.path.join(results_dir,'bands.tex'),'a')
                    fout.write("\\includegraphics[width=\\textwidth]{%s}\n\\newpage\n"%str('bs-'+uid+'.ps'))
                    fout.close()
            if silent==False:
                print statusout
        if (outtype=='FPLO') and (results_dir!=None):
            fout=open(os.path.join(results_dir,'bands.tex'),'a')
            fout.write("\\end{document}\n")
            fout.close()
            
                    
        if False:
            calc_bs=self.setup_calculator(uid, calc_scheme, bandstructure=True)
            if (calc_bs!=None) and (os.path.isdir(bspath)):
                if ((sympoints=='FPLO') or (outtype=='FPLO')) and (calc.name.lower()!='fplo'):
                    # create FPLO input files for band structure plot
                    fplodir=os.path.join(bspath,'tmp-fplo')
                    if os.path.isdir(fplodir)==False:
                        os.mkdir(fplodir)
                        os.chdir(fplodir)
                        FPLO(job_settings=self.get_job_commands('fplo')).initialize(dbentry.atoms)
                        os.chdir(parentdir)
                    rec_latt=get_reciprocal_fplo(pathname=fplodir, outfilename='+symminfo')
                #
                if dbentry.check_convergency(bspath,calc_bs):
                    #we have a bandstructure calculation and just need to plot it
                    if calc.name.lower()=='vasp':
                        scdir=os.path.join(dbentry.calcdir,calc_scheme,'bsini')
                        bs=get_vasp_bandstructure(pathname=bsdir, filename_eigenval='EIGENVAL', pathname_sc=scdir)
                        write_bandstructure_xny(bs, xscale=1.0, unit='eV', E_Fermi_zero=True, shifte=0.0, ktrafo=None, filename=os.path.join(bspath,'bands.dat'))
                        if outtype=='FPLO':
                            
                            trafo=inv(rec_latt)
                            btrafo=rec_latt
                        
            calc=self.setup_calculator(uid, calc_scheme)
            #check if parent calculation is converged
            if (calc!=None) and (dbentry.check_convergency(calc_scheme,calc)):
                calc_bs=self.setup_calculator(uid, calc_scheme, bandstructure=True)
                if ((sympoints=='FPLO') or (outtype=='FPLO')) and (calc.name.lower()!='fplo'):
                    # create FPLO input files for band structure plot
                    fplodir=os.path.join(dbentry.calcdir,calc_scheme,'tmp-fplo')
                    if os.path.isdir(fplodir)==False:
                        os.mkdir(fplodir)
                        os.chdir(fplodir)
                        FPLO(job_settings=self.get_job_commands('fplo')).initialize(dbentry.atoms)
                        os.chdir(parentdir)
                    rec_latt=get_reciprocal_fplo(pathname=fplodir, outfilename='+symminfo')
                #check if bs calculation is converged
                subdir=os.path.join(calc_scheme,'bs')
                if dbentry.check_convergency(subdir,calc):
                    scdir=os.path.join(dbentry.calcdir,calc_scheme)
                    bsdir=os.path.join(dbentry.calcdir,subdir)
                    if calc.name.lower()=='vasp':
                        bs=get_vasp_bandstructure(pathname=bsdir, filename_eigenval='EIGENVAL', pathname_sc=scdir, filename_outcar_sc='OUTCAR')
                        if outtype=='FPLO':
                            trafo=inv(rec_latt)
                            btrafo=rec_latt
                            xnyfile=os.path.join(bsdir,'+band')
                        else:
                            xnyfile=os.path.join(bsdir,'bands.dat')
                        write_bandstructure_xny(bs, xscale=1.0, unit='eV', E_Fermi_zero=True, shifte=0.0, ktrafo=btrafo, filename=xnyfile)
                    if outtype=='FPLO':
                        title=self.structureDB[uid].get_composition(reduce=True)+'('+str(scdir)+')'
                        run_fplo_bandplot(calcdir=bsdir,title=title,job_commands=self.get_job_commands('fplo'))
                elif update==True:
                    if (sympoints=='FPLO') and (calc.name.lower()!='fplo'):
                        # try to get sympoints from fplo output
                        if rec_latt!=None:
                            trafo=inv(rec_latt)
                            btrafo=rec_latt
                            symp=get_sympoints_fplo(pathname=fplodir, outfilename='=.in')
                        else:
                            symp=None
                    else:
                        symp=sympoints
                    if symp==None:
                        print 'get_bandstructure(): Could not get symmetry points for ',uid,', nothing done!'
                    else:
                        bsdir=os.path.join(dbentry.calcdir,calc_scheme,'bs')
                        if os.path.isdir(bsdir)==False:
                            os.mkdir(bsdir)
                            os.chdir(bsdir)
                            if calc.name.lower()=='vasp':
                                shutil.copy('../POSCAR','./')
                                shutil.copy('../POTCAR','./')
                                commands.getstatusoutput('cat ../INCAR |grep -v ISMEAR > INCAR')
                                shutil.copy('../CHGCAR','./')
                                shutil.copy('../CHG','./')
                                export_sympoints_to_vasp(symp,pathname='./',ktrafo=trafo)
                                export_points_file(symp, pathname='./')
                                commands.getstatusoutput('echo \' ICHARGE = 11 \' >> INCAR')
                            elif calc.name.lower()=='fplo':
                                shutil.copy('../=.dens','./')
                                calc.initialize(dbentry.atoms, bs_calc=True)
                            os.chdir(parentdir)
                        dbentry.run_calculation(subdir, calc, job_commands=self.get_job_commands(calc.name.lower()))
                        print 'bs not ready, updating'
                else:
                    print 'bs not ready, nothing done'
            else:
                if update==True:
                    print 'Running calculation for ',uid
                    dbentry.run_calculation(calc_scheme, calc, job_commands=self.get_job_commands(calc.name.lower()))
                else:
                    print 'Nothing done for ',uid

    def get_DOS(self, uid_list, calc_scheme, update=False, results_dir=None, dosfile_name=None, pDOS=[], comment="!", DOS_per_atom=False):
        DOS={}
        if (results_dir!=None) and (os.path.isdir(results_dir)==False):
            os.mkdir(results_dir)
        for uid in uid_list:
            if uid in self.structureDB:
                dbentry=self.structureDB[uid]
            else:
                print 'get_DOS(): unknown uid %s, abort'%uid
                return False
            factor=1.0
            if DOS_per_atom:
                factor=1.0/float(len(dbentry.atoms_initial))
            dosdir=os.path.join(dbentry.calcdir,calc_scheme,'bs-symp-fplo')
            calc=self.setup_calculator(uid, calc_scheme)
            dos=None
            if calc!=None:
                if calc.name.lower()=='fplo':
                    for sp in self.get_searchpaths():
                        dosdirx=os.path.join(sp,dbentry.calcdir,calc_scheme,'bs-symp-fplo')
                        dos=calc.get_DOS(calcdir=dosdirx,pDOS=pDOS,normalization=factor)
                        if dos!=None:
                            break
                else:
                    print "DOS for ",calc.name,' not implemented, abort!'
                    return None
                if update and (dos==None):
                    self.get_bandstructure([uid],calc_scheme, update=update)
                #output gle files:
                if (dos!=None) and (results_dir!=None):
                    if dosfile_name==None:
                        dosfile=os.path.join(results_dir,"dos-%s-%s-%s.dat"%(self.structureDB[uid].get_composition(reduce=True),uid,calc_scheme))
                    else:
                        dosfile=os.path.join(results_dir,dosfile_name)
                    fout=open(dosfile,"w")
                    fout.write("%s energy total"%comment)
                    for orbits in pDOS:
                        fout.write(" %s"%orbits)
                    fout.write("\n")
                    for i in range(len(dos['energy'])):
                        fout.write("%8.5f %8.5f"%(dos['energy'][i],dos['total'][i]))
                        for orbits in pDOS:
                            fout.write(" %8.5f"%dos[orbits][i])
                        fout.write("\n")
                        #todo: spinpolarization
                    fout.close()
            DOS[uid]=dos
        return DOS
                                     
                 
    def get_transport_properties(self, uid, calc_scheme, properties=['|S|','zT','PF'], T=300, tau=2E-14, kappa_l=2, dopings=['p','n'], dopingrange=[0,1.E20],ksdmin=80,interpolate=True):
        """Get transport properties (as requested in 'properties') as function of carrier
        concentration for database entry with unique identifier uid, using the calculational
        scheme(s) calc_scheme (may be also list).
        Scattering time tau and lattice thermal conductivity kappa_l are used as
        parameters to estimate zT and PF (see New.J.Phys. 15 (2013) 105010
        for details).
        Return value: dictionary tp with transport properties as function of 
                      carrier concentration for dopings if requested data is
                      available, otherwise None.
                      Contents of tp:
                       'p'/'n': sorted list with transport properties [(N_i,properties)]
                       'info': dictionary with information about parameters and
                               calculational settings
                      database entry with unique identifier uid
        Pararmeters:
         properties: list with transport properties to be calculated,
                     possible properties are:
                     S: Seebeck coeffient in [\mu V/K]
                     |S|: absolute value of Seebeck coeffient in [\mu V/K]
                     zT: estimated figure of merit
                     PF: power factor in [\mu W/K cm^2]
         T: temperature in K
         tau: scattering time in [s]
         kappa_l: lattice thermal conductivity in [W/(m K)]
         dopingrange: doping range in [e/cm^3]
         ksdmin: minimal value for k-space density
         interpolate: linear interpolation of values at boundary of doping range
                      (may be used to obtain values at fixed doping)
        """
        L=2.44E-8
        tp=None
        transpdata={}
        transpdata['p']=[]
        transpdata['n']=[]
        #allow for multiple calc_schemes
        if isinstance(calc_scheme,list):
            calc_schemes=calc_scheme
        else:
            calc_schemes=[calc_scheme]
        tracefile=None
        tpaths=[]
        for cs in calc_schemes:
            tpaths=self.get_converged_calculation_paths(uid, cs, property='transport',ksdmin=ksdmin)
            ao=self.get_atoms_object(uid,calc_scheme=cs)
            if (tpaths!=[]):
                tracefile=os.path.join(tpaths[0],'boltzgen','hte.trace')
                if ao!=None:
                    vol=ao.get_volume()*1.E-24
                else:
                    print "Strange: no atoms object, but trace file: %s"%tracefile
                    tracefile=None
            if tracefile!=None:
                break
        if tracefile!=None:
            tr=get_transport_boltztrap(filename=tracefile, Tmin=T-0.001,Tmax=T+0.001)
            for i in range(len(tr['N'])):
                if (tr['N'][i]>0.0):
                    doping='p'
                else:
                    doping='n'
                N=abs(tr['N'][i])
                conc=N/vol
                S=tr['Seebeck'][i]
                sigotau=tr['sigotau'][i]
                p=[conc]
                for prop in properties:
                    if prop=='S':
                        p.append(S)
                    elif prop=='|S|':
                        p.append(abs(S))
                    elif prop=='zT':
                        ZT=S*S*sigotau*tau*float(T)/(kappa_l+sigotau*tau*L*float(T))
                        p.append(ZT)
                    elif prop=='PF':
                        PF=S*S*sigotau*tau*1E4
                        p.append(PF)
		    elif prop=='sigotau':
		   	p.append(sigotau)
                transpdata[doping].append(p)
            #sort and interpolate
            tp={}
            for doping in dopings:
                tp[doping]=[]
                N0=None
                for p1 in sorted(transpdata[doping]):
                    N1=p1[0]
                    if interpolate and (N0!=None):
                        if (N0<dopingrange[0]) and (N1>dopingrange[0]):
                            ip=[]
                            for i in range(len(p)):
                                y=p0[i]+(p1[i]-p0[i])/(N1-N0)*(dopingrange[0]-N0)
                                ip.append(y)
                            tp[doping].append(ip)
                    if (N1<dopingrange[1]) and (N1>dopingrange[0]):
                        tp[doping].append(p1)
                    if interpolate and (N0!=None) and (dopingrange[0]<dopingrange[1]):
                        if (N0<dopingrange[1]) and (N1>dopingrange[1]):
                            ip=[]
                            for i in range(len(p)):
                                y=p0[i]+(p1[i]-p0[i])/(N1-N0)*(dopingrange[1]-N0)
                                ip.append(y)
                            tp[doping].append(ip)
                    p0=p1
                    N0=p0[0]
            info={}
            info['structure']="%s(%s), uid=%s"%(self.structureDB[uid].get_composition(reduce=True)
                                                ,self.get_structure_type(uid),uid)
            info['tracefile']=tracefile
            info['volume']=vol
            tp['info']=info
        return tp

    def get_zT_max(self, uid, calc_scheme, T=300, tau=2E-14, kappa_l=2, doping=['p','n'], dopingrange=[0,1.E20],ksdmin=80,returns=[],interpolate=True):
        """Estimate figure of merit zT from a Boltztrap calculation with scattering time tau
        and lattice thermal conductivity kappa_l as parameters (see New.J.Phys. 15 (2013) 105010
        for details).
        """
        if isinstance(doping,list):
            dopings=doping
        else:
            dopings=[doping]
        zTmax=None
        doplevel=None
        retval=[]
        zTdata=self.get_transport_properties(uid, calc_scheme, properties=['zT'], T=T, tau=tau, kappa_l=kappa_l, dopings=dopings, dopingrange=dopingrange,ksdmin=ksdmin,interpolate=interpolate)
        if zTdata!=None:
            for dop in dopings:
                for (N,zT) in zTdata[dop]:
                    if (zTmax==None) or (zTmax<zT):
                        zTmax=zT
                        doplevel=N
            for name in returns:
                if name=='dopinglevel':
                    retval.append(doplevel)
                elif name=='tracefile':
                    retval.append(zTdata['info']['tracefile'])
                elif name=='volume':
                    retval.append(zTdata['info']['volume'])
        if returns!=[]:
            return zTmax,retval
        return zTmax


    def get_PF_max(self, uid, calc_scheme, T=300, tau=2E-14, doping=['p','n'], dopingrange=[0,1E20],ksdmin=80,returns=[],interpolate=True):
        """Estimate power factor PF from a Boltztrap calculation with scattering time tau
        and as parameter (see New.J.Phys. 15 (2013) 105010 for details).
        """
        if isinstance(doping,list):
            dopings=doping
        else:
            dopings=[doping]
        PFmax=None
        doplevel=None
        retval=[]
        PFdata=self.get_transport_properties(uid, calc_scheme, properties=['PF'], T=T, tau=tau, dopings=dopings, dopingrange=dopingrange,ksdmin=ksdmin,interpolate=interpolate)
        
        if PFdata!=None:
            for dop in dopings:
                for (N,PF) in PFdata[dop]:
                    if (PFmax==None) or (PFmax<PF):
                        PFmax=PF
                        doplevel=N
            for name in returns:
                if name=='dopinglevel':
                    retval.append(doplevel)
                elif name=='tracefile':
                    retval.append(PFdata['info']['tracefile'])
                elif name=='volume':
                    retval.append(PFdata['info']['volume'])
        if returns!=[]:
            return PFmax,retval
        return PFmax

    
    def get_zT_max_old(self, uid, calc_scheme, T=300, tau=2E-14, kappa_l=2, doping='all', dopingrange=[0,1.E20],ksdmin=80,returns=[]):
        """Estimate figure of merit zT from a Boltztrap calculation with scattering time tau
        and lattice thermal conductivity kappa_l as parameters (see New.J.Phys. 15 (2013) 105010
        for details).
        """
        ZTmax=None
        doplevel=None
        vol=None
        L=2.44E-8
        #allow for multiple calc_schemes
        if isinstance(calc_scheme,list):
            calc_schemes=calc_scheme
        else:
            calc_schemes=[calc_scheme]
        tracefile=None
        tpaths=[]
        for cs in calc_schemes:
            ao=self.get_atoms_object(uid,calc_scheme=cs)
            if ao==None:
              ao=self.get_atoms_object(uid)
            tpaths=self.get_converged_calculation_paths(uid, cs, property='transport')
            #print cs,tpaths,ao
            if tpaths!=[]:
                ksdmax=0
                for tp in tpaths:
                    ksd=int(tp.split('transp-k')[1])
                    if ksd>ksdmax:
                        tpath=tp
                        ksdmax=ksd
                if (ksdmax>=ksdmin):
                    tracefile=os.path.join(tpath,'boltzgen','hte.trace')
            if tracefile!=None:
                vol=ao.get_volume()*1.E-24
                break
        if tracefile!=None:
            #found requested transport data 
            ZTmax=0.0
            doplevel=0.0
            tr=get_transport_boltztrap(filename=tracefile, Tmin=T-0.001,Tmax=T+0.001)
            for i in range(len(tr['N'])):
                N=tr['N'][i]
                conc=abs(N/vol)
                if (doping=='all') or ((doping=='p') and (N>0.0)) or ((doping=='n') and (N<0.0)):
                    if (dopingrange==None) or ((conc>dopingrange[0]) and (conc<dopingrange[1])):
                        S=abs(tr['Seebeck'][i])
                        sigotau=tr['sigotau'][i]
                        ZT=S*S*sigotau*tau*float(T)/(kappa_l+sigotau*tau*L*float(T))
                        if (ZT>ZTmax):
                            ZTmax=ZT
                            doplevel=conc
        retval=[]
        for name in returns:
            if name=='dopinglevel':
                retval.append(doplevel)
            elif name=='tracefile':
                retval.append(tracefile)
            elif name=='volume':
                retval.append(vol)
        if retval!=[]:
            return ZTmax,retval
        return ZTmax
            
    def get_PF_max_old(self, uid, calc_scheme, T=300, tau=2E-14, doping='all', dopingrange=[0,1E20],ksdmin=80,returns=[]):
        """Estimate power factor PF from a Boltztrap calculation with scattering time tau
        and as parameter (see New.J.Phys. 15 (2013) 105010 for details).
        """
        PFmax=None
        doplevel=None
        vol=None
        L=2.44E-8
        #allow for multiple calc_schemes
        if isinstance(calc_scheme,list):
            calc_schemes=calc_scheme
        else:
            calc_schemes=[calc_scheme]
        tracefile=None
        for cs in calc_schemes:
            ao=self.get_atoms_object(uid,calc_scheme=cs)
            if ao==None:
              ao=self.get_atoms_object(uid)
            tpaths=self.get_converged_calculation_paths(uid, cs, property='transport')
            if (ao!=None) and (tpaths!=[]):
                ksdmax=0
                for tp in tpaths:
                    ksd=int(tp.split('transp-k')[1])
                    if ksd>ksdmax:
                        tpath=tp
                        ksdmax=ksd
                if (ksdmax>=ksdmin):
                    tracefile=os.path.join(tpath,'boltzgen','hte.trace')
            if tracefile!=None:
                vol=ao.get_volume()*1.E-24
                break
        if tracefile!=None:
            #found requested transport data 
            PFmax=0.0
            doplevel=0.0
            tr=get_transport_boltztrap(filename=tracefile, Tmin=T-0.001,Tmax=T+0.001)
            for i in range(len(tr['N'])):
                N=tr['N'][i]
                conc=abs(N/vol)
                if (doping=='all') or ((doping=='p') and (N>0.0)) or ((doping=='n') and (N<0.0)):
                    if (dopingrange==None) or ((conc>dopingrange[0]) and (conc<dopingrange[1])):
                        S=abs(tr['Seebeck'][i])
                        sigotau=tr['sigotau'][i]
                        PF=S*S*sigotau*tau
                        if (PF>PFmax):
                            PFmax=PF
                            doplevel=conc
        retval=[]
        for name in returns:
            if name=='dopinglevel':
                retval.append(doplevel)
            elif name=='tracefile':
                retval.append(tracefile)
            elif name=='volume':
                retval.append(vol)
        if retval!=[]:
            return PFmax,retval
        return PFmax
    
    def get_transport(self, uid_list, calc_scheme, lab=None, update=False, plotit=False, kspace_density=None, boltz_generic=True, results_dir=None, cmp_schemes=[], xaxis='E_Fermi', props=['Seebeck', 'DOS', 'specific_heat'], debug=False):
        # run transport calculations
        # debug: export wien struct file
        uids_ready=[]
        parentdir=os.getcwd()
        for uid in uid_list:
            print 'Transport for',uid
            dbentry=self.structureDB[uid]
            if kspace_density==None:
                transpdir="transp-ini"
            else:
                transpdir="transp-k%d"%kspace_density
            if boltz_generic==True:
                boltzdir="boltzgen"
            else:
                boltzdir="boltzw2k"
            # check if Boltztrap calculation is ready..somewhere
            for sp in self.get_searchpaths():
                tracefile=os.path.join(sp,dbentry.calcdir, calc_scheme, transpdir,boltzdir,'hte.trace')
                condtensfile=os.path.join(sp,dbentry.calcdir, calc_scheme, transpdir,boltzdir,'hte.condtens')
                if (os.path.isfile(tracefile)) and (os.path.isfile(condtensfile)):
                    print "... Boltztrap calculation ready ",uid
                    uids_ready.append(uid)
                    if plotit==True:
                        # props=['Seebeck', 'DOS', 'specific_heat']
                        for prop in props:
                            transp=get_transport_boltztrap(pathname=os.path.join(dbentry.calcdir, calc_scheme, transpdir,boltzdir),filename='hte.trace')
                            if lab==None:
                                label=calc_scheme+transpdir
                            else:
                                label=lab
                            plot(transp[xaxis],transp[prop],label=label,linewidth=2.0)
                            suptitle(prop+' '+dbentry.get_composition(reduce=True))
                            for i in range(len(cmp_schemes)):
                                if 'dir' in cmp_schemes[i]:
                                    dir=os.path.join(sp,dbentry.calcdir, cmp_schemes[i]['dir'])
                                    transp=get_transport_boltztrap(pathname=dir,filename='hte.trace')
                                else:
                                    if 'calc_scheme' in cmp_schemes[i]:
                                        cs=cmp_schemes[i]['calc_scheme']
                                    else:
                                        cs=calc_scheme
                                    if 'kspace_density' in cmp_schemes[i]:
                                        ksd="transp-k%d"%cmp_schemes[i]['kspace_density']
                                    else:
                                        ksd="transp-k%d"%kspace_density
                                    dir=os.path.join(dbentry.calcdir, cs, ksd, boltzdir)
                                    transp=get_transport_boltztrap(pathname=os.path.join(dbentry.calcdir, cs, ksd, boltzdir),filename='hte.trace')
                                if 'label' in  cmp_schemes[i]:
                                    label=cmp_schemes[i]['label']
                                else:
                                    label=dir
                                if transp!=None:
                                    plot(transp[xaxis],transp[prop],label=label)
                            legend()
                            xlabel(xaxis)
                        #axis([-0.2,0.2,-0.0005,0.0005])
                        show()
                    elif (results_dir!=None) and (os.path.isdir(results_dir)):
                        fname="%s-%s-kd%d.trace"%(uid,calc_scheme,kspace_density)
                        shutil.copy(tracefile,os.path.join(results_dir,fname))
                        fname="%s-%s-kd%d.condtens"%(uid,calc_scheme,kspace_density)
                        shutil.copy(condtensfile,os.path.join(results_dir,fname))
                    break
            if (update==True) and (uid not in uids_ready):
                # run set of calculations to obtain transport properties
                calc_parent=self.setup_calculator(uid, calc_scheme)
                calc_tr=self.setup_calculator(uid, calc_scheme, transport=True, transport_kspace_density=kspace_density)
                #calc_initr=self.setup_calculator(uid, calc_scheme, transport=True)
                calc_initr=self.setup_calculator(uid, calc_scheme, bandstructure_init=True)
                # check if transport calculation is ready
                if dbentry.check_convergency(os.path.join(calc_scheme, transpdir),calc_tr):
                    # run Boltztrap (copy boltzini, write bs, run)
                    bspath=os.path.join(dbentry.calcdir,calc_scheme, transpdir)
                    bs=get_vasp_bandstructure(pathname=bspath)
                    boltzpath=os.path.join(bspath,boltzdir)
                    write_intrans_boltztrap(boltzdir=boltzpath, boltz_generic=boltz_generic, n_electrons=bs['n_electrons']) 
                    shutil.copy('BoltzTraP.def',boltzpath)
                    if boltz_generic==True:
                        tmpatoms=self.get_atoms_object(uid,calc_scheme=calc_scheme)
                        #write_structure_boltztrap(dbentry.atoms_initial,pathname=boltzpath,filename="hte.struct")
                        write_structure_boltztrap(tmpatoms,pathname=boltzpath,filename="hte.struct")
                        write_bandstructure_boltztrap(bs,pathname=boltzpath, filename="hte.energy", ktrafo=None, runBoltz=True)
                    else:
                        print 'boltzw2k still to do, nothing done'
                # check if transport calculation was already set up
                elif os.path.isdir(os.path.join(dbentry.calcdir, calc_scheme, transpdir)):
                    print 'transport calculation not converged, updating',uid
                    dbentry.run_calculation(os.path.join(calc_scheme, transpdir),calc_tr, job_commands=self.get_job_commands(calc_scheme))
                    # check if ini transport calculation is converged
                elif dbentry.check_convergency(os.path.join(calc_scheme, 'transp-ini'),calc_initr):
                    # initialize transport calculation
                    print 'transport ini calculation converged, running vasp transport'
                    dbentry.setup_calculation(os.path.join(calc_scheme, transpdir),calc_tr)
                    shutil.copy(os.path.join(dbentry.calcdir,calc_scheme,'transp-ini','CHGCAR'),os.path.join(dbentry.calcdir,calc_scheme, transpdir))
                    dbentry.run_calculation(os.path.join(calc_scheme, transpdir),calc_tr, job_commands=self.get_job_commands(calc_scheme))
                    # check if ini transport was already set up and resubmit 
                elif os.path.isdir(os.path.join(dbentry.calcdir, calc_scheme, 'transp-ini')):
                    print 'transport ini calculation not converged, updating'
                    dbentry.run_calculation(os.path.join(calc_scheme, 'transp-ini'), calc_initr, job_commands=self.get_job_commands(calc_initr.name.lower()))
                    # check if parent calculation is converged and initialize ini transport calculation
                elif dbentry.check_convergency(calc_scheme,calc_parent):
                    print 'parent calculation converged, running vasp ini transport'
                    dbentry.setup_calculation(os.path.join(calc_scheme, 'transp-ini'),calc_initr)
                    dbentry.run_calculation(os.path.join(calc_scheme, 'transp-ini'),calc_initr, job_commands=self.get_job_commands(calc_initr.name.lower()))
                else:
                    print 'parent not converged, updating'
                    calc_parent=self.setup_calculator(uid, calc_scheme)
                    dbentry.run_calculation(calc_scheme, calc_parent, job_commands=self.get_job_commands(calc_parent.name.lower()))
            if debug:
                #use with care, double submission if job is still in queue
                wien_dir=os.path.join(dbentry.calcdir,'wien2k')
                if not os.path.isdir(wien_dir):
                    os.mkdir(wien_dir)
                    wien_dir_prim=os.path.join(wien_dir,'prim')
                    os.mkdir(wien_dir_prim)
                    wien_structfile=os.path.join(wien_dir_prim,'prim.struct')
                    write_struct(wien_structfile,atoms2=dbentry.atoms_initial)
                    os.chdir(wien_dir_prim)
                    commands.getstatusoutput("init_lapw -b")
                    shutil.copy('prim.struct_sgroup','../wien2k.struct')
                    os.chdir(parentdir)
                    if os.path.isdir('w2k_jobfiles/'):
                        flist=glob.glob('w2k_jobfiles/*')
                        for fname in flist:
                            shutil.copy(fname,wien_dir)
                        os.chdir(wien_dir)
                        commands.getstatusoutput("init_lapw -b")
                        exitcode, out = commands.getstatusoutput("qsub job-wien.sh")
                        if exitcode == 0:
                            outspl=out.split()
                            jobid=outspl[2]
                            cmd="echo %s > jobid"%str(jobid)
                            print 'Submitted wien ini job for uid ',uid,' jobid:',jobid #,cmd
                            commands.getstatusoutput(cmd)
                        os.chdir(parentdir)
                        print 'Submitted wien ini job for uid ',uid
                elif not os.path.isdir(os.path.join(wien_dir,'ini')):
                    #check if job is still running
                    jobid_file=os.path.join(wien_dir,'jobid')
                    running=False
                    if os.path.isfile(jobid_file):
                        infile=open(jobid_file,'r')
                        jobid=infile.readline().strip()
                        infile.close()
                        qstatcmd='qstat -j '+jobid #todo
                        exitcode, out = commands.getstatusoutput(qstatcmd)
                        if exitcode==0:
                            running=True
                            print 'wien ini calculation running, jobid is',jobid
                    if running==False:
                        #shutil.copy('w2k_jobfiles/*',wien_dir)
                        os.chdir(wien_dir)
                        commands.getstatusoutput("save_lapw -d ini")
                        commands.getstatusoutput("clean_lapw -s")
                        commands.getstatusoutput("init_lapw -b -s kgen -e kgen -numk 12000")
                        exitcode, out = commands.getstatusoutput("qsub job-wien.sh")
                        if exitcode == 0:
                            outspl=out.split()
                            jobid=outspl[2]
                            cmd="echo %s > jobid"%str(jobid)
                            print 'Submitted wien transp-k12000 job for uid ',uid,' jobid:',jobid #,cmd
                            commands.getstatusoutput(cmd)
                        os.chdir(parentdir)
                elif not os.path.isdir(os.path.join(wien_dir,'transp-k12000')):
                    #check if job is still running
                    jobid_file=os.path.join(wien_dir,'jobid')
                    running=False
                    if os.path.isfile(jobid_file):
                        infile=open(jobid_file,'r')
                        jobid=infile.readline().strip()
                        infile.close()
                        qstatcmd='qstat -j '+jobid #todo
                        exitcode, out = commands.getstatusoutput(qstatcmd)
                        if exitcode==0:
                            running=True
                            print 'wien transp-k12000 calculation running, jobid is',jobid
                    if running==False:
                        #save and clean up
                        os.chdir(wien_dir)
                        commands.getstatusoutput("save_lapw -d transp-k12000")
                        shutil.copy('wien2k.energy','transp-k12000/')
                        commands.getstatusoutput("clean_lapw -s")
                        commands.getstatusoutput("init_lapw -b -s kgen -e kgen -numk 24000")
                        exitcode, out = commands.getstatusoutput("qsub job-wien.sh")
                        if exitcode == 0:
                            outspl=out.split()
                            jobid=outspl[2]
                            cmd="echo %s > jobid"%str(jobid)
                            print 'Submitted wien transp-k24000 job for uid ',uid,' jobid:',jobid
                            commands.getstatusoutput(cmd)
                        os.chdir('transp-k12000/')
                        bs=get_wien_bandstructure(pathname='./', filename='wien2k.energy', scf_file='wien2k.scf')
                        if bs!=None:
                            print 'Running Boltztrap for uid ',uid
                            boltzdir='boltzw2k'
                            write_intrans_boltztrap(boltzdir=boltzdir,filename='hte.intrans', boltz_generic=False, E_Fermi=0.0, n_electrons=bs['n_electrons'])
                            shutil.copy(os.path.join(parentdir,'BoltzTraP.def'),boltzdir)
                            shutil.copy('wien2k.struct',os.path.join(boltzdir,'hte.struct'))
                            write_bandstructure_wien(bs, pathname=boltzdir, filename="hte.energy", runBoltz=True) 
                        os.chdir(parentdir)
                elif not os.path.isdir(os.path.join(wien_dir,'transp-k24000')):
                    #check if job is still running
                    jobid_file=os.path.join(wien_dir,'jobid')
                    running=False
                    if os.path.isfile(jobid_file):
                        infile=open(jobid_file,'r')
                        jobid=infile.readline().strip()
                        infile.close()
                        qstatcmd='qstat -j '+jobid
                        exitcode, out = commands.getstatusoutput(qstatcmd)
                        if exitcode==0:
                            running=True
                            print 'wien transp-k24000 calculation running, jobid is',jobid
                    if running==False:
                        #save and clean up
                        os.chdir(wien_dir)
                        commands.getstatusoutput("save_lapw -d transp-k24000")
                        shutil.copy('wien2k.energy','transp-k24000/')
                        commands.getstatusoutput("clean_lapw -s")
                        os.chdir('transp-k24000/')
                        bs=get_wien_bandstructure(pathname='./', filename='wien2k.energy', scf_file='wien2k.scf')
                        if bs!=None:
                            print 'Running Boltztrap for uid ',uid
                            boltzdir='boltzw2k'
                            write_intrans_boltztrap(boltzdir=boltzdir,filename='hte.intrans', boltz_generic=False, E_Fermi=0.0, n_electrons=bs['n_electrons'])
                            shutil.copy(os.path.join(parentdir,'BoltzTraP.def'),boltzdir)
                            shutil.copy('wien2k.struct',os.path.join(boltzdir,'hte.struct'))
                            write_bandstructure_wien(bs, pathname=boltzdir, filename="hte.energy", runBoltz=True) 
                        os.chdir(parentdir)
        return uids_ready

                   

    def get_source_info_old_version(self,uid,structureDB):
        """try to reconstruct source information from old version
        """
        source_info_failed=[('non_HTE',{'comment':'upgrade from version 0.x to %s failed, original comment: %s'%(self.get_version(),structureDB[uid].comment)})]
        substitutions=uid.split('subst_')
        uid_org=substitutions[0]
        del substitutions[0]
        if uid_org in structureDB:
            dbentry=structureDB[uid_org]
            if 'reference structure from ase data' in dbentry.comment:
                source_info=[('get_reference_state',
                              {'comment':'reference structure from ase data'})]
            elif 'imported from ciffile <' in dbentry.comment:
                try:
                    ciffile=dbentry.comment.split('<')[1].split('>')[0]
                    index=int(dbentry.comment.split('>, entry ')[1].split()[0])
                    cif_blocks=split_cif(ciffile)
                    (cif_uid,block)=cif_blocks[index]
                    cif_sources={'cif_file_external':os.path.abspath(ciffile),
                                 'cif_index_external':index,
                                 'cif_id':cif_uid}
                    if self.get_cif_repository()!=None:
                        outfile=open(os.path.join(self.get_cif_repository(),
                                                  uid+'.cif'),"w")
                        outfile.write(block)
                        outfile.close()
                        cif_sources['cif_file']=os.path.join(self.get_cif_repository(),
                                                             uid+'.cif')
                    source_info=[('import_cif',cif_sources)]
                except:
                    return source_info_failed
            else:
                return source_info_failed
        uid_history=uid_org
        for subst in substitutions:
            atlist={}
            i=0
            el1=None
            el2=None
            n=len(subst)
            while i<n:
                #print i,subst[i:]
                if (el1==None):
                    if subst[i+1]=='-':
                        el1=subst[i]
                        i=i+2
                    else:
                        el1=subst[i:i+2]
                        i=i+3
                else:
                    if (i+1<n) and (subst[i+1].isalpha) and (subst[i+1]==subst[i+1].lower()):
                        el2=subst[i:i+2]
                        i=i+2
                    else:
                         el2=subst[i]
                         i=i+1
                    atlist[el1]=el2
                    el1=None
                    el2=None
                #print atlist,el1,el2
            source_info.append(('substitute_atoms',{'atom_list':[atlist],
                                                        'initial_uid':uid_history}))
            uid_history=uid_history+'subst_'+subst
        source_info.append([('upgraded from version 0.x'),{'comment':'original comment: '+structureDB[uid].comment}])
        return source_info

    def get_distance_from_convex_hull(self,uid,calc_scheme, magsettings={'submitted':True}, sub_directories={}, reference_potential='auto',sys=[], chull_object=None):
        """get distance from calculated convex hull for given uid
        reference_potential: 'energy' or 'enthalpy', default: 'auto' uses enthalpy if pstress tag is in calc_scheme, otherwise energy
        chull_object: default: None (determine chull), otherwise chull obeject is used
        """
        dEh=None
        if (not (uid in self.structureDB)):
            return dEh
        elif (not (calc_scheme in self.calc_schemes)) and (not (isinstance(calc_scheme,tuple))):
            return dEh
        if reference_potential=='auto':
            if isinstance(calc_scheme,tuple):
                calculator_name, settings=self.calc_schemes[calc_scheme[0]]
            else:
                calculator_name, settings=self.calc_schemes[calc_scheme]
            if 'pstress' in settings:
                ref_potential='enthalpy'
            else:
                ref_potential='energy'
        elif (reference_potential=='energy') or (reference_potential=='enthalpy'):
            ref_potential=reference_potential
        else:
            return dEh
        comp=self.structureDB[uid].get_composition(reduce=False)
        system=[]
        for el in comp:
            system.append(el)
        system=sorted(system)
        if sys!=[]:
            system=sorted(sys) #workaround to get chull if no compounds are in DB (e.g. to compare with exp. reported phases)
        strsys=str(system)
        strcs=str(calc_scheme)
        if (chull_object!=None):
            chull=chull_object
        elif (strsys in self.tmpdata['chulls']) and (strcs in self.tmpdata['chulls'][strsys]):
            chull=self.tmpdata['chulls'][strsys][strcs]
        else:
            sel=self.select(elements=system, nelements_max=len(system))
            dataset=[]
            for uidx in sel:
                datapoint=[]
                for el in system:
                    datapoint.append(self.structureDB[uidx].get_atfraction(el))
                if ref_potential=='enthalpy':
                    ef=self.get_formation_enthalpy(uidx,calc_scheme)
                else:
                    ef=self.get_formation_energy(uidx,calc_scheme)
                if ef!=None:
                    p=tuple(datapoint),ef,uidx
                    dataset.append(p)
            chull=CHULL(len(system))
            chull.get_convex_hull(dataset)
            if not (strsys in self.tmpdata['chulls']):
                    self.tmpdata['chulls'][strsys]={}
            self.tmpdata['chulls'][strsys][strcs]=chull
        datapoint=[]
        for el in system:
            datapoint.append(self.structureDB[uid].get_atfraction(el))
        x=tuple(datapoint)
        if ref_potential=='enthalpy':
            ef=self.get_formation_enthalpy(uid,calc_scheme)
        else:
            ef=self.get_formation_energy(uid,calc_scheme, magsettings=magsettings, sub_directories=sub_directories)
        if ef!=None:
            yhull,supp=chull.get_value(x)
            return ef-yhull
        return None
   
    def get_chull(self,sel,calc_scheme, system, magsettings={'submitted':True}, sub_directories={}, reference_potential='auto'):
        """returns convex hull (chull object) calculated from uids in sel
        can be used to get convex hull of limited set of structures (e.g. exp or FM or ambient pressure phases)
        system: list of elements in system (order matters!)
        reference_potential: 'energy' or 'enthalpy', default: 'auto' uses enthalpy if pstress tag is in calc_scheme, otherwise energy
        """
        if reference_potential=='auto':
            if isinstance(calc_scheme,tuple):
                calculator_name, settings=self.calc_schemes[calc_scheme[0]]
            else:
                calculator_name, settings=self.calc_schemes[calc_scheme]
            if 'pstress' in settings:
                ref_potential='enthalpy'
            else:
                ref_potential='energy'
        elif (reference_potential=='energy') or (reference_potential=='enthalpy'):
            ref_potential=reference_potential
        else:
            return None
        dataset=[]
        for uidx in sel:
            datapoint=[]
            for el in system:
                    datapoint.append(self.structureDB[uidx].get_atfraction(el))
            if ref_potential=='enthalpy':
                ef=self.get_formation_enthalpy(uidx,calc_scheme)
            else:
                ef=self.get_formation_energy(uidx,calc_scheme, magsettings=magsettings, sub_directories=sub_directories)
            if ef!=None:
                p=tuple(datapoint),ef,uidx
                dataset.append(p)
        chull=CHULL(len(system))
        chull.get_convex_hull(dataset)
        return chull

    def run_calculation(self,uid,calc_scheme, rerun=False, options={}, magsettings={'submitted':True}, sub_directories={}, nsub_max=-1):
        """run HTE calculation for structure database entry with unique identifier 'uid' and calculational scheme 'calc_scheme'
        optional arguments:
            rerun: (not ready)
            options: (not ready)
            magsettings: define magnetic structures, see setup_magnetic_structures() for details
            sub_directories: (not ready, remove?) 
        """
        print "check_point115, entering hte.run_calculation"
        if not ('qstat' in self.tmpdata):
            exitcode, out = commands.getstatusoutput("squeue|wc -l")
            self.tmpdata['qstat']=(int(out),0)
            print "check_point118,queue status:",self.tmpdata['qstat']
        nqueue,nsub=self.tmpdata['qstat']
        #print self.tmpdata['qstat'],nsub
        if nsub>50:
            exitcode, out = commands.getstatusoutput("squeue|wc -l")
            self.tmpdata['qstat']=(int(out),0)
            print "check_point119,queue status:",self.tmpdata['qstat']
            nqueue,nsub=self.tmpdata['qstat']
        self.tmpdata['qstat']=(nqueue,nsub+1)
        if ((nqueue+nsub)>self.max_jobs_in_queue['total']):
            print "check_point120,return True:",
            self.add_logmessage("run_calculation(): %s(%s) not submitted because queue is full"%(uid,calc_scheme))
            return True
        if nsub_max<0:
            nsub_max=self.nsub_max
        parentdir=os.getcwd()
        if sub_directories!={}:
            subdirs=sub_directories
            print "check_point121,subdirs:",subdirs
        else:
            subdirs={}
            magconfigs=self.setup_magnetic_structures(uid,calc_scheme, magsettings=magsettings)
            for magconf in magconfigs:
                subdir=os.path.join(calc_scheme,magconf)
                subdirs[subdir]=magconfigs[magconf]
        sp=self.get_scratch_directory()
        print "check_point124,sp=:",sp
        if sp==None:
            print "check_point122,sp==None:"
            return False
        elif (rerun==True):
            sp=os.path.join(sp,"rerun")
            print "check_point123,sp=:",sp
        if not (os.path.isdir(sp)):
            os.mkdir(sp)
        os.chdir(sp)
        calc_name,settings=self.calc_schemes[calc_scheme]
        print "check_point125, calc_name,settings are:",calc_name,settings
        for subdir in subdirs:
            print "check_point126, entering subdir"
            calc,settings=self.setup_calculator(uid,calc_scheme,afm=subdirs[subdir],update=True, options=options,return_settings=True)
            print "check_point127, calc,settings in subdir are:", calc,settings
            if (calc!=None) and (not (os.path.isfile(os.path.join(uid,subdir,"POSCAR")))) and (calc_name.lower()=='vasp'):
                print "check_point128,calc!=None"
                ao=self.structureDB[uid].atoms
                submit=True
                if ('natoms_max' in self.max_jobs_in_queue):
                    for natmax in self.max_jobs_in_queue['natoms_max']:
                        if (len(ao)>natmax) and ((nqueue+nsub)>self.max_jobs_in_queue['natoms_max'][natmax]):
                           submit=False
                           break
                if submit==False:
                    self.add_logmessage("run_calculation(): %s(%s) not submitted because queue is full (nat=%d)"%(uid,calc_scheme,len(ao)))
                    continue
                structure_info={'cell':ao.get_cell(),'chemical_symbols':ao.get_chemical_symbols(),'scaled_positions':ao.get_scaled_positions()}
                print "check_point117, structure_info:",structure_info
                setup_vasp_calculation(structure_info,settings,pathname=os.path.join(uid,subdir))
            try:
                energy=self.structureDB[uid].get_energy(subdir,calc, sloppy_mode=False, update=True, job_commands=self.get_job_commands(calc_scheme),nsub_max=nsub_max, settings=subdirs[subdir])
            except:
                self.add_logmessage("WARNING(run_calculation): Failed to run calculation in %s"%os.path.join(sp,uid,subdir))
            if (rerun==True) and (energy!=None):
                calc_name,settings=self.calc_schemes[calc_scheme]
                calcpath=os.path.join(sp,uid,subdir)
                if calc_name.lower()=='vasp':
                    pdict=get_properties_vasp(calcdir=calcpath, convert_strings=True)
                else:
                    pdict={}
                if pdict!={}:
                    self.tmpdata['prop_dict'][uid][subdir]=pdict
        os.chdir(parentdir)
        return True
        

    def submit_array_jobs(self):
        parentdir=os.getcwd()
        for d in glob.glob("tmp-array-*/"):
            os.chdir(d)
            if (os.path.isfile("njobs")):
                infile=open("njobs", 'r')
                njobs=int(infile.readline().strip("\n"))
                infile.close()
            else:
                os.chdir(parentdir)
                continue
            exitcode, out = commands.getstatusoutput("sbatch job-array.sh")
            print d,exitcode, out
            #for the moment slurm only
            outspl=out.split('Submitted batch job')
            jobid=outspl[1].split()[0]
            for i in range(1,njobs+1):
                fname="subdir_%d.txt"%i
                task_id="%s_%d"%(jobid,i)
                infile=open(fname, 'r')
                calcdirabs=infile.readline()
                uid=infile.readline().strip("\n")
                calcdir=infile.readline().strip("\n")
                infile.close()
                if not (calcdir in self.structureDB[uid].submitted_jobs):
                    self.structureDB[uid].submitted_jobs[calcdir]={'nsubmit':1}
                else:
                    self.structureDB[uid].submitted_jobs[calcdir]['nsubmit']=self.structureDB[uid].submitted_jobs[calcdir]['nsubmit']+1
                self.structureDB[uid].submitted_jobs[calcdir]['jobid']=task_id
                print uid,task_id,calcdirabs,calcdir
            exitcode, out = commands.getstatusoutput("rm njobs")
            os.chdir(parentdir)


    def write_mcif(self,uid, calc_scheme, magsettings={'submitted':True},qaxis=[0,0,1], rotation=[], file="",pd={},initial_magnetic_moments=False):
        if file=="":
            file="%s-%s-%s.mcif"%(self.structureDB[uid].get_composition(reduce=True),uid,calc_scheme)
        if pd=={}:
            pdict=self.get_properties(uid, calc_scheme, magsettings=magsettings)
        else:
            pdict=pd
        #lines=['data_%s'%uid,'_cell_angle_alpha              90','_cell_angle_beta               90','_cell_angle_gamma              90','loop_','_space_group_symop_magn_operation.id','_space_group_symop_magn_operation.xyz','1 x,y,z,+1']
        if True:
            cell=pdict['cell']
            a = norm(cell[0])
            b = norm(cell[1])
            c = norm(cell[2])
            alpha = arccos(dot(cell[1], cell[2])/(b*c))*180./pi
            beta = arccos(dot(cell[0], cell[2])/(a*c))*180./pi
            gamma = arccos(dot(cell[0], cell[1])/(a*b))*180./pi
            lines=['data_%s'%uid,'_cell_angle_alpha              %5.2f'%alpha,'_cell_angle_beta               %5.2f'%beta,'_cell_angle_gamma              %5.2f'%gamma,'loop_','_space_group_symop_magn_operation.id','_space_group_symop_magn_operation.xyz','1 x,y,z,+1','']
            abc=['a','b','c']
            for i in range(3):
                lines.append("_cell_length_%s\t %.3f"%(abc[i],norm(pdict['cell'][i])))
            lines=lines+['loop_','_atom_site_label','_atom_site_type_symbol','_atom_site_fract_x','_atom_site_fract_y','_atom_site_fract_z']
            for i in range(len(pdict['chemical_symbols'])):
                lines.append("%s%d %s %.8f %.8f %.8f"%(pdict['chemical_symbols'][i],i+1,pdict['chemical_symbols'][i],pdict['scaled_positions'][i][0],pdict['scaled_positions'][i][1],pdict['scaled_positions'][i][2]))
            lines=lines+['','loop_','_atom_site_moment.label','_atom_site_moment.crystalaxis_x','_atom_site_moment.crystalaxis_y','_atom_site_moment.crystalaxis_z']
            if initial_magnetic_moments==True:
                magnetic_moments=pdict['initial_magnetic_moments']
            else:
                magnetic_moments=pdict['magnetic_moments']
            for i in range(len(pdict['chemical_symbols'])):
                if isinstance(magnetic_moments[i],float):
                    moms=np.dot(magnetic_moments[i],qaxis/norm(qaxis))
                else:
                  moms=magnetic_moments[i]
                if rotation!=[]:
                    moms=np.dot(rotation,moms)
                lines.append("%s%d %.4f %.4f %.4f"%(pdict['chemical_symbols'][i],i+1,moms[0],moms[1],moms[2]))
            outfile=open(file,"w")
            for line in lines:
                outfile.write("%s\n"%line)
            outfile.close()
        else: #except:
            print "write_mcif: failed %s"%uid

    # def add_m_constr