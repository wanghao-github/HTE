import os
from transport import *
import numpy as np


###############################################################
########### routines for VASP without use of ASE ##############
###############################################################

def setup_vasp_calculation(structure_info, settings, pathname="./"
                            , kpoints_tags=['kpts','gamma']
                            ,potcar_tags=['xc','setups']
                            ,incar_tags=['magmom', 'encut', 'lcharg', 'prec', 'ispin', 'ismear', 'lwave', 'lorbit', 'nsw']):
    lines={}
    success=True
    pos_ok,lines['POSCAR']=write_vasp_poscar(structure_info)
    incar_settings={}
    kpoints_settings={}
    potcar_settings={}
    for tag in settings:
        if tag.lower() in kpoints_tags:
            kpoints_settings[tag.lower()]=settings[tag]
        elif tag.lower() in potcar_tags:
            potcar_settings[tag.lower()]=settings[tag]
        else: #if tag.lower() in incar_tags:
            incar_settings[tag.lower()]=settings[tag]
    inc_ok,lines["INCAR"]=write_vasp_incar(incar_settings)
    kpo_ok,lines["KPOINTS"]=write_vasp_kpoints(kpoints_settings)
    pot_ok,pot_command=write_vasp_potcar(structure_info, potcar_settings, pathname=pathname,pp_paths=[])
    if (pos_ok==True) and (inc_ok==True) and (kpo_ok==True) and (pot_ok):
        try:
            if (not os.path.isdir(pathname)):
                os.makedirs(pathname)
            for fname in lines:
                out=open(os.path.join(pathname,fname),"w")
                for line in lines[fname]:
                    out.write("%s\n"%line)
                out.close()
            os.system(pot_command)
            shutil.copy(os.path.join(pathname,"POSCAR"),os.path.join(pathname,"POSCAR.INI"))
        except:
            success=False
    else:
        success=False
    return success



def write_vasp_poscar(structure_info, pathname=""):
    """write POSCAR file for Vasp calculation
    structure_info: dictionary with structure information, must contain:
        'cell': array with lattice vectors
        'chemical_symbols': array with atom names
        'scaled_positions': array with atomic positions (relativ coordinates)
        or alternatively:
        'positions': array with atomic positions (absolute coordinates, only used if 'scaled_positions' not given)
    """
    success=True
    linesposcar=[]
    try:
        symbsposcar=[]
        multposcar=[]
        if 'chemical_symbols' in structure_info:
            for el in structure_info['chemical_symbols']:
                if (not (el in symbsposcar)) or (el!=symbsposcar[-1]):
                    symbsposcar.append(el)
                    multposcar.append(1)
                else:
                    multposcar[-1]=multposcar[-1]+1
                #todo: check resort if atoms are not ordered
            lsym=""
            lmult=""
            for i in range(len(symbsposcar)):
                lsym=lsym+"%s "%symbsposcar[i]
                lmult=lmult+" %d"%multposcar[i]
        else:
            success=False
        if 'cell' in structure_info:
            linesposcar.append(lsym)
            linesposcar.append("%19.14f"%1.0)
            for i in range(3):
                ai=structure_info['cell'][i]
                linesposcar.append(" %21.16f %21.16f %21.16f"%(ai[0],ai[1],ai[2]))
            linesposcar.append(lmult)
        else:
            success=False
        if ('scaled_positions' in structure_info) and (len(structure_info['scaled_positions'])==len(structure_info['chemical_symbols'])):
            linesposcar.append("Direct")
            for pos in structure_info['scaled_positions']:
                linesposcar.append(" %19.16f %19.16f %19.16f"%(pos[0],pos[1],pos[2]))
        elif ('positions' in structure_info) and (len(structure_info['positions'])==len(structure_info['chemical_symbols'])):
            linesposcar.append("Cartesian")
            for pos in structure_info['positions']:
                linesposcar.append(" %19.16f %19.16f %19.16f"%(pos[0],pos[1],pos[2]))
        else:
            success=False
        if (success==True) and (pathname!=""):
            #write POSCAR
            if (not os.path.isdir(pathname)):
                os.makedirs(pathname)
            out=open(os.path.join(pathname,"POSCAR"),"w")
            for line in linesposcar:
                out.write("%s\n"%line)
            out.close()
            shutil.copy(os.path.join(pathname,"POSCAR"),os.path.join(pathname,"POSCAR.INI"))
    except:
        success=False
    return success,linesposcar


def write_vasp_incar(settings, pathname=""):
    """write INCAR file for Vasp calculation
    settings: dictionary with calculation settings
        'cell': array with lattice vectors
        'chemical_symbols': array with atom names
        'scaled_positions': array with atomic positions (relativ coordinates)
        or alternatively:
        'positions': array with atomic positions (absolute coordinates, only used if 'scaled_positions' not given)
    """
    success=True
    linesincar=["INCAR created by HTE"]
    if True:
        for tag in sorted(settings):
            tup=tag.upper()
            if tup=="MAGMOM":
                line=" MAGMOM ="
                mcoll=[0,""]
                for mom in settings[tag]:
                    if (isinstance(mom,list)) or (isinstance(mom,np.ndarray)):
                        moms=mom
                    else:
                        moms=[mom]
                    for m in moms:
                        #print m,mcoll
                        if abs(m)<0.01:
                            mstr="0"
                        else:
                            mstr="%.2f"%m
                        if mstr!=mcoll[1]:
                            if mcoll[0]>1:
                                line=line+" %d*%s"%(mcoll[0],mcoll[1])
                            elif mcoll[0]==1:
                                line=line+" %s"%mcoll[1]
                            mcoll=[1,mstr]
                        else:
                            mcoll[0]=mcoll[0]+1
                        #print line,mstr,mcoll
                if mcoll[0]>1:
                    line=line+" %d*%s"%(mcoll[0],mcoll[1])
                elif mcoll[0]==1:
                    line=line+" %s"%mcoll[1]
                    #for m in moms:
                    #    if abs(m)<0.01:
                    #        line=line+" 0"
                    #    else:
                    #        line=line+" %.2f"%m
                linesincar.append(line)
                
            if tup=="M_CONSTR":
                line=" M_CONSTR ="
                mcoll=[0,""]
                for mom in settings[tag]:
                    if (isinstance(mom,list)) or (isinstance(mom,np.ndarray)):
                        moms=mom
                    else:
                        moms=[mom]
                    for m in moms:
                        #print m,mcoll
                        if abs(m)<0.01:
                            mstr="0"
                        else:
                            mstr="%.2f"%m
                        if mstr!=mcoll[1]:
                            if mcoll[0]>1:
                                line=line+" %d*%s"%(mcoll[0],mcoll[1])
                            elif mcoll[0]==1:
                                line=line+" %s"%mcoll[1]
                            mcoll=[1,mstr]
                        else:
                            mcoll[0]=mcoll[0]+1
                        #print line,mstr,mcoll
                if mcoll[0]>1:
                    line=line+" %d*%s"%(mcoll[0],mcoll[1])
                elif mcoll[0]==1:
                    line=line+" %s"%mcoll[1]
                    #for m in moms:
                    #    if abs(m)<0.01:
                    #        line=line+" 0"
                    #    else:
                    #        line=line+" %.2f"%m
                linesincar.append(line)
                
            elif isinstance(settings[tag],list):
                line=" %s ="%tup
                for val in settings[tag]:
                    if isinstance(val,int):
                        line=line+" %d"%val
                    elif isinstance(val,float):
                        line=line+" %.2f"%val
                    else:
                        line=line+" %s"%str(val)
                linesincar.append(line)
            elif isinstance(settings[tag],bool):
                if settings[tag]==True:
                    linesincar.append(" %s = .TRUE."%(tag.upper()))
                else:
                    linesincar.append(" %s = .FALSE."%(tag.upper()))
            elif isinstance(settings[tag],float):
                linesincar.append(" %s = %.3f"%(tag.upper(),settings[tag]))
            elif isinstance(settings[tag],int):
                linesincar.append(" %s = %d"%(tag.upper(),settings[tag]))
            elif isinstance(settings[tag],str):
                linesincar.append(" %s = %s"%(tag.upper(),settings[tag]))
            #else:
            #    success=False
        if (success==True) and (pathname!=""):
            #write INCAR
            if (not os.path.isdir(pathname)):
                os.makedirs(pathname)
            out=open(os.path.join(pathname,"INCAR"),"w")
            for line in linesincar:
                out.write("%s\n"%line)
            out.close()
    else:
        success=False
    return success,linesincar
        
def vasp_check_pp_files(element, settings, pppaths):
        #TODO: make this more transparent, no need for old ASE conventions
        if 'xc' in settings:
            if settings['xc'] == 'PW91':
                xc = '_gga/'
            elif settings['xc'] == 'PBE':
                xc = '_pbe/'
            else:
                success=False
        if 'setups' in settings:
            setups=settings['setups']
        else:
            setups={}
        pp=element
        if element in setups:
            pp=element+setups[element]
        for path in pppaths:
            ppfile=os.path.join(path,'potpaw'+xc.upper()+str(pp),'POTCAR')
            if isfile(ppfile) or islink(ppfile):
                return True,ppfile
            elif isfile(ppfile+'.Z') or islink(ppfile+'.Z'):
                return True,ppfile+'.Z'
        return False,""

def write_vasp_potcar(structure_info, settings, pathname="",pp_paths=[]):
    """return command to write POTCAR file for Vasp calculation
    """
    if (pp_paths==[]) and ('VASP_PP_PATH' in os.environ):
        pppaths = os.environ['VASP_PP_PATH'].split(':')
    else:
        pppaths=pppaths
    success=True
    elpot=""
    commandline=""
    pipe=" > "
    for el in structure_info['chemical_symbols']:
        if el!=elpot:
            elpot=el
            isok,filename=vasp_check_pp_files(el,settings,pppaths)
            if isok==True:
                if filename.endswith(".Z"):
                    command="gunzip -c %s"%filename
                else:
                    command="cat %s"%filename
                if commandline=="":
                    commandline="%s > %s/POTCAR "%(command,pathname)
                else:
                    commandline=commandline+"; %s >> %s/POTCAR "%(command,pathname)
            else:
                success=False
    return success,commandline


def write_vasp_kpoints(settings, pathname=""):
    """write INCAR file for Vasp calculation
    settings: dictionary with calculation settings
        'cell': array with lattice vectors
        'chemical_symbols': array with atom names
        'scaled_positions': array with atomic positions (relativ coordinates)
        or alternatively:
        'positions': array with atomic positions (absolute coordinates, only used if 'scaled_positions' not given)
    """
    success=True
    lines=["KPOINTS created by HTE","0"]
    if True:
        lines=["KPOINTS created by HTE","0"]
        if ('gamma' in settings) and (settings['gamma']==True):
            lines.append("Gamma")
        else:
            lines.append("Monkhorst-Pack")
        if 'kpts' in settings:
            N=settings['kpts']
            lines.append("%d %d %d"%(N[0],N[1],N[2]))
            lines.append("0 0 0")
        else:
            success=False
        if (success==True) and (pathname!=""):
            #write KPOINTS
            if (not os.path.isdir(pathname)):
                os.makedirs(pathname)
            out=open(os.path.join(pathname,"KPOINTS"),"w")
            for line in lines:
                out.write("%s\n"%line)
            out.close()
    else:
        success=False
    return success,lines

        
###############################################################
########### routines for VASP (may use ASE)      ##############
###############################################################

def get_vasp_number_of_valence_electrons(outcar='OUTCAR'):
    """returns the total number of valence electrons as given
    in the OUTCAR file of a vasp pseudopotential calculation"""
    num_eval=None
    if (os.path.isfile(outcar)):
        neval=0
        zval=[]
        mult=[]
        fin=open(outcar,"r")
        line=fin.readline() #??/??/??/nspin
        while line:
            if ("POMASS" in line) and ("ZVAL   =" in line):
                zval.append(line.split("ZVAL   =")[1].split()[0])
            if ("ions per type =" in line):
                mult=line.split("ions per type =")[1].split()
            line=fin.readline()
        fin.close()
        if (zval!=[]) and (len(zval)==len(mult)):
            num_eval=0
            for i in range(len(zval)):
                num_eval=num_eval+float(zval[i])*float(mult[i])
    return num_eval
                        
def get_vasp_bandgap(pathname='./',outcar='OUTCAR',eigenval='EIGENVAL',pathname_sc=None,methods=['nval','E_fermi']):
    """returns the band gap of a vasp pseudopotential calculation"""
    gap=None
    if pathname_sc==None:
        pathname_sc=pathname
    bandstructure=get_vasp_bandstructure(pathname=pathname,filename_eigenval=eigenval,filename_outcar_sc=outcar,pathname_sc=pathname_sc)
    nval=get_vasp_number_of_valence_electrons(outcar=os.path.join(pathname_sc,outcar))
    if (bandstructure==None):
        return gap
    for method in methods:
        if (method=='E_fermi'):
            #determine the valence band maximum and conduction band minimum
            e_cbmin=None
            e_vbmax=None
            efermi=bandstructure['E_Fermi']
            eb=bandstructure['e_n_k']
            for ikp in range(len(eb)):
                enk=eb[ikp]
                if (bandstructure['nspin']==1):
                    spindirs=['energies_up']
                else:
                    spindirs=['energies_dn']
                for ispin in spindirs:
                    ex=enk[ispin]
                    for iband in range(len(ex)):
                        e=ex[iband]
                        #valence band maximum
                        if (e<efermi) and ((e_vbmax==None) or (e_vbmax<e)):
                            e_vbmax=e
                        if (e>efermi) and ((e_cbmin==None) or (e_cbmin>e)):
                            e_cbmin=e                     
            gap=e_cbmin-e_vbmax
            break
        elif (method=='nval') and (bandstructure['nspin']==1) and (nval!=None):
            e_cbmin=None
            e_vbmax=None
            ib=int(round(0.5*nval))
            if (abs(nval-2.0*ib)>0.001):
                return 0.0
            eb=bandstructure['e_n_k']
            for ikp in range(len(eb)):
                enk=eb[ikp]
                ex=enk['energies_up']
                if (e_vbmax==None) or (ex[ib-1]>e_vbmax):
                    e_vbmax=ex[ib-1]
                if (e_cbmin==None) or (ex[ib]<e_cbmin):
                    e_cbmin=ex[ib]
            gap=e_cbmin-e_vbmax
            if (gap<0.0):
                gap=0.0
            break
    return gap

def get_properties_vasp(calcdir='./',outcar='OUTCAR', contcar='CONTCAR', incar='INCAR', convert_strings=False, resort=True, eps=1e-5):
    """Evaluate vasp calculation and return properties in dictionary (as string values)"""
    prop_dict={}
    magmoms={}
    ispin=1
    errors=False
    outcar_file=os.path.join(calcdir,outcar)
    contcar_file=os.path.join(calcdir,contcar)
    incar_file=os.path.join(calcdir,incar)
    #use CONTCAR to get chemical species and structure
    if os.path.isfile(contcar_file):
        try:
            chem_symb=[]
            infile=open(contcar_file,'r')
            line=infile.readline() #comment line
            scale=infile.readline().strip() #possible scaling of lattice vectors
            cell=[]
            for i in range(3):
                cell.append(infile.readline().split())
            fscale=float(scale)
            scell=[]
            fcell=np.zeros((3,3))
            for i in range(3):
                slatt=[]
                for j in range(3):
                    fcell[i][j]=fscale*float(cell[i][j])
                    slatt.append(str(fcell[i][j]))
                scell.append(slatt)
            if abs(fscale-1.0)<eps:
                prop_dict['cell']=cell
            else:
                prop_dict['cell']=scell
                print "***",scell,cell,fscale
            if convert_strings==True:
                prop_dict['cell']=fcell
            #patch for CONTCARS where atom types/multiplicities extend over several lines
            atom_types=infile.readline().split()
            lsplit=infile.readline().split()
            while (lsplit[0].isdigit()==False):
                atom_types=atom_types+lsplit
                lsplit=infile.readline().split()
            atom_mult=lsplit
            lsplit=infile.readline().split()
            while (lsplit[0].isdigit()==True):
                atom_mult=atom_mult+lsplit
                lsplit=infile.readline().split()
            for i in range(len(atom_mult)):
                for j in range(int(atom_mult[i])):
                    chem_symb.append(atom_types[i])
            prop_dict['chemical_symbols']=chem_symb    
            atpos_type=lsplit
            #atom_types=infile.readline().split()
            #atom_mult=infile.readline().split()
            #for i in range(len(atom_mult)):
            #    for j in range(int(atom_mult[i])):
            #        chem_symb.append(atom_types[i])
            #prop_dict['chemical_symbols']=chem_symb    
            #atpos_type=infile.readline().strip()
            if atpos_type=='Direct':
                scaled_positions=[]
                for i in range(len(chem_symb)):
                    scpos=infile.readline().split()
                    fscpos=[]
                    for j in range(len(scpos)):
                        fscpos.append(float(scpos[j]))
                    if convert_strings==True:
                        scaled_positions.append(fscpos)
                    else:
                        scaled_positions.append(scpos)
                prop_dict['scaled_positions']=scaled_positions
                print "check_point107, print scaled_positions",  scaled_positions               
            infile.close()
        except:
            errors=True
    #use INCAR to get initial magnetization
    if os.path.isfile(incar_file):
        noncol=False
        magmoms_ini=[]
        try:
            infile=open(incar_file,'r')
            line=infile.readline()
            while line:
                if 'MAGMOM' in line:
                    moms=line.split('=')[1].split()
                    magmoms_ini=[]
                    for mom in moms:
                        if '*' in mom:
                            for i in range(int(mom.split('*')[0])):
                                magmoms_ini.append(mom.split('*')[1])
                        else:
                            magmoms_ini.append(mom)
                if ('LNONCOLLINEAR' in line) and ('TRUE' in line.split('=')[1]):
                    noncol=True
                line=infile.readline()
            infile.close()
            if magmoms_ini!=[]:
                if (noncol==False) and (len(magmoms_ini)==len(prop_dict['chemical_symbols'])):
                    prop_dict['initial_magnetic_moments']=magmoms_ini
                elif (noncol==True) and (len(magmoms_ini)==3*len(prop_dict['chemical_symbols'])):
                    mvec=[]
                    mom3=[]
                    for mom in magmoms_ini:
                        mom3.append(mom)
                        if len(mom3)==3:
                            mvec.append(mom3)
                            mom3=[]
                    prop_dict['initial_magnetic_moments']=mvec
                else:
                    errors=True
        except:
            errors=True
    #use OUTCAR to get energy etc.
    if os.path.isfile(outcar_file):
        infile=open(outcar_file,'r')
        line=infile.readline()
        while line:
            if (line.startswith('   ISPIN  =')):
                try:
                    ispin=int(line.split("=")[1].split()[0])
                except:
                    print "ISPIN"
                    errors=True
            if (line.startswith('  energy  without entropy')):
                val=line.split()[-1]
                try:
                    fval=float(val)
                    if convert_strings==True:
                        prop_dict['energy']=fval
                    else:
                        prop_dict['energy']=val
                except:
                    print "energy"
                    errors=True
            if (line.startswith('  enthalpy is  TOTEN    =')):
                val=line.split()[4]
                try:
                    fval=float(val)
                    if convert_strings==True:
                        prop_dict['enthalpy']=fval
                    else:
                        prop_dict['enthalpy']=val
                except:
                    errors=True
            if (line.startswith(' number of electron')):
                try:
                    val=line.split()[-1]
                    fval=float(val)
                    if noncol==True:
                        val=line.split('magnetization')[1].split()
                        fval=norm(np.array((float(val[0]),float(val[1]),float(val[2]))))
                        val="%.4f"%fval
                    if convert_strings==True:
                        prop_dict['magnetic_moment']=fval
                    else:
                        prop_dict['magnetic_moment']=val
                except:
                    if ispin==2:
                        print 'magnetic_moment'
                        errors=True
            if (line.startswith(' magnetization (')):
                dir=line.split('(')[-1][0]
                if dir in ['x','y','z']:
                    magmoms[dir]=[]
                    while line and not (line.startswith('-----------------')):
                        line=infile.readline()
                    line=infile.readline()
                    while (line) and (len(line.split())>0) and (not (line.startswith('-----------------'))):
                        val=line.split()[-1]
                        try:
                            fval=float(val)
                            if convert_strings==True:
                                magmoms[dir].append(fval)
                            else:
                                magmoms[dir].append(val)
                        except:
                            errors=True
                        line=infile.readline()
                    prop_dict['magnetic_moments']=magmoms
                    print "check_point177,prop_dict['magnetic_moments'] is", prop_dict['magnetic_moments']
            if (line.startswith(' TOTAL ELASTIC MODULI (kBar)')):
                #adopted from D. Ohmer
                try:
                    line=infile.readline()
                    line=infile.readline()
                    EMtens=[]
                    for i in range(6):
                        line=infile.readline()
                        lsplit=line.split()
                        EM=[]
                        for j in range(1,7):
                            val=lsplit[j]
                            fval=float(val)
                            if convert_strings==True:
                                EM.append(fval)
                            else:
                                EM.append(val)
                        EMtens.append(EM)
                    prop_dict['total_elastic_moduli']=EMtens
                except:
                    errors=True
            line=infile.readline()
        infile.close()
        if 'x' in magmoms:
            dirs=sorted(magmoms)
            if len(dirs)==1:
                prop_dict['magnetic_moments']=magmoms['x']
            else:
                prop_dict['magnetic_moments']=[]
                for i in range(len(magmoms['x'])):
                    mm=[]
                    for dir in dirs:
                        mm.append(magmoms[dir][i])
                    prop_dict['magnetic_moments'].append(mm)
    prop_dict['errors']=str(errors)
    #restore intitial structure if it was resorted by ase
    sort_file=os.path.join(calcdir,"ase-sort.dat")
    if (resort==True) and (os.path.isfile(sort_file)):
        try:
            sortini=[]
            backsort=[]
            infile=open(sort_file,'r')
            for i in range(len(prop_dict['chemical_symbols'])):
                line=infile.readline()
                sortini.append(int(line.split()[0]))
                backsort.append(int(line.split()[1]))
            infile.close()
            for arg in ['chemical_symbols','scaled_positions','initial_magnetic_moments','magnetic_moments']:
                pdsorted=[]
                for i in range(len(backsort)):
                    pdsorted.append(prop_dict[arg][backsort[i]])
                prop_dict[arg]=pdsorted
            prop_dict['sort']=sortini
        except:
            errors=True
            
    return prop_dict
