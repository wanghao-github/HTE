# simple version of an fplo calculator

from ase.calculators.general import Calculator
import os
import commands
import shutil
from ase import Atoms, Atom
from hte_dbentry import *


def get_properties_fplo(calcdir='./', outfile='out'):
    """Evaluate FPLO calculation and return properties in dictionary (as string values)
       Units are in eV and Angstroem
    """
    prop_dict={}
    out=os.path.join(calcdir, outfile)
    finished=False
    if os.path.isfile(out):
        exitcode, mes=commands.getstatusoutput("tail -n1 %s"%out)
        if 'TERMINATION: Finished' in mes:
            finished=True
            if 'SCF calculation' in mes:
                scf=True
            elif 'single step calculation' in mes:
                scf=False
            else:
                prop_dict={'errors':"Check type of calculation"}
                return prop_dict
    if finished==False:
        return prop_dict
    au2AA=get_conversion('a.u.','AA')
    #Evaluate outfile
    infile=open(out,'r')
    line=infile.readline()
    while line:
        if line.startswith('lattice vectors'):
            cell=[]
            try:
                for i in range(3):
                    line=infile.readline().split(':')[1].split()
                    cell.append([str(float(line[0])*au2AA),str(float(line[1])*au2AA),str(float(line[2])*au2AA)])
            except:
                prop_dict={'errors':"Check lattice vectors"}
                return prop_dict
            prop_dict['cell']=cell
        elif line.startswith('Number of sites :'):
            nsite=int(line.split(':')[1])
        elif line.startswith('No.  Element WPS CPA-Block    X                      Y                      Z'):
            chem_symb=[]
            atpos=[]
            try:
                for i in range(nsite):
                    line=infile.readline().split()
                    chem_symb.append(line[1])
                    atpos.append([str(float(line[4])*au2AA),str(float(line[5])*au2AA),str(float(line[6])*au2AA)])
            except:
                prop_dict={'errors':"Check atomic positions"}
                return prop_dict
            prop_dict['chemical_symbols']=chem_symb
            prop_dict['positions']=atpos
        elif ("MAG.MOMENT" in line):
            magmoms=[]
            try:
                line=infile.readline()
                for i in range(nsite):
                    line=infile.readline().split('|')
                    fmom=float(line[3])
                    magmoms.append(line[3])
            except:
                prop_dict={'errors':"Check magnetic moments"}
                return prop_dict
            prop_dict['magnetic_moments']=magmoms
        elif  (scf==True) and ("EE" in line):
            try:
                energy=float(line.split()[1])*get_conversion('Hartree','eV')
            except:
                prop_dict={'errors':"Check total energy"}
                return prop_dict
            prop_dict['energy']=str(energy)
        elif 'mag. moment' in line:
            line=infile.readline()
            try:
                magmom=line.split()[2]
                fmom=float(magmom)
            except:
                prop_dict={'errors':"Check total magnetic moment"}
                return prop_dict
            prop_dict['magnetic_moment']=magmom
        elif (scf==False) and ('Band energy:' in line):
            try:
                band_energy=float(line.split(":")[1].split(',')[0])*get_conversion('Hartree','eV')
            except:
                prop_dict={'errors':"Check band energy"}
                return prop_dict
            prop_dict['band_energy']=str(band_energy)
        elif 'BANDWT: total     gap' in line:
            try:
                gap=line.split()[4]
                fgap=float(gap)
            except:
                prop_dict={'errors':"Check band gap"}
                return prop_dict
            prop_dict['band_gap']=gap
        line=infile.readline()
    infile.close()
    return prop_dict    

def FPLO_densfile_is_spinpolarized(self, calcdir='./'):
    densfile=os.path.join(calcdir,"=.dens")
    if os.path.isfile(densfile):
        exitcode, mes=commands.getstatusoutput('grep "mspin=" %s | tail -n1'%densfile)
        if '2;' in mes:
            return True
    return False

class FPLO(Calculator):
    def __init__(self, **kwargs):
        self.name='FPLO'
        self.fedit=None
        if ('job_settings' in kwargs):
            job_settings=kwargs['job_settings']
            if 'fedit' in job_settings:
                self.fedit=job_settings['fedit']
        if 'fedit' in kwargs:
            self.fedit=kwargs['fedit']
            del kwargs['fedit']
        #else:
        #    self.fedit='fedit9.01-35-i386'
        self.settings=kwargs


    def initialize(self, atoms_object, **kwargs):
        # set up FPLO calculation for atoms_object
        for x in self.settings:
            if not (x in kwargs):
                kwargs[x]=self.settings[x]
        use_str=True
        if 'update_symmetry' in kwargs:
            update_symmetry=kwargs['update_symmetry']
            del kwargs['update_symmetry']
        else:
            update_symmetry=True
        if 'use_str' in kwargs:
            use_str=kwargs['use_str']
            del kwargs['use_str']
        if 'calcdir' in kwargs:
            calc_dir=kwargs['calcdir']
        else:
            calc_dir='./'
        if 'TOL' in kwargs:
            TOL=kwargs['TOL']
            del kwargs['TOL']
        else:
            TOL=1.e-6
        if use_str:
            update_symmetry=False
            if 'magmom' in kwargs:
                status,kw=self.update_symmetry_str(atoms_object, magmoms=kwargs['magmom'], calcdir=calc_dir,TOL=TOL)
                del kwargs['magmom']
                for x in kw:
                    kwargs[x]=kw[x]
            else:
                status=self.export_str_file(atoms_object, calcdir=calc_dir,TOL=TOL)
                if status==True:
                    #check this
                    exitcode, out = commands.getstatusoutput('timeout 60s str')
                if exitcode!=0:
                    status=False
            if status==True:
                exitcode, out = commands.getstatusoutput('%s -pipe < +pipe_wyck > xxx'%self.fedit)
                if 'kspace_density' in kwargs:
                    symminfo=self.read_symminfo() #is not up to date at this time
                    N=[12,12,12]
                    if 'reciprocal_lattice' in symminfo:
                        rc=symminfo['reciprocal_lattice']/0.5291772083
                    else:
                        rc=atoms_object.get_reciprocal_cell() #TODO: may not match actual BZ
                    Nx=1
                    for i in range(3):
                        x2=0.0
                        for j in range(3):
                            x2=x2+rc[i][j]*rc[i][j]
                        N[i]=int(ceil(kwargs['kspace_density']*sqrt(x2)))
                        Nx=max(N[i],Nx)
                    kwargs['k-mesh']=[Nx,Nx,Nx] #TODO
            else:
                return False
        if self.export_fplo_pipefile(atoms_object, update_symmetry=update_symmetry, **kwargs):
            exitcode, out = commands.getstatusoutput('%s -pipe < +pipe > xxx'%self.fedit)
            return True
        # set up failed for some reason
        print 'WARNING: FPLO.initialize(): setup failed'
        return False
            
    def export_fplo_pipefile(self, atoms_object, calcdir='./', pipefile='+pipe', bs_calc=False, xc_vers=5, update_symmetry=True,k_mesh=[12,12,12],ispin=1, **kwargs):
        # create a pipe file with structural information for FPLO
        # returns False upon failure, True upon success
        #
        # check if atoms_object.info has sufficient information: TODO: remove this, no longer needed!
        if update_symmetry:
            info_ok=HTEdbentry(atoms_obj=atoms_object).update_symmetry_information()
        else:
            info_ok=True
        if info_ok==False:
            print 'export_fplo_pipefile(): insufficient structure information, nothing done!'
            return False
        if os.path.isdir(calcdir):
            outfile=os.path.join(calcdir,pipefile)
            fout=open(outfile,'w')
            fout.write("### pipe file created by HTE\n")
            if update_symmetry==True:
                nwyckoff=len(atoms_object.info['_atom_site_type_symbol'])
                # enter symmetry menu:
                fout.write("# enter symmetry menu:\n@+@\n")
                # set unit to Angstroem:
                fout.write("# set unit to Angstroem:\n@u@\n@a@\n@x@\n")
                # set spacegroup
                fout.write("# set spacegroup\n@s@\n@%d@\n@x@\n"%atoms_object.info['spacegroup'].no)
                # set lattice parameters
                fout.write("# set lattice parameters\n")
                fout.write("@l@%8.4f %8.4f %8.4f\n"%(atoms_object.info['_cell_length_a'], atoms_object.info['_cell_length_b'], atoms_object.info['_cell_length_c']))
                # set angles
                fout.write("# set angles\n")
                fout.write("@a@%8.4f %8.4f %8.4f\n"%(atoms_object.info['_cell_angle_alpha'], atoms_object.info['_cell_angle_beta'], atoms_object.info['_cell_angle_gamma']))
                # set number of atoms
                fout.write("# set number of atoms\n")
                fout.write("@n@%d\n"%nwyckoff)
                # set Wyckoff positions
                fout.write("# set Wyckoff positions\n")
                for i in range(nwyckoff):
                    atom_symbol=''
                    atom_site_symbol=atoms_object.info['_atom_site_type_symbol'][i].strip()
                    for j in range(len(atom_site_symbol)+1):               
                        if atom_site_symbol[0:j].isalpha():
                            atom_symbol=atom_site_symbol[0:j]
                    fout.write("@%d@%s@%12.8f %12.8f %12.8f\n"%(i+1,atom_symbol,atoms_object.info['_atom_site_fract_x'][i],atoms_object.info['_atom_site_fract_y'][i],atoms_object.info['_atom_site_fract_z'][i]))
                # update and leave symmetry menu
                fout.write("# update and leave symmetry menu\n")
                fout.write("@+@\n")
                fout.write("@x@\n")
            # set Vxc
            if xc_vers!=None:
                fout.write("# set Vxc\n")
                fout.write("@v@\n@%d@\n"%xc_vers)        
                fout.write("@x@\n")
            if bs_calc==True:
                # calculate bandstructure and IDOS
                fout.write("# calculate bandstructure\n")
                fout.write("@ b@\n@b@+\n@p@+\n@x@\n")
            # spinpolarized?
            fout.write("# spinpolarization\n")
            fout.write("@s@%d\n"%ispin)
            if ('initial_spin_split' in kwargs) and (kwargs['initial_spin_split']) in ['t','f']:
                fout.write("@i@%s\n"%kwargs['initial_spin_split'])
            if ('initial_spin_split_sorts' in kwargs):
                fout.write("@ i@\n")
                for i in range(len(kwargs['initial_spin_split_sorts'])):
                    fout.write("@%d@%.1f\n"%(i+1,kwargs['initial_spin_split_sorts'][i]))
                fout.write("@x@\n")
            # relativistic
            if ('relativistic' in kwargs):
                fout.write("# relativistic settings\n")
                relat=kwargs['relativistic'][0]
                if relat.lower() in ['n','s','k','f']:
                    fout.write("@r@\n@%s@\n@x@\n"%relat)     
            if ('quantization_axis' in kwargs):
                fout.write("# quantization_axis\n")
                qax=kwargs['quantization_axis']
                fout.write("@u@%d %d %d\n"%(qax[0],qax[1],qax[2]))
            if ('niter' in kwargs):
                fout.write("@n@%d\n"%kwargs['niter'])
            if ('iterat_version' in kwargs):
                fout.write("@ t@\n@v@\n@%d@\n@x@\n@x@\n"%kwargs['iterat_version'])
            if ('occupied_bands' in kwargs):
                fout.write("@o@%d\n"%kwargs['occupied_bands'])
            # set k-mesh
            if ('kmesh' in kwargs): #remove elif part (check if used!)
                fout.write("# set k-mesh\n")
                k_mesh=kwargs['kmesh']
                fout.write("@k@%d %d %d\n"%(k_mesh[0],k_mesh[1],k_mesh[2]))        
            elif ('k-mesh' in kwargs):
                fout.write("# set k-mesh\n")
                k_mesh=kwargs['k-mesh']
                fout.write("@k@%d %d %d\n"%(k_mesh[0],k_mesh[1],k_mesh[2]))        
            # leave fedit
            fout.write("# leave fedit\n")
            fout.write("@q@\n")
            fout.write("\n")
            fout.close()
            return True
        print 'export_fplo_pipefile(): ',calcdir, ' not found, nothing done!'
        return False
        
    def read_convergence(self, outfile='out'):
        """Returns True if FPLO calculation terminated with status 'TERMINATION: Finished',
        otherwise False
        """
        if os.path.isfile(outfile):
            exitcode, mes=commands.getstatusoutput("tail -n1 %s"%outfile)
            if 'TERMINATION: Finished' in mes:
                return True
        return False

    def is_spinpolarized_densfile(self, calcdir='./'):
        densfile=os.path.join(calcdir,"=.dens")
        if os.path.isfile(densfile):
            exitcode, mes=commands.getstatusoutput('grep "mspin=" %s | tail -n1'%densfile)
            if '2;' in mes:
                return True
        return False

    def read_convergence_level(self, calcdir='./', outfile='out'):
        grit=None
        try:
            outf=os.path.join(calcdir,outfile)
            if os.path.isfile(outf):
                exitcode, mes=commands.getstatusoutput('grep "last deviation" %s | tail -n1'%outf)
                if 'last deviation u=' in mes:
                    grit=float(mes.split("u=")[1].split()[0])
                elif 'last deviation=' in mes:
                    grit=float(mes.split("deviation=")[1].split()[0])
        except:
            grit=None
        return grit

    def estimate_occupied_bands(self, calcdir='./', outfile='out', offset=[10,1.1]):
        loccu=-1
        outf=os.path.join(calcdir,outfile)
        if os.path.isfile(outf):
            try:
                exitcode, mes=commands.getstatusoutput('grep "BANDWT: lowest fully unoccupied band:" %s | tail -n5'%outf)
                for line in mes.split("\n"):
                    if "BANDWT: lowest fully unoccupied band:" in line:
                        lo=int(line.split("lowest fully unoccupied band:")[1])
                        if lo>loccu:
                            loccu=lo
            except:
                loccu=-1
        if loccu==-1:
            return loccu
        return max(loccu+offset[0],int(loccu*offset[1]))
    
    
    def read_energy(self, outfile='out', converged_only=True):
        energy=None
        if (converged_only) and (self.read_convergence(outfile=outfile)==False):
            return None
        if os.path.isfile(outfile):
            exitcode, mes=commands.getstatusoutput("grep EE %s | tail -n1"%outfile)
            if 'EE:' in mes:
                energy=float(mes.split()[1])
        return energy


    def read_bandgap(self, outfile='out',method='outfile'):
        # estimate band gap (for the moment only info from outfile)
        gap=None
        if self.read_convergence(outfile=outfile):
            exitcode, mes=commands.getstatusoutput("grep \"BANDWT: total     gap\" %s | tail -n1"%outfile)
            if 'gap' in mes:
                gap=float(mes.split()[4])
        return gap

    def read_DOS_at_E_F(self, outfile='out'):
        """read DOS at Fermi level from FPLO outfile
        """
        DOS_at_E_F=None
        if self.read_convergence(outfile=outfile):
            exitcode, mes=commands.getstatusoutput("grep \"BANDWT: Density of states at E_f:   N(E_f) =\" %s | tail -n1"%outfile)
            if 'Density of states at E_f:   N(E_f)' in mes:
                DOS_at_E_F=float(mes.split("=")[1].split()[0])
        return DOS_at_E_F

    def read_band_energy(self, outfile='out'):
        """read band energy from FPLO outfile
        status: experimental, may need stricter convergency settings and a large number of k-points to yield a reasonable number
        """
        eband=None
        if self.read_convergence(outfile=outfile):
            exitcode, mes=commands.getstatusoutput("grep \"Band energy: \" %s | tail -n1"%outfile)
            if 'Band energy:' in mes:
                eband=float(mes.split(":")[1].split(',')[0])
        return eband

    def read_input_file(self, calcdir='./'):
        """returns a dictionary with (some) settings of FPLO =.in file
        """
        sym_info={}
        infile_name=os.path.join(calcdir,"=.in")
        if os.path.isfile(infile_name):
            infile=open(infile_name,'r')
            line=infile.readline()
            while line:
                if 'struct {char[2] element;real tau[3];} wyckoff_positions[nsort]' in line:
                    sym_info['wyckoff_elements']=[]
                    while (line) and (not ("};" in line)):
                        if '{"' in line:
                            el=line.split('"')[1]
                            sym_info['wyckoff_elements'].append(el)
                        line=infile.readline()
                if 'struct {int number;char[*] symbol;} spacegroup' in line:
                    if not('={' in line):
                        line=line+infile.readline()
                    sym_info['spacegroup_number']=int(line.split('={')[1].split(',')[0])
                    sym_info['spacegroup_symbol']=line.split(',"')[1].split('"')[0]
                elif 'lengthunit' in line:
                    if not('={' in line):
                        line=line+infile.readline()
                    sym_info['lengthunit']=line.split('"')[1]
                elif 'lattice_constants' in line:
                    if not('={' in line):
                        line=line+infile.readline()
                    latstr=line.split('={')[1]
                    if "};" in latstr:
                        latstr=latstr.split("};")[0]
                    sym_info['lattice_constants']=latstr.split(',')
                line=infile.readline()
        return sym_info
    
    def read_symminfo(self, calcdir='./',file='+symmetry'):
        """reads (some) information from +symminfo created by FPLO
        """
        sym_info={}
        infile_name=os.path.join(calcdir,file)
        if os.path.isfile(infile_name):
            infile=open(infile_name,'r')
            line=infile.readline()
            while line:
                if 'reciprocial lattice vectors' in line:
                    rec_latt=np.eye(3)
                    for i in range(3):
                        g=infile.readline().split(':')[1].split()
                        for j in range(3):
                            rec_latt[i][j]=float(g[j])
                    sym_info['reciprocal_lattice']=rec_latt
                line=infile.readline()
        return sym_info
        
    def get_symmetry_settings(self, calcdir='./'):
        """read symmetry information from =.in
        returns a dictionary (obsolete, use read_input_file instead)
        """
        sym_info={}
        infile_name=os.path.join(calcdir,"=.in")
        if os.path.isfile(infile_name):
            infile=open(infile_name,'r')
            line=infile.readline()
            while line:
                if 'struct {int number;char[*] symbol;} spacegroup' in line:
                    if not('={' in line):
                        line=line+infile.readline()
                    sym_info['spacegroup_number']=int(line.split('={')[1].split(',')[0])
                    sym_info['spacegroup_symbol']=line.split(',"')[1].split('"')[0]
                elif 'lengthunit' in line:
                    if not('={' in line):
                        line=line+infile.readline()
                    sym_info['lengthunit']=line.split('"')[1]
                elif 'lattice_constants' in line:
                    if not('={' in line):
                        line=line+infile.readline()
                    latstr=line.split('={')[1]
                    if "};" in latstr:
                        latstr=latstr.split("};")[0]
                    sym_info['lattice_constants']=latstr.split(',')
                line=infile.readline()
        
    def get_spacegroup(self, calcdir='./'):
        """read the spacegroup from =.in
        returns a tuple (number,symbol)
        """
        number,symbol=None,None
        infile_name=os.path.join(calcdir,"=.in")
        if os.path.isfile(infile_name):
            infile=open(infile_name,'r')
            line=infile.readline()
            while line:
                if 'struct {int number;char[*] symbol;} spacegroup' in line:
                    line=line+infile.readline()
                    number=int(line.split('={')[1].split(',')[0])
                    symbol=line.split(',"')[1].split('"')[0]
                    break
                line=infile.readline()
            infile.close()
        return (number,symbol)


    def read_basis_states(self, calcdir='./', outfile='out', core=False):
        out=os.path.join(calcdir, outfile)
        if os.path.isfile(out):
            infile=open(out,'r')
            line=infile.readline()
            inbasis=False
            nwyck=0
            basis_states=[]
            while line:
                if inbasis:
                    if line.startswith('---------------'):
                        inbasis=False
                        break
                    sline=line.split()
                    if core or (sline[4].startswith('val')):
                        state={'sort':int(sline[0]),'element':sline[1],'nl':int(sline[2]),'orbital':sline[3],'type':sline[4]}
                        #dosfile='+dos.sort%3d.nl%3d'%(state['sort'],state['nl'])
                        #print state,dosfile.replace(' ', '0')
                        basis_states.append(state)
                elif line.startswith('   sort     state'): #commom for fp9/14               vat-class   Q        N       p        px'):
                    line=infile.readline()
                    inbasis=True
                line=infile.readline()
            #print basis_states
            return basis_states
        else:
            print 'Could not open ',out
            return None

        
    def export_str_file(self, atoms_object, calcdir='./',TOL=1.e-6):
        au2ang=0.5291772083
        if os.path.isdir(calcdir):
            outfile=os.path.join(calcdir,'=.str')
            fout=open(outfile,'w')
            fout.write("%s\n 1\n"%str(TOL))
            cell=atoms_object.get_cell()
            for i in range(3):
                fout.write("%30.14f %30.14f %30.14f\n"%(cell[i][0]/au2ang,cell[i][1]/au2ang,cell[i][2]/au2ang))
            atpos=atoms_object.get_positions()
            symb=atoms_object.get_chemical_symbols()
            symb_unique={}
            nat=0
            for el in symb:
                if (not (el in symb_unique)):
                    nat=nat+1
                    symb_unique[el]=nat 
            fout.write("%d\n F\n"%len(symb))
            for i in range(len(symb)):
                fout.write("%s   %d %30.14f %30.14f %30.14f\n"%(symb[i],symb_unique[symb[i]],atpos[i][0]/au2ang,atpos[i][1]/au2ang,atpos[i][2]/au2ang))
            fout.write(" T\n   1     0     0\n   0     1     0\n   0     0     1\n")
            fout.close()
            return True
        return False


    def update_symmetry_str(self, atoms_object, magmoms=[], calcdir='./',TOL=1.e-6):
        au2ang=0.5291772083
        status=True
        kw={}
        if os.path.isdir(calcdir):
            outfile=os.path.join(calcdir,'=.str')
            fout=open(outfile,'w')
            fout.write("%s\n 1\n"%str(TOL))
            cell=atoms_object.get_cell()
            for i in range(3):
                fout.write("%30.14f %30.14f %30.14f\n"%(cell[i][0]/au2ang,cell[i][1]/au2ang,cell[i][2]/au2ang))
            atpos=atoms_object.get_positions()
            symb=atoms_object.get_chemical_symbols()
            dummy_symb=[]
            dummy_tab={}
            dumno=0
            if len(symb)==len(magmoms):
                for el,mom in zip(symb,magmoms):
                    newdum=True
                    for dum_el in dummy_tab:
                        if (el,mom) == dummy_tab[dum_el]:
                            dummy_symb.append(dum_el)
                            newdum=False
                            break
                    if newdum==True:
                        dum_el="X%d"%dumno
                        dummy_tab[dum_el]=(el,mom)
                        dummy_symb.append(dum_el)
                        dumno=dumno+1
            else:
                dummy_symb=symb
            #print symb,dummy_symb,dummy_tab
            symb_unique={}
            nat=0
            for el in dummy_symb:
                if (not (el in symb_unique)):
                    nat=nat+1
                    symb_unique[el]=nat 
            fout.write("%d\n F\n"%len(dummy_symb))
            for i in range(len(dummy_symb)):
                fout.write("%s   %d %30.14f %30.14f %30.14f\n"%(dummy_symb[i],symb_unique[dummy_symb[i]],atpos[i][0]/au2ang,atpos[i][1]/au2ang,atpos[i][2]/au2ang))
            fout.write(" T\n   1     0     0\n   0     1     0\n   0     0     1\n")
            fout.close()
            #run str and modify pipe file
            parentdir=os.getcwd()
            os.chdir(calcdir)
            exitcode, out = commands.getstatusoutput('timeout 60s str')
            if (exitcode==0) and (os.path.isfile("+pipe_wyck")):
                infile=open("+pipe_wyck","r")
                lines=[]
                line=infile.readline()
                init_spin=[]
                while line:
                    for dum_el in dummy_tab:
                        if dum_el in line:
                            el,mom=dummy_tab[dum_el]
                            lspl=line.split("@")
                            init_spin.append(mom)
                            lspl[2]=" %s "%el
                            line=""
                            for x in lspl:
                                line=line+x
                                if x!=lspl[-1]:
                                    line=line+"@"
                            break
                    lines.append(line)
                    line=infile.readline()
                infile.close()
                outfile=open("+pipe_wyck","w")
                for line in lines:
                    outfile.write(line)
                outfile.close()
                #print init_spin
                kw['initial_spin_split_sorts']=init_spin
            else:
                status=False
            os.chdir(parentdir)
        else:
            status=False
        return status,kw


    def get_DOS(self, pDOS=[], calcdir='./', outfile='out', normalization=1.0, converged_only=True,iDOS=False):
        """return dictionary with DOS data:
        'energy', 'total', partial DOS.
        partial contributions are summed over the nonequivalent atoms, e.g. pDOS=['Fe 3d'] will return
        the partial contributions of all Fe 3d orbitals; default is pDOS=[] (only total DOS),
        pDOS=['all'] will return all partial contributions.
        iDOS=True will return the integrated DOS instead (default: iDOS=False)
        """
        DOS={}
        outfilep=os.path.join(calcdir,outfile)
        if (converged_only) and (self.read_convergence(outfile=outfilep)==False):
            return None
        # get DOS filenames
        dosstates={'total':['+dos.total']}
        if iDOS==True:
            dosstates={'total':['+idos.total']}
        if pDOS!=[]:
            basis_states=self.read_basis_states(outfile=outfilep)
            if 'all' in pDOS:
                pDOS.remove('all')
                for state in basis_states:
                    orb="%s %s"%(state['element'],state['orbital'])
                    if not (orb in pDOS):
                        pDOS.append(orb)
            #print basis_states
            for orbitals in pDOS:
                el,orb=orbitals.split()
                dosfiles=[]
                for state in basis_states:
                    if (el==state['element']) and (orb==state['orbital']):
                        fname='+dos.sort%3d.nl%3d'%(state['sort'],state['nl'])
                        if iDOS==True:
                            fname='+idos.sort%3d.nl%3d'%(state['sort'],state['nl'])
                        dosfiles.append(fname.replace(' ', '0'))
                dosstates[orbitals]=dosfiles
        energy=[]
        do_energy=True
        for orbitals in dosstates:
            init_orbit=True
            print orbitals
            dos=[]
            dos2=[]
            ispin=1
            for fname in dosstates[orbitals]:
                dosfile=os.path.join(calcdir,fname)
                print dosfile
                if os.path.isfile(dosfile):
                    infile=open(dosfile,'r')
                    line=infile.readline()
                    i=0
                    while line:
                        if (line.startswith('#')) and ('spin =' in line):
                            ispin=int(line.split('spin =')[1])
                            print "ispin=%d"%ispin
                            if ispin==2:
                                do_energy=False
                                i=0
                        elif (line.startswith('#')) or (line.isspace()):
                            # ispin
                            print 'skip comment/empty line (todo!):',line
                        else:
                            #print line
                            if do_energy:
                                x=float(line.split()[0])
                                energy.append(x)
                            y=float(line.split()[1])*normalization
                            if ispin==1:
                                if init_orbit:
                                    dos.append(y)
                                else:
                                    dos[i]=dos[i]+y
                            else:
                                if init_orbit:
                                    dos2.append(y)
                                else:
                                    dos2[i]=dos2[i]+y
                            i=i+1
                        line=infile.readline()
                    do_energy=False
                    init_orbit=False
            if dos2==[]:
                DOS[orbitals]=dos
            else:
                DOS[orbitals]=[dos,dos2]
        DOS['energy']=energy
        return DOS
                
    def get_iDOS(self, pDOS=[], calcdir='./', outfile='out', normalization=1.0, converged_only=True):
        """return dictionary with integrated DOS data:
        'energy', 'total', partial DOS.
        partial contributions are summed over the nonequivalent atoms, e.g. pDOS=['Fe 3d'] will return
        the partial contributions of all Fe 3d orbitals; default is pDOS=[] (only total DOS),
        pDOS=['all'] will return all partial contributions.
        """
        return self.get_DOS(pDOS=pDOS, calcdir=calcdir, outfile=outfile, normalization=normalization, converged_only=converged_only, iDOS=True)
    
    def get_lDOS(self, ldosstates=[], calcdir='./', outfile='out', normalization=1.0, converged_only=True):
        """read local DOS: 'Fe 3d' -> lm rsolved DOS for 3d states of all Fe sites
        """
        ao=self.read_atoms(calcdir=calcdir, outfile=outfile)
        basis_states=self.read_basis_states(calcdir=calcdir, outfile=outfile) #outfile=os.path.join(calcdir,outfile))
        sites=ao.get_chemical_symbols()
        DOS={}
        energy=[]
        do_energy=True
        for ldosstate in ldosstates:
            sym,orb=ldosstate.split()
            for state in basis_states:
                if (sym==state['element']) and (orb==state['orbital']):
                    nl=state['nl']
                    break
            for isite in range(len(sites)):
                if (sites[isite]==sym):
                    ldosfile=os.path.join(calcdir,str("+ldos.site%3d.nl%3d"%(isite+1,nl)).replace(' ', '0'))
                    print ldosstate,ldosfile
                    infile=open(ldosfile,'r')
                    line=infile.readline()
                    i=0
                    while line:
                        if (line.startswith('#')):
                            m=int(line.split("=")[3].split()[0])
                            s=int(line.split("=")[4])
                            label="%s_%d %s"%(sym,isite,orb)
                            if not label in DOS:
                                DOS[label]={}
                            print ldosstate,m,s,label
                            dos=[]
                        elif (line.isspace()):
                            DOS[label][(m,s)]=dos
                            do_energy=False
                        else:
                            if do_energy:
                                x=float(line.split()[0])
                                energy.append(x)
                            y=float(line.split()[1])*normalization
                            dos.append(abs(y))
                        line=infile.readline()
        DOS['energy']=energy
        return DOS
        
    def read_atoms(self, calcdir='./', outfile='out', core=False):
        au2AA=0.5291772083
        out=os.path.join(calcdir, outfile)
        if os.path.isfile(out):
            infile=open(out,'r')
            line=infile.readline()
            while line:
                if line.startswith('lattice vectors'):
                    a=[]
                    for i in range(3):
                        line=infile.readline().split(':')[1].split()
                        a.append((float(line[0])*au2AA,float(line[1])*au2AA,float(line[2])*au2AA))
                elif line.startswith('Number of sites :'):
                    nsite=int(line.split(':')[1])
                elif line.startswith('No.  Element WPS CPA-Block    X                      Y                      Z'):
                    symbols=[]
                    atpos=[]
                    for i in range(nsite):
                        line=infile.readline().split()
                        symbols.append(line[1])
                        atpos.append((float(line[4])*au2AA,float(line[5])*au2AA,float(line[6])*au2AA))
                line=infile.readline()
            infile.close()
            ao=Atoms(symbols=symbols,positions=atpos,pbc=True,cell=a)
            return ao
        else:
            print 'Could not open ',out
            return None

      
    def get_magnetic_moments(self, calcdir='./', outfile='out', converged_only=True):
        """
        """
        magmoms=None
        outfilep=os.path.join(calcdir,outfile)
        if (converged_only) and (self.read_convergence(outfile=outfilep)==False):
            return None
        if os.path.isfile(outfilep):
            infile=open(outfilep,'r')
            line=infile.readline()
            nsite=0
            while line:
                if line.startswith('Number of sites :'):
                    nsite=int(line.split(':')[1])
                elif ("MAG.MOMENT" in line) and (nsite>0):
                    magmoms=[]
                    line=infile.readline()
                    for i in range(nsite):
                        line=infile.readline().split('|')
                        magmoms+=[float(line[3])]
                line=infile.readline()
            infile.close()
        else:
            print 'Could not open ',outfilep
        return magmoms

    def get_orbital_moments(self, calcdir='./', outfile='out', converged_only=True):
        """
        """
        magmoms=None
        outfilep=os.path.join(calcdir,outfile)
        if (converged_only) and (self.read_convergence(outfile=outfilep)==False):
            return None
        if os.path.isfile(outfilep):
            infile=open(outfilep,'r')
            line=infile.readline()
            nsite=0
            while line:
                if line.startswith('Number of sites :'):
                    nsite=int(line.split(':')[1])
                elif ("ONSITE ORBITAL MOMENTS" in line) and (nsite>0):
                    magmoms=[]
                    line=infile.readline()
                    line=infile.readline()
                    for i in range(nsite):
                        line=infile.readline().split()
                        magmoms+=[float(line[3])]
                line=infile.readline()
            infile.close()
        else:
            print 'Could not open ',outfilep
        return magmoms

    def get_orbital_moment(self, calcdir='./', outfile='out', converged_only=True):
        """
        """
        magmoms=None
        outfilep=os.path.join(calcdir,outfile)
        if (converged_only) and (self.read_convergence(outfile=outfilep)==False):
            return None
        if os.path.isfile(outfilep):
            infile=open(outfilep,'r')
            line=infile.readline()
            nsite=0
            while line:
                if line.startswith('Number of sites :'):
                    nsite=int(line.split(':')[1])
                elif ("ONSITE ORBITAL MOMENTS" in line) and (nsite>0):
                    magmoms=0.
                    line=infile.readline()
                    line=infile.readline()
                    for i in range(nsite):
                        line=infile.readline().split()
                        magmoms=magmoms+float(line[3])
                line=infile.readline()
            infile.close()
        else:
            print 'Could not open ',outfilep
        return magmoms

    def read_magnetic_moment(self, calcdir='./', outfile='out', converged_only=True):
        """
        """
        magmom=None
        outfilep=os.path.join(calcdir,outfile)
        if (converged_only) and (self.read_convergence(outfile=outfilep)==False):
            return None
        if os.path.isfile(outfilep):
            infile=open(outfilep,'r')
            line=infile.readline()
            while line:
                if 'mag. moment' in line:
                    line=infile.readline()
                    magmom=float(line.split()[2])
                line=infile.readline()
            infile.close()
        else:
            print 'Could not open ',outfilep
        return magmom

    def get_number_of_spins(self, calcdir='./', outfile='out'):
        """
        """
        mspin=None
        outfilep=os.path.join(calcdir,outfile)
        if os.path.isfile(outfilep):
            infile=open(outfilep,'r')
            line=infile.readline()
            while line:
                if "int mspin" in line:
                    if (not '=' in line) or (not ';' in line.split('=')[1]):
                        line=line+infile.readline()
                    mspin=int(line.split("={")[1].split(",")[0])
                line=infile.readline()
            infile.close()
        return mspin
