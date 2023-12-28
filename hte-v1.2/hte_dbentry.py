#
import sys
import os
import commands
import shutil
from copy import deepcopy
from ase import Atoms, Atom
from ase.calculators.vasp import Vasp
from ase.io.cif import *
from ase.io import read
from ase.io import write
from ase.io.vasp import read_vasp
from ase.io.vasp import write_vasp
from ase.visualize import view
from ase.data import atomic_numbers, chemical_symbols
from ase.lattice.spacegroup import *
from transport import *
from tools import *

try:
    import spglib
    has_spglib=True
except:
    has_spglib=False

class HTEdbentry(object):
    def __init__(self, atoms_obj, calcdir, comment=None,_source_info=None):
        """internal class for database entries of HTE not intended for users
        should only be accessed via routines of HTE
        """
        self.calcdir=calcdir # typically uid (no abs. path!)
        self.atoms_initial=atoms_obj # always keep the initial structure
        #info about structure history: tuples (hte_method,{arguments})
        if (_source_info==None):
            self.source_info=[('non_HTE',{'comment':comment})]
        else:
            self.source_info=_source_info
        self.cif_info={} # store information from cif file
        self.user_info={} # store information from user, e.g. if some data has to be corrected
        self.spglib_info=None
        self.submitted_jobs={}
        self.stored_calc_results={}
        self.do_not_calc=False
        self.comment=comment
        self.similar_to=[]
        self.atoms=None # 

    def get_cif_source(self,methods=['internal','external']):
        """returns a tuple (cif file, index, cif_id)
        if the structure is based on cif data
        """
        cif_file,cif_id,cif_index=None,None,None
        method,info=self.source_info[0]
        if (method=='import_cif'):
            if 'cif_id' in info:
                cif_id=info['cif_id']
            for meth in methods:
                if meth=='internal':
                    if ('cif_file' in info) and (os.path.isfile(info['cif_file'])):
                        cif_file,cif_index=info['cif_file'],0
                        break
                elif meth=='external':
                    if ('cif_file_external' in info) and (os.path.isfile(info['cif_file_external'])):
                        cif_file,cif_index=info['cif_file_external'],info['cif_index_external']
                        break
        return cif_file,cif_id,cif_index
    
    def get_source_info(self):
        return self.source_info

    def get_structure_type(self,methods=['user','cif','atoms_info']):
        """get the structure type
        """
        for method in methods:
            if (method=='user') and ('_chemical_name_structure_type' in self.user_info):
                return self.user_info['_chemical_name_structure_type']
            if (method=='cif') and ('_chemical_name_structure_type' in self.cif_info):
                return self.cif_info['_chemical_name_structure_type']
            if (method=='atoms_info') and ('structure_type' in self.atoms_initial.info):
                #if old version could not be converted
                return self.atoms_initial.info['structure_type']
        return None

    def get_wyckoff(self,sort=None, store=False):
        """get the wyckoff positions, multiplicities and element names.
           if sort=True: returns a list of tuples, where each tuple contains all information regarding one
           position,  e.g.: At '4a' sits 'Fe':sort:[('a','Fe',4)], unsorted:[['a'],['Fe'],[4]]
           if sort=None: return dictionary.
           if store=True: store the wyckoff information in structureDB, return nothing.
        """
        wyckoff_info={}
        if ('_atom_site_wyckoff_symbol' in self.cif_info) and ('wyckoff_symbol' not in wyckoff_info):
            wyckoff_info['wyckoff_symbol']=self.cif_info['_atom_site_wyckoff_symbol']

        if ('wyckoff_elem' in self.cif_info) and ('_atom_site_type_symbol' in self.cif_info):
            wyckoff_info['wyckoff_elem']=self.cif_info['wyckoff_elem']
        elif (not 'wyckoff_elem' in self.cif_info) and ('_atom_site_type_symbol' in self.cif_info):
            wyckoff_info['wyckoff_elem']=self.cif_info['_atom_site_type_symbol']

        if ('_atom_site_symmetry_multiplicity' in self.cif_info) and ('multiplicity' not in wyckoff_info):
            wyckoff_info['multiplicity']=self.cif_info['_atom_site_symmetry_multiplicity']
        if (wyckoff_info) and (sort==True):
            return zip(*wyckoff_info.values())
        if (wyckoff_info) and (sort==False):
            return wyckoff_info.values()
        if (wyckoff_info) and (sort==None):
            return wyckoff_info
        if (wyckoff_info) and (sort=='string'):
            sort_ready=zip(*wyckoff_info.values())
            sort_ready.sort(key=lambda x:x[1])
            return ''.join(zip(*sort_ready)[2])
        return None

    def is_runnable(self):
        if self.do_not_calc:
            return False
        if self.atoms_initial==None:
            return False
        if self.is_disordered():
            return False
        return True

    
    def is_disordered(self):
        """ todo
        """
        if '_atom_site_occupancy' in self.cif_info:
            occ=self.cif_info['_atom_site_occupancy']
            for i in range(len(occ)):
                if float(occ[i])<0.99:
                    return True
        else:
            print '_atom_site_occupancy not given'
        return False
                
    def get_data_from_cif_source(self,cif_file=None,index=0):
        """extract data from cif source
        """
        warnings=[]
        if cif_file==None:
            cif_file,cif_id,cif_index=self.get_cif_source()
        if (cif_file) and (os.path.isfile(cif_file)):
            try:
                x,self.cif_info=parse_cif(cif_file)[index]
            except:
                warnings.append('failed to read cif file %s, index %d'%(cif_file,index))
        return warnings
        
    def get_spacegroup_number(self, methods=['spglib','atoms.info','cif','structure_type']):
        spg_no=None
        for method in methods:
            if (method=='atoms.info'):
                if 'spacegroup' in self.atoms_initial.info:
                    spg_no=self.atoms_initial.info['spacegroup'].no
                    break
            elif (method=='spglib'):
                if (self.spglib_info!=None):
                    spg_no=self.spglib_info['number']
                    break
            elif (method=='structure_type'):
                if (self.get_structure_type()!=None) and (len(self.get_structure_type().split(','))==3):
                    try:
                        spg_no=int(self.get_structure_type().split(',')[2])
                        break
                    except:
                        continue
            elif (method=='cif'):
                if '_symmetry_int_tables_number' in self.cif_info:
                    spg_no=int(self.cif_info['_symmetry_int_tables_number'])
                    break
                elif '_space_group.it_number' in self.cif_info:
                    spg_no=int(self.cif_info['_space_group.it_number'])
                    break
        return spg_no

        
    def get_composition(self, reduce=False, methods=['atoms','cif'], latex=False, gle=False, sort='alpha'):
        """returns the chemical composition as formula sum (reduce=True)
        or as dictionary (reduce=False); possible methods:
        'atoms': get composition from atoms object
        'cif': get composition from information in cif file
        default settings: reduce=False, methods=['atoms','cif']
        """
        comp={}
        formula=''
        for method in methods:
            if (method=='atoms') and (self.atoms_initial!=None):
                #return composition of atoms object
                for sym in self.atoms_initial.get_chemical_symbols():
                    if sym in comp:
                        comp[sym]=comp[sym]+1
                    else:
                        comp[sym]=1
                formula=composition2formula(comp,latex=latex, gle=gle, sort=sort)
                #for sym in sorted(comp.keys()):
                #    formula=formula+sym+str(comp[sym])
                break
            elif method=='cif_org':
                if '_chemical_formula_sum' in self.cif_info:
                    formula=self.cif_info['_chemical_formula_sum']
                    comp=formula2composition(formula)
            elif method=='cif':
                if '_chemical_formula_sum' in self.cif_info:
                    formula=self.cif_info['_chemical_formula_sum']
                    comp=formula2composition(formula)
                    source_info=self.get_source_info()
                    for i in range(len(source_info)):
                        method,args=source_info[i]
                        if method=='substitute_atoms':
                            for substitution in args['atom_list']:
                                chem_comp={}
                                for atom_name in comp:
                                    if atom_name in substitution:
                                        newel=substitution[atom_name]
                                    else:
                                        newel=atom_name
                                    if newel in chem_comp:
                                        chem_comp[newel]=chem_comp[newel]
                                        +comp[atom_name]
                                    else:
                                        chem_comp[newel]=comp[atom_name]
                                comp=chem_comp
                        elif method=='create_test_set':
                            chem_comp={}
                            for i in range(len(args['elements'])):
                                chem_comp[args['elements'][i]]=comp[args['cif_elements'][i]]
                            comp=chem_comp
                    formula=composition2formula(comp,latex=latex, gle=gle, sort=sort)
        if reduce:
            return formula
        return comp
           
    def get_atfraction(self,atname):
        n=0
        for Z in self.atoms_initial.arrays['numbers']:
            if chemical_symbols[Z]==atname:
                n=n+1
        return 1.0*n/self.atoms_initial.get_number_of_atoms()


    def self_check(self, warn_list=['non_HTE_source','incomplete_cifdata','missing_cif_source','missing_spglib_symmetry'],error_list=['spacegroup','composition'],force_checks=False, symprec=1e-3):
        """check consistency of HTEdbentry and update information
        if necessary
        """
        warnings=[]
        errors=[]
        if (self.spglib_info==None) or (force_checks):
            if has_spglib:
                try:
                    self.spglib_info=spglib.get_symmetry_dataset(self.atoms_initial, symprec=symprec)
                except:
                    self.spglib_info=None
        method,args=self.get_source_info()[0]
        if method=='non_HTE':
            if 'non_HTE_source' in error_list:
                errors.append('non_HTE_source')
            elif 'non_HTE_source' in warn_list:
                warnings.append('non_HTE_source')
        if method=='import_cif':
            cif_file,cif_id,cif_index=self.get_cif_source()
            if cif_file==None:
                if  'missing_cif_source' in error_list:
                    errors.append('missing_cif_source')
                elif 'missing_cif_source' in warn_list:
                    warnings.append('missing_cif_source')
            # check if information on composition is consistent:
            comp_atoms=self.get_composition(reduce=False, methods=['atoms'])
            comp_cif=self.get_composition(reduce=False, methods=['cif'])
            comp_ok=True
            if ((len(comp_cif)==0) and (len(comp_atoms)==1)):
                if 'incomplete_cifdata' in warn_list:
                    warnings.append('composition not given in cif source')
            elif len(comp_atoms)!=len(comp_cif):
                comp_ok=False
            else:
                fac=None
                for el in comp_atoms:
                    if el in comp_cif:
                        if fac==None:
                            fac=1.0*comp_atoms[el]/comp_cif[el]
                        if fabs(1.0*comp_atoms[el]-comp_cif[el]*fac)>0.001:
                            comp_ok=False
                            break
                    else:
                        comp_ok=False
                        break
            if comp_ok==False:
                if 'composition' in error_list:
                    errors.append('composition of atoms object=%s inconsistent with cif source=%s!'%(str(comp_atoms),str(comp_cif)))
                elif 'composition' in warn_list:
                    warnings.append('composition of atoms object=%s inconsistent with cif source=%s!'%(str(comp_atoms),str(comp_cif)))
        # check if symmetry information is consistent
        spg_atoms=self.get_spacegroup_number(methods=['atoms.info'])
        spg_cif=self.get_spacegroup_number(methods=['cif'])
        spg_type=self.get_spacegroup_number(methods=['structure_type'])
        spg_spglib=self.get_spacegroup_number(methods=['spglib'])
        spg_ok=True
        spgs=[spg_atoms,spg_cif,spg_type,spg_spglib]
        for sg1 in spgs:
            for sg2 in spgs:
                if (sg1!=None) and (sg2!=None) and (sg1!=sg2):
                    spg_ok=False
        if spg_ok==False:
            if 'spacegroup' in error_list:
                errors.append('inconsistent space group numbers: %s (atoms object) - %s (cif) %s (structure type) %s (spglib)!'%(str(spg_atoms),str(spg_cif),str(spg_type),str(spg_spglib)))
            elif 'spacegroup' in warn_list:
                warnings.append('inconsistent space group numbers: %s (atoms object) - %s (cif) %s (structure type) %s (spglib)!'%(str(spg_atoms),str(spg_cif),str(spg_type),str(spg_spglib)))
        if spg_spglib==None:
            if 'missing_spglib_symmetry' in error_list:
                errors.append('no symmetry information from spglib')
            elif 'missing_spglib_symmetry' in warn_list:
                warnings.append('no symmetry information from spglib')
        #
        if len(warnings)==0:
            warnings=None
        if len(errors)==0:
            errors=None
        return warnings,errors
    
    def update_symmetry_information(self, initial=False, symprec=1e-3):
        # check if atoms.info has necessary symmetry information/update if possible
        tags_for_fplo=['_atom_site_type_symbol','spacegroup','_cell_length_a','_cell_length_b','_cell_length_c','_cell_angle_alpha','_cell_angle_beta','_cell_angle_gamma','_atom_site_fract_x','_atom_site_fract_y','_atom_site_fract_z']
        info_ok=True
        if initial:
            atoms_object=self.atoms_initial
        else:
            atoms_object=self.atoms
        for tag in tags_for_fplo:
            if not (tag in atoms_object.info):
                info_ok=False
        if info_ok:
            return True
        if has_spglib:
            # try to get information from spglib
            cell, scaled_positions, numbers = spglib.refine_cell(atoms_object, symprec=symprec)
            a=norm(cell[0])
            b=norm(cell[1])
            c=norm(cell[2])
            atoms_object.info['_cell_length_a']=a
            #print atoms_object.info['_cell_length_a'],self.atoms_initial.info['_cell_length_a']
            atoms_object.info['_cell_length_b']=b
            #print atoms_object.info['_cell_length_b'],self.atoms_initial.info['_cell_length_b']
            atoms_object.info['_cell_length_c']=c
            #print atoms_object.info['_cell_length_c'],self.atoms_initial.info['_cell_length_c']
            atoms_object.info['_cell_angle_alpha']=arccos(dot(cell[1], cell[2])/(b*c))*180./pi
            #print atoms_object.info['_cell_angle_alpha'],self.atoms_initial.info['_cell_angle_alpha']
            atoms_object.info['_cell_angle_beta']=arccos(dot(cell[0], cell[2])/(a*c))*180./pi
            #print atoms_object.info['_cell_angle_beta'],self.atoms_initial.info['_cell_angle_beta']
            atoms_object.info['_cell_angle_gamma']=arccos(dot(cell[0], cell[1])/(a*b))*180./pi
            #print atoms_object.info['_cell_angle_gamma'],self.atoms_initial.info['_cell_angle_gamma']
            sym_data=spglib.get_symmetry_dataset(atoms_object, symprec=symprec)
            if ('spacegroup' in atoms_object.info) and (atoms_object.info['spacegroup'].no!=sym_data['number']):
                print '*** FATAL ERROR in update_symmetry_information():  spglib spacegroup differs from atoms.info !!!!'
            atoms_object.info['spacegroup']=Spacegroup(sym_data['number'])
            atoms_object.info['_atom_site_type_symbol']=[]
            atoms_object.info['_atom_site_fract_x']=[]
            atoms_object.info['_atom_site_fract_y']=[]
            atoms_object.info['_atom_site_fract_z']=[]
            for i in np.unique(sym_data['equivalent_atoms']):
                atoms_object.info['_atom_site_type_symbol'].append(chemical_symbols[numbers[i]])
                atoms_object.info['_atom_site_fract_x'].append(scaled_positions[i][0])
                atoms_object.info['_atom_site_fract_y'].append(scaled_positions[i][1])
                atoms_object.info['_atom_site_fract_z'].append(scaled_positions[i][2])
            info_ok=True
        return info_ok
    
    def setup_calculation(self, subdir, calc, bs_calc=False):
        if (self.atoms==None) or (calc==None):
            return False
        parentdir=os.getcwd()
        calcdir=os.path.join(self.calcdir,subdir)
        if not os.path.isdir(calcdir):
            print 'Creating directory ',calcdir
            os.makedirs(calcdir)
        if (calc.name.lower()=='fplo') and ('initdir' in calc.settings):
            for fname in glob.glob(os.path.join(calc.settings['initdir'],"=.*")):
                shutil.copy(fname,calcdir)
            for fname in glob.glob(os.path.join(calc.settings['initdir'],"+*")):
                shutil.copy(fname,calcdir)
        os.chdir(calcdir)
        if calc.name.lower()=='fplo':
                if 'kmesh' in calc.settings:
                    kmesh=calc.settings['kmesh']
                    print kmesh
                else:
                    kmesh=[12,12,12]
                for TOL in [1.e-6,1.e-5,1.e-4]:
                    status=calc.initialize(self.atoms, bs_calc=bs_calc,k_mesh=kmesh,TOL=TOL)
                    no,sym=calc.get_spacegroup()
                    if no==self.get_spacegroup_number():
                        break
                    #else: #TODO: cross check against spglib
                    #    status=False
        elif calc.name.lower()=='vasp':
                calc.initialize(self.atoms)
                # Write input (may be done twice for older ase versions)
                write_vasp('POSCAR', calc.atoms_sorted, symbol_count = calc.symbol_count)
		shutil.copy('POSCAR','POSCAR.INI')
                calc.write_incar(self.atoms)
                calc.write_potcar()
                calc.write_kpoints()
                calc.write_sort_file()
                status=True #hopefully
        else:
            status=False
        os.chdir(parentdir)
        return status

    def get_bandgap(self, subdir, calc, update=False, silent=False, job_commands=None, unit='eV', sloppy_mode=True, nsub_max=2):
        gap=None
        if (calc==None):
            return gap
        if calc.name.lower()=='fplo':
            gap=calc.read_bandgap(outfile=os.path.join(self.calcdir,subdir,'out'))
        elif calc.name.lower()=='vasp':
            gap=get_vasp_bandgap(pathname=os.path.join(self.calcdir,subdir))
        else:
            return gap
        if (update==True) and (gap==None):
            self.run_calculation(subdir, calc, silent=silent, job_commands=job_commands, nsub_max=nsub_max)
        return gap
        
    def get_DOS_at_E_F(self, subdir, calc, update=False, silent=False, job_commands=None):
        if calc==None:
            return None
        if calc.name.lower()=='fplo':
            return calc.read_DOS_at_E_F(outfile=os.path.join(self.calcdir,subdir,'out'))
        else:
            return None
        
    def get_energy(self, subdir, calc, update=False, silent=True, job_commands=None, unit='eV', sloppy_mode=True, nsub_max=2, settings={}):
        if (sloppy_mode==True) and ('energy' in self.stored_calc_results):
            edict=self.stored_calc_results['energy']
            if subdir in edict:
                #print 'using stored energy',self.calcdir,subdir,edict[subdir]
                return edict[subdir]*get_conversion('eV',unit) 
        if calc==None:
            return None
        parentdir=os.getcwd()
        energy=None
        calcdir=os.path.join(self.calcdir,subdir)
        if calc.name.lower()=='fplo':
            energy=calc.read_energy(outfile=os.path.join(calcdir,'out'))
            if energy!=None:
                energy=energy*get_conversion('Hartree',unit)
        elif calc.name.lower()=='vasp':
            if self.check_convergency(subdir,calc):
                outcar=os.path.join(calcdir,'OUTCAR')
                if os.path.isfile(outcar):
                    try:
                        for line in open(outcar, 'r'):
                            if line.startswith('  energy  without entropy'):
                                energy=float(line.split()[-1])
                    except:
                        energy=None
                        print "WARNING(get_energy): check %s"%outcar
            if energy!=None:
                energy=energy*get_conversion('eV',unit) 
        if (update==True) and (energy==None):
            self.run_calculation(subdir, calc, silent=silent, check_convergency=False, job_commands=job_commands, nsub_max=nsub_max, settings=settings)
        if (energy!=None):
            if 'energy' in self.stored_calc_results:
               edict=self.stored_calc_results['energy']
            else:
                edict={}
            edict[subdir]=energy*get_conversion(unit,'eV')
            self.stored_calc_results['energy']=edict
            #print 'storing energy',subdir,edict[subdir]
        return energy

    def check_convergency(self, subdir, calc, maxdV=0.1):
        print "check_point77, entering the check_convergency subroutine"
        if calc==None:
            return False
        parentdir=os.getcwd()
        converged=False
        calcdir=os.path.join(self.calcdir,subdir)
        print "check_point78, calcdir is:",calcdir
        if calc.name.lower()=='vasp':
            outcar=os.path.join(calcdir,'OUTCAR')
            if os.path.isfile(outcar):
                os.chdir(calcdir)
                try:
                    converged=calc.read_convergence()
                    print "check_point79, here is converged", converged
                except:
                    print "WARNING: Check file ",outcar
                if converged==True:
                    try:
                        #check if volume changed more than maxdV
                        ao_out=read_vasp()
                        ao_in=read_vasp(filename='POSCAR')
                        if (abs(1.0-(ao_in.get_volume()/ao_out.get_volume()))>maxdV):
                            converged=False
                            print "Restart in ",calcdir," because volume change larger than maxdV"
                    except:
                        print "WARNING: Check files ",outcar
                os.chdir(parentdir)
        elif calc.name.lower()=='fplo':
            outfilename=os.path.join(calcdir,'out')
            converged=calc.read_convergence(outfile=outfilename)
            
            print "check_point80, final return converged is:",converged
        return converged
    
    def get_converged_structure(self, subdir, calc, silent=False):
        atoms=None
        parentdir=os.getcwd()
        if self.check_convergency(subdir, calc):
            calcdir=os.path.join(self.calcdir,subdir)
            try:
                os.chdir(calcdir)
                if calc.name.lower()=='fplo':
                    atoms=calc.read_atoms()
                elif calc.name.lower()=='vasp':
                    atoms=read_vasp()
            except:
                if silent==False:
                    print "WARNING(get_converged_structure): failed to read structure in ",calcdir
        os.chdir(parentdir)
        return atoms

    def get_jobstatus(self, subdir, calc, job_environment='SGE'):
        #returns a dictionary with status information about the calculation
        #keys: cputime, queue_status, converged
        status={}
        status['converged']=False
        status['energy']=None
        status['cputime']=None
        status['qstat']=None
        status['nsubmit']=None
        #
        parentdir=os.getcwd()
        calcdir=os.path.join(self.calcdir,subdir,"")
        print "check_point59 , calcdir is:",calcdir
        #check if job is still in the queue:
        if calcdir in self.submitted_jobs:
            jobid=self.submitted_jobs[calcdir]['jobid']
            print "check_point55 , jobid is:",jobid
            print "check_point56 , job_environment is:",job_environment
            if job_environment=='SLURM':
                qstatcmd='squeue -j '+str(jobid)+'| grep '+str(jobid) #squeue seems to have some delay 
                print "check_point57, qstatcmd is:",qstatcmd
            else: # SGE
                qstatcmd='qstat -j '+str(jobid)
            exitcode, out = commands.getstatusoutput(qstatcmd)
            if exitcode==0:
                status['qstat']='in queue(%s)'%(jobid)
            else:
                print "check_point58 ,status['qstat']= is done"
                status['qstat']='done(%s)'%(jobid)
            status['nsubmit']=self.submitted_jobs[calcdir]['nsubmit']
        if calc==None:
            print "check_point60 , calc is None, return status"
            return status
        #check if calculation is converged and get energy etc.
        if self.check_convergency(subdir, calc):
            print "check_point46, jobstatus check_convergency"
            status['converged']=True
            if calc.name.lower()=='vasp':
                outcar=os.path.join(calcdir,'OUTCAR')
                infile=open(outcar,'r')
                line=infile.readline()
                while line:
                    if 'Total CPU time used' in line:
                        cputag=line.split(':')
                        status['cputime']=float(cputag[1])
                    line=infile.readline()
                infile.close()
            status['energy']=self.get_energy(subdir, calc)
        return status
    
    def get_energy_per_atom(self, subdir, calc, update=False, silent=False, job_commands=None, unit='eV', sloppy_mode=True, nsub_max=2):
        E=self.get_energy(subdir, calc, update=update, silent=silent, job_commands=job_commands, unit=unit, sloppy_mode=sloppy_mode, nsub_max=nsub_max)
        if E!=None:
            return E/self.atoms.get_number_of_atoms()
        else:
            return None
        
    def run_calculation(self, subdirx, calc, silent=False, check_convergency=True, nsub_max=2, job_commands=None, settings={}):
        update_input=False
        if calc.name.lower()=='fplo': #TODO (for the moment o.k.)
            update_input=True
        # check status of calculation and run if necessary
        if (calc==None) or (job_commands==None):
            if (silent==False):
                print "WARNING(run_calculation): Check calculator=%s/job_commands=%s"%(str(calc),str(job_commands))
            return False
        if ('job_environment' in job_commands) and (job_commands['job_environment'] in ['SGE','SLURM']):
            job_environment=job_commands['job_environment']
        else:
            print "WARNING(run_calculation): job_environment %s, no job control!"
            job_environment=None
        parentdir=os.getcwd()
        subdir=subdirx
        if (calc.name.lower()=='fplo') and ('subdir' in calc.settings):
            subdir=os.path.join(subdirx,calc.settings['subdir'])
        calcdir=os.path.join(self.calcdir,subdir)
        if (check_convergency==True) and (self.check_convergency(subdir, calc)):
            #Nothing to do, calculation is converged, but check if calculations is in self.submitted_jobs:
            if not (calcdir in self.submitted_jobs):
                self.submitted_jobs[calcdir]={'jobid':'unknown','nsubmit':-1}
                if settings!={}:
                    self.submitted_jobs[calcdir]['settings']=settings
            return False
        #check if job is still in the queue:
        jobid=None
        if (job_environment!=None) and (calcdir in self.submitted_jobs):
            jobid=self.submitted_jobs[calcdir]['jobid']
            if job_environment=='SLURM':
                qstatcmd='squeue -j '+str(jobid)+'| grep '+str(jobid)
            else: # SGE
                qstatcmd='qstat -j '+str(jobid)
            exitcode, out = commands.getstatusoutput(qstatcmd)
            if exitcode==0:
                if silent==False:
                    print 'calculation running, jobid is',self.submitted_jobs[calcdir]['jobid']
                return False
            if silent==False:
                print 'jobid ',self.submitted_jobs[calcdir]['jobid'],' has finished'
            # do not resubmit a job more than nsub_max
            if self.submitted_jobs[calcdir]['nsubmit']>nsub_max:
                print 'WARNING: check calculation in ',calcdir
                return False
        jobid1=jobid
        fname=os.path.join(calcdir,'jobid')
        if os.path.isfile(fname):
            infile=open(fname,'r')
            jobid=infile.readline().strip()
            infile.close()
            if jobid!=jobid1:
                if job_environment=='SLURM':
                    qstatcmd='squeue -j '+str(jobid)+'| grep '+str(jobid)
                else: # SGE
                    qstatcmd='qstat -j '+str(jobid)
                exitcode, out = commands.getstatusoutput(qstatcmd)
                if exitcode==0:
                    if silent==False:
                        print 'calculation running, jobid is',jobid
                    return False
        if (not (os.path.isdir(calcdir))) or (update_input==True):
            if silent==False:
                print 'set up calculation in ',subdir
            cont=self.setup_calculation(subdir, calc)
            if cont==False:
                #setup_calculation failed for some reason(e.g. we are still waiting for some other calculation to finish)
                return cont
        #copy jobfile(s)?
        if ('job_file' in job_commands) and (job_commands['job_file']!=None):
            jobfiles=glob.glob(job_commands['job_file'])
            for jobfile in jobfiles:
                if silent==False:
                    print "run_calculation: Copying %s in directory %s"%(jobfile,calcdir)
                shutil.copy(jobfile,calcdir)
        os.chdir(calcdir)
        # create a file with estimated resources (for the moment just number of atoms, second parameter number of submissions)
        outfile=open("hte_estimate_resource.txt","w")
        if (job_environment!=None) and (calcdir in self.submitted_jobs):
            nsubmit=self.submitted_jobs[calcdir]['nsubmit']+1
        else:
            nsubmit=1
        if 'array_job' in job_commands:
            outfile.write("%d %d %s %s"%(len(self.atoms),nsubmit,self.calcdir,calcdir))
        else:
            outfile.write("%d %d"%(len(self.atoms),nsubmit))
        outfile.close()
        if calc.name.lower()=='vasp':
            # continue calculation from CONTCAR if present
            if os.path.isfile('CONTCAR'):
                try:
                    if silent==False:
                        print "Trying to continue from CONTCAR"
                    atoms=read_vasp()
                    if not os.path.isfile('POSCAR.INI'):
                        shutil.copy('POSCAR','POSCAR.INI')
                    shutil.copy('CONTCAR','POSCAR')
                except:
                    #restart from scratch, no valid CONTCAR
                    if silent==False:
                        print "failed, restarting from scratch"
                    shutil.copy('POSCAR.INI','CONTCAR')
		    #os.chdir(parentdir)
		    #cont=self.setup_calculation(subdir, calc)
		    #if cont==False:
	            #    #setup_calculation failed for some reason(e.g. we are still waiting for some other calculation to finish)
		    #	print "restart in %s failed, abort"%calcdir
        	    #    return cont
		    #os.chdir(calcdir)
        if ('job_command' in job_commands):
            exitcode, out = commands.getstatusoutput(job_commands['job_command'])
        else:
            if silent==False:
                print "WARNING(run_calculation): No job command"
            return False
        if (job_environment==None) or ('array_job' in job_commands):
            os.chdir(parentdir)
            print exitcode, out
            return True
        if not (calcdir in self.submitted_jobs):
            self.submitted_jobs[calcdir]={'nsubmit':1}
        else:
            self.submitted_jobs[calcdir]['nsubmit']=self.submitted_jobs[calcdir]['nsubmit']+1
        if settings!={}:
            self.submitted_jobs[calcdir]['settings']=settings
        jobid=-1
        if exitcode == 0:
            if job_environment=='SLURM':
                outspl=out.split('Submitted batch job')
                jobid=outspl[1].split()[0]
            else: # job_environment=='SGE':
                outspl=out.split()
                jobid=outspl[2]
            commands.getstatusoutput('echo %s > jobid'%jobid)
            if silent==False:
                print 'Job submitted, id=',jobid
        elif silent==False:
            print 'Job submission failed'
        self.submitted_jobs[calcdir]['jobid']=jobid
        os.chdir(parentdir)
        return True
    
    def get_number_of_elements(self):
        return len(self.get_composition())
    
        
    def get_spacegroup(self):
        if (hasattr(self.atoms_initial.info, 'spacegroup')):
            return self.atoms_initial.info['spacegroup']
        else:
            return 'None'
