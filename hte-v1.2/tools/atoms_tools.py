from ase import Atoms, Atom
#from math import *
#from numpy import *
from hte import *

def get_minimal_bond_distance(ao, element1, element2):
    dist=None
    if ao==None:
        return dist
    for i in range(len(ao.arrays['numbers'])):
        for j in range(len(ao.arrays['numbers'])):
            if (ao.get_chemical_symbols()[i]==element1) and (ao.get_chemical_symbols()[j]==element2):
                if i==j:
                    cell=ao.get_cell()
                    d=min(norm(cell[0]),norm(cell[1]),norm(cell[2]))
                else:
                    d=ao.get_distance(i,j,mic=True)
                if (dist==None) or (d<dist):
                    dist=d
    return dist


def get_latticeparameters(ao):
    cell=ao.get_cell()
    lattice_parameters=(norm(cell[0]),norm(cell[1]),norm(cell[2]))
    angles=(arccos(dot(cell[1], cell[2])/(b*c))*180./pi,arccos(dot(cell[0], cell[2])/(a*c))*180./pi,arccos(dot(cell[0], cell[1])/(a*b))*180./pi)
    return lattice_parameters,angles

def scale_lattice(ao, constraints):
    #start with a volume scaling
    vol=ao.get_volume()
    scale_coa=False
    if ('spacegroup' in ao.info) and (ao.info['spacegroup'].no<195):
        scale_coa=True
        cell=ao.get_cell()
        coa=norm(cell[2])/norm(cell[0])
    sum_atvol=0.0
    if 'atomic_volumes' in constraints:
        volconstr=constraints['atomic_volumes']
        for el in ao.get_chemical_symbols():
            sum_atvol=sum_atvol+volconstr[el]
    if 'bond_distances' in constraints:
        bonddis=constraints['bond_distances']
    #check how initial structure matches constraints
    sigma=0.5*sum_atvol
    match=exp(-0.5*pow((vol-sum_atvol)/sigma,2))#/(sigma*sqrt(2*pi))
    if 'bond_distances' in constraints:
        bonddis=constraints['bond_distances']
        for (el1,el2) in bonddis:
            dopt,sigma=bonddis[(el1,el2)]
            #sigma=sigma*dopt
            dao=get_minimal_bond_distance(ao, el1, el2)
            match=match*exp(-0.5*pow((dopt-dao)/sigma,2))#/(sigma*sqrt(2*pi))
    print ao.get_cell()
    print '...match of initial structure:',match
    best_match=match
    for vscale in linspace(0.8,1.2,10):
        for scalec in linspace(0.5,1.5,10):
          for scaleb in linspace(0.5,1.5,10):
            print "**",vscale,scalec,scaleb
            match=1.0
            aonew=ao.copy()
            cell=aonew.get_cell()
            if scale_coa:
                cell[2]=cell[1]*scalea
                cell[1]=cell[1]*scaleb
                aonew.set_cell(cell, scale_atoms=True)
            vol=aonew.get_volume()
            aonew.set_cell(pow(vscale*sum_atvol/vol,1./3.)*cell, scale_atoms=True)
            vol=aonew.get_volume()
            #check how scaled lattice matches
            sigma=0.5*sum_atvol
            match=exp(-0.5*pow((vol-sum_atvol)/sigma,2))#/(sigma*sqrt(2*pi))
            print 'volumes:',vol,sum_atvol,match,sigma
            if 'bond_distances' in constraints:
                bonddis=constraints['bond_distances']
                for (el1,el2) in bonddis:
                    dopt,sigma=bonddis[(el1,el2)]
                    dao=get_minimal_bond_distance(aonew, el1, el2)
                    #sigma=sigma*dopt
                    matchbd=exp(-0.5*pow((dopt-dao)/sigma,2))#/(sigma*sqrt(2*pi))
                    print 'bd',el1,el2,dopt,dao,matchbd
                    match=match*matchbd
            #print aonew.get_cell()
            print vscale,scalec,'... match=',match
            if match>best_match:
                best_match=match
                ao_best=aonew.copy()
                print 'new best match:',best_match
                print ao_best.get_cell()
    print 'best match:',best_match
    print ao_best.get_cell()
            

def substitute_atoms(ao,atom_list):
    new_atoms_object=ao.copy()
    new_chem_symbols=new_atoms_object.get_chemical_symbols()
    #print "SUBSTIn",new_chem_symbols,atom_list
    try:
        wyckoffs=spglib.get_symmetry_dataset(new_atoms_object)['wyckoffs']
    except:
        wyckoffs=len(new_chem_symbols)*['']
    uidxts={}
    for i in range(len(new_chem_symbols)):
        el_w="%s_%s"%(new_chem_symbols[i],wyckoffs[i])
        #print new_chem_symbols[i], new_chem_symbols[i] in atom_list
        if new_chem_symbols[i] in atom_list:
            uidxts[new_chem_symbols[i]]=atom_list[new_chem_symbols[i]]
            new_chem_symbols[i]=atom_list[new_chem_symbols[i]]
        elif el_w in atom_list:
            uidxts[el_w]=atom_list[el_w]
            new_chem_symbols[i]=atom_list[el_w]
    new_atoms_object.set_chemical_symbols(new_chem_symbols)
    #print "SUBST",new_chem_symbols
    return new_atoms_object,uidxts
