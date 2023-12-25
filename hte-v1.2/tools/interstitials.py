
from hte import *
from ase.data import *
try:
    import spglib
    has_spglib=True
except:
    has_spglib=False

def CreateInterstitialDefect(atoms, interstitial_atom, rmin='auto', Nmax=5, nscan=24, symprec=1e-3):
    """ Identify potential interstitial positions and create a defect cells
    with one interstiatial atom at this position
    """
    interst_pos={}
    newatoms=deepcopy(atoms)
    spglib_info=spglib.get_symmetry_dataset(atoms, symprec=symprec)
    if (rmin=='auto'):
        Z=atomic_numbers[interstitial_atom]
        rmin=covalent_radii[Z]*0.7
        print "rmin=",rmin
    N=1
    while (N<=Nmax):
        newpos=Identify_Interstitial_Position(newatoms,rmin, nscan=nscan, symprec=symprec,interstitial_atom=interstitial_atom)
        N=N+1
        if newpos==[]:
            break
        # create new atoms object
        new_symbols=newatoms.get_chemical_symbols()
        new_scaled_positions=[]
        for x in newatoms.get_scaled_positions():
            new_scaled_positions.append(x)
        multiplicity=0
        for rot,trans in zip(spglib_info['rotations'],spglib_info['translations']):
            gpos=(np.dot(rot,newpos)+trans)%1.
            is_new=True
            for pos in new_scaled_positions:
                if np.sum(np.fabs(gpos-pos))<0.0001:
                    is_new=False
                    break
            if (is_new==True):
                new_symbols.append(interstitial_atom)
                new_scaled_positions.append(gpos)
                multiplicity=multiplicity+1
        newatoms=Atoms(symbols=new_symbols,scaled_positions=new_scaled_positions,cell=atoms.get_cell(),pbc=True)
        spglib_new=spglib.get_symmetry_dataset(newatoms, symprec=symprec)
        label="%s_%d%s"%(interstitial_atom,multiplicity,spglib_new['wyckoffs'][-1])
        interst_pos[label]=newpos
    #create atoms objects with one interstitial at the positions identified
    aolist={}
    for label in interst_pos:
        new_symbols=atoms.get_chemical_symbols()
        new_scaled_positions=[]
        for x in atoms.get_scaled_positions():
            new_scaled_positions.append(x)
        new_symbols.append(interstitial_atom)
        new_scaled_positions.append(interst_pos[label])
        newatoms=Atoms(symbols=new_symbols,scaled_positions=new_scaled_positions,cell=atoms.get_cell(),pbc=True)
        #spglib_new=spglib.get_symmetry_dataset(newatoms, symprec=symprec)
        #newatoms.info={}
        #newatoms.info['spacegroup'].no=spglib_new['number']
        #newatoms.info['spacegroup'].symbol=spglib_new['international']
        aolist[label]=newatoms
    return aolist
                                    
def Identify_Interstitial_Position(atoms, rmin, nscan=24, radii='auto',N_radii=0, symprec=1e-3,interstitial_atom='C', scale_radii=0.7,use_sym=True, prefer_high_sym=False):
    """ Find position with largest empty sphere
    """
    if use_sym==True:
        is_done=np.zeros((nscan,nscan,nscan),int)
        spglib_info=spglib.get_symmetry_dataset(atoms, symprec=symprec)
    scaled_pos=[]
    scaled_pos_mult=1000
    if radii=='auto':
        radii=[]
        ir=0
        for Z in atoms.get_atomic_numbers():
            ir=ir+1
            if (N_radii>0) and (ir>N_radii):
                radii.append(0.0)
            else:
                radii.append(covalent_radii[Z]*scale_radii)
    elif radii=='None':
        radii=len(atoms)*[0.0]
    print radii
    #create supercell
    cell=atoms.get_cell()
    N=[1,1,1]
    #for i in range(3):
    #    d=0.0
    #    for j in range(3):
    #        d=d+cell[i][j]*cell[i][j]
    #    d=np.sqrt(d)
    #    N.append(int(cutoff/d)+1)
    #subdivide original cell and scan for neighbor distances
    dmax=0.
    for ix in range(nscan):
        for iy in range(nscan):
            for iz in range(nscan):
                if (use_sym==True) and (is_done[ix][iy][iz]!=0):
                    break
                mult=1
                ri=cell[0]*ix/nscan+cell[1]*iy/nscan+cell[2]*iz/nscan
                if (use_sym):
                    sci=np.array([1.*ix/nscan,1.*iy/nscan,1.*iz/nscan])
                    genpos=[[ix,iy,iz]]
                    for rot,trans in zip(spglib_info['rotations'],spglib_info['translations']):
                        gpos=(np.dot(rot,sci)+trans)%1.
                        gind=[ix,iy,iz]
                        for j in range(3):
                            gind[j]=int(gpos[j]*nscan)%nscan
                        #print "XXX",sci,gpos,gind
                        if not (gind in genpos):
                            mult=mult+1
                            genpos.append(gind)
                        is_done[gind[0]][gind[1]][gind[2]]=1
                        #print ix,iy,iz,gind
                is_interst=True
                if (prefer_high_sym==True) and (mult>scaled_pos_mult):
                    is_interst=False
                #print mult
                dmin=1000.
                for iatom in range(len(atoms)):
                    if (is_interst==False):
                        break
                    for n1 in range(-N[0], N[0] + 1):
                        if (is_interst==False):
                            break
                        for n2 in range(-N[1], N[1] + 1):
                            if (is_interst==False):
                                break
                            for n3 in range(-N[2], N[2] + 1):
                                Rat=np.dot((n1, n2, n3), cell)+atoms.get_positions()[iatom]
                                d=0.0
                                for j in range(3):
                                    d=d+(Rat[j]-ri[j])*(Rat[j]-ri[j])
                                d=np.sqrt(d)-radii[iatom]
                                if d<rmin:
                                    is_interst=False
                                    break
                                elif (d<dmin):
                                    dmin=d
                if (is_interst==True):
                    if (dmin>dmax) or ((prefer_high_sym==True) and (mult<scaled_pos_mult)):
                        dmax=dmin
                        rmax=ri
                        scaled_pos=np.array([ix/float(nscan),iy/float(nscan),iz/float(nscan)])
                        scaled_pos_mult=mult
    if (scaled_pos!=[]):
        print "interstitial at ",scaled_pos," dist=",dmax
    return scaled_pos


def Create_Interstitial_Supercell(ao, interstial_atoms=[], interstitial_blocks=[[0,0,0]], N=[2,2,2], silent=False):
    """Create a N_1xN_2xN_3 supercell with interstitials in interstitial_blocks"""
    sc_cell=ao.get_cell()
    for i in range(3):
        sc_cell[i]=ao.get_cell()[i]*N[i]
    pos=ao.get_scaled_positions()
    elements=ao.get_chemical_symbols()
    sc_positions=[]
    sc_elements=[]
    for N_0 in range(N[0]):
        for N_1 in range(N[1]):
            for N_2 in range(N[2]):
                N_i=[N_0,N_1,N_2]
                for isite in range(len(ao.get_scaled_positions())):
                    if (not (isite in interstial_atoms)) or ([N_0,N_1,N_2] in interstitial_blocks):
                        sc_elements.append(elements[isite])
                        npos=[]
                        for i in range(3):
                            npos.append(1.0*(pos[isite][i]+N_i[i])/N[i])
                        sc_positions.append(npos)
    supercell=Atoms(symbols=sc_elements,cell=sc_cell,scaled_positions=sc_positions,pbc=True)
    if (silent==False):
        print "Create_Interstitial_Supercell: original cell:"
        print ao.get_cell()
        for el,pos in zip(ao.get_chemical_symbols(),ao.get_scaled_positions()):
            print el,pos
        print "Create_Interstitial_Supercell: supercell:"
        print supercell.get_cell()
        for el,pos in zip(supercell.get_chemical_symbols(),supercell.get_scaled_positions()):
            print el,pos
    return supercell



def Create_Magnetic_Supercell(prop_dict, interstial_atoms=[], distortion=[], N=[2,2,2], silent=False):
    """Create supercell with interstitial atoms and same magnetic structure as parent
    prop_dict: dictionary with properties of original cell ('cell','scaled_positions' and (optionally) 'magnetic_moments')
    """
    pd={}
    pd_sc={}
    for prop in ['cell','chemical_symbols','scaled_positions','magnetic_moments']:
        if prop in prop_dict:
            pd[prop]=prop_dict[prop]
        elif (prop=='magnetic_moments'):
            pd['magnetic_moments']=np.zeros(len(pd['chemical_symbols']))
        else:
            return pd_sc
    chem_symb=[]
    magmoms=[]
    sc_pos=[]
    for (el,pos,mom) in zip(pd['chemical_symbols'],pd['scaled_positions'],pd['magnetic_moments']):
        for N_0 in range(N[0]):
            for N_1 in range(N[1]):
                for N_2 in range(N[2]):
                    N_i=[N_0,N_1,N_2]
                    chem_symb.append(el)
                    magmoms.append(mom)
                    npos=[]
                    for i in range(3):
                        npos.append(1.0*(pos[i]+N_i[i])/N[i])
                    sc_pos.append(npos)
    for el,pos,ib,mom in interstial_atoms:
        chem_symb.append(el)
        npos=[]
        for i in range(3):
            npos.append(1.0*(pos[i]+ib[i])/N[i])
        sc_pos.append(npos)
        if mom!=None:
            magmoms.append(mom)
        elif isinstance(pd['magnetic_moments'][0],list):
            magmoms.append([0,0,0])
        else:
            magmoms.append(0.0)
    new_cell=np.eye(3)
    for i in range(3):
        new_cell[i]=np.array(pd['cell'][i])*N[i]
    if distortion!=[]:
        new_cell=np.dot(distortion,new_cell)
    pd_sc['cell']=new_cell
    pd_sc['chemical_symbols']=chem_symb
    pd_sc['scaled_positions']=sc_pos
    pd_sc['magnetic_moments']=magmoms
    return pd_sc
