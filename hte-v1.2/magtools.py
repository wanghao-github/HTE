from hte import *
from ase.lattice.spacegroup.cell import *
import re
import scipy
from scipy import linalg, matrix

try:
    import spglib
    has_spglib=True
except:
    has_spglib=False

class MSG(object):
    def __init__(self,mcif=None,symelem=None,eps=0.01):
        self.eps=eps
        #print "MSG_init(): eps=",self.eps
        self.elements=[]
        self.symbols=[]
        self.unity=None
        self.name=None
        self.multiplication_table=None
        if mcif!=None:
            self.read_mcif(mcif)
        if symelem!=None:
            for g in symelem:
                #print g
                sym=self.matrix2symbol(g)
                #print sym
                self.symbols.append(sym)
                gnew=self.symbol2matrix("1 "+sym)
                #print gnew
                if self.is_equal_element(g,gnew,ignore_lattice_translations=True):
                    self.elements.append(gnew)
                else:
                    print "MSG __init__(): Failed conversion: ",g,sym,gnew
                    bla #normalize translations first
            #self.elements=symelem
            
    def get_elements(self):
        return self.elements
    
    def get_symbols(self):
        return self.symbols

    def get_order(self):
        return len(self.elements)
    
    def get_multiplication_table(self):
        if isinstance(self.multiplication_table, type(None))==False:
            return self.multiplication_table
        order=self.get_order()
        multtab=np.zeros((order,order),dtype=int)
        for i in range(order):
            g=self.elements[i]
            for j in range(order):
                h=self.elements[j]
                gh=self.mult(g,h)
                for k in range(order):
                    if self.is_equal_element(gh,self.elements[k])==True:
                        multtab[i,j]=k
                        break
        self.multiplication_table=multtab
        return self.multiplication_table
    
    
    def fraction_to_decimal(fraction):
        sign = -1 if fraction.startswith('-') else 1
        numerator, denominator = map(int, fraction.lstrip('+-').split('/'))
        return sign * numerator / denominator
    
    
    def symbol2matrix(self,symb):
        #syntax from Bilbao: '1 x,y,z,+1 mx,my,mz'
        #'11 -y+1/2,x+1/2,-z+3/4,+1 my,-mx,mz'
        # try to use integers to avoid rounding errors
	#print symb
        sign={'+':1.0,'-':-1.0}
        m={'x':np.array([1,0,0]),'y':np.array([0,1,0]),'z':np.array([0,0,1]),'-x':np.array([-1,0,0]),'-y':np.array([0,-1,0]),'-z':np.array([0,0,-1])}
        symb_split=symb.split(' ')[1].split(',')
        if len(symb_split)==4:
            [x,y,z,t]=symb_split
        elif len(symb_split)==3:
            [x,y,z]=symb_split
        rot=np.zeros((3,3),dtype=int)
        trans=np.zeros(3)
        # new version, can handle more expressions
        for i in range(len([x,y,z])):
            sym=[x,y,z][i]+"+"
            prefak=""
            for a in sym:
                if a in ['+','-']:
                    # evaluate xpr
                    if prefak.strip()!="":
                        #print prefak
                        if '/' in prefak:
                            tss=prefak.split('/')
                            trans[i]=trans[i]+float(tss[0])/float(tss[1])
                        else:
                            trans[i]=trans[i]+float(prefak)
                    prefak=a
                elif a in ['x','y','z']:
                    if prefak=="-":
                        rot[i]=rot[i]-m[a]
                    elif (prefak=="+") or (prefak==""):
                        rot[i]=rot[i]+m[a]
                    else:
                        rot[i]=rot[i]+int(prefak)*m[a]
                    prefak=""
                else:
                    prefak=prefak+a
        #print symb,rot,trans
        if len(symb_split)==3:
            return rot,trans
        return rot,trans,int(t)
        # old version
        for i in range(len([x,y,z])):
            sym=[x,y,z][i]
            sgn=1
            for a in sym:
                if a in ['+','-']:
                    sgn=sign[a]
                elif a in ['x','y','z']:
                    rot[i]=rot[i]+sgn*m[a]
                elif a.isdigit() or a=='.':
                    ts=[x,y,z][i].split('+')[-1]
                    if '/' in ts[1]:
                        tss=ts.split('/')
                        trans[i]=float(tss[0])/float(tss[1])
                    else:
                        trans[i]=float(ts)
        if len(symb_split)==3:
            return rot,trans
        return rot,trans,int(t)
        #old
        rot=np.array([m[x.split('+')[0]],m[y.split('+')[0]],m[z.split('+')[0]]])
        trans=np.array([0.,0.,0.])
        for i in range(len([x,y,z])):
            ts=[x,y,z][i].split('+')
            if len(ts)==1:
                trans[i]=0
            elif '/' in ts[1]:
                tss=ts[1].split('/')
                trans[i]=float(tss[0])/float(tss[1])
            else:
                trans[i]=float(ts[1])
        return rot,trans,int(t)

    def matrix2symbol(self,mat,print_mag=True, eps=None):
        if eps==None:
            eps=self.eps
        # use integer numbers if possible to avoid rounding errors in multiplications
        m=['x','y','z']
        tstr={'1/2':0.5,'1/3':1./3.,'2/3':2./3.,'1/6':1./6.,'5/6':5./6.,'1/4':0.25,'3/4':0.75,'1/8':0.125,'3/8':0.375,'5/8':0.625,'7/8':0.875}
        if len(mat)==3:
            rot,trans,t=mat
        else:
            rot,trans=mat
            t=''
        symb=''
        for i in range(3):
            sym=''
            for j in range(3):
                if (abs(rot[i][j])<eps):
                    continue
                else:
                    if (rot[i][j]>0):
                        if sym!='':
                            sym=sym+'+'
                    if (abs(rot[i][j]-1)<eps):
                        sym=sym+m[j]
                    elif (abs(rot[i][j]+1)<eps):
                        sym=sym+'-'+m[j]
                    else:
                        symx=str(rot[i][j])
                        foundint=False
                        for ival in range(-4,5):
                            if (abs(rot[i][j]-ival)<eps):
                                symx=str(ival)
                                foundint=True
                                break
                        sym=sym+symx+m[j]
                        if (foundint==False):
                            print "matrix2symbol(): No interger representation found for rotation matrix ",rot, self.name
                #
                #if (rot[i][j]>0):
                #    if sym!='':
                #        sym=sym+'+'
                #    sym=sym+m[j]
                #elif (rot[i][j]<0):
                #    sym=sym+'-'+m[j]
            if abs(trans[i])>eps: #!=0:
                done=False
                for x in tstr:
                    if fabs(tstr[x]-trans[i])<eps:
                        sym=sym+'+'+x
                        done=True
                        break
                for x in range(-5,5):
                    if (abs(trans[i]-x)<eps):
                        done=True
                        break                        
                if not done:
                    print "matrix2symbol(): No interger representation found for translation ",trans, self.name
                    #sym=sym+'+'+str(trans[j])
            if symb!='':
                symb=symb+','
            symb=symb+sym
        if len(mat)==3:
            symb=symb+',%+d'%t
            symb_m=' '
            RM=np.linalg.det(rot)*t*rot
            for i in range(3):
                sym=''
                for j in range(3):
                    if (RM[i][j]>0):
                        if sym!='':
                            sym=sym+'+'
                        sym=sym+'m'+m[j]
                    elif (RM[i][j]<0):
                        sym=sym+'-'+'m'+m[j]
                if symb_m!=' ':
                    symb_m=symb_m+','
                symb_m=symb_m+sym
            if print_mag==True:
                symb=symb+symb_m
        return symb

    def get_null_space(self,A, eps=1e-12):
        u, s, vh = scipy.linalg.svd(A)
        padding = max(0,np.shape(A)[1]-np.shape(s)[0])
        null_mask = np.concatenate(((s <= eps), np.ones((padding,),dtype=bool)),axis=0)
        null_space = scipy.compress(null_mask, vh, axis=0)
        return null_space
    
    def read_mcif(self,fname):
        #for the moment read in x,y,z,t operations, expand to mcif reader later
        fin=open(fname,'r')
        line=fin.readline()
        while not (line.startswith("1 ")):
            #print "Skipping",line
            #tmp solution to read group name
            if line.startswith("# "):
                self.name=line.split("# ")[1].strip("\n")
            line=fin.readline()
        while line.strip():
            #print line
            self.symbols.append(line.strip('\n'))
            self.elements.append((self.symbol2matrix(line)))
            line=fin.readline()
        fin.close()
        return self.elements

    def read_full_mcif(self,fname):
        #
        ciftags={}
        in_loop=False
        fin=open(fname,'r')
        line=fin.readline()
        print "check_point184, line is: ", line
        while line:
            # print line
            line=line.strip()
            print line
            lspl=line.split()
            if (line=="") or (line=="\0"):
                in_loop=False
                loop_tags=[]
            if (in_loop==True):
                if line.startswith("_"):
                    loop_tags.append(line)
                elif (len(lspl)==len(loop_tags)):
                    for x,y in zip(loop_tags,lspl):
                      if x in ciftags:
                          ciftags[x].append(y)
                      else:
                          ciftags[x]=[y]
                #print "INLOOP",zip(loop_tags,lspl)
            elif line.startswith("data_"):
                uid=line.strip()
            elif line.startswith("loop_"):
                in_loop=True
                loop_tags=[]
                #print "LOOP START"
            elif line.startswith("_"):
                ciftags[lspl[0]]=lspl[-1]
            line=fin.readline()
        #print uid,":", ciftags
        pd={}
        #for the moment: use ASE to construct cell
        cellpar=[]
        for x in ['_cell_length_a','_cell_length_b','_cell_length_c','_cell_angle_alpha','_cell_angle_beta','_cell_angle_gamma']:
            # cellpar.append(float(ciftags[x]))
            cellpar.append(float(re.sub(r'\(.*?\)', '', ciftags[x])))
            print " ciftags[x] is ", ciftags[x]
            print "cellpar is ", cellpar
        pd['cell']=cellpar_to_cell(cellpar)
        print "check_point185, pd['cell'] is: ", pd['cell']
        msg=MSG()
        print "check_point186, msg=MSG() is: ", msg
        
        # Assuming `ciftags` is a dictionary containing the CIF data
        if '_space_group_symop_magn_operation.id' in ciftags and '_space_group_symop_magn_operation.xyz' in ciftags:
            id_field = '_space_group_symop_magn_operation.id'
            xyz_field = '_space_group_symop_magn_operation.xyz'
        elif '_space_group_symop.magn_id' in ciftags and '_space_group_symop.magn_operation_xyz' in ciftags:
            id_field = '_space_group_symop.magn_id'
            xyz_field = '_space_group_symop.magn_operation_xyz'
        else:
            raise ValueError("Required CIF fields for magnetic symmetry operations not found.")

        
        # init_op = []
        # Now use the determined field names in your loop
        for i, g in zip(ciftags[id_field], ciftags[xyz_field]):
            # symb = i + " " + g        
            # print "check_point204, self.symbol2matrix(symb) is ", self.symbol2matrix(symb)
            # init_op.append(self.symbol2matrix(symb))
            # # msg.symbols.append(symb)
            # # msg.elements.append(self.symbol2matrix(symb))
            
            
            symb=i+" "+g
            msg.symbols.append(symb)
            msg.elements.append(self.symbol2matrix(symb))
        # print "check_point206 init_op is ", init_op

        # for i,g in zip(ciftags['_space_group_symop_magn_operation.id'],ciftags['_space_group_symop_magn_operation.xyz']):
        #     symb=i+" "+g
        #     msg.symbols.append(symb)
        #     msg.elements.append(self.symbol2matrix(symb))
            
        # print "check_point187, msg.symbols. is: ", msg.symbols
        # print "check_point188, msg.elements is: ", msg.elements
        
        # if '_space_group_symop_magn_centering.xyz' in ciftags:
        #     for i,g in  zip(ciftags['_space_group_symop_magn_centering.id'],ciftags['_space_group_symop_magn_centering.xyz']):
        #         symb=i+" "+g
        #         if (not (symb in msg.symbols)):
        #             msg.symbols.append(symb)
        #             msg.elements.append(self.symbol2matrix(symb))
        #             #print "added _space_group_symop_magn_centering.xyz:",symb
        # msg.is_group(complete=True)
        # print "is group:",msg.is_group()

        if '_space_group_symop_magn_centering.xyz' in ciftags and '_space_group_symop_magn_centering.id' in ciftags:
            id_center_field = '_space_group_symop_magn_centering.id'
            xyz_center_field = '_space_group_symop_magn_centering.xyz'
        elif '_space_group_symop.magn_centering_xyz' in ciftags and '_space_group_symop.magn_centering_id' in ciftags:
            id_center_field = '_space_group_symop.magn_centering_id'
            xyz_center_field = '_space_group_symop.magn_centering_xyz'
        else:
            raise ValueError("Required CIF fields for magnetic centering operations not found.")

        
        print "check_point212, ciftags[id_center_field] is", ciftags[id_center_field]
        # # Process the data using the determined field names
        for i, g in zip(ciftags[id_center_field], ciftags[xyz_center_field]):
            symb=i+" "+g
            if (not (symb in msg.symbols)):
                msg.symbols.append(symb)
                msg.elements.append(self.symbol2matrix(symb))
        #     symb = i + " " + g
        #     print "check_point205, g is" , g
        #     # if not symb in msg.symbols:
        #     #     msg.symbols.append(symb)
        #     #     msg.elements.append(self.symbol2matrix(symb))
        #     # for j in range(len(init_op)):
        #     #     init_trans = init_op[j][1]
        #     #     print "check_point208, init_trans is: ",init_trans
        #     if len(self.symbol2matrix(symb)) == 2:

        #         lat_rot, lat_trans = self.symbol2matrix(symb)
        #         print "check_point206, lat_rot and lat_trans is, ", lat_rot,lat_trans
                    
        #     elif len(self.symbol2matrix(symb)) == 3:          
        #         lat_rot, lat_trans, timeinv = self.symbol2matrix(symb)
        #         print "check_point207, lat_rot, lat_trans and time inv is, ", lat_rot,lat_trans,timeinv
                
        #     else:
        #         print "check_point213, stoped here "
        #     # print "Before copying:", type(init_op)
        #     # new_ops = init_op.copy()
        #     # print "After copying, new_ops:", new_ops
        #     # new_ops = init_op.copy()
        #     new_ops = list(init_op)
        #     print "check_point211, new_ops is ", new_ops
        #     print "check_point210, ciftags[id_center_field] is ",ciftags[id_center_field]
            
        #     for rotation, translation, timeinv in init_op:
        #         new_trans = translation + lat_trans
        #         new_ops.append((rotation, new_trans, timeinv))
                    
        #     print "check_point209, new_ops is ", new_ops
        #     # for i, g in zip(ciftags[id_field], ciftags[xyz_field]):
        #     #     symb = i + " " + g        
        #     #     print "check_point204, self.symbol2matrix(symb) is ", self.symbol2matrix(symb)
        #     #     init_op.append(self.symbol2matrix(symb))
        #     # msg.symbols.append(symb)
        # if len(ciftags[id_center_field]) == 1:
        #     msg.elements = init_op
        # elif len(ciftags[id_center_field]) >= 2:
        #     msg.elements = new_ops
        
        # print "check_point214, msg.elements is " ,   msg.elements     
                    
        # # Uncomment the following line to see logs of added symbols
        #         print "check_point189, symb added:", symb

        # After processing, check if the group is complete
        # complete_group, = msg.is_group(complete=True)
        # print "is group:", msg.is_group()
        
        pd['chemical_symbols']=[]
        pd['scaled_positions']=[]
        pd['initial_magnetic_moments']=[]
        
        
        pd_temp_no_center = {}
            
        pd_temp_no_center['chemical_symbols']=[]
        pd_temp_no_center['scaled_positions']=[]
        pd_temp_no_center['initial_magnetic_moments']=[]
        
        for (el, pos, lab) in zip(ciftags['_atom_site_type_symbol'],
                          zip(ciftags['_atom_site_fract_x'], ciftags['_atom_site_fract_y'], ciftags['_atom_site_fract_z']),
                          ciftags['_atom_site_label']):
            mom = np.zeros(3)  # Initialize magnetic moment as zero vector

            # Determine the correct field names
            label_field = '_atom_site_moment.label' if '_atom_site_moment.label' in ciftags else '_atom_site_moment_label'
            x_field = '_atom_site_moment.crystalaxis_x' if '_atom_site_moment.crystalaxis_x' in ciftags else '_atom_site_moment_crystalaxis_x'
            y_field = '_atom_site_moment.crystalaxis_y' if '_atom_site_moment.crystalaxis_y' in ciftags else '_atom_site_moment_crystalaxis_y'
            z_field = '_atom_site_moment.crystalaxis_z' if '_atom_site_moment.crystalaxis_z' in ciftags else '_atom_site_moment_crystalaxis_z'
    
            if lab in ciftags[label_field]:
                label_index = ciftags[label_field].index(lab)
                # mom[0] = float(ciftags[x_field][label_index])
                # mom[1] = float(ciftags[y_field][label_index])
                # mom[2] = float(ciftags[z_field][label_index])
                print "ciftags[x_field][label_index] is ",ciftags[x_field][label_index] 
                mom[0] = float(re.sub(r'\(.*?\)', '', ciftags[x_field][label_index]))
                mom[1] = float(re.sub(r'\(.*?\)', '', ciftags[y_field][label_index]))
                mom[2] = float(re.sub(r'\(.*?\)', '', ciftags[z_field][label_index]))
                print "mom is: ", mom               
                        
            for g in msg.get_elements()[:len(ciftags[id_field])]:
                # print "g is, ", g
                # print "pos is, ", pos
                               
                clear_pos = [
                re.sub(r'\(.*?\)', '', pos[0]),
                re.sub(r'\(.*?\)', '', pos[1]),
                re.sub(r'\(.*?\)', '', pos[2])
                ]
                
                print "clear_pos is ",clear_pos
                # npos = msg.symop_pos(g, (float(re.sub(r'\(.*?\)', '', pos[0]))), (float(re.sub(r'\(.*?\)', '', pos[1]))), (float(re.sub(r'\(.*?\)', '', pos[2]))))
                npos = msg.symop_pos(g, (float(clear_pos[0]), float(clear_pos[1]), float(clear_pos[2])))
                print "npos is ",npos
                nmom = msg.symop_mag(g, mom)
                print "mom is ",mom
                is_new = True
                
                for px in pd['scaled_positions']:
                    if msg.is_equal_site(npos, px):
                        print "msg.is_equal_site "
                        is_new = False
                        break
                if is_new:
                    pd['chemical_symbols'].append(el)
                    pd['scaled_positions'].append(npos)
                    pd['initial_magnetic_moments'].append(nmom)         
                    pd_temp_no_center['chemical_symbols'].append(el)
                    pd_temp_no_center['scaled_positions'].append(npos)
                    pd_temp_no_center['initial_magnetic_moments'].append(nmom)
        
        
        print "check_point221, no problem here pdpd_temp_no_center is ",pd_temp_no_center    
        # pd = pd_temp_no_center
            
            
        for j in range(len(pd_temp_no_center['scaled_positions'])):
            for g in msg.get_elements()[-(len(msg.get_elements()) - len(ciftags[id_field])):]:
                print "check_point220 ,-(len(msg.get_elements()) - len(ciftags[id_field])): is ",g
                npos2 = msg.symop_pos(g, (float(pd_temp_no_center['scaled_positions'][j][0]), float(pd_temp_no_center['scaled_positions'][j][1]), 
                                          float(pd_temp_no_center['scaled_positions'][j][2])))
                nmom2 = msg.symop_mag(g, pd_temp_no_center['initial_magnetic_moments'][j])
                
                print "npos2 is ",npos2
                print "nmom2 is ",nmom2
                
                is_new = True
                
                for px in pd['scaled_positions']:
                    if msg.is_equal_site(npos2, px):
                        print "msg.is_equal_site "
                        is_new = False
                        break            
                if is_new:
                    pd['chemical_symbols'].append(pd_temp_no_center['chemical_symbols'][j])
                    pd['scaled_positions'].append(npos2)
                    pd['initial_magnetic_moments'].append(nmom2)

        #             pd_temp_no_center['chemical_symbols'].append(el)
        #             pd_temp_no_center['scaled_positions'].append(npos)
        #             pd_temp_no_center['initial_magnetic_moments'].append(nmom)
                    
                    
        pd['initial_magnetic_moments'] = np.array(pd['initial_magnetic_moments'])
        # print "check_point219 ,initial_pos is ",initial_pos        
        print "uid,pd,ciftags are", uid,pd,ciftags
        
        return uid,pd,ciftags

    def write_mcif(self,fname,print_mag=True):
        #for the moment write only x,y,z,t operations, expand to mcif later
        fout=open(fname,'w')
        if self.name!=None:
            fout.write("# %s\n"%self.name)
        i=0
        for g in self.get_elements():
            i=i+1
            sym=self.matrix2symbol(g,print_mag=print_mag)
            fout.write("%d %s\n"%(i,sym))
        fout.close()
        return i

    def symop_pos(self,g,pos):
        if len(g)==2:
            (rot,trans)=g
        elif len(g)==3:
            (rot,trans,t)=g
        else:
            return None
        # print "in symop_pos, g is,", g
        newpos=(np.dot(rot,pos)+trans)%1.
        return newpos

    def is_equal_site(self,a,b,eps=-1):
        """returns True if (a-b)_i in [0,-1,+1] for all i
        """
        if eps==-1:
            eps=self.eps
        for i in range(len(a)):
            is_equal=False
            for d in [0.0,1.0,-1.0]:
                if abs(a[i]-b[i]-d)<eps:
                    is_equal=True
                    break
            if is_equal==False:
                break
        return is_equal
    
            
    def symop_mag(self,g,mag):
        if len(g)==3:
            (rot,trans,t)=g
        else:
            return None
        newmag=np.linalg.det(rot)*t*(np.dot(rot,mag))
        return newmag

    def gauss_elim(self,M):
        eps=0.001
        irank=0
        for j in range(len(M[0])):
            #swap if ith element is zero
            nonzero=False
            for i in range(irank,len(M)):
                if fabs(M[i][j])>eps:
                    nonzero=True
                    break
            if nonzero==False:
                continue
            if irank!=i:
                x=copy(M[irank])
                M[irank]=M[i]/M[i][j]
                M[i]=copy(x)
            for i in range(irank+1,len(M)):
                M[i]=M[i]-M[i][j]/M[irank][j]*M[irank]
            irank=irank+1
        G=[]
        for i in range(irank):
            G.append(M[i])
        return np.array(G)
    
    def get_magnetic_configurations(self,ao_in,name='xx',reduce_collinear=True,exclude_ferromagnetic=True, max_configs=10000, silent=False, AFatoms=['Mn'], forceMAG=[], return_status=False):
        """get magnetic configurations for atoms object ao which are allowed by the symmetry of the present magnetic subgroup
        AFatoms: list of chemical elements which are expected to couple AF; alternatively ['Mn_0','Mn_2'] can be used to
                 couple atoms on selected sites
        forceMAG: reject solutions where atoms in list have no magnetic moment
        """
        status={}
        #print AFatoms
        #if silent==False:
        #    print "get_magnetic_configurations:", name
        magconfigs={}
        ao=ao_in
        if AFatoms!=[]:
            elsAF=[]
            posAF=[]
            mappingAF=[]
            # create tmp ao with AF atoms only
            for i in range(len(ao.get_chemical_symbols())):
                if (ao.get_chemical_symbols()[i] in AFatoms):
                    elsAF.append(ao.get_chemical_symbols()[i])
                    posAF.append(ao.get_scaled_positions()[i])
                    mappingAF.append(i)
                elif ((ao.get_chemical_symbols()[i]+"%d") in AFatoms):
                    elsAF.append(ao.get_chemical_symbols()[i])
                    posAF.append(ao.get_scaled_positions()[i])
                    mappingAF.append(i)
            print "check_point194, before aoAF, elsAF and posAF are ",elsAF,posAF
            aoAF=Atoms(symbols=elsAF,scaled_positions=posAF,cell=ao.get_cell(),pbc=True)
            print "check_point195, aoAF is",aoAF
            symmetry_constraints_r=self.get_symmetry_equivalent_positions(aoAF)
            symmetry_constraints={'mag_constraints':[]}
            for sol_r in symmetry_constraints_r['mag_constraints']:
                #print "sol_r",sol_r
                sol=[]
                for i in range(len(ao.get_chemical_symbols())):
                    sol.append(np.zeros(3))
                for i in range(len(sol_r)):
                    sol[mappingAF[i]]=sol_r[i]
                #print sol
                symmetry_constraints['mag_constraints'].append(np.array(sol))
        else:
            symmetry_constraints=self.get_symmetry_equivalent_positions(ao)
        #print "Number of independent solutions:",len(symmetry_constraints['mag_constraints']) #,symmetry_constraints['mag_constraints']
        status['dimension']=len(symmetry_constraints['mag_constraints'])
        if len(symmetry_constraints['mag_constraints'])==0:
            if return_status==True:
                return magconfigs,status
            return magconfigs
        #build +- linear combinations of magnetic sublattices compatible with msg (only for magnetic elements)
        updn=[[1]]
        #if silent==False:
        #    print len(symmetry_constraints['mag_constraints'])
        #print symmetry_constraints #['mag_constraints']
        updnchoice= [-1,1]
        for i in range(1,len(symmetry_constraints['mag_constraints'])):
            updnn=[]
            for config in updn:
                for spin in updnchoice: #[-1,1]:
                    confign=deepcopy(config)
                    confign.append(spin)
                    updnn.append(confign)
            updn=updnn
            if len(updn)>max_configs:
               updnchoice= [-1]
        ini_mags=[]
        for config in updn:
            add=False
            magmom=0*symmetry_constraints['mag_constraints'][0]
            for i in range(len(config)):
                magmom=magmom+config[i]*symmetry_constraints['mag_constraints'][i]
                ini_mag=[]
                for el,vec in zip(ao.get_chemical_symbols(),magmom):
                    if (ismagnetic(el)) and (np.sum(np.fabs(vec))>0.001):
                        muvec=np.dot(vec,ao.get_cell())
                        #check!!
                        muvec=vec
                        muvec=muvec/norm(muvec)
                        mom=get_ini_magmom(el)*muvec
                        add=True
                    else:
                        mom=np.zeros(3)
                    ini_mag.append(mom)
            if add==True:
                ini_mags.append(np.array(ini_mag))
            #    print "added:",ini_mag
            #else:
            #    print "not added:",ini_mag
            max_configs=10000
            if len(ini_mags)>=max_configs:
                print "*** More configs possible, stopping here..."
                break 
        if len(ini_mags)==0:
            if return_status==True:
                return magconfigs,status
            return magconfigs
        mid=0
        for magmom in ini_mags:
            noncol=True
            if reduce_collinear:
                noncol=False
                momcol=[]
                momref=[]
                for mom1 in magmom:
                    if (norm(mom1)>0.01) and (momref==[]):
                        momref=mom1
                    if (momref!=[]) and (abs(np.dot(mom1,momref))>0.01):
                        momcol.append(norm(mom1)*np.dot(mom1,momref)/abs(np.dot(mom1,momref)))
                    else:
                        momcol.append(0.)
                    for mom2 in magmom:
                        if norm(np.cross(mom1,mom2))>0.01:
                            noncol=True
                            break
                if noncol==False:
                    magmom=np.array(momcol)
            add=True
            for conf in magconfigs:
                if (noncol==magconfigs[conf]['lnoncollinear']) and (norm(magmom-magconfigs[conf]['magmom'])<0.01): #CHECK
                    add=False
            if (exclude_ferromagnetic==True) and (noncol==False):
                ferro=True
                for mom in magmom:
                    if mom<0.0:
                        ferro=False
                        break
                if ferro==True:
                    add=False
            if forceMAG!=[]:
                for i in range(len(magmom)):
                    el=ao.get_chemical_symbols()[i]
                    if ((el in forceMAG) or ((el+"_%d"%(i+1)) in forceMAG)) and (norm(magmom[i])<0.05):
                        add=False
                        print "rejected:","forceMAG=",forceMAG,zip(ao.get_chemical_symbols(),magmom)
                    #print el+"_%d"%(i+1), norm(magmom[i]),magmom[i], (el+"%d"%(i+1)) in forceMAG, norm(magmom[i])<0.05
                    #if ((el+"_%d"%(i+1)) in forceMAG) and (norm(magmom[i])<0.05):
                    #    add=False
                    #    print "CCC rejected",zip(ao.get_chemical_symbols(),magmom)
            if add:
                i=0
                nameconf=name
                nameconfsh=name
                for el,mom in zip(ao.get_chemical_symbols(),magmom):
                    i=i+1
                    if (noncol==True) and (norm(mom)>0.05):
                        nameconf=nameconf+"%s%d_%.1f_%.1f_%.1f"%(el,i,mom[0],mom[1],mom[2])
                        nameconfsh=nameconfsh+"%d_%.0fx%.0fx%.0f"%(i,mom[0],mom[1],mom[2])
                    elif (noncol==False) and (abs(mom)>0.05):
                        nameconf=nameconf+"%s%d_%.1f"%(el,i,mom)
                        nameconfsh=nameconfsh+"%s%d_%.0f"%(el,i,mom)
                if len(nameconf)>250:
                    nameconf=nameconfsh
                magconfigs[nameconf]={'magmom':magmom,'atoms_object':ao,'lnoncollinear':noncol}
                mid=mid+1
                if mid>max_configs:
                    print "WARNING(get_magnetic_configurations): too many spin configurations (%d), increase max_configs=%d"%(mid,max_configs)
                    status['failed']=True
                    if return_status==True:
                        return magconfigs,status
                    return {}
        if return_status==True:
            status['success']=True
            return magconfigs,status
        return magconfigs
    
    def get_magnetic_configurations_old(self,ao,name='xx',reduce_collinear=True):
        magconfigs={}
        symmetry_constraints=self.get_symmetry_equivalent_positions(ao)
        #build +- linear combinations of magnetic sublattices compatible with msg (only for magnetic elements)
        magvecs=[]
        for mvec in symmetry_constraints['mag_constraints']:
            ini_mag=[]
            add=False
            for el,vec in zip(ao.get_chemical_symbols(),mvec):
                if (ismagnetic(el)) and (np.sum(np.fabs(vec))>0.001):
                    muvec=np.dot(vec,ao.get_cell())
                    muvec=muvec/norm(muvec)
                    mom=get_ini_magmom(el)*muvec
                    add=True
                else:
                    mom=np.zeros(3)
                ini_mag.append(mom)
            if add==True:
                magvecs.append(np.array(ini_mag))
        if len(magvecs)==0:
            return magconfigs
        updn=[[1]]
        for i in range(1,len(magvecs)):
            updnn=[]
            for config in updn:
                for spin in [-1,1]:
                    confign=deepcopy(config)
                    confign.append(spin)
                    updnn.append(confign)
            updn=updnn
        mid=0
        for config in updn:
            magmom=0*magvecs[0]
            for i in range(len(config)):
                magmom=magmom+config[i]*magvecs[i]
            if reduce_collinear:
                noncol=False
                momcol=[]
                for mom1 in magmom:
                    momcol.append(norm(mom1))
                    for mom2 in magmom:
                        if norm(np.cross(mom1,mom2))>0.01:
                            noncol=True
                            break
                if noncol==False:
                    magmom=np.array(momcol)
            magconfigs["%s_%d"%(name,mid)]={'magmom':magmom,'atoms_object':ao,'lnoncollinear':True}
            mid=mid+1
        return magconfigs
        mid=0
        for mvec in symmetry_constraints['mag_constraints']:
            ini_mag=[]
            for el,vec in zip(ao.get_chemical_symbols(),mvec):
                if (ismagnetic(el)) and (np.sum(np.fabs(mvec))>0.001):
                    mom=get_ini_magmom(el)*vec
                else:
                    mom=np.zeros(3)
                ini_mag.append(mom)
                print el,mom
            if (np.sum(np.fabs(np.array(ini_mag)))>0.001):
                magconfigs["%s_%d"%(name,mid)]=ini_mag
                mid=mid+1
        return magconfigs
    
    def get_symmetry_equivalent_positions(self,ao,positions=[], cell=[], magmoms=[], complete=False, silent=True):
        symmetry_constraints={}
        if ao!=None:
            positions=zip(ao.get_chemical_symbols(),ao.get_scaled_positions())
        equiv_set=range(len(positions))
        mag_constraints={}
        for i in range(len(positions)):
            el,pos=positions[i]
            #for (g,symb) in zip(self.get_elements(),self.get_symbols()):
            for g in self.get_elements():
                gpos=self.symop_pos(g,pos)
                gpos_in_pos=False
                #print pos,'-->',gpos
                for j in range(len(positions)):
                    el2,pos2=positions[j]
                    #if np.sum(np.fabs(gpos-pos2))<0.0001:
                    if self.is_equal_site(gpos,pos2):
                        gpos_in_pos=True
                        if (equiv_set[j]!=equiv_set[i]):
                            if (equiv_set[j]!=j):
                                print "symmetry not consistent with positions"
                                print self.symbols
                                print pos,g,gpos
                                print positions
                                print equiv_set,i,j
                                symmetry_constraints['isconsistent']=False
                                return symmetry_constraints
                            else:
                                equiv_set[j]=equiv_set[i]
                        break
                if gpos_in_pos==False:
                    print "get_symmetry_equivalent_positions:",el,gpos," not in ",positions
                    print pos,'-->',gpos,self.matrix2symbol(g)
                    if (complete==True):
                        positions.append((el,gpos))
                        equiv_set.append(i)
                        if magmoms!=[]:
                            magmoms.append(self.symop_mag(g,magmoms[i]))
                            print magmoms
                    else:
                        bla
                if (len(g)==3) and (i==j):
                    rot,trans,t=g
                    if (i in mag_constraints):
                        mag_constraints[i]=np.append(mag_constraints[i],np.linalg.det(rot)*t*rot-np.eye(3),axis=0)
                    else:
                        mag_constraints[i]=np.linalg.det(rot)*t*rot-np.eye(3)
        if silent==False:
            print "equivalent positions: ",equiv_set,len(equiv_set)
        if complete==True:
		els=[]
		spos=[]
		for i in range(len(positions)):
			el,pos=positions[i]
			els.append(el)
			spos.append(pos)
		for pos,mom in zip(positions,magmoms):
			print pos,mom 
		ao=Atoms(els,cell=cell, scaled_positions=spos,pbc=True) #pd_sc['chemical_symbols'],cell=pd_sc['cell'],scaled_positions=pd_sc['scaled_positions'],pbc=True)
		return ao,magmoms
        #constraints on magnetic moments
        magcon=[]
        for g in self.get_elements():
            if (len(g)==3):
                rot,trans,t=g
            else:
                break
            #print self.matrix2symbol(g),"***"
            magc=np.zeros((3*len(positions),3*len(positions)))
            for i in range(len(positions)):
                el,pos=positions[i]
                gpos=self.symop_pos(g,pos)
                for j in range(len(positions)):
                    el2,pos2=positions[j]
                    #if np.sum(np.fabs(gpos-pos2))<0.0001:
                    if self.is_equal_site(gpos,pos2):
                        break
                #print "#",i,pos,'-->',j,gpos
                RM=np.linalg.det(rot)*t*rot
                for i0 in range(3):
                    for j0 in range(3):
                        magc[3*i+i0][3*i+j0]=magc[3*i+i0][3*i+j0]+RM[i0][j0]
                    magc[3*i+i0][3*j+i0]=magc[3*i+i0][3*j+i0]-1.
            #print "***",self.matrix2symbol(g)
            if magcon==[]:
                magcon=magc
            else:
                magcon=np.append(magcon,magc,axis=0)
        #print "X"
        #print magcon
        nulsp=self.get_null_space(matrix(magcon))
        #print "********"
        #print np.round(self.get_null_space(matrix(magcon)),decimals=3)
        #print "********"
        #print self.gauss_elim(magcon)
        #print "********"
        symmetry_constraints['equivalent_atoms']=equiv_set
        symmetry_constraints['mag_constraints']=[]
        #print "Y"
        for m in nulsp:
            moms=m.reshape(len(positions),3)
            for i in range(len(moms)):
                norm=np.sqrt(np.dot(moms[i],moms[i]))
                if norm<0.001:
                    moms[i]=np.zeros(3)
                else:
                    moms[i]=moms[i]/norm
            symmetry_constraints['mag_constraints'].append(m.reshape(len(positions),3))
        return symmetry_constraints

    def mult(self,a,b):
        if len(a)==2:
            r_a,t_a=a
            r_b,t_b=b
            r_ab=np.dot(r_a,r_b)
            t_ab=np.round(np.dot(r_a,t_b)+t_a,decimals=8)%1.
            return (r_ab,t_ab)
        if len(a)==3:
            r_a,t_a,ti_a=a
            r_b,t_b,ti_b=b
            r_ab=np.dot(r_a,r_b)
            #check this
            t_ab=np.dot(r_a,t_b)+t_a
            for i in range(3):
                if np.round(t_ab[i],decimals=5)%1.<self.eps:
                    t_ab[i]=0.
                else:
                    t_ab[i]=t_ab[i]%1.
            #t_ab=np.round(np.dot(r_a,t_b)+t_a,decimals=8)%1.
            ti_ab=ti_a*ti_b
            return (r_ab,t_ab,ti_ab)
        return None
    
    def get_pure_translations(self):
        trans=[]
        for g in self.get_elements():
            if len(g)==2:
                r,t=g
            elif len(g)==3:
                r,t,ti=g
                if ti==-1:
                    continue
            R_u,t_u=self.symbol2matrix('1 x,y,z') # unity
            #print R_u,t_u
            if (np.sum(np.fabs(r-R_u))==0) and  (np.sum(np.fabs(t-t_u))>0.00001):
                trans.append(self.matrix2symbol(g))
        return trans
            
            
    def is_group(self, complete=False):
        """check if present elements form a group
        complete=True: complete if necessary"""
        # print "check_point190, is group" 
        G=self.get_elements()
        # print "check_point191, G is" , G
        isgroup=True
        #check completeness:
        iscomplete=False
        while (iscomplete==False) and (len(G)<1000):
            iscomplete=True
            for a in G:
                # print "a is ", a
                for b in G:
                    # print "b is  ", b
                    ab=self.mult(a,b)
                    if not self.is_element(ab):
                        isgroup=False
                        iscomplete=False
                        if complete==False:
                            print ab," not in G",self.matrix2symbol(a),"*",self.matrix2symbol(b),"=",self.matrix2symbol(ab)
                            return False
                        else:
                            G.append(ab)
                            #print "completing: ",len(G),self.matrix2symbol(a),"*",self.matrix2symbol(b),"=",self.matrix2symbol(ab)
        self.elements=G
        hasunity=False
        for a in G:
            if self.is_unity(a,G):
                hasunity=True
                unity=a
                #print "unity: ",unity
        if not hasunity:
            print "1 not found in G"
            return False
        for a in G:
            hasinverse=False
            for b in G:
                if self.is_equal_element(unity,self.mult(a,b)):
                    hasinverse=True
                    break
            if not hasinverse:
                print "No inverse found for ",a
                return False
        if complete==True:
            return isgroup,G
        return isgroup

    def is_subgroup(self,H,check_group=True):
        """returns True if H is a (proper) subgroup of the present group G"""
        if (len(H.get_elements())>=len(self.get_elements())):
            return False
        if (check_group==True) and (not H.is_group()):
            return False
        for g in H.get_elements():
            if not self.is_element(g):
                #print "Not in group:",g
                return False
        return True

    def get_left_coset(self,g):
        gG=MSG()
        for h in self.get_elements():
            gh=self.mult(g,h)
            gG.elements.append(gh)
        return gG
    
    def is_maximal_subgroup(self,H,check_group=True, silent=False):
        """brute force: returns True if H is a maximal subgroup of the present group G"""
        if not self.is_subgroup(H,check_group=check_group):
            if silent==False:
                print "Not subgroup"
            return False
        ordG=len(self.get_elements())
        #get cosets with respect to primitive translations and tr:
        tG={}
        symbs=['1 x,y,z+1/2,+1 mx,my,mz','1 x,y+1/2,z,+1 mx,my,mz','1 x+1/2,y,z,+1 mx,my,mz','1 x,y,z,-1 -mx,-my,-mz']
        for symb in symbs:
            g=self.symbol2matrix(symb)
            if H.is_element(g):
                if silent==False:
                    print "Excluded, because of element ",g
                return False #exclude grey subgroups etc.
            if self.is_element(g):
                tG[symb]=H.get_left_coset(g)
                #print "left coset:",symb,len(tG[symb].get_elements()),len(self.get_elements())
        for g in self.get_elements():
            notincoset=True
            for symb in tG:
                if tG[symb].is_element(g):
                    notincoset=False
            #print H.is_element(g),not notincoset
            if (not H.is_element(g)) and (notincoset):
                Hp=deepcopy(H)
                Hp.elements.append(g)
                Hp.is_group(complete=True)
                if len(Hp.get_elements())!=ordG:
                    if silent==False:
                        print "not maximal: found subgroup < G",g
                        print "ord(G)=",ordG,"ord(<H,g>)=",len(Hp.get_elements()),Hp.is_group()
                    return False
                #else:
                #    print "ord(G)=ord(<G,g>) for g=",g
        return True

    
    def get_unity_i(self):
        if self.unity==None:
            ordG=self.get_order()
            multtab=self.get_multiplication_table()
            for i in range(ordG):
                unity=i
                for j in range(ordG):
                    if multtab[i][j]!=j:
                        unity=-1
                        break
                if unity==i:
                    break
            if unity==i:
                self.unity=unity
        return self.unity

    
    def get_inverse_element_i(self, i):
        unity=self.get_unity_i()
        ordG=self.get_order()
        inverse=None
        multtab=self.get_multiplication_table()
        for j in range(ordG):
            if multtab[i,j]==unity:
                inverse=j
                break
        return inverse

    def get_element_indices(self, elements):
        """returns indices of list of elements in this group"""
        elements_i=[]
        ordG=self.get_order()
        for g in elements:
            for i in range(ordG):
                if self.is_equal_element(g, self.get_elements()[i]):
                    elements_i.append(i)
                    break
        if len(elements_i)==len(elements):
            return elements_i
        return None #not in group
                    
    def complete_subgroup_i(self, subg_i):
        """returns indeces of completed subgroup"""
        multtab=self.get_multiplication_table()
        is_complete=False
        subg=deepcopy(subg_i)
        while is_complete==False:
            complete_set=deepcopy(subg)
            for g in subg:
                for h in subg:
                    k=multtab[g][h]
                    if (not (k in complete_set)):
                        complete_set.append(k)
            if len(complete_set)==len(subg):
                is_complete=True
            subg=complete_set
        return complete_set
    
    def get_conjugate_subgroups_i(self, subg_i):
        """get conjugate subgroups"""
        #print "SUBG_I",subg_i
        ordG=self.get_order()
        multtab=self.get_multiplication_table()
        #print multtab
        csub_i=[subg_i]
        for i in range(ordG):
            im=self.get_inverse_element_i(i)
            #print i,"inv",im,"1",self.get_unity_i()
            cg=[]
            for j in subg_i:
                #print "ZZZ",j,multtab[i,j],im,subg_i
                ghgm=multtab[multtab[i,j],im]
                #print "XXX",ghgm
                cg.append(ghgm)
            #print "CG",cg
            cg_sorted=sorted(cg)
            if not cg_sorted in csub_i:
                csub_i.append(cg_sorted)
        return csub_i
    

    def is_maximal_subgroup_from_multtab(self,subg_i,excl,silent=True,check_group=True,elspr=[]):
        """brute force: returns True if H is a maximal subgroup of the present group G"""
        elsprobed=deepcopy(elspr)
        ordG=self.get_order()
        multtab=self.get_multiplication_table()
        if (check_group==False):
            if (len(subg_i)*2==ordG):
                return True
            if (len(subg_i)*3==ordG):
                return True
            elif (len(subg_i)==ordG):
                return False
        else:
            for g in subg_i:
                for h in subg_i:
                    if not (multtab[g][h] in subg_i):
                        return False
            #more checks needed?
        #get cosets with respect to primitive translations and tr:
        lcoset={}
        for i in excl:
            lcoset[i]=[]
            for j in subg_i:
                lcoset[i].append(multtab[i][j])
        #print "XX5: cosets"
        for i in range(ordG): #range(self.get_order()):
            notincoset=True
            for g in lcoset:
                if i in lcoset[g]:
                    notincoset=False
            if (not (i in subg_i)) and (notincoset):
                #print "XX6: completing",len(subg_i),i,len(elsprobed)
                subg=deepcopy(subg_i)
                for j in subg_i:
                    subg.append(multtab[i][j])
                #subg.append(i)
                is_complete=False
                skip=False
                while (is_complete==False):
                    complete_set=deepcopy(subg)
                    for g in subg:
                        if skip==True:
                            break
                        for h in subg[len(subg_i):]:
                            k=multtab[g][h]
                            if (k in excl):
                                #print "XX7",k,"in excl"
                                skip=True
                                break
                            elif (not (k in complete_set)):
                                complete_set.append(k)
                    if len(complete_set)==len(subg):
                        is_complete=True
                    subg=complete_set
                if (skip==False) and (len(subg)!=ordG):
                    #if silent==False:
                    #    print self.matrix2symbol(self.elements[i]),len(subg),subg
                    return False
        return True


    def get_maximal_subgroups_old(self,subgr=None,step=0,elementset='auto',excl=[],maxsubs=[],subgrlist=[],maxsteps=7,symdir=None,sym_parent=None,q=[0,0,0]):
        """brute force search for maximal subgroups: returns dictionary/list of maximal subgroups"""
        if maxsubs==None:
            maxsubs=[]
        if subgrlist==None:
            subgrlist=[]
        if excl==None:
            excl=[]
        #check if maximal subgroups are in symdir:
        if (step==0) and (symdir!=None) and (sym_parent!=None) and ('number' in sym_parent):
            msg_parent=MSG(symelem=zip(sym_parent['rotations'],sym_parent['translations']))
            symdir_parent=os.path.join(symdir, str(sym_parent['number']))
            if (not os.path.isdir(symdir_parent)):
                os.makedirs(symdir_parent)
            symfile_parent=os.path.join(symdir_parent,'sym_parent.txt')
            qstr="q_%5.3f_%5.3f_%5.3f"%(q[0],q[1],q[2])
            symdirq=os.path.join(symdir_parent,qstr)
            if (os.path.isfile(symfile_parent)):
                msg1=MSG(mcif=symfile_parent)
                if not msg_parent.is_equal(msg1):
                    print "parent symmetry inconsistent with %s, call programmer!"%symfile_parent
                    return maxsubs
                else:
                    if os.path.isdir(symdirq):
                        maxsubs={}
                        mcifs=glob.glob(os.path.join(symdirq,"*.mcif"))
                        for mcif in mcifs:
                            msg=MSG(mcif=mcif)
                            name=qstr+mcif.split('/')[-1]
                            maxsubs[name]=msg
                        return maxsubs
            else:
                msg_parent.write_mcif(symfile_parent)
        #determine max subgroups
        if (step==0):
            print "Determining maximal subgroups, this will take some time..."
        if step>=maxsteps:
            print "Reached step limit"
            return maxsubs
        if elementset=='auto':
            elementset=[]
            symbs=['x,y,z+1/2,+1','x,y+1/2,z,+1','x+1/2,y,z,+1','x+1/2,y+1/2,z,+1','x+1/2,y,z+1/2,+1','x,y+1/2,z+1/2,+1','x+1/2,y+1/2,z+1/2,+1','x,y,z,-1']
            excl=[]
            for g in self.get_elements():
                if self.matrix2symbol(g,print_mag=False) in symbs:
                    print "not using ",self.matrix2symbol(g)
                    excl.append(g)
                else:
                    elementset.append(g)
        for i in range(len(elementset)):
            if (step==0):
                print "%d of %d.."%(i+1,len(elementset))
            g=elementset[i]
            if subgr==None:
                subg=MSG()
                subg.elements.append(self.get_unity())
            else:
                subg=deepcopy(subgr)
            if not subg.is_element(g):
                subg.elements.append(g)
                subg.is_group(complete=True)
                subg_i=[]
                for j in range(len(elementset)):
                    if subg.is_element(elementset[j]):
                        subg_i.append(j)
                subg_i=sorted(subg_i)
                skip=False
                for h in excl:
                    if subg.is_element(h):
                        skip=True
                        break
                for G in subgrlist:
                    if G==subg_i:
                        skip=True
                for G in maxsubs:
                    if G==subg_i:
                        skip=True
                if skip==True:
                    continue
                subgrlist.append(subg_i)
                print len(subgrlist)," subgroups found"
                if self.is_maximal_subgroup(subg):
                    maxsubs.append(subg)
                    print "Found max sub of order ",len(subg.get_elements()),"-- #",len(maxsubs)
                else:
                    maxsubs=self.get_maximal_subgroups(subgr=subg,elementset=elementset,excl=excl,subgrlist=subgrlist,maxsubs=maxsubs,step=step+1)
        #ouput of symmetry cards in symdir
        if (step==0) and (symdir!=None) and (sym_parent!=None) and ('number' in sym_parent):
            os.mkdir(symdirq)
            i=0
            maxs={}
            for msg in maxsubs:
                fname=os.path.join(symdirq,"msg_%d.mcif"%i)
                msg.write_mcif(fname)
                name=qstr+"msg_%d.mcif"%i
                maxs[name]=msg
                i=i+1
            maxsubs=maxs
        return maxsubs


    def get_maximal_subgroups(self,subgr=None,step=0,elementset='autox',excl=None,maxsubs=None,subgrlist=None,maxsteps=10,symdir=None,symdirhte=None,sym_parent=None,q=[0,0,0],exclude_conjugate_subgroups=True,engine_new=False,silent=True):
        """brute force search for maximal subgroups: returns dictionary/list of maximal subgroups"""
        if maxsubs==None:
            maxsubs=[]
        if subgrlist==None:
            subgrlist=[]
        if excl==None:
            excl=[]
        #print "MSG:",step,subgr
        #check if maximal subgroups are in symdir:
        if (step==0) and (symdir!=None) and (sym_parent!=None) and ('number' in sym_parent):
            msg_parent=MSG(symelem=zip(sym_parent['rotations'],sym_parent['translations']),eps=self.eps)
            if  (msg_parent.is_group()==False):
                print "Check parent group:"
                for g in msg_parent.elements:
                    print g,self.matrix2symbol(g)
                return maxsubs
            symdir_parent=os.path.join(symdir, str(sym_parent['number']))
            if (not os.path.isdir(symdir_parent)):
                os.makedirs(symdir_parent)
            symfile_parent=os.path.join(symdir_parent,'sym_parent.txt')
            qstr="q"
            for i in range(3):
                if fabs(q[i])<0.01:
                    qstr=qstr+"_0"
                else:
                    qstr=qstr+"_%5.3f"%q[i]
            symdirq=os.path.join(symdir_parent,qstr)
            #print symdirhte,str(sym_parent['number']),qstr
            if (symdirhte!=None) and (os.path.isdir(os.path.join(symdirhte,str(sym_parent['number']),qstr))):
                #todo: checks below
                symdirq=os.path.join(symdirhte,str(sym_parent['number']),qstr)
                maxsubs={}
                mcifs=glob.glob(os.path.join(symdirq,"msg*.mcif"))
                for mcif in mcifs:
                    msg=MSG(mcif=mcif,eps=self.eps)
                    name=qstr+mcif.split('/')[-1]
                    maxsubs[name]=msg
                print "Using MSG tables in %s"%symdirq,symfile_parent
                #print maxsubs
                return maxsubs
            print "Using MSG tables in %s"%symdirq,symfile_parent
            if (os.path.isfile(symfile_parent)):
                msg1=MSG(mcif=symfile_parent,eps=self.eps)
                print msg1.is_group()
                if not msg_parent.is_equal(msg1):
                    print "parent symmetry inconsistent with %s, call programmer!"%symfile_parent
                    #msg_parent.write_mcif("%s-debug"%symfile_parent)
                    return maxsubs
                else:
                    if os.path.isdir(symdirq):
                        maxsubs={}
                        mcifs=glob.glob(os.path.join(symdirq,"msg*.mcif"))
                        for mcif in mcifs:
                            msg=MSG(mcif=mcif, eps=self.eps)
                            name=qstr+mcif.split('/')[-1]
                            maxsubs[name]=msg
                        return maxsubs
            else:
                msg_parent.write_mcif(symfile_parent)
        #determine max subgroups
        if (step==0):
            print "Determining maximal subgroups, this will take some time...",self.get_order()
            #tmp
            multfile=os.path.join(symdir_parent,"%s_multtab.dat"%qstr)
            order=self.get_order()
            if os.path.isfile(multfile):
                fin=open(multfile,"r")
                multtab=np.zeros((order,order),dtype=int)
                print "multtab from",multfile
                for i in range(order):
                    line=fin.readline()
                    lsp=line.split()
                    for j in range(order):
                        multtab[i][j]=int(lsp[j])
                fin.close()
                self.multiplication_table=multtab
                print "multtab done"
            else:
                multtab=self.get_multiplication_table()
                print "multtab done"
                fout=open(multfile,"w")
                for i in range(order):
                    for j in range(order):
                        fout.write("%d "%multtab[i][j])
                    fout.write("\n")
                fout.close()
                
        if step>=maxsteps:
            print "Reached step limit",subgr
            return maxsubs
        if elementset=='auto':
            elementset=[]
            #symbs=['x,y,z+1/2,+1','x,y+1/2,z,+1','x+1/2,y,z,+1','x+1/2,y+1/2,z,+1','x+1/2,y,z+1/2,+1','x,y+1/2,z+1/2,+1','x+1/2,y+1/2,z+1/2,+1','x,y,z,-1'] #todo: exclude all primitive translations
            symbs=['x,y,z+1/2,+1','x,y+1/2,z,+1','x+1/2,y,z,+1'] #todo: exclude all primitive translations
            symbsinc=['x,y,z+1/2,-1','x,y+1/2,z,-1','x+1/2,y,z,-1'] #todo: exclude all primitive translations
	    if qstr=='q_0_0_1.000':
             #print "centered group"
             symbsinc=['x,y,z+1/2,-1','x,y+1/2,z,-1','x+1/2,y,z,-1','x+1/2,y+1/2,z+1/2,-1','x+1/2,y+1/2,z,-1','x+1/2,y,z+1/2,-1','x,y+1/2,z+1/2,-1']
            symb_ti='x,y,z,-1'
            excl=[]
            subgr_0=[]
            for i in range(self.get_order()):
                g=self.elements[i]
                symb=self.matrix2symbol(g,print_mag=False)
                if symb==symb_ti:
                    print "not using ",i,symb
                    excl.append(i)
                elif symb in symbs:
                    print "not using ",self.matrix2symbol(g)
                else:
                    if symb in symbsinc: #that is what is done on the Bilbao server?
                        print "forcing inclusion of ",i,symb
                        subgr_0.append(i)
                    elementset.append(i)
            subgr=self.complete_subgroup_i(subgr_0)
            print "Starting from ",subgr
        multtab=self.get_multiplication_table()
	#print "multtab done"
        elsprobed=[]
        if elementset=='autox': #make tests
            elementset=[]
            #symbs=['x,y,z+1/2,+1','x,y+1/2,z,+1','x+1/2,y,z,+1','x+1/2,y+1/2,z,+1','x+1/2,y,z+1/2,+1','x,y+1/2,z+1/2,+1','x+1/2,y+1/2,z+1/2,+1','x,y,z,-1'] #todo: exclude all primitive translations
            symbs=['x,y,z+1/2,+1','x,y+1/2,z,+1','x+1/2,y,z,+1','x,y,z+1/4,+1','x,y+1/4,z,+1','x+1/4,y,z,+1','x,y,z+1/4,-1','x,y+1/4,z,-1','x+1/4,y,z,-1'] #todo: exclude all primitive translations
            symbsinc=['x,y,z,+1','x,y,z+1/2,-1','x,y+1/2,z,-1','x+1/2,y,z,-1'] #
	    if qstr=='q_0_0_1.000':
             print "centered group"
             symbsinc=['x,y,z,+1','x,y,z+1/2,-1','x,y+1/2,z,-1','x+1/2,y,z,-1','x+1/2,y+1/2,z+1/2,-1','x+1/2,y+1/2,z,-1','x+1/2,y,z+1/2,-1','x,y+1/2,z+1/2,-1']
             symbsinc=['x,y,z,+1','x,y,z+1/2,-1','x,y+1/2,z,-1','x+1/2,y,z,-1','x+1/2,y+1/2,z+1/2,-1','x+1/2,y+1/2,z,-1','x+1/2,y,z+1/2,-1']
            if (self.elements!=[]) and (len(self.elements[0])==2):
                symbsinc=['x,y,z']
            symb_ti='x,y,z,-1'
            excl=[]
            subgr_0=[]
            for i in range(self.get_order()):
                g=self.elements[i]
                symb=self.matrix2symbol(g,print_mag=False)
                if symb==symb_ti:
                    print "not using ",i,symb
                    excl.append(i)
                elif symb in symbs:
                    print "not using ",self.matrix2symbol(g)
                else:
                    if symb in symbsinc: #that is what is done on the Bilbao server?
                        print "forcing inclusion of ",i,symb
                        subgr_0.append(i)
                    elementset.append(i)
            subgr=self.complete_subgroup_i(subgr_0)
            print "Starting from ",subgr
        if (engine_new==True):
            maxsubs=self.get_maximal_subgroups_engine_new(step=0,subgr=subgr,elementset_to_be_probed=elementset,excl=excl,maxsubs=[],subgrlist=[])
            print 'get_maximal_subgroups_engine_new() done'
            elementset=[]
        for i in elementset:
            if subgr==None:
                subg=[]
                #subg=deepcopy(subgr_0) # complete it here?
            else:
                subg=deepcopy(subgr)
            ordsg=len(subg)
            if not (i in subg):
                subg.append(i)
                is_complete=False
                skip=False
                #print "compl: step, el, ord:",step,i,len(subg)
                while (is_complete==False) and (skip==False):
                    complete_set=deepcopy(subg)
                    for g in subg:
                        if skip==True:
                            break
                        for h in subg[ordsg:]:
                            k=multtab[g][h]
                            if (k in elsprobed) or (not (k in elementset)):
                                skip=True
                                break
                            if (not (k in complete_set)):
                                complete_set.append(k)
                    if len(complete_set)==len(subg):
                        is_complete=True
                    subg=complete_set
                subg_i=sorted(subg)
                elsprobed.append(i)
                if skip==True:
                    continue
                #print "XX1: step, el, ord:",step,i,len(subg_i)
                #print "subgroup:",subg_i
                skip=False
                for h in subg:
                    if not (h in elementset):
                        skip=True
                if skip==True:
                    continue
                #for h in excl:
                #    if h in subg:
                #        skip=True
                #        break
                for G in subgrlist:
                    if G==subg_i:
                        skip=True
                        break
                if skip==True:
                    continue
                for G in maxsubs:
                    if (skip==True) or (G==subg_i):
                        skip=True
                        break
                if skip==True:
                    continue
                subgrlist.append(subg_i)
                #print "XX2: step, el, ord:",step,i,len(subgrlist)
                for subgc in self.get_conjugate_subgroups_i(subg_i):
                    inlist=False
                    for G in subgrlist:
                        if G==subgc:
                            inlist=True
                            break
                    if inlist==False:
                        subgrlist.append(subgc)
                #print "XX3: step, el, nsub:",step,i,len(subgrlist)
                if self.is_maximal_subgroup_from_multtab(subg_i,excl,silent=silent,check_group=False,elspr=elsprobed):
                    maxsubs.append(subg_i)
                    print "Found maximal subgroup of order ",len(subg_i),"-- #",len(maxsubs)
                    print len(subgrlist)," subgroups found"
                    line="-- "
                    for i in subg_i:
                        line=line+self.matrix2symbol(self.elements[i])+";"
                    #print line
                else:
                    #print "XX4"
                    #print "next step"
                    maxsubs=self.get_maximal_subgroups(subgr=subg_i,elementset=elementset,excl=excl,subgrlist=subgrlist,maxsubs=maxsubs,step=step+1)
        if (step==0):
            #print "ZZZ"
            mxs=[]
            #check for conjugate subgroups
            cggroups=[]
            for i in range(len(self.elements)):
                unity=i
                for j in range(len(self.elements)):
                    if multtab[i][j]!=j:
                        unity=-1
                        break
                if unity==i:
                    break
            for sub in maxsubs:
                if sub in cggroups:
                    continue
                for i in range(len(self.elements)):
                    for j in range(len(self.elements)):
                        if multtab[i][j]==unity:
                            im1=j
                            break 
                    cg=[]
                    for j in sub:
                        k=multtab[i][j]
                        cg.append(multtab[k][im1])
                    cg=sorted(cg)
                    for sub2 in maxsubs:
                        if (sub!=sub2) and (cg==sub2) and (not (sub2 in cggroups)):
                            if (exclude_conjugate_subgroups==True):
                                cggroups.append(sub2)
                                #print "conjugate:",sub,sub2
            cmxsubs=[]
            for sub in maxsubs:
                if sub in cggroups:
                    continue
                elements=[]
                for i in sub:
                    elements.append(self.elements[i])
                mxs.append(MSG(symelem=elements,eps=self.eps))
                #tmp for debug
                cmsg=[]
                for cg in self.get_conjugate_subgroups_i(sub):
                    elements=[]
                    for i in cg:
                        elements.append(self.elements[i])
                    cmsg.append(MSG(symelem=elements,eps=self.eps))
                cmxsubs.append(cmsg)
            maxsubs=mxs
        #ouput of symmetry cards in symdir
        if (step==0) and (symdir!=None) and (sym_parent!=None) and ('number' in sym_parent):
            os.mkdir(symdirq)
            fname=os.path.join(symdirq,"grey.mcif")
            self.write_mcif(fname)
            i=0
            maxs={}
            for msg in maxsubs:
                fname=os.path.join(symdirq,"msg_%d.mcif"%i)
                msg.write_mcif(fname)
                name=qstr+"msg_%d.mcif"%i
                maxs[name]=msg
                #tmp debug
                j=0
                for msg in cmxsubs[i]:
                    msg.write_mcif(fname+"_%d"%j)
                    j=j+1
                i=i+1
            maxsubs=maxs
	    #bla
        return maxsubs


    def get_maximal_subgroups_engine_new(self,step=0,subgr=[],elementset_to_be_probed='auto',
                                         excl=None,maxsubs=[],subgrlist=[],maxsteps=10,symdir="msg"
                                         ,exclude_conjugate_subgroups=True,sym_parent=None,silent=True):
                                         #msgpaths=['msg/'],msgnames=[],sym_files_parent=["grey.mcif"],
                                         #elementset='autox',excl=None,maxsubs=None,subgrlist=None,maxsteps=10,symdir=None,symdirhte=None,sym_parent=None,q=[0,0,0]
                                         #,exclude_conjugate_subgroups=True,silent=True):
        """brute force search for maximal subgroups: returns dictionary/list of maximal subgroups"""
        # todo: step=0: initialize/check for tabulated maximal subgroups
        if step>=maxsteps:
            print "Reached step limit",subgr
            return maxsubs
        multtab=self.get_multiplication_table()
        for i in elementset_to_be_probed: # set of elements which hasn't been probed for maximal subgroup before
            subg=deepcopy(subgr)
            ordsg=len(subg)
            if not (i in subg):
                subg.append(i)
                is_complete=False
                skip=False
                #print "compl: step, el, ord:",step,i,len(subg)
                while (is_complete==False) and (skip==False):
                    complete_set=deepcopy(subg)
                    for g in subg:
                        if skip==True:
                            break
                        for h in subg[ordsg:]:
                            k=multtab[g][h]
                            if (not (k in subgr+elementset_to_be_probed)):
                                skip=True
                                break
                            if (not (k in complete_set)):
                                complete_set.append(k)
                    if len(complete_set)==len(subg):
                        is_complete=True
                    subg=complete_set
                subg_i=sorted(subg)
                if skip==True:
                    continue
                #print "XX1: step, el, ord:",step,i,len(subg_i)
                for G in subgrlist:
                    if G==subg_i:
                        #print "Should not happen G==subg_i"
                        skip=True
                        break
                if skip==True:
                    continue
                for G in maxsubs:
                    if (skip==True) or (G==subg_i):
                        #print "Should not happen G==maxsubs"
                        skip=True
                        break
                if skip==True:
                    continue
                subgrlist.append(subg_i)
                #print "XX2: step, el, ord:",step,i,len(subgrlist)
                for subgc in self.get_conjugate_subgroups_i(subg_i):
                    inlist=False
                    for G in subgrlist:
                        if G==subgc:
                            inlist=True
                            break
                    if inlist==False:
                        subgrlist.append(subgc)
                #print "XX3: step, el, nsub:",step,i,len(subgrlist)
                if self.is_maximal_subgroup_from_multtab(subg_i,excl,silent=silent,check_group=False):
                    maxsubs.append(subg_i)
                    print "Found maximal subgroup of order ",len(subg_i),"-- #",len(maxsubs)
                    print len(subgrlist)," subgroups found"
                    line="-- "
                    for i in subg_i:
                        line=line+self.matrix2symbol(self.elements[i])+";"
                    #print line
                else:
                    #print "XX4"
                    #print "next step"
                    elsprobn=deepcopy(elementset_to_be_probed[elementset_to_be_probed.index(i)+1:])
                    maxsubs=self.get_maximal_subgroups_engine_new(step=step+1,subgr=subg_i,elementset_to_be_probed=elsprobn,excl=excl,subgrlist=subgrlist,maxsubs=maxsubs)
            #else:
            #    print "Should not happen: i in subg"
        if (step==0):
            return maxsubs
            #print "ZZZ"
            mxs=[]
            #check for conjugate subgroups
            cggroups=[]
            for i in range(len(self.elements)):
                unity=i
                for j in range(len(self.elements)):
                    if multtab[i][j]!=j:
                        unity=-1
                        break
                if unity==i:
                    break
            for sub in maxsubs:
                if sub in cggroups:
                    continue
                for i in range(len(self.elements)):
                    for j in range(len(self.elements)):
                        if multtab[i][j]==unity:
                            im1=j
                            break 
                    cg=[]
                    for j in sub:
                        k=multtab[i][j]
                        cg.append(multtab[k][im1])
                    cg=sorted(cg)
                    for sub2 in maxsubs:
                        if (sub!=sub2) and (cg==sub2) and (not (sub2 in cggroups)):
                            if (exclude_conjugate_subgroups==True):
                                cggroups.append(sub2)
                                #print "conjugate:",sub,sub2
            cmxsubs=[]
            for sub in maxsubs:
                if sub in cggroups:
                    continue
                elements=[]
                for i in sub:
                    elements.append(self.elements[i])
                mxs.append(MSG(symelem=elements,eps=self.eps))
                #tmp for debug
                cmsg=[]
                for cg in self.get_conjugate_subgroups_i(sub):
                    elements=[]
                    for i in cg:
                        elements.append(self.elements[i])
                    cmsg.append(MSG(symelem=elements,eps=self.eps))
                cmxsubs.append(cmsg)
            maxsubs=mxs
        #ouput of symmetry cards in symdir
        if (step==0) and (symdir!=None) and (sym_parent!=None) and ('number' in sym_parent):
            os.mkdir(symdirq)
            fname=os.path.join(symdirq,"grey.mcif")
            self.write_mcif(fname)
            i=0
            maxs={}
            for msg in maxsubs:
                fname=os.path.join(symdirq,"msg_%d.mcif"%i)
                msg.write_mcif(fname)
                name=qstr+"msg_%d.mcif"%i
                maxs[name]=msg
                #tmp debug
                j=0
                for msg in cmxsubs[i]:
                    msg.write_mcif(fname+"_%d"%j)
                    j=j+1
                i=i+1
            maxsubs=maxs
	    #bla
        return maxsubs

    def is_element(self,a,G=None):
        if G==None:
            g=self.get_elements()
        else:
            g=G.get_elements()
        if len(a)==2:
            r,t=a
            for (r_g,t_g) in g:
                if (np.sum(np.fabs(r-r_g))==0) and (np.sum(np.fabs(t-t_g))<0.00001):
                    return True
        if len(a)==3:
            r,t,ti=a
            for (r_g,t_g,ti_g) in g:
                if (np.sum(np.fabs(r-r_g))==0) and (np.sum(np.fabs(t-t_g))<0.00001) and (ti==ti_g):
                    return True
        return False

    def is_equal(self, G):
        """returns true if the (magnetic) space group G has the same elements as the present MSG object"""
        if (len(G.get_elements())!=len(self.get_elements())):
            return False
        for g in self.get_elements():
            ginG=False
            for h in G.get_elements():
                if len(g)!=len(h):
                    return False
                if self.is_equal_element(g,h)==True:
                    ginG=True
                    break
            if ginG==False:
		#print "not in G:",self.matrix2symbol(g)
                return False
        return True

    
    def is_equal_element(self,a,b,ignore_lattice_translations=False):
        if len(a)==2:
            r,t=a
            r_g,t_g=b
            if (np.sum(np.fabs(r-r_g))==0) and (np.sum(np.fabs(t-t_g))<0.00001):
                return True
            if (ignore_lattice_translations==True) and (np.sum(np.fabs(r-r_g))==0):
                for i in range(3):
                    is_ok=False
                    for j in range(-3,3):
                        #print 'ZZZ',(abs(t[i]-t_g[i]+j),self.eps)
                        if (abs(t[i]-t_g[i]+j)<self.eps):
                            is_ok=True
                    if (is_ok==False):
                        return False
                    #if ((np.round(abs(t[i]-t_g[i]),decimals=5)%1.)>self.eps):
                    #    print t[i],t_g[i],"%.10f"%abs(t[i]-t_g[i]),np.round(abs(t[i]-t_g[i]),decimals=5)%1.,self.eps
                    #    return False
                return True
        if len(a)==3:
            r,t,ti=a
            r_g,t_g,ti_g=b
            if (np.sum(np.fabs(r-r_g))==0) and (np.sum(np.fabs(t-t_g))<0.00001) and (ti==ti_g):
                return True
            if (ignore_lattice_translations==True) and (np.sum(np.fabs(r-r_g))==0) and (ti==ti_g):
                for i in range(3):
                    if ((abs(t[i]-t_g[i])%1)>self.eps):
                        return False
                return True                
        return False

    def is_unity(self,a,G):
        for g in G:
            ag=self.mult(a,g)
            if not self.is_equal_element(ag,g):
                return False
        return True

    def get_unity(self):
        for g in self.get_elements():
            isunity=True
            for h in self.get_elements():
                gh=self.mult(g,h)
                if not self.is_equal_element(gh,h):
                    isunity=False
                    break
            if isunity==True:
                return g
        return None
        
    def grey_msg(self,sg=None):
        spg=self
        if sg!=None:
            spg=sg
        else:
            spg=self
        msg=MSG(eps=self.eps)
        for (r,t) in spg.get_elements():
            for ti in [+1,-1]:
                msg.elements.append((r,t,ti))
        return msg


def ismagnetic(el):
    # magatoms={'Mn':3.5,'Fe':2.5,'Co':1.5,'Ni':0.8}
    magatoms={
        'Mn': 4.0, 'Fe': 3.5, 'Co': 2.0, 'Ni': 1.5,
        'Cr': 3.0, 'Os': 1.5, 'Cu': 0.5, 'Eu': 6.5,
        'Gd': 6.0, 'Tb': 5.0, 'Dy': 5.0, 'Ho': 6.5,
        'Er': 6.0, 'Tm': 6.0, 'Yb': 4.5, 'Pr': 3.2,
        'Nd': 3.6, 'Sm': 0.8, 'V':2.0,'U':1.0,'Nd':2.0,'Ru':2.0,'Ir':1.0
        }
    if el in magatoms:
        return True
    return False

def get_ini_magmom(el):
#  magatoms={'Mn':3.5,'Fe':2.5,'Co':1.5,'Ni':0.8}
    magatoms={
        'Mn': 4.0, 'Fe': 3.5, 'Co': 2.0, 'Ni': 1.5,
        'Cr': 3.0, 'Os': 1.5, 'Cu': 0.5, 'Eu': 6.5,
        'Gd': 6.0, 'Tb': 5.0, 'Dy': 5.0, 'Ho': 6.5,
        'Er': 6.0, 'Tm': 6.0, 'Yb': 4.5, 'Pr': 3.2,
        'Nd': 3.6, 'Sm': 0.8, 'V':2.0,'U':1.0,'Nd':2.0,'Ru':2.0,'Ir':1.0
        }
    if el in magatoms:
        return magatoms[el]
    return 0

def get_magnetic_sublattices(ao,spglib_info=None,return_afm=False, symprec=1e-3, Nmax=-1):
    """returns the number of magnetic sublattices which are not connected by symmetry (default: return_afm=False)
    or a dictionary with magnetic configurations arising from AF coupling of sublattices (return_afm=True).
    parameters:
    symprec: precision to determine space group symmetry with spglib (default: symprec=1e-3)
    Nmax: maximal number of magnetic configurations to be returned (default: Nmax=-1, all configurations)
    """
    #ase resorts atoms, so it is better to run spglib on current ao
    if spglib_info==None:
        try:
            spglib_info=spglib.get_symmetry_dataset(ao, symprec=symprec)
        except:
            print "Could not get symmetry information, check spglib!"
            return False
    equiv=spglib_info['equivalent_atoms']
    elements=ao.get_chemical_symbols()
    magsublatt=[]
    for (el,sort) in zip(elements,equiv):
        if (ismagnetic(el)) and (not (sort in magsublatt)):
            magsublatt.append(sort)
    if return_afm==False:
        return magsublatt
    afmconfigs={}
    if (return_afm==True) and (len(magsublatt)>1):
        updn=[{magsublatt[0]:get_ini_magmom(elements[magsublatt[0]])}]
        xspins=[-1,1]
        for i in range(1,len(magsublatt)):
            updnn=[]
            for config in updn:
                for spin in xspins:
                    confign=deepcopy(config)
                    confign[magsublatt[i]]=get_ini_magmom(elements[magsublatt[i]])*spin
                    updnn.append(confign)
            updn=updnn
            if (Nmax>0) and (len(updn)>Nmax):
                xspins=[-1]
        for config in updn:
            if (Nmax>0) and (len(afmconfigs)>=Nmax):
                break
            magmom=[]
            name="afm_q0_"
            isafm=False
            for i in sorted(config):
                if config[i]>=0:
                    name=name+elements[i]+"%d"%i+"+"
                else:
                    isafm=True
                    name=name+elements[i]+"%d"%i+"-"
            if isafm==True:
                for (el,sort) in zip(elements,equiv):
                    if sort in magsublatt:
                        magmom.append(config[sort])
                    else:
                        magmom.append(0)
                afmconfigs[name]={'magmom':magmom,'atoms_object':ao}
    return afmconfigs



def setup_supercell(ao, q=[0,0,1]):
    newcell=ao.get_cell()
    newpos=ao.get_scaled_positions()
    newelements=ao.get_chemical_symbols()
    for i in range(3):
        if q[i]!=0:
            n=int(1/q[i])
            newcell[i]=newcell[i]*n
            pos=[]
            elem=[]
            for j in range(n):
                t=1.0*j/n
                for (el,p) in zip(newelements,newpos):
                    np=copy(p)
                    np[i]=p[i]/n+t
                    pos.append(np)
                    elem.append(el)
            newpos=pos
            newelements=elem
    newao=Atoms(symbols=newelements,cell=newcell,scaled_positions=newpos,pbc=True)
    return newao


def Read_Magnetic_Space_Groups(fname, check_group=False):
    MagSpaceGroups={}
    symbsSC={}
    symbsSChex={}
    fin=open(fname,'r')
    line=fin.readline()
    while line:
        if line.startswith('Definition of point operator symbols (from Litvin)'):
            break
        line=fin.readline()
    line=fin.readline()
    print "* Non-hexagonal symbols"
    while line.startswith('Hexagonal groups:')==False:
        lspl=line.split()
        if lspl[0].strip() in symbsSC:
            print "Check %s, already  in..."%lspl
        symbsSC[lspl[0].strip()]=lspl[1].strip()
        line=fin.readline()
    line=fin.readline()
    print len(symbsSC)
    print symbsSC
    print "* Hexagonal symbols"
    while line.startswith('----------')==False:
        lspl=line.split()
        if lspl[0].strip() in symbsSChex:
            print "Check %s, already  in..."%lspl
        symbsSChex[lspl[0].strip()]=lspl[1].strip()
        line=fin.readline()
    print len(symbsSChex)
    print symbsSChex
    inOP=False
    while line:
        if line.startswith('BNS:'):
            bla,BNSno,BNSsym,bla,OGnum,OGsym=line.split()
            print "* ",BNSno,BNSsym
        if (line.startswith('Operators:')) or (line.startswith('Operators (BNS):')):
            line=line.split(':')[1]
            OP=[]
            num=0
            inOP=True
        if (inOP==True) and (line.startswith('Wyckoff positions')):
            centering=[]
            if line.split(':')[1].strip()!='':
                #check this
                ts=line.split(':')[1].strip().split(' ')
                for trans in ts:
                    symb='1 '
                    t=trans.split(',')
                    #print "xxx:",ts,trans,t
                    for i in range(3):
                        symb=symb+['x','y','z'][i]
                        c=t[i].strip('(').strip(')+')
                        if c.strip()!='0':
                            symb=symb+'+'+c.strip()
                        symb=symb+','
                    centering.append(symb+'+1')
                    #print centering
            msg=MSG()
            for gsym in OP:
                msg.symbols.append(gsym)
                g=msg.symbol2matrix(gsym)
                msg.elements.append(msg.symbol2matrix(gsym))
                for gc in centering:
                    msg.elements.append(msg.mult(msg.symbol2matrix(gsym),msg.symbol2matrix(gc)))
                    msg.symbols.append(msg.matrix2symbol(msg.mult(msg.symbol2matrix(gsym),msg.symbol2matrix(gc))))
            if (check_group==True) and (msg.is_group()==False):
                print "Check %s:"%BNSno
                print msg.symbols
                print msg.elements
                return {}
            MagSpaceGroups[BNSno]=msg
            inOP=False
        if inOP==True:
            symbs=line.split()
            for symb in symbs:
                num=num+1
                if "'" in symb:
                    ti=-1
                else:
                    ti=1
                rotSC=symb.split('|')[0].strip('(')
                if int(BNSno.split('.')[0]) in range(143,195):
                    rot=symbsSChex[rotSC].split(',')
                else:
                    rot=symbsSC[rotSC].split(',')
                trans=symb.split('|')[1].split(')')[0].split(',')
                #print rot,trans,ti
                symbol='%d '%num
                for i in range(3):
                    if trans[i]=='0':
                        symbol=symbol+rot[i]+','
                    else:
                        symbol=symbol+rot[i]+'+'+trans[i]+','
                symbol=symbol+'%+d'%ti
                #print symbol
                OP.append(symbol)
        line=fin.readline()
    return MagSpaceGroups


def get_spin_rotation(M1,M2ini,eps=0.3):
    """get_spin_rotation(): returns a rotation matrix R which rotates the magnetic moments of spin configuration
    M2 into spin configuration M1 (if R*M2-M1<eps)
    test it!!!
    """
    M2=M2ini
    fixedvec=[]
    print "initial vectors:",M1,M2
    if len(M1)!=len(M2):
        return True
    #find magnetization vectors with large m1 x m2 and align them /todo: mirror symm?
    maxcr=0.0
    imaxcr=-1
    for i  in range(len(M1)):
        m1=M1[i]
        m2=M2[i]
        dp=np.dot(m1,m2)
        n=np.cross(m1,m2)
        if norm(n)>maxcr:
            maxcr=norm(n)
            imaxcr=i
            print imaxcr,maxcr
        elif (dp>1) and (norm(n)<eps):
            fixedvec.append(m1/norm(m1))
    if maxcr>eps:
        m1=M1[imaxcr]
        m2=M2[imaxcr]
        n=np.cross(m1,m2)/maxcr
        (ux,uy,uz)=n
        theta0=-math.acos(np.dot(m1,m2)/norm(m1)/norm(m2))
        found=False
        for x in [0,math.pi]: # acos gives angle +-pi
            theta=theta0+x
            Rn=np.array([[math.cos(theta)+ux*ux*(1-math.cos(theta)), ux*uy*(1-math.cos(theta))-uz*math.sin(theta), ux*uz*(1-math.cos(theta))+uy*math.sin(theta)],
                         [uy*ux*(1-math.cos(theta))+uz*math.sin(theta), math.cos(theta)+uy*uy*(1-math.cos(theta)), uy*uz*(1-math.cos(theta))-ux*math.sin(theta)],
                         [uz*ux*(1-math.cos(theta))-uy*math.sin(theta), uz*uy*(1-math.cos(theta))+ux*math.sin(theta), math.cos(theta)+uz*uz*(1-math.cos(theta))]])
            m2rot=np.dot(Rn,m2)
            if (norm(np.cross(m2rot,m1))<eps) and (np.dot(m2rot,m1)>0):
                print "THETA1",theta/math.pi,"pi",n," axis",m2rot,np.cross(m2rot,m1),np.dot(m2rot,m1)
                found=True
                fixedvec.append(m1/norm(m1))
                break
        if found==False:
            print "Ooops - check algorithm"
            bla
        #rotate all vectors in M2
        M2n=[]
        for j in range(len(M2)):
            m2n=np.dot(Rn,M2[j])
            M2n.append(m2n)
        M2=M2n
        print "after first rotation:",M2,M1
    # second rotation around fixed vector
    if fixedvec==[]:
        print "Ooops - check algorithm"
        bla
    #find magnetization vectors with large m1 x m2rot and rotate around fixed axis /todo: mirror symm?
    maxcr=0.0
    imaxcr=-1
    for i  in range(len(M1)):
        m1=M1[i]
        m2=M2[i]
        n=np.cross(m1,m2)
        if norm(n)>maxcr:
            maxcr=norm(n)
            imaxcr=i
            print "2nd rotation",imaxcr,maxcr
    if maxcr>eps:
        n=fixedvec[0]
        (ux,uy,uz)=n
        m1=M1[imaxcr]
        m2=M2[imaxcr]
        v1=m1-n*np.dot(n,m1)
        v2=m2-n*np.dot(n,m2)
        theta0=-math.acos(np.dot(v1,v2)/norm(v1)/norm(v2)) #check sign!!
        for theta in [theta0,-theta0,theta0+math.pi,-theta0+math.pi]: # acos gives angle +-pi, check if all make sense
            #theta=theta0+x
            Rn=np.array([[math.cos(theta)+ux*ux*(1-math.cos(theta)), ux*uy*(1-math.cos(theta))-uz*math.sin(theta), ux*uz*(1-math.cos(theta))+uy*math.sin(theta)],
                         [uy*ux*(1-math.cos(theta))+uz*math.sin(theta), math.cos(theta)+uy*uy*(1-math.cos(theta)), uy*uz*(1-math.cos(theta))-ux*math.sin(theta)],
                         [uz*ux*(1-math.cos(theta))-uy*math.sin(theta), uz*uy*(1-math.cos(theta))+ux*math.sin(theta), math.cos(theta)+uz*uz*(1-math.cos(theta))]])
            m2rot=np.dot(Rn,m2)
            print "THETA2",theta/math.pi,"pi",n," axis",m2rot,np.cross(m2rot,m1),np.dot(m2rot,m1)
            if (norm(np.cross(m2rot,m1))<eps) and (np.dot(m2rot,m1)>0):
                print "THETA2",theta/math.pi,"pi",n," axis",m2rot,np.cross(m2rot,m1),np.dot(m2rot,m1)
                found=True
                fixedvec.append(n)
                break
        #rotate all vectors in M2
        M2n=[]
        for j in range(len(M2)):
            m2n=np.dot(Rn,M2[j])
            M2n.append(m2n)
        M2=M2n
        print "after second rotation:",M2,M1
    # check if rotated spin structure is same as M1
    differ=False
    for j in range(len(M1)):
        if norm(np.cross(M1[j],M2[j]))>eps:
            differ=True
    print "spin configurations differ:",differ
    return differ
    
    R=np.array([[1,0,0],[0,1,0],[0,0,1]])
    trafos=[]
    for i in range(len(M1)):
        m1=M1[i]
        m2=M2[i]
        if trafos==[]:
            n=np.cross(m1,m2)
        else:
            n=M1[0] #np.cross(trafos[0],m2)
        print "XXX", n,trafos
        if norm(n)<eps:
                print "id:",R,m1,m2
                continue
        (ux,uy,uz)=n/norm(n)
        print "ccc",n/norm(n),m1,m2,np.dot(n,m2)
        theta=-math.acos(np.dot(m1,m2)/norm(m1)/norm(m2))
        print "BBBB",theta
        Rn=np.array([[math.cos(theta)+ux*ux*(1-math.cos(theta)), ux*uy*(1-math.cos(theta))-uz*math.sin(theta), ux*uz*(1-math.cos(theta))+uy*math.sin(theta)],
                     [uy*ux*(1-math.cos(theta))+uz*math.sin(theta), math.cos(theta)+uy*uy*(1-math.cos(theta)), uy*uz*(1-math.cos(theta))-ux*math.sin(theta)],
                     [uz*ux*(1-math.cos(theta))-uy*math.sin(theta), uz*uy*(1-math.cos(theta))+ux*math.sin(theta), math.cos(theta)+uz*uz*(1-math.cos(theta))]])
        if trafos!=[]:
            n=n/norm(n)
            for k in range(-10,10):
                theta=3.14*k/20
                Rn=np.array([[math.cos(theta)+ux*ux*(1-math.cos(theta)), ux*uy*(1-math.cos(theta))-uz*math.sin(theta), ux*uz*(1-math.cos(theta))+uy*math.sin(theta)],
                     [uy*ux*(1-math.cos(theta))+uz*math.sin(theta), math.cos(theta)+uy*uy*(1-math.cos(theta)), uy*uz*(1-math.cos(theta))-ux*math.sin(theta)],
                     [uz*ux*(1-math.cos(theta))-uy*math.sin(theta), uz*uy*(1-math.cos(theta))+ux*math.sin(theta), math.cos(theta)+uz*uz*(1-math.cos(theta))]])
                #print theta,np.dot(Rn,m2)
            v1=m1-n*np.dot(n,m1)
            v2=m2-n*np.dot(n,m2)
            theta=-math.acos(np.dot(v1,v2)/norm(v1)/norm(v2)) #check sign!!
            print "BBBB",v1,v2,theta,n
            Rn=np.array([[math.cos(theta)+ux*ux*(1-math.cos(theta)), ux*uy*(1-math.cos(theta))-uz*math.sin(theta), ux*uz*(1-math.cos(theta))+uy*math.sin(theta)],
                     [uy*ux*(1-math.cos(theta))+uz*math.sin(theta), math.cos(theta)+uy*uy*(1-math.cos(theta)), uy*uz*(1-math.cos(theta))-ux*math.sin(theta)],
                     [uz*ux*(1-math.cos(theta))-uy*math.sin(theta), uz*uy*(1-math.cos(theta))+ux*math.sin(theta), math.cos(theta)+uz*uz*(1-math.cos(theta))]])
        #Rn=np.array([[-0.208686 ,  -0.78937749 , 0.57735027], [ 0.78796396 , 0.21396136 , 0.57735027], [-0.57927796 , 0.57541612 , 0.57735027]])

        trafos.append(n/norm(n))
        M2n=[]
        for j in range(len(M2)):
            m2n=np.dot(Rn,M2[j])
            M2n.append(m2n)
        print "ZZZ",M2n
        M2=M2n
        if i==1:
            differ=False
            for j in range(len(M1)):
                if norm(np.cross(M1[j],M2[j]))>eps:
                    differ=True
            print "spin configurations differ:",differ
            return differ
    for m1,m2 in zip(M1,M2):
        n=np.cross(m1,m2)
        if norm(n)<eps:
            print "id:",R,m1,m2
            continue
        (ux,uy,uz)=n/norm(n)
        theta=-math.acos(np.dot(m1,m2)/norm(m1)/norm(m2))
        Rn=np.array([[math.cos(theta)+ux*ux*(1-math.cos(theta)), ux*uy*(1-math.cos(theta))-uz*math.sin(theta), ux*uz*(1-math.cos(theta))+uy*math.sin(theta)],
                     [uy*ux*(1-math.cos(theta))+uz*math.sin(theta), math.cos(theta)+uy*uy*(1-math.cos(theta)), uy*uz*(1-math.cos(theta))-ux*math.sin(theta)],
                     [uz*ux*(1-math.cos(theta))-uy*math.sin(theta), uz*uy*(1-math.cos(theta))+ux*math.sin(theta), math.cos(theta)+uz*uz*(1-math.cos(theta))]])
        R=np.dot(R,Rn)
        for i in range(len(M2)):
            print (ux,uy,uz),np.dot(Rn,M2[i])
        #print np.dot(M2,R)
        #print n,theta,R
        #print M1,M2


def check_ao_equivalence(ao1, ao2, symprec=5e-2):
    lattice, scaled_positions, numbers=spglib.find_primitive(ao1, symprec=symprec)
    aop1=Atoms(numbers=numbers,cell=lattice,scaled_positions=scaled_positions,pbc=True)
    lattice, scaled_positions, numbers=spglib.find_primitive(ao2, symprec=symprec)
    aop2=Atoms(numbers=numbers,cell=lattice,scaled_positions=scaled_positions,pbc=True)
    if len(aop1)!=len(aop2):
        print "check_ao_equivalence: primitive cells diifer in size:",len(aop1),"!=",len(aop2)
        return False
    else:
        print "check_ao_equivalence: atoms objects seem to be equivalent"
    return True
