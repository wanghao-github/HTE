import os
import commands

def formula2composition(formula_sum):
    """returns dictionary with elements and occupation numbers
    """
    comp={}
    cfsum=''
    for i in range(len(formula_sum)):
        if (formula_sum[i].isalpha()) or (formula_sum[i].isdigit()):
            cfsum=cfsum+formula_sum[i]
    #print cfsum
    i=0
    n=len(cfsum)
    while (i<n):
        if (i+1<n) and (cfsum[i+1].isalpha()) and (cfsum[i+1]==cfsum[i+1].lower()):
            el=cfsum[i:i+2]
            i=i+2
        else:
            el=cfsum[i:i+1]
            i=i+1
        #print el
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
    
def composition2formula(chem_composition, latex=False, gle=False, sort='alpha'):
    formula=''
    if sort=='num':
        for sym,comp in sorted(chem_composition.items(), key=lambda x: x[1], reverse=True):
            if chem_composition[sym]==1:
                formula=formula+sym
            elif latex==True:
                formula=formula+sym+'$_{'+str(chem_composition[sym])+'}$'
            elif gle==True:
                formula=formula+sym+'_{'+str(chem_composition[sym])+'}'
            else:
                formula=formula+sym+str(chem_composition[sym])
        return formula
    elif sort=='alpha_pretty':
        for sym in sorted(chem_composition.keys()):
            if chem_composition[sym]==1:
                formula=formula+sym            
            elif latex==True:
                formula=formula+sym+'$_{'+str(chem_composition[sym])+'}$'
            elif gle==True:
                formula=formula+sym+'_{'+str(chem_composition[sym])+'}'
            else:
                formula=formula+sym+str(chem_composition[sym])
        return formula
    elif isinstance(sort,list):
        comp=sorted(chem_composition.keys())
        for sym in sort:
            if sym in comp:
                if chem_composition[sym]==1:
                    formula=formula+sym            
                elif latex==True:
                    formula=formula+sym+'$_{'+str(chem_composition[sym])+'}$'
                elif gle==True:
                    formula=formula+sym+'_{'+str(chem_composition[sym])+'}'
                else:
                    formula=formula+sym+str(chem_composition[sym])
        for sym in comp:
            if not (sym in sort):
                if chem_composition[sym]==1:
                    formula=formula+sym            
                elif latex==True:
                    formula=formula+sym+'$_{'+str(chem_composition[sym])+'}$'
                else:
                    formula=formula+sym+str(chem_composition[sym])
        return formula
    for sym in sorted(chem_composition.keys()):
        if latex==True:
            formula=formula+sym+'$_{'+str(chem_composition[sym])+'}$'
        elif gle==True:
            formula=formula+sym+'_{'+str(chem_composition[sym])+'}'
        else:
            formula=formula+sym+str(chem_composition[sym])
    return formula

    
def store_files(source_dir,target_dir, files={'bzip2':['OUTCAR','out'],'cp':['POSCAR.INI','POSCAR','CONTCAR','KPOINTS','INCAR','hte_propdict.txt']}):
    print "store_files:",source_dir,target_dir
    if not os.path.isdir(target_dir):
        exitcode, out = commands.getstatusoutput("mkdir -p %s"%target_dir)
    for method in files:
        if method=='bzip2':
            for sfname in files['bzip2']:
                fname=os.path.join(source_dir,sfname)
                if os.path.isfile(fname):
                    target_fname=os.path.join(target_dir,sfname)+".bz2"
                    exitcode, out = commands.getstatusoutput("bzip2 -c %s > %s"%(fname,target_fname))
        if method=='gzip':
            for sfname in files['gzip']:
                fname=os.path.join(source_dir,sfname)
                if os.path.isfile(fname):
                    target_fname=os.path.join(target_dir,sfname)+".gz"
                    exitcode, out = commands.getstatusoutput("gzip -c %s > %s"%(fname,target_fname))
        if method=='cp':
            for sfname in files['cp']:
                fname=os.path.join(source_dir,sfname)
                if os.path.isfile(fname):
                    target_fname=os.path.join(target_dir,sfname)
                    exitcode, out = commands.getstatusoutput("cp %s %s"%(fname,target_fname))
    #TODO: include some checks



def value2string(value,latex=False, format_string="%.3f"):
    if (value==None) and (latex==True):
        string="--"
    elif (latex==True):
        string=format_string%value
    elif  isinstance(value,float):
        string="%.8f"%value
    else:
        string=str(value)
    return string+" "


def string2prop_dict(line, float_props=['magnetic_moments','magnetic_moment','energy','enthalpy','cell','scaled_positions','initial_magnetic_moments','band_energy','band_gap','positions'],list_props=['chemical_symbols','cell','scaled_positions','positions','magnetic_moments']):
    prop,val=line.strip('\n').split('=')
    prop=prop.strip()
    if (',' in val) or (prop in list_props):
        val2=val.split(',')
        val=[]
        for i in range(len(val2)):
            if ';' in val2[i]:
                val3=val2[i].split(';')
                valx=[]
                print"check_point225 ,val3",val3
                for j in range(len(val3)):
                    if prop in float_props:
                        valx.append(float(val3[j]))
                    else:
                        valx.append(val3[j].strip())
                val.append(valx)
            else:
                if prop in float_props:
                    val.append(float(val2[i]))
                else:
                    val.append(val2[i].strip())
    else:
        if prop in float_props:
            val=float(val)
        else:
            val=val.strip()
    return prop,val

def split_cif(ciffile):
    cif_blocks=[]
    if os.path.isfile(ciffile):
        infile=open(ciffile,'r')
        cif_block=''
        uid=None
        line=infile.readline()
        while line:
            if (line.startswith('# End of data set')) and (uid!=None):
                cif_block=cif_block+line
                cif_blocks.append((uid,cif_block))
                cif_block=''
                uid=None
            elif (line.startswith('data_')):
                if uid!=None:
                    cif_blocks.append((uid,cif_block))
                    cif_block=line
                else:
                    cif_block=cif_block+line
                    uid=line.strip()
            else:
                cif_block=cif_block+line
            line=infile.readline()
        if uid!=None:
            cif_blocks.append((uid,cif_block))
    return cif_blocks
