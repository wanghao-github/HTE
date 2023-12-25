from numpy  import *
from numpy.linalg import *
from ase.visualize import view
from ase.calculators.vasp import *
from hte import *
try:
 from pyspglib import spglib
 has_spglib=True
except:
 has_spglib=False

######################################################################          
def get_vasp_bandstructure(pathname='./', filename_eigenval='EIGENVAL', pathname_sc=None, filename_outcar_sc='OUTCAR', silent=True):
    # return a dictionary with band structure information, keywords are:
    # 'unit_E': energy unit ('eV', 'Hartree', 'Rydberg')
    # 'E_Fermi': Fermi energy in above unit
    # 'type': 'lines' or 'BZ'
    # 'nkp': number of k-points
    # 'nband': number of bands
    # 'nspin': number of spins
    # 'e_n_k': list with actual band structure;
    # e_n_k[ikp] is a dictionary with 'kpoint' (k-vector) and 'energies' (e_n(k)) for k-point #ikp

    bandstructure={}
    parentdir=os.getcwd()
    # get the Fermi energy from a selfconsistent OUTCAR file (not from bs calculation!)
    if pathname_sc==None:
        pathname_sc=pathname
    if filename_outcar_sc!=None:
        outcarfile=os.path.join(pathname_sc,filename_outcar_sc)
        if os.path.isfile(outcarfile):
            os.chdir(pathname_sc)
            calc=Vasp()
            ef=calc.read_fermi()
            nelect=calc.read_number_of_electrons()
            os.chdir(parentdir)
        else:
            if not silent:
                print 'get_vasp_bandstructure(): could not find OUTCAR file <%s>' % outcarfile
            return None
    else:
        ef=None
    bandstructure['unit_E']='eV'
    bandstructure['E_Fermi']=ef    
    bandstructure['n_electrons']=nelect    
    # evaluate eigenval file
    fin=None
    if (pathname!=None) and (filename_eigenval!=None):
        eigenvalfile=os.path.join(pathname,filename_eigenval)
        if os.path.isfile(eigenvalfile):
            fin=open(eigenvalfile,"r")
    else:
        eigenvalfile=None
    if fin==None:
        print 'get_vasp_bandstructure(): could not find EIGENVAL file <%s>' % eigenvalfile
        return None
    #read header
    line=fin.readline() #??/??/??/nspin
    sline=line.split()
    nspin=int(sline[3])
    bandstructure['nspin']=nspin
    line=fin.readline() #unknown data
    line=fin.readline() #unknown data
    line=fin.readline() #unknown data
    line=fin.readline() #unknown data
    line=fin.readline() #??/nkp/nband
    sline=line.split()
    nband=int(sline[2])
    nkp=int(sline[1])
    bandstructure['nkp']=nkp
    bandstructure['nband']=nband
    bandstructure['e_n_k']=[]
    for ikp in range(nkp):
        ek={}
        line=fin.readline() #empty line
        line=fin.readline() #kx ky kz wk
        sline=line.split()
        ek['kpoint']=[float(sline[0]),float(sline[1]),float(sline[2])]
        energies_up=[]
        energies_dn=[]
        for iband in range(nband):
            line=fin.readline() #no. energy
            sline=line.split()
            energies_up.append(float(sline[1]))
            if nspin==2:
                energies_dn.append(float(sline[2]))
        ek['energies_up']=energies_up
        if nspin==2:
            ek['energies_dn']=energies_dn
        bandstructure['e_n_k'].append(ek)
    fin.close()
    return bandstructure


def get_wien_bandstructure(pathname=None, filename='wien2k.energy', scf_file='wien2k.scf', skiplines=4):
    # read bandstructure from w2k calculation
    bandstructure={}
    bandfile=filename
    nband_xny=None
    ef=None
    nelect=None
    if (scf_file!=None) and (os.path.isfile(scf_file)):
        #get Fermi energy and number of electrons
        fin=open(scf_file,"r")
        line=fin.readline()
        while line:
            if line.startswith(':FER'):
                ef=float(line[38:])
                #print 'E_f=',ef
            if line.startswith(':NOE'):
                nelect=float(line[38:])
                #print 'nel=',nelect
            line=fin.readline()
        fin.close()
    bandstructure['E_Fermi']=ef    
    bandstructure['n_electrons']=nelect    
    if (pathname!=None) and (filename!=None):
        bandfile=os.path.join(pathname,filename)
    if (bandfile!=None) and (os.path.isfile(bandfile)):
        fin=open(bandfile,"r")
        #skip first header lines
        for i in range(skiplines):
            line=fin.readline()
        nkp=0
        bandstructure['e_n_k']=[]
        # line with k-point and nband
        line=fin.readline()
        while (line):
            ek={}
            ek['kpoint']=[float(line[0:19]),float(line[19:38]),float(line[38:57])]
            kname=line[57:67]
            nband=int(line[73:79])
            if (nband_xny==None) or (nband<nband_xny):
                nband_xny=nband
            #print ek,kname,nband
            #sline=line.split()
            #print sline
            #nkp=int(sline[4])
            #nband=int(sline[5])
            #ek['kpoint']=[float(sline[0]),float(sline[1]),float(sline[2])]
            energiestmp=[]
            #read bandenergies
            for iband in range(nband):
                line=fin.readline()
                sline=line.lower().replace('d','E').split()
                #print sline
                energiestmp.append(float(sline[1]))
            ek['energies_up']=energiestmp
            ek['energies_dn']=None
            bandstructure['e_n_k'].append(ek)
            #print 'XXX',nkp, bandstructure['e_n_k'][nkp]
            nkp=nkp+1
            # line with k-point and nband
            line=fin.readline()
        fin.close()
        bandstructure['nkp']=nkp
        bandstructure['nband']=nband_xny #may differ between k-points, take max value with energies for all k-points
        bandstructure['nspin']=1 #for the moment only nonspinpolarized
        bandstructure['unit_E']='Ry'
        #bandstructure['E_Fermi']=None
        return bandstructure
    else:
        print 'get_wien_bandstructure: Could not open <%s> in <%s>, abort!' %(filename,pathname)
        return None


        
def get_fplo_bandstructure(pathname=None, filename='+band'):
    bandstructure={}
    # change to directory if pathname is given, otherwise look in current directory
    parentdir=os.getcwd()
    if pathname!=None:
        bandfile=os.path.join(pathname,filename)
    else:
        bandfile=filename
    if os.path.isfile(bandfile):
        print 'Reading band structure from file <%s>...' %bandfile
        bandstructure['unit_E']='eV'
        bandstructure['E_Fermi']=0.0
        fin=open(bandfile,"r")
        #read header line
        line=fin.readline()
        sline=line.split()
        nsite=int(sline[1])
        nkp=int(sline[3])
        nband=int(sline[5])
        nspin=int(sline[6])
        first_band=int(sline[7])
        last_band=int(sline[8])
        bandstructure['nkp']=nkp
        bandstructure['nband']=nband
        bandstructure['e_n_k']=[]
        bandstructure['nspin']=nspin
        print nkp,nband,nspin
        for ikp in range(nkp):
            #line with k-point
            ek={}
            line=fin.readline()
            sline=line.split()
            ek['kpoint']=[float(sline[1]),float(sline[2]),float(sline[3])]
            #line(s) with band energies for spin up and down
            for ispin in range(nspin):
                energiestmp=[]
                line=fin.readline()
                sline=line.split()
                for iband in range(nband):
                    energiestmp.append(float(sline[iband+1]))
                if ispin==0:
                    ek['energies_up']=energiestmp
                    print ikp
                else:
                    ek['energies_dn']=energiestmp
            bandstructure['e_n_k'].append(ek)
        fin.close()
    else:
       print 'Could not open %s, nothing done!',filename
    os.chdir(parentdir)
    return bandstructure

def get_conversion(unitin, unitout):
    # return conversion factor between unitin and unitout
    # from physics.nist.gov/constants
    # (units must be of same type (energy, length,..)
    
    # convert unitin to eV (energies) or \AA (length)
    if unitin=='eV':
        factorin=1.0
    elif unitin=='Hartree':
        factorin=27.2113834
    elif (unitin=='Rydberg') or (unitin=='Ry'):
        factorin=27.2113834/2.0
    elif unitin=='a.u.':
        factorin=0.5291772083
    elif unitin=='AA':
        factorin=1.0
    else:
        print 'Unknown unit ',unitin
        return None
    # convert unitout to eV (energies) or \AA (length)
    if unitout=='eV':
        factorout=1.0
    elif unitout=='Hartree':
        factorout=27.2113834
    elif (unitout=='Rydberg') or (unitout=='Ry'):
        factorout=27.2113834/2.0
    elif unitout=='a.u.':
        factorout=0.5291772083
    elif unitout=='AA':
        factorout=1.0
    else:
        print 'Unknown unit ',unitout
        return None
    return factorin/factorout

def write_bandstructure_xny(bandstructure, xscale=1.0, unit=None, E_Fermi_zero=True, shifte=0.0, filename='bands.dat', first_band=None, last_band=None, ktrafo=None):
    # output of bandstructure in a format which can be used by xmgrace and other plotting tools
    # input: bandstructure as dictionary like in get_vasp_bandstructure

    eb=bandstructure['e_n_k']
    yscale=1.0
    if (unit!=None) and (unit!=bandstructure['unit_E']):
        yscale=get_conversion(bandstructure['unit_E'],unit)
        if yscale==None:
            print 'write_bandstructure_xny: unit conversion failed, nothing done!'
            return False
    if (first_band==None):
        first_band=1
    if (last_band==None):
        last_band=bandstructure['nband']

    fout=open(filename,"w")
    nsite=1 #for the moment
    line='# %3d%12.7f%9d%9d%9d%9d%9d%9d' % (nsite,bandstructure['E_Fermi']*yscale,bandstructure['nkp'],0,last_band-first_band+1,bandstructure['nspin'],first_band,last_band)
    fout.write(line+'\n')
    x=0.0
    for ikp in range(len(eb)):
        # print comment line with k-point
        line='#'
        enk=eb[ikp]
        # determine x (distance from previous k-points)
        if ktrafo==None:
            kp=enk['kpoint']
        else:
            kp=[0.0,0.0,0.0]
            for i in range(3):
                for j in range(3):
                    kp[i]=kp[i]+ktrafo[j][i]*enk['kpoint'][j]
                    #kp[i]=kp[i]+ktrafo[i][j]*enk['kpoint'][j]
        dx=0.0
        for i in range(3):
            line=line+' '+str(kp[i])
            if ikp>0:
                dx=dx+(kp[i]-kpold[i])*(kp[i]-kpold[i])
        fout.write(line+'\n')
        x=x+sqrt(dx)*xscale
        # print line with x e1 e2 ... for spin up
        line=str(x)
        #for i in range(len(enk['energies_up'])):
        for i in range(first_band-1,last_band):
            e_up=enk['energies_up'][i]
            if (E_Fermi_zero==True):
                e_up=e_up-bandstructure['E_Fermi']
            e_up=e_up*yscale+shifte
            line=line+' %15.7f' % e_up
        fout.write(line+'\n')
        if bandstructure['nspin']==2:
            # print line with x e1 e2 ... for spin up
            line=str(x)
            for i in range(len(enk['energies_dn'])):
                e_dn=enk['energies_dn'][i]
                if (E_Fermi_zero==True):
                    e_dn=e_dn-bandstructure['E_Fermi']
                e_dn=e_dn*yscale+shifte
                line=line+' %15.7f' % e_dn
                fout.write(line+'\n')
        kpold=kp
    fout.close()
    return True

def write_bandstructure_boltztrap(bandstructure, unit='Ry', yscale=None, E_Fermi_zero=True, pathname='./', filename="energies.boltztrap", ktrafo=None, runBoltz=False):
    # output of bandstructure in a format which can be used by boltztrap code
    # input: bandstrcuture as dictionary like in get_vasp_bandstructure
    # default is unit conversion to Ry; if yscale is given this is taken as conversion factor
    if bandstructure['nspin']!=1:
        print 'write_bandstructure_boltztrap: No idea what to do for nspin=%d' % bandstructure['nspin']
        return False
    if  (E_Fermi_zero==True) and (bandstructure['E_Fermi']==None):
        print 'write_bandstructure_boltztrap: Fermi energy not known, nothing done'
        return False
    #    
    eb=bandstructure['e_n_k']
    if yscale==None:
        if (unit!=None) and (bandstructure['unit_E']!=None):
            yscale=get_conversion(bandstructure['unit_E'],unit)
        if yscale==None:
            print 'write_bandstructure_boltztrap: unit conversion failed, nothing done!'
            return False
    outfile=os.path.join(pathname,filename)
    fout=open(outfile,"w")
    fout.write('HTE output'+'\n') # title
    fout.write(str(len(eb))+'\n') # no. of k-points
    print '***************'
    print ktrafo
    for ikp in range(len(eb)):
        enk=eb[ikp]
        #print ikp,enk['kpoint']
        if ktrafo==None:
            kp=enk['kpoint']
        else:
            kp=[0.0,0.0,0.0]
            for i in range(3):
                for j in range(3):
                    kp[i]=kp[i]+ktrafo[j][i]*enk['kpoint'][j]
                    #kp[i]=kp[i]+ktrafo[i][j]*enk['kpoint'][j]
        #print ikp,enk['kpoint'],' -> ',kp
        #line=''
        #for i in range(3):
        #    line=line+' '+str(kp[i])
        #line=line+' '+str(len(enk['energies_up']))+'\n'
        #fout.write(line)
        fout.write("%12.8f %12.8f %12.8f %d\n" %(kp[0],kp[1],kp[2],len(enk['energies_up'])))
        line=''
        for i in range(len(enk['energies_up'])):
            e_up=enk['energies_up'][i]
            if (E_Fermi_zero==True):
                e_up=e_up-bandstructure['E_Fermi']
            e_up=e_up*yscale
            fout.write("%18.8f\n" %e_up)
    fout.close()
    if runBoltz==True:
        parentdir=os.getcwd()
        os.chdir(pathname)
        exitcode, out=commands.getstatusoutput('/home/users/opahlivs/bin/BoltzTraP BoltzTraP.def')
        os.chdir(parentdir)
        print exitcode, out
    return 0
    
def write_structure_boltztrap(ao, pathname='./', filename="hte.struct"):
    if not os.path.isdir(pathname):
        os.path.mkdir(pathname)
    fname=os.path.join(pathname,filename)
    fout=open(fname,"w")
    fout.write('HTE output'+'\n') # title
    latt=ao.get_cell()/0.5291772083
    for i in range(3):
        line=''
        for j in range(3):
            line=line+"%12.5f"%latt[i][j]
        fout.write(line+'\n')
    krot=get_kspace_operations(ao)
    fout.write(str(len(krot))+'\n')
    for iop in range(len(krot)):
        for i in range(3):
            for j in range(3):
                fout.write(str(krot[iop][i][j])+' ')
            fout.write('\n')
    fout.close()
    return 0


def export_kpoints_to_fplo(bandstructure, ktrafo=None, filename='=.kp'):
    fout=open(filename,"w")
    fout.write(str(bandstructure['nkp'])+' f 0 0\n')
    eb=bandstructure['e_n_k']
    for ikp in range(len(eb)):
        enk=eb[ikp]
        if ktrafo==None:
            kp=enk['kpoint']
        else:
            kp=[0.0,0.0,0.0]
            #k=k_W*ktrafo
            for i in range(3):
                for j in range(3):
                    kp[i]=kp[i]+ktrafo[j][i]*enk['kpoint'][j]
        fout.write(str(kp[0])+' '+str(kp[1])+' '+str(kp[2])+'\n')
    fout.close()
    return True

def export_kpoints_to_vasp(bandstructure, ktrafo=None, filename='KPOINTS'):
    fout=open(filename,"w")
    fout.write('HTE export_kpoints_to_vasp\n'+str(bandstructure['nkp'])+'\nreciprocal\n')
    eb=bandstructure['e_n_k']
    for ikp in range(len(eb)):
        enk=eb[ikp]
        if ktrafo==None:
            kp=enk['kpoint']
        else:
            kp=[0.0,0.0,0.0]
            #k=k_W*ktrafo
            for i in range(3):
                for j in range(3):
                    kp[i]=kp[i]+ktrafo[j][i]*enk['kpoint'][j]
        fout.write(str(kp[0])+' '+str(kp[1])+' '+str(kp[2])+' 1\n')
    fout.close()
    return True

def write_bandstructure_wien(bandstructure, unit='Ry', yscale=None, E_Fermi_zero=True, pathname='./', filename="energy.wien", ktrafo=None, headerfile=None, runBoltz=False):
    eb=bandstructure['e_n_k']
    if yscale==None:
        if (unit!=None) and (bandstructure['unit_E']!=None):
            yscale=get_conversion(bandstructure['unit_E'],unit)
        if yscale==None:
            print 'write_bandstructure_boltztrap: unit conversion failed, nothing done!'
            return False
    outfile=os.path.join(pathname,filename)
    fout=open(outfile,"w")
    #header-to be done
    if (headerfile!=None) and (os.path.isfile(headerfile)):
        fin=open(headerfile,"r")
        lines=fin.readlines()
        fin.close()
        for line in lines:
            print 'YYY',line
            fout.write(line)
    else:
        #to be done
        fout.write('200.30000200.30000200.44500  0.30000  0.30000  0.30000  0.30000  0.30000  0.30000  0.30000  0.30000  0.30000  0.30000\n')
        fout.write('0.30000  0.30000  0.44500999.00000997.00000 -3.73750997.00000999.00000999.00000999.00000999.00000999.00000\n')
        fout.write('200.30000200.30000198.41000  0.30000  0.30000  0.30000  0.30000  0.30000  0.30000  0.30000  0.30000  0.30000  0.30000\n')
        fout.write('0.30000  0.30000 -1.59000999.00000997.00000997.00000  0.30000999.00000999.00000999.00000999.00000999.00000\n')
    for ikp in range(len(eb)):
        enk=eb[ikp]
        if ktrafo==None:
            kp=enk['kpoint']
        else:
            kp=[0.0,0.0,0.0]
            for i in range(3):
                for j in range(3):
                    kp[i]=kp[i]+ktrafo[j][i]*enk['kpoint'][j]
                    #kp[i]=kp[i]+ktrafo[i][j]*enk['kpoint'][j]
        fout.write('%19.12E%19.12E%19.12E%10s%6d%6d  1.0\n' %(kp[0],kp[1],kp[2],str(ikp+1),819,len(enk['energies_up'])))
        for i in range(len(enk['energies_up'])):
            e_up=enk['energies_up'][i]
            if (E_Fermi_zero==True):
                e_up=e_up-bandstructure['E_Fermi']
            e_up=e_up*yscale
            fout.write(str(i)+' '+str(e_up)+'\n')
    fout.close()
    if runBoltz==True:
        parentdir=os.getcwd()
        os.chdir(pathname)
        exitcode, out=commands.getstatusoutput('/home/users/opahlivs/bin/BoltzTraP BoltzTraP.def')
        os.chdir(parentdir)
        print exitcode, out
    return True

def get_transport_boltztrap(pathname='./', filename='*.trace', Tmin=299.0, Tmax=301.0):
    tracefile=os.path.join(pathname,filename)
    if os.path.isfile(tracefile):
        fin=open(tracefile,"r")
        #read header line
        line=fin.readline()
        line=fin.readline()
        transport={}
        transport['E_Fermi']=[]
        transport['T']=[]
        transport['N']=[]
        transport['DOS']=[]
        transport['Seebeck']=[]
        transport['sigotau']=[]
        transport['specific_heat']=[]
        transport['R_H']=[]
        while line:
            property={}
            sline=line.split()
            Ef=float(sline[0]) #*27.2113834/2.0 #use eV
            T=float(sline[1])
            N=float(sline[2])
            DOS=float(sline[3])
            S=float(sline[4])
            sigotau=float(sline[5]) 
            R_H=float(sline[6])
            kap_0=float(sline[7])
            c=float(sline[8])
            chi=float(sline[9])
            if (T>Tmin) and (T<Tmax):
                transport['E_Fermi'].append(Ef)
                transport['T'].append(T)
                transport['N'].append(N)
                transport['DOS'].append(DOS)
                transport['Seebeck'].append(S)
                transport['sigotau'].append(sigotau)
                transport['specific_heat'].append(c)
                transport['R_H'].append(R_H)
            line=fin.readline()
        fin.close()
        return transport
    else:
        print 'Could not open ',tracefile
    return None


def get_sympoints_fplo(pathname='./', outfilename='out'):
    outfile=os.path.join(pathname,outfilename)
    nsymp=0
    kp=[]
    k=[0.0,0.0,0.0]
    if os.path.isfile(outfile):
        fin=open(outfile,"r")
        line=fin.readline()
        while line:
            if 'special_sympoints' in line:
                line=fin.readline()
                #print line
                if '={' in line:
                    line=fin.readline()
                while ((line) and (not ('};' in line))):
                    #print line
                    sline=line.strip().replace('{','').replace('}','').replace('\"','').lstrip(',').split(',')
                    for i in range(3):
                        if '/' in sline[i+1]:
                            kdum=sline[i+1].split('/')
                            k[i]=float(kdum[0])/float(kdum[1])
                        else:
                            k[i]=float(sline[i+1])
                    kpoint={}
                    kpoint['name']=sline[0]
                    kpoint['vector']=array(k)
                    kp.append(kpoint)
                    line=fin.readline()
            line=fin.readline()
        fin.close()
    if len(kp)==0:
        kp=None
    return kp
                  
def get_reciprocal_fplo(pathname='./', outfilename='out',scale=True):
    outfile=os.path.join(pathname,outfilename)
    #print 'checking ',outfile
    rec=[]
    rec_latt=None
    if os.path.isfile(outfile):
        fin=open(outfile,"r")
        line=fin.readline()
        while line:
            if 'lattice constants:' in line:
                sline=line.split()
                alatt=float(sline[2])
                #print alatt
            if 'reciprocial lattice vectors' in line:
                for i in range(3):
                    line=fin.readline()
                    sline=line.split()
                    rv=[]
                    for j in range(3):
                        rv.append(float(sline[j+2]))
                    rec.append(rv)
                #print rec
            line=fin.readline()            
        if len(rec)==3:
            rec_latt=array(rec)
            if (scale):
                rec_latt=alatt*rec_latt
    return rec_latt
                                  

def export_points_file(sympoints, ktrafo=None, pathname='./', filename='+points'):
    fout=open(os.path.join(pathname,filename),"w")
    fout.write("#    %d\n"%len(sympoints))
    x=0.0
    ymax=544.22769
    for ikp in range(len(sympoints)):
        kpoint=sympoints[ikp]
        if ikp>0:
            dk=kpoint['vector']-kp_prev['vector']
            x=x+sqrt(dot(dk.T,dk))
        fout.write("# \' %s \'\n"%kpoint['name'])
        fout.write("%10.5f%12.5f\n"%(x,-ymax))
        fout.write("%10.5f%12.5f\n\n"%(x,ymax))
        kp_prev=kpoint
    fout.close()
    return True
        

def export_sympoints_to_vasp(sympoints, ktrafo=None, pathname='./', filename='KPOINTS'):
    kpointsfile=os.path.join(pathname,filename)
    fout=open(kpointsfile,"w")
    fout.write('HTE export_kpoints_to_vasp\n')
    fout.write('10\nline\nreciprocal\n')
    for ikp in range(len(sympoints)-1):
        kpoint=sympoints[ikp]
        if ktrafo==None:
            kp=kpoint['vector']
        else:
            kp=dot(ktrafo,kpoint['vector'])
        fout.write("%8.5f %8.5f %8.4f %s\n" %(kp[0],kp[1],kp[2],kpoint['name']))
        kpoint=sympoints[ikp+1]
        if ktrafo==None:
            kp=kpoint['vector']
        else:
            kp=dot(ktrafo,kpoint['vector'])
        fout.write("%8.5f %8.5f %8.4f %s\n\n" %(kp[0],kp[1],kp[2],kpoint['name']))
    fout.close()
    return True

def export_sympoints_to_wien(sympoints, ktrafo=None, pathname='./', filename='hte.klist', ndiv=10):
    # status: not ready
    kpointsfile=os.path.join(pathname,filename)
    fout=open(kpointsfile,"w")
    for ikp in range(len(sympoints)-2):
        kpoint=sympoints[ikp]
        if ktrafo==None:
            kp=kpoint['vector']
        else:
            kp=dot(ktrafo,kpoint['vector'])
        fout.write("%10s%10d%10d%10d%5.2f%5.2f8.5f %8.5f %8.4f %s\n" %(kp[0],kp[1],kp[2],kpoint['name']))
        kpoint=sympoints[ikp+1]
        if ktrafo==None:
            kp=kpoint['vector']
        else:
            kp=dot(ktrafo,kpoint['vector'])
        fout.write("%8.5f %8.5f %8.4f %s\n\n" %(kp[0],kp[1],kp[2],kpoint['name']))
    fout.close()
    return True


def export_fplo_pipefile(atoms_object, calcdir='./', pipefile='+pipe', bs_calc=False, xc_vers=5):
    # create a pipe file with structural information for FPLO
    # returns False upon failure, True upon success
    #
    # check if atoms_object.info has sufficient information:
    tags_for_fplo=['_atom_site_type_symbol','spacegroup','_cell_length_a','_cell_length_b','_cell_length_c','_cell_angle_alpha','_cell_angle_beta','_cell_angle_gamma']
    info_ok=True
    for tag in tags_for_fplo:
        if tag in atoms_object.info==False:
            info_ok=False
    if info_ok==False:
        print 'export_fplo_pipefile(): insufficient structure information, nothing done!'
        return False
    nwyckoff=len(atoms_object.info['_atom_site_type_symbol'])
    if os.path.isdir(calcdir):
        outfile=os.path.join(calcdir,pipefile)
        fout=open(outfile,'w')
        fout.write("### pipe file created by HTE\n")
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
        fout.write("@a@%6.2f %6.2f %6.2f\n"%(atoms_object.info['_cell_angle_alpha'], atoms_object.info['_cell_angle_beta'], atoms_object.info['_cell_angle_gamma']))
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
                    print 'ZZZ',j,atom_symbol,atom_site_symbol,' xx '
            fout.write("@%d@%s@%8.4f %8.4f %8.4f\n"%(i+1,atom_symbol,atoms_object.info['_atom_site_fract_x'][i],atoms_object.info['_atom_site_fract_y'][i],atoms_object.info['_atom_site_fract_z'][i]))
        # update and leave symmetry menu
        fout.write("# update and leave symmetry menu\n")
        fout.write("@+@\n")
        fout.write("@x@\n")
        # set Vxc
        fout.write("# set Vxc\n")
        fout.write("@v@\n@%d@\n"%xc_vers)        
        fout.write("@x@\n")
        if bs_calc==True:
            # calculate bandstructure
            fout.write("# calculate bandstructure\n")
            fout.write("@ b@\n@b@+\n@x@\n")
        # leave fedit
        fout.write("# leave fedit\n")
        fout.write("@q@\n")
        fout.write("\n")
        fout.close()
        return True
    print 'export_fplo_pipefile(): ',calcdir, ' not found, nothing done!'
    return False

def check_convergency_fplo(path='./', outfile_name='out'):
    outfile=os.path.join(path,outfile_name)
    if os.path.isfile(outfile):
        exitcode, mes=commands.getstatusoutput("tail -n1 %s"%outfile)
        print mes
        if 'TERMINATION: Finished : SCF calculation' in mes:
            return True
    return False
        
def run_fplo_bandplot(calcdir='./', bandfiles=['+band'], lower_energy_bound=-10.0, upper_energy_bound=5.0, title=None, job_commands=None):
    fedit="fedit9.01-35-x86_64"
    if (job_commands!=None) and ('fedit' in job_commands):
        fedit=job_commands['fedit']
    parentdir=os.getcwd()
    if os.path.isdir(calcdir):
        os.chdir(calcdir)
        fout=open('+pipe_bs','w')
        # set energy bounds
        fout.write("# set energy bounds\n")
        fout.write("@l@%5.2f\n@u@%5.2f\n"%(lower_energy_bound,upper_energy_bound))
        # set band files to plot
        fout.write("# set band files to plot\n")
        numfiles=len(bandfiles)
        fout.write("@n@%d\n"%numfiles)
        for i in range(numfiles):
            fout.write("@%d@%s@%s\n"%(i+1,bandfiles[i],bandfiles[i]))
        if title!=None:
            fout.write("# set title\n")
            fout.write("@t@%s\n"%title)
        # run band plot and quit
        fout.write("# run band plot\n")
        fout.write("@+@\n@q@\n")
        fout.close()
        commands.getstatusoutput("%s -bandplot -pipe < +pipe_bs"%fedit)
        os.chdir(parentdir)
        return True
    print "run_fplo_bandplot() failed"
    return False


def read_genstruct_boltztrap(filename):
    latt=None
    kops=None
    if os.path.isfile(filename):
        fin=open(filename,"r")
        line=fin.readline()
        #lattice vectors
        latt=[]
        for i in range(3):
            sline=fin.readline().split()
            vec=[]
            for j in range(3):
                vec.append(float(sline[j]))
            latt.append(vec)
        print latt
        #kops
        nkops=int(fin.readline())
        kops=[]
        for iops in range(nkops):
            mat=[]
            for i in range(3):
                sline=fin.readline().split()
                vec=[]
                for j in range(3):
                    vec.append(float(sline[j]))
                mat.append(vec)
            kops.append(mat)
        print kops
    return latt,kops

def cmp_mat(A,B, Tol=0.001):
    #returns True if two 3x3 matrices agree within Tol
    for i in range(3):
        for j in range(3):
            d=A[i][j]-B[i][j]
            if (sqrt(d*d)>Tol):
                return False
    return True

def cmp_matrix_set(A,B, no_multiple_matches=False):
    nfound=0
    if (A==None) or (B==None):
        print 'cmp_matrix_set(): empty set, nothing done'
        return nfound
    B_range=range(len(B))
    for i in range(len(A)):
        found=False
        for j in B_range:
            #print B_range
            if cmp_mat(A[i],B[j]):
                print 'A_',i,' = B_',j
                if no_multiple_matches:
                    B_range.remove(j)
                #print B_range
                found=True
                break
        if found==False:
            print 'A_',i,' not found in set B'
        else:
            nfound=nfound+1
    print 'Found ',nfound,'of ',len(A),len(B)
    return nfound


def write_intrans_boltztrap(boltzdir='./',filename='hte.intrans', boltz_generic=True, E_Fermi=0.0, n_electrons=1.0):
    parentdir=os.getcwd()
    if not os.path.isdir(boltzdir):
        os.mkdir(boltzdir)
    os.chdir(boltzdir)
    fout=open(filename,'w')
    if boltz_generic:
        fout.write("GENE          # use generic interface\n")
    else:
        fout.write("WIEN          # use wien interface\n")
    fout.write("0 0 0 0.0         # iskip (not presently used) idebug setgap shiftgap \n")
    fout.write("%7.5f 0.0005 0.4 %6.1f     # Fermilevel (Ry), energygrid, energy span around Fermilevel, number of electrons\n"%(E_Fermi,n_electrons))
    fout.write("CALC                    # CALC (calculate expansion coeff), NOCALC read from file\n")
    fout.write("5                         # lpfac, number of latt-points per k-point\n")
    fout.write("BOLTZ                     # run mode (only BOLTZ is supported)\n")
    fout.write(".15                       # (efcut) energy range of chemical potential\n")
    fout.write("800. 50.                  # Tmax, temperature grid\n")
    fout.write("-1.                       # energyrange of bands given individual DOS output sig_xxx and dos_xxx (xxx is band number)\n")
    fout.write("HISTO\n")
    fout.close()
    os.chdir(parentdir)
    return True


def analyze_bandgap(bs, set_E_Fermi='mid_gap', debug=False):
    # analysis of bandstructure: gap, valence band maximum and conduction band minimum
    # bs: dict with banstructure information
    # set_E_Fermi: 'mid_gap', 'vb_max', 'cb_min'
    # check if all necessary information is contained in bs:
    if (bs!=None) and ('n_electrons' in bs) and (bs['n_electrons']!=None) and ('nspin' in bs) and (bs['nspin']==1):
        soc=2.0
        if ('fullrel' in bs) and (bs['fullrel']==True):
            soc=1.0
        nel=int(round(bs['n_electrons']/soc))
        if abs(nel*soc-bs['n_electrons'])>0.0001:
            bs['gap']=0.0
            if debug==True:
                print 'analyze_bandgap(): electron count -> metallic band structure'
            return True
        # determine valence band maximum and conduction band minimum
        ikp_vbmax=0
        ikp_cbmin=0
        eb=bs['e_n_k']
        for ikp in range(len(eb)):
            enk=eb[ikp]['energies_up']
            if (eb[ikp]['energies_up'][nel-1]>eb[ikp_vbmax]['energies_up'][nel-1]):
                ikp_vbmax=ikp
            if (eb[ikp]['energies_up'][nel]<eb[ikp_cbmin]['energies_up'][nel]):
                ikp_cbmin=ikp
        if debug==True:
            print 'vbmax=%6.2f %s at k=(%5.2f,%5.2f,%5.2f)' %(eb[ikp_vbmax]['energies_up'][nel-1],bs['unit_E'],eb[ikp_vbmax]['kpoint'][0],eb[ikp_vbmax]['kpoint'][1],eb[ikp_vbmax]['kpoint'][2])
            print 'cbmin=%6.2f %s at k=(%5.2f,%5.2f,%5.2f)' %(eb[ikp_cbmin]['energies_up'][nel],bs['unit_E'],eb[ikp_cbmin]['kpoint'][0],eb[ikp_cbmin]['kpoint'][1],eb[ikp_cbmin]['kpoint'][2])
        if (eb[ikp_vbmax]['energies_up'][nel-1]<eb[ikp_cbmin]['energies_up'][nel]):
            gap=eb[ikp_cbmin]['energies_up'][nel]-eb[ikp_vbmax]['energies_up'][nel-1]
            new_E_Fermi=bs['E_Fermi']
            if set_E_Fermi=='mid_gap':
                new_E_Fermi=eb[ikp_vbmax]['energies_up'][nel-1]+0.5*gap
            elif set_E_Fermi=='vb_max':
                new_E_Fermi=eb[ikp_vbmax]['energies_up'][nel-1]
            elif set_E_Fermi=='cb_min':
                new_E_Fermi=eb[ikp_cbmin]['energies_up'][nel]
            if debug==True:
                print 'New Fermi energy (%s): E_F=%6.3f %s (old value %6.3f)' %(set_E_Fermi,new_E_Fermi,bs['unit_E'],bs['E_Fermi'])
                print 'band gap: %6.3f %s (%6.3f eV)'  %(gap,bs['unit_E'],gap*get_conversion(bs['unit_E'], 'eV'))
            bs['E_Fermi']=new_E_Fermi
            bs['gap']=gap
        else:
            bs['gap']=0.0
            if debug==True:
                print 'analyze_bandgap(): metallic bandstructure'
        return True
    else:
        print 'analyze_bandgap(): insufficient information, nothing done!'
        return False
                                           

def get_kspace_operations(ao, methods=['atoms_info','spglib'], symprec_spglib=1e-5):
    # returns k-space operations for the atoms object ao
    kops=None
    for method in methods:
        if method=='atoms_info':
            # get operations from ao.info
            if ('spacegroup' in ao.info) and (ao.info['spacegroup']!=None):
                if ('unit_cell' in ao.info):
                    if (ao.info['unit_cell']=='conventional'):
                        primitive_cell=False
                    else:
                        primitive_cell=True
                else:
                    print 'get_kspace_operations(): Warning, assuming primitive cell'
                rot=ao.info['spacegroup'].rotations
                p2c_dir=ao.info['spacegroup'].scaled_primitive_cell
                p2c_rec=ao.info['spacegroup'].reciprocal_cell 
                kops=[]
                for iop in range(len(rot)):
                    if primitive_cell:
                        mat=dot(p2c_rec.T, dot(rot[iop],p2c_dir))
                    else:
                        mat=rot[iop]
                    newop=True
                    for i in range(len(kops)):
                        if cmp_mat(kops[i],mat):
                            newop=False
                            break
                    if newop:
                        kops.append(mat)
                return kops
            else:
                print 'get_kspace_operations(): atoms object has no space group information'
        if method=='spglib':
            if has_spglib:
                rot=spglib.get_symmetry_dataset(ao, symprec=symprec_spglib)['rotations']
                kops=[]
                for iop in range(len(rot)):
                    mat=rot[iop]
                    newop=True
                    for i in range(len(kops)):
                        if cmp_mat(kops[i],mat):
                            newop=False
                            break
                    if newop:
                        kops.append(mat)
                return kops
            else:
                print 'get_kspace_operations(): spglib not found'
        
