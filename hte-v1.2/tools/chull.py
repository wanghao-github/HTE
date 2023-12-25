from math import *
import numpy as np

class CHULL(object):
    """Calculation of convex hulls
    STATUS: untested, use at own risk!
    author: IO
    """
    def __init__(self, dim=2):
        self.dimension=dim
        # list with datapints on convex hull
        self.chull_points=[]
        for i in range(dim):
            p=[]
            for j in range(dim):
                if (i==j):
                    p.append(1.0)
                else:
                    p.append(0.0)
            px=tuple(p),0.0
            self.chull_points.append(px)
        # regions with linear dependence, border defined by indices of chull_points
        self.linsupport=[]
        p=range(dim)
        self.linsupport.append(p)
        
    def print_points(self):
        print "data points on convex hull:"
        for points in self.chull_points:
            print points
            
    def get_value(self,x,Debug=True):
        """get_value(self,x,Debug=True)
        returns the value of the convex hull for x
        input: x tuple with mole fractions
        """
        dim=self.dimension
        ymin=None
        supp_min=None
        for supp in self.linsupport:
            #print supp
            A=np.zeros([dim,dim])
            for i in range(dim):
                dp=self.chull_points[supp[i]]
                p=dp[0]
                for j in range(dim):
                    #print p[j]
                    A[j,i]=p[j]
            #print 'A=',A
            B=np.matrix(x).T
            #print 'B=',B
            try:
                alph=np.linalg.solve(A,B)
                y_x=0.0
                for i in range(dim):
                    dp=self.chull_points[supp[i]]
                    y=dp[1]
                    y_x=y_x+float(alph[i])*y
                    if (alph[i]<0.0) or (alph[i]>1.0):
                        supp=None
                        break
                #print 'XXX',x,alph,supp,y_x
                if (supp!=None) and ((ymin==None) or (y_x<ymin)):
                    ymin=y_x
                    supp_min=supp
                    #print 'min:',supp_min,ymin
            except np.linalg.linalg.LinAlgError as err:
                if 'Singular matrix' in err.message:
                    z=0
                    #print "Singular matrix for ",supp
                else:
                    raise
        if supp_min==None:
            print "Warning: Could not find support for ",x
            return None,None
        return ymin,supp_min
    
    def get_convex_hull(self,dataset,TOL_E=0.0000001):
        pmin=None
        ymin=None
        dataset_new=[]
        for p in dataset:
            x=p[0]
            y=p[1]
            yhull,supp=self.get_value(x)
            if (y-yhull)<=TOL_E:
                if p not in dataset_new:
                    dataset_new.append(p)
                if (pmin==None) or (y-yhull<=ymin):
                    pmin=p
                    ymin=y-yhull
        #for the moment: connect all points with triangle
        if (pmin!=None) and (ymin<=TOL_E) and (pmin not in self.chull_points):
            np=len(self.chull_points)
            for i in range(np):
                if self.dimension==2:
                    self.linsupport.append([i,np])
                elif self.dimension==3:
                    for j in range(np):
                        if i<j:
                            self.linsupport.append([i,j,np])
            self.chull_points.append(pmin)
            dataset_new.remove(pmin)
        if len(dataset_new)!=0:
            #print 'New candlist:',dataset_new
            #print (pmin!=None), (ymin<=TOL_E), (pmin not in self.chull_points)
            #for p in self.chull_points:
            #    print  p==pmin,p,' xxx ',pmin
            self.get_convex_hull(dataset_new)
        return self
    
def get_convex_hull(dataset,dim):
    chull=CHULL(dim=dim)
    pmin=None
    ymin=None
    candlist={}
    for supp in chull.linsupport:
        candlist[supp]=[]
    for p in dataset:
        x=p[0]
        y=p[1]
        yhull,supp=chull.get_value(x)
        if (y<yhull) and ((pmin==None) or (y-yhull<ymin)):
            pmin=p
            ymin=y-yhull
                
    #print pmin,ymin
    return None
