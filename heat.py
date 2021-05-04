import torch
import torch.nn as nn
import torch.nn.functional as F
import matplotlib.pylab as plt
import pickle
import scipy.io as sio
import numpy as np 
from numpy.linalg import matrix_power
import cvxpy as cp
import random 

def l1solv(Fw1,meas,norm):
# Create variable.
   x_l1 = cp.Variable(shape=(meas.shape[0],1))
   constraints = [x_l1>=0]
# Form objective.
   obj = cp.Minimize(cp.norm(x_l1, norm)+cp.norm(Fw1.numpy()@x_l1-meas, 2)**2*1e2)
# Form and solve problem.
   prob = cp.Problem(obj, constraints)
   prob.solve()
   x2 = torch.from_numpy(x_l1.value[:,0]).cpu().type(dtype=torch.float32)
   return torch.norm(x2,0.5)
   #return net1(x2)

def l1rec(Fw1,meas,norm):
# Create variable.
   x_l1 = cp.Variable(shape=(meas.shape[0],1))
   constraints = [x_l1>=0]
# Form objective.
   obj = cp.Minimize(cp.norm(x_l1, norm)+cp.norm(Fw1.numpy()@x_l1-meas, 2)**2*1e2)
# Form and solve problem.
   prob = cp.Problem(obj, constraints)
   prob.solve()
   x2 = torch.from_numpy(x_l1.value[:,0]).cpu().type(dtype=torch.float32)
   return x2


def slope(N,A,meas):
    F_Nr=l1solv(torch.matrix_power(A, int(N+20)),meas,2)
    F_Nl=l1solv(torch.matrix_power(A, int(N-20)),meas,2)
    return (F_Nr-F_Nl)/40        


def search2(lf=100,rb=2500,walk=100,iter=100):
    meas, measn, u, A,Fw = gensingledata(n=500,alpha=20,m=5,Ncount=1000,sigma=1)
    for j in range(iter):
        if j==10:
            walk = 50
        slf =slope(lf,A,measn)
        srb = slope(rb,A,measn)
        if slf<0 and srb>0:
            lf=lf+walk
            rb=rb-walk
        elif slf<0 and srb<0:
            rb=rb+walk
        elif slf>0 and srb>0:
            lf=lf-walk  
        elif slf>0 and srb<0:
            lf=lf-walk
            rb=rb+walk    
        print(lf,rb,slf,srb)
            
        
        


def genData(n=500,alpha=20,k=5,N=1000,Ncount=1000):
    source = torch.zeros(n,N)
    wsource = torch.zeros(n,N)
    A=makeA1d(n,alpha)
    Fw = torch.matrix_power(A, Ncount)
    for i in range(0,N):
        Unit=torch.zeros(n,1)
        Unit[torch.randint(0,n-2,(k,))+1,0]=1000*torch.rand(k)
        source[:,i]=Unit.T
        meas = Fw @ Unit
        offset = torch.randint(-900,900,(1,))
        wsource[:,i] = l1rec(torch.matrix_power(A, Ncount+int(offset)),meas,2)
    return source, wsource    

def makeA1d(n=500,alpha=20):
    A=torch.zeros((n,n))
    s=1/(2*alpha)
    for i in range(0,n-1):
        A[i,i]=1-2*s
        A[i,i-1]=s
        A[i,i+1]=s
    A[:,0]=0
    A[0,:]=0
    A[:,n-1]=0
    A[n-1,:]=0
    return A
    

def makeA2d(n=10,alpha=20):
    s=torch.zeros((n,n))
    s[0:n-1,0:n-1]=1/(2*alpha)
    sx=torch.reshape(s,(n**2,1))
    A=torch.zeros((n**2,n**2))
    sy=sx
    ##interior points
    for i in range(1,(n**2)-2):
        A[i,i]=(1-2*sx[i]-2*sy[i])
        A[i,i-1]=sx[i]
        A[i,i+1]=sx[i]
        A[i,((i+n)%(n**2))]=sy[i]
        A[i,((i-n)%(n**2))]=sy[i]
    
    ##x-axis boundaries
    A[0:n-1,:]=0
    A[:,0:n-1]=0
    A[((n**2)-n):((n**2)-1),:]=0
    A[:,((n**2)-n):((n**2)-1)]=0
        
    ##y-axis boundaries
    for i in range(1,n):
        A[(i*n)-1,:]=0
        A[(i*n)-2,:]=0
    return A


            
def gensingledata(n=500,alpha=20,m=5,Ncount=1000,sigma=0.1):
    torch.random.seed()
    A=makeA1d(n,alpha)
    u=torch.zeros((n,1))
    np.random.seed(0)
    k = np.random.randint(0,n-50,m)
    torch.manual_seed(0)
    u[k+25] = 1000*torch.rand(m,1)
    Fw=torch.matrix_power(A, Ncount)
    meas=Fw@u
    measn= meas+sigma*torch.rand(meas.shape)
    return meas, measn,u,A,Fw



def gensingledata2d(n=50,alpha=50,m=4,Ncount=1000,sigma=0.1):
    A2d=makeA2d(n,alpha)
    u=torch.zeros(n,n)
    np.random.seed(0)
    l = [(random.randrange(10,n-10), random.randrange(10,n-10)) for i in range(m)]
    for j in l:
        u[j]=1000*torch.rand(1)
    ures=torch.reshape(u,(n**2,1))
    Fw=torch.matrix_power(A2d, Ncount)
    mres=Fw@ures
    meas=torch.reshape(mres,(n,n))
    measn=meas + sigma*torch.rand(meas.shape)
    #x=l1rec(Fw,mres,1)
    #rec = torch.reshape(x,(n,n))
    return u,ures,meas,mres,A2d
        
    
    
def searchalgo(lf=400,rb=3000,walk=300,iter=100):
    mid = (lf+rb)/2
    meas, measn, u, A,Fw = gensingledata(n=500,alpha=20,m=5,Ncount=1000,sigma=0.04)
    Fl = l1solv(torch.matrix_power(A, int(lf)),meas,2)
    Fr = l1solv(torch.matrix_power(A, int(rb)),meas,2)
    
    #thr=1e-5*np.linalg.norm(meas)
    thr=20e4
    for j in range(iter):
        if j==40:
            walk = 300
        if (Fl>thr) & (Fr>thr):
          mid = lf + np.random.randint(-walk,walk)
          Fm = l1solv(torch.matrix_power(A, int(mid)),meas,2)
          if Fm>thr:
              if mid > lf:
                 lf = mid
          else:
            if mid <rb:
              rb = mid

        elif Fl<thr:
            lf = lf/2
        else:
            rb = rb+1000
        est=max(lf,rb)
        Fl= l1solv(torch.matrix_power(A, int(lf)),meas,2)
        Fr= l1solv(torch.matrix_power(A, int(rb)),meas,2)
        #sl=slope(est)
        #print(sl)
        #print(Fm)
        print('left:{:.1f},right:{:.1f}, estimation:{:.1f}'.format(lf,rb,est))
    return est        
            
def plotnormvsN(lw=1000,rb=9000,points=200):
    N1 = torch.linspace(lw,rb,points)
    F=torch.zeros_like(N1)
    meas, measn, u, A,Fw = gensingledata(n=500,alpha=20,m=15,Ncount=5000,sigma=0.5)
    #u,ures,meas,mres,A2d= gensingledata2d(n=30,alpha=20,m=4,Ncount=500,sigma=0.1)
    for i in range(0,points):
        F[i]= l1solv(torch.matrix_power(A, int(N1[i])),measn,2)
    plt.plot(N1.detach().numpy(),F.detach().numpy())
    plt.xlabel('Ncount')
    plt.ylabel('Defined norm')
    plt.title('Norm vs ncount for 1D with Ncount=5000 but with 5 spikes with noise amplitude 1.')
    plt.savefig('figures/Ncount5000.svg', format='svg', dpi=1200)
    return F,N1

def perreco(noise=False,sigma=0.04):
   meas, measn, u, A,Fw = gensingledata(n=500,alpha=20,m=5,Ncount=2000,sigma=0.04)
   x_l1 = cp.Variable(shape=(meas.shape[0],1))
   if (noise==False): 
       constraints = [Fw.numpy()@x_l1==meas]
   else:
       eps=5     
       constraints = [cp.norm(Fw.numpy()@x_l1-measn, 1) <=eps ]     
   obj = cp.Minimize(cp.norm(x_l1, 1))
# Form and solve problem.
   prob = cp.Problem(obj, constraints)
   prob.solve()
   x2 = torch.from_numpy(x_l1.value).cpu().type(dtype=torch.float32)
    
   plt.plot(x2,label='Reconl1')
   print(torch.norm(x2,1))
   plt.plot(u,label='Original')
   print(torch.norm(u,1))
   plt.legend() 
   return u,x2 