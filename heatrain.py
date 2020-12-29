

import torch
import torch.nn as nn
import torch.nn.functional as F



class TransferNet(nn.Module):

    def __init__(self):
        super(TransferNet, self).__init__()
        self.conv1 = nn.Conv2d(nChannels, 60, 3)
        self.conv2 = nn.Conv2d(60, 90, 3)
        self.conv3 = nn.ConvTranspose2d(90, 60, 3)
        self.conv4 = nn.ConvTranspose2d(60, nChannels, 3)
        self.thr = nn.ReLU()
        
        
    def forward(self,x):
        # x0=x
        #x = (self.bn1(x))
        x = (self.thr(self.conv1(x)))
        x = (self.thr(self.conv2(x)))
        x = (self.thr(self.conv3(x)))
        x = self.conv4(x) 
        return x  




import hdf5storage
mat = hdf5storage.loadmat('Heatdata.mat') 
import torch
import h5py
import numpy as np
import torch.nn as nn
import torch.nn.functional as F
import matplotlib.pylab as plt
import pickle
import scipy.io as sio
dtype = torch.double
device = torch.device("cpu")
# torch.set_printoptions(precision=2)
# device = torch.device("cuda:0") # Uncomment this to run on GPU


Z_in = mat['meas']
Z_tgt= mat['source']

#Convert into torch arrays
Z_in=torch.from_numpy(Z_in)
Z_tgt=torch.from_numpy(Z_tgt)


Z_tgt = Z_tgt.reshape([Z_tgt.shape[0],Z_tgt.shape[1],Z_tgt.shape[2]]) 
Z_in = Z_in.reshape([Z_in.shape[0],Z_in.shape[1],Z_in.shape[2]]) 

nValidation = 1000
nData = Z_in.shape[2]-nValidation


print(nData)

# initial scrambling 
perm = torch.randperm(nData+nValidation)
Z_in = Z_in[:,:,perm]
Z_tgt = Z_tgt[:,:,perm]

# validation dataset
val_in = Z_in[:,:,-nValidation:]
val_tgt = Z_tgt[:,:,-nValidation:]



# N is batch size
N = 50
nChannels = 1
nEpochs =100
beta=3

tl=torch.zeros(nEpochs,1)
vl=torch.zeros(nEpochs,1)


net = TransferNet()
net = net.float()

# Loss and optimizer
learning_rate = 5e-1
criterion = nn.MSELoss()  
optimizer = torch.optim.Adagrad(net.parameters(), lr=learning_rate)

for t in range(nEpochs):
    tloss=0
    val_loss = 0
    #perm = torch.randperm(nData)
    for b_ix in np.arange(0,nData,N):
        
        
        ##Pre-processing sparse codes for input 
        x1 =  Z_in[:,:,perm[b_ix:b_ix+N]].view(N,nChannels,Z_in.shape[0],Z_in.shape[1])
        y1 =  Z_tgt[:,:,perm[b_ix:b_ix+N]].view(N,nChannels,Z_tgt.shape[0],Z_tgt.shape[1])
        

        ##zero parameter gradients
        optimizer.zero_grad()


        ##Passing the sparse code 
        
        x_out1 = net(x1.float())
        
        
        ##calculating the loss for passed input and target
        
        loss1 = criterion(x_out1,y1.float()).float() #+ beta*torch.sum(torch.abs(x_out1))
        
        
        # Backward and optimize
        
        loss1.backward()
        optimizer.step()

        with torch.no_grad():
            tloss+=loss1.item()
            if b_ix % nData == 0:
                val_pred = net(val_in.reshape(nValidation,1,val_in.shape[0],val_in.shape[1]).float())
                val_loss = criterion(val_pred, val_tgt.reshape(nValidation,1,val_in.shape[0],val_in.shape[1])).float() 
                print('Epoch:{:d}| T Loss:{:.6f}| V Loss:{:.6f}'.format(t+1, tloss/nData, val_loss.item()/nValidation))
                tl[t]=tloss/nData
                vl[t]=val_loss.item()/nValidation
                tloss = 0
                val_loss = 0




##Testing 



Z_in = mat['meas']
Z_tgt= mat['source']

#Convert into torch arrays
Z_in=torch.from_numpy(Z_in)
Z_tgt=torch.from_numpy(Z_tgt)


tn=15

zin=Z_in[:,:,0:tn]
ztgt=Z_tgt[:,:,0:tn]

x=zin.reshape([tn,nChannels,zin.shape[0],zin.shape[1]])
y=ztgt.reshape([tn,nChannels,ztgt.shape[0],ztgt.shape[1]])

xout= net(x.float())

zrec=xout.reshape([xout.shape[2],xout.shape[3],tn])


zin=zin.detach().numpy()
ztgt=ztgt.detach().numpy()
zrec=zrec.detach().numpy()

sio.savemat('source.mat',{'zin': zin})  
sio.savemat('meas.mat',{'ztgt': ztgt}) 
sio.savemat('rec.mat',{'zrec': zrec}) 
sio.savemat('tloss.mat',{'tl': tl})
sio.savemat('vloss.mat',{'vl': vl})
