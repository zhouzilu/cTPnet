# source /opt/rh/rh-python36/enable
# source /home/stat/zhouzilu/project/zhouzilu/expressionGAN/.env/bin/activate
import numpy as np
import torch.nn as nn
import pandas as pd
import matplotlib.pyplot as plt
import torch.nn.functional as F
import torch
import torch.optim as optim
import sys
from scipy import stats
from matplotlib.backends.backend_pdf import PdfPages
import re
from sklearn.metrics import pairwise_distances
from matplotlib.ticker import FormatStrFormatter

#################
### Figure 1b ###
#################

loc='../cTP_net_figures_MB/'
# Take a look at the PBMC data itself
X = pd.read_csv('../../data/GSE100866_PBMC_vs_flow_10X-RNA_umi.csv') # Raw RNA data
y = pd.read_csv('../../data/GSE100866_PBMC_vs_flow_10X-ADT_umi.csv') # Raw cite-seq protein data
gene = X.iloc[:,0]
X=X.drop(columns=X.columns[0])
protein = y.iloc[:,0]
y=y.drop(columns=y.columns[0])
y.index=protein
# sanity check
(X.columns==y.columns).all()
# load celltype classification
celltype=(pd.read_csv('/home/stat/zhouzilu/project/Zilu_Chengzhong/CITEseq/CellType/cite_pbmc.cluster.csv',index_col='cell')).transpose()
# sanity check
(X.columns==celltype.columns).all()
# filtering
human_perc=X.loc[gene.str.contains('HUMAN').tolist(),:].sum(axis=0)/X.sum(axis=0)
X=X.loc[gene.str.contains('HUMAN').tolist(),human_perc>0.9]
y=y.loc[:,human_perc>0.9]
gene=gene[gene.str.contains('HUMAN').tolist()]
celltype=celltype.loc[:, human_perc>0.9]
# sanity check
(X.columns==celltype.columns).all()
protein0=['CD3', 'CD4', 'CD8', 'CD2', 'CD45RA', 'CD57', 'CD16', 'CD14','CD11c', 'CD19']
y=y.loc[protein0]
y=y.apply(lambda x: np.log((x+1.0)/stats.gmean(x+1.0)), axis=0)
# Remove low expressed genes
gene = gene[(X.sum(axis=1)>10).tolist()]
X=X.loc[X.sum(axis=1)>10,:]
# Use direct UMI count is hard to transfer across experiment
# Let's try normalize the data with seurat like method
X=X.apply(lambda x: np.log((x*10000.0/sum(x))+1.0), axis=0)
# random cell order
X=X.T
X=X.sample(frac=1,random_state=4905)
celltype_tab=['CD4 T','CD14+CD16+ Mono','NK','CD8 T 1','B','CD8 T 2','CD14-CD16+ Mono','DC']
# separate test data from train data
X_test=X.sample(n=667,random_state=4905)
y_test=y[X_test.index]
celltype_test=celltype[X_test.index]
raw_expression_testindex = X_test.index
# Plot correlation between cell expression and protein levels in the test set
gene0_expand=['HUMAN_CD3','HUMAN_CD4','HUMAN_CD8']
f,axs=plt.subplots(3,3,figsize=(7.5,7.5))
axs = axs.ravel(order='F')
protein0_expland=['CD3', 'CD4', 'CD8']
#gene[gene == 'HUMAN_ITGAX']]
for i,p in enumerate(protein0_expland):
	if p=='CD3':
		axs[i].scatter(X_test.loc[:,[x in ['HUMAN_CD3D','HUMAN_CD3E','HUMAN_CD3G','HUMAN_CD247'] for x in gene]].sum(1),y_test.loc[p,:],alpha=0.5,marker ='.')
		axs[i].set_xlabel('RNA Count\nCD3s')
		axs[i].set_ylabel('Protein\n'+p)
		#axs[i].plot([-2,5],[-2,5],ls='--',c='.3')
		axs[i].set_title('cor: '+str(np.round_(np.corrcoef(X_test.loc[:,[x in ['HUMAN_CD3D','HUMAN_CD3E','HUMAN_CD3G','HUMAN_CD247'] for x in gene]].sum(1),y_test.loc[p,:])[1,0],2)))
		corDF.loc['rawRNA',p]=np.corrcoef(X_test.loc[:,[x in ['HUMAN_CD3D','HUMAN_CD3E','HUMAN_CD3G','HUMAN_CD247'] for x in gene]].sum(1),y_test.loc[p,:])[1,0]
	elif p=='CD8':
		axs[i].scatter(X_test.loc[:,[x in ['HUMAN_CD8A','HUMAN_CD8B'] for x in gene]].sum(1),y_test.loc[p,:],alpha=0.5,marker ='.')
		axs[i].set_xlabel('RNA Count\nCD8s')
		axs[i].set_ylabel('Protein\n'+p)
		#axs[i].plot([-2,5],[-2,5],ls='--',c='.3')
		axs[i].set_title('cor: '+str(np.round_(np.corrcoef(X_test.loc[:,[x in ['HUMAN_CD8A','HUMAN_CD8B'] for x in gene]].sum(1),y_test.loc[p,:])[1,0],2)))
		corDF.loc['rawRNA',p]=np.corrcoef(X_test.loc[:,[x in ['HUMAN_CD8A','HUMAN_CD8B'] for x in gene]].sum(1),y_test.loc[p,:])[1,0]
	else:
		axs[i].scatter(X_test.loc[:,gene==gene0_expand[i]],y_test.loc[p,:],alpha=0.5,marker ='.')
		axs[i].set_xlabel('RNA Count\n'+gene0_expand[i][6:])
		axs[i].set_ylabel('Protein\n'+p)
		#axs[i].plot([-2,5],[-2,5],ls='--',c='.3')
		axs[i].set_title('cor: '+str(np.round_(np.corrcoef(X_test.loc[:,gene==gene0_expand[i]].values.flatten(),y_test.loc[p,:])[1,0],2)))
		corDF.loc['rawRNA',p]=np.corrcoef(X_test.loc[:,gene==gene0_expand[i]].values.flatten(),y_test.loc[p,:])[1,0]

f.subplots_adjust(left=0.12, bottom=0.11, right=0.9, top=0.88, wspace=0.8, hspace=0.8)

#################################
# With denoised expression data #
#################################
X = pd.read_csv('/home/stat/zhouzilu/project/Zilu_Chengzhong/CITEseq/CBMC_PBMC_svx_denoised/PBMC_tl.csv') # Denoised RNA data
y = pd.read_csv('../../data/GSE100866_PBMC_vs_flow_10X-ADT_umi.csv') # Raw cite-seq protein data
gene = X.index
protein = y.iloc[:,0].values
y=y.drop(columns=y.columns[0])
y.index=protein
y=y[X.columns]
# sanity check
(X.columns==y.columns).all()
# filtering
protein0=['CD3', 'CD4', 'CD8', 'CD2', 'CD45RA', 'CD57', 'CD16', 'CD14','CD11c', 'CD19']
y=y.loc[protein0]
y=y.apply(lambda x: np.log((x+1.0)/stats.gmean(x+1.0)), axis=0)
# Remove low expressed genes
X=X.loc[X.sum(axis=1)>10,:]
# Use direct UMI count is hard to transfer across experiment
# Let's try normalize the data with seurat like method
X=X.apply(lambda x: np.log((x*10000.0/sum(x))+1.0), axis=0)
# random cell order
X=X.T
X=X.sample(frac=1,random_state=4905)
# separate test data from train data
X_test=X.sample(n=667,random_state=4905)
y_test=y[X_test.index]
X=X.drop(X_test.index)
y=y.drop(columns=y_test.columns)
y=y[X.index]
# sanity check
(X_test.index==raw_expression_testindex).all()
# Sanity check have passed. They are the same set of cells
# Plot correlation between cell expression and protein levels in the test set
gene0_expand=['CD3','CD4','CD8']
# gene for CD3: CD3D, CD3E, CD3G, CD247
# gene for CD8: CD8A, CD8B
# By protein-type diagnosis
#with PdfPages(loc+'figure_rawRNA.pdf') as pdf:
protein0_expland=['CD3', 'CD4', 'CD8']
for i,p in enumerate(protein0_expland):
	if p=='CD3':
		axs[i+3].scatter(X_test.loc[:,['CD3D','CD3E','CD3G','CD247']].sum(1),y_test.loc[p,:],alpha=0.5,marker ='.')
		axs[i+3].set_xlabel('Denoised RNA\nCD3s')
		axs[i+3].set_title('cor: '+str(np.round_(np.corrcoef(X_test.loc[:,['CD3D','CD3E','CD3G','CD247']].sum(1),y_test.loc[p,:])[1,0],2)))
		corDF.loc['dnRNA',p]=np.corrcoef(X_test.loc[:,['CD3D','CD3E','CD3G','CD247']].sum(1),y_test.loc[p,:])[1,0]
	elif p=='CD8':
		print(1)
		axs[i+3].scatter(X_test.loc[:,['CD8A','CD8B']].sum(1),y_test.loc[p,:],alpha=0.5,marker ='.')
		axs[i+3].set_xlabel('Denoised RNA\nCD8s')
		axs[i+3].set_title('cor: '+str(np.round_(np.corrcoef(X_test.loc[:,['CD8A','CD8B']].sum(1),y_test.loc[p,:])[1,0],2)))
		corDF.loc['dnRNA',p]=np.corrcoef(X_test.loc[:,['CD8A','CD8B']].sum(1),y_test.loc[p,:])[1,0]
	else:
		axs[i+3].scatter(X_test.loc[:,gene0_expand[i]],y_test.loc[p,:],alpha=0.5,marker ='.')
		axs[i+3].set_xlabel('Denoised RNA\n'+gene0_expand[i])
		axs[i+3].set_title('cor: '+str(np.round_(np.corrcoef(X_test.loc[:,gene0_expand[i]].values.flatten(),y_test.loc[p,:])[1,0],2)))
		corDF.loc['dnRNA',p]=np.corrcoef(X_test.loc[:,gene0_expand[i]].values.flatten(),y_test.loc[p,:])[1,0]

#################################
# With predicted protein data #
#################################

# covert to tenor
X=torch.tensor(X.values)
X=X.type(torch.FloatTensor)
y=torch.tensor(y.values)
y=y.type(torch.FloatTensor)
y=torch.t(y)
X_test=torch.tensor(X_test.values)
X_test=X_test.type(torch.FloatTensor)
y_test=torch.tensor(y_test.values)
y_test=y_test.type(torch.FloatTensor)
y_test=torch.t(y_test)
n_batches=10


class Net(nn.Module):
    def __init__(self):
        super(Net, self).__init__()
        #self.fc1 = nn.Linear(14505, 1000)
        self.fc1 = nn.Linear(X.shape[1], 1000)
        self.fc3 = nn.Linear(1000, 128)
        self.fc40 = nn.Linear(128, 64)
        self.fc50 = nn.Linear(64, 1)
        self.fc41 = nn.Linear(128, 64)
        self.fc51 = nn.Linear(64, 1)
        self.fc42 = nn.Linear(128, 64)
        self.fc52 = nn.Linear(64, 1)
        self.fc43 = nn.Linear(128, 64)
        self.fc53 = nn.Linear(64, 1)
        self.fc44 = nn.Linear(128, 64)
        self.fc54 = nn.Linear(64, 1)
        self.fc45 = nn.Linear(128, 64)
        self.fc55 = nn.Linear(64, 1)
        self.fc46 = nn.Linear(128, 64)
        self.fc56 = nn.Linear(64, 1)
        self.fc47 = nn.Linear(128, 64)
        self.fc57 = nn.Linear(64, 1)
        self.fc48 = nn.Linear(128, 64)
        self.fc58 = nn.Linear(64, 1)
        self.fc49 = nn.Linear(128, 64)
        self.fc59 = nn.Linear(64, 1)
    def forward(self, x):
        x = F.relu(self.fc1(x))
        x = F.relu(self.fc3(x))
        o0 = self.fc50(F.relu(self.fc40(x)))
        o1 = self.fc51(F.relu(self.fc41(x)))
        o2 = self.fc52(F.relu(self.fc42(x)))
        o3 = self.fc53(F.relu(self.fc43(x)))
        o4 = self.fc54(F.relu(self.fc44(x)))
        o5 = self.fc55(F.relu(self.fc45(x)))
        o6 = self.fc56(F.relu(self.fc46(x)))
        o7 = self.fc57(F.relu(self.fc47(x)))
        o8 = self.fc58(F.relu(self.fc48(x)))
        o9 = self.fc59(F.relu(self.fc49(x)))
        return o0,o1,o2,o3,o4,o5,o6,o7,o8,o9

net = Net()
net.load_state_dict(torch.load('../../PbyGE/model_CITE_PBMC_norm_protein_denoised_allprotein_MB_L1/model_rep1_ep119')) # Best optimal NN model, available in Github

# optimizer = optim.Adagrad(net.parameters(), lr=lr_vec[lri], lr_decay=0.001)
test_outputs = list(net(X_test))
test_outputs=torch.transpose(torch.stack(test_outputs),0,1).view(X_test.shape[0],-1) # This is the protein prediction data
# By protein-type diagnosis
protein_tmp=['CD3', 'CD4', 'CD8']
for i,p in enumerate(protein_tmp):
	# ax.scatter(y_truth.drop(['CD4','CD8']),y_pred.drop(['CD4','CD8']))
	axs[i+6].scatter(test_outputs.detach().numpy()[:,i],y_test.detach().numpy()[:,i],alpha=0.5,marker ='.')
	axs[i+6].plot([-2,5],[-2,5],ls='--',c='.3')
	axs[i+6].set_title('cor: '+str(np.round_(np.corrcoef(y_test.detach().numpy()[:,i],test_outputs.detach().numpy()[:,i])[1,0],2)))
	axs[i+6].set_xlabel('Predicted Protein\n'+p)
	corDF.loc['noholdoutdn',p]=np.corrcoef(y_test.detach().numpy()[:,i],test_outputs.detach().numpy()[:,i])[1,0]

f.subplots_adjust(left=0.12, bottom=0.11, right=0.9, top=0.88, wspace=0.8, hspace=0.8)
f.savefig(loc+'figure1b.pdf')
plt.close()

