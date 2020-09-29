# Last update time 05/15/2020
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

interactive=True

protein_list=pd.read_csv('../data/protein_list/protein_list_01072020',header=None)[0].values.tolist()
X_list=['../data/PBMC_Stoeckius/PBMC_tl.csv','../data/CBMC_Stoeckius/CBMC_tl.csv','../data/BMMC_Stuart/GSE128639_RNA_counts_denoised_filt.csv','../data/PBMC_Ye_Pre/p1_illumina_raw_filtered_RNA_denoised.csv','../data/PBMC_Ye_Pre/p2_illumina_raw_filtered_RNA_denoised.csv']
y_list=['../data/PBMC_Stoeckius/GSE100866_PBMC_vs_flow_10X-ADT_umi.csv','../data/CBMC_Stoeckius/GSE100866_CBMC_8K_13AB_10X-ADT_umi.csv','../data/BMMC_Stuart/GSE128639_citeseq_adt_counts.tsv','../data/PBMC_Ye_Pre/p1_illumina_raw_filtered_Antibody.csv','../data/PBMC_Ye_Pre/p2_illumina_raw_filtered_Antibody.csv']
header_list=['PBMC_Stoeckius','CBMC_Stoeckius','BMMC_Stuart','PBMC_Ye_Pre_p1','PBMC_Ye_Pre_p2']

loc='../model_optimal/training_05152020'

repi=0

gene_list=None
X_final=None
y_final=None
for i,X_file in enumerate(X_list):
	print(X_file)
	X = pd.read_csv(X_file)
	if i==2:
		y = pd.read_csv(y_list[i],sep='\t')
	else:
		y = pd.read_csv(y_list[i])
	if interactive:
		X=X.transpose().sample(frac=0.1,random_state=4905).transpose()
	# Dealing with X's dimensionality
	if i==3 or i==4:
		gene=X.index.tolist()
		gene=[x[0:(len(x)-2)] for x in gene]
		X.index=gene
		gene=set(gene)
	else:
		gene = set(X.index.tolist())
	if gene_list is None:
		gene_list=gene
	else:
		gene_list=set(gene_list).intersection(gene)
	gene_list=list(gene_list)
	gene_list.sort()
	X=X.loc[gene_list,]
	if not X_final is None:
		X_final=X_final.loc[gene_list,]
	# Dealing with Y's dimensionality
	if i==2:
		protein=y.index.tolist()
		protein[protein.index('CD8a')]='CD8'
		protein[protein.index('CD127-IL7Ra')]='CD127'
		protein[protein.index('HLA.DR')]='HLA-DR'
		protein[protein.index('CD197-CCR7')]='CD197'
		protein[protein.index('CD278-ICOS')]='CD278'
	elif i==3 or i==4:
		protein=y.index.tolist()
		protein=[x.split("__")[0] for x in protein]
	else:
		protein = y.iloc[:,0].values.tolist()
		y=y.drop(columns=y.columns[0])
		if i==1:
			protein[protein.index('CCR5')]='CD195'
			protein[protein.index('CCR7')]='CD197'
	print(protein)
	y.index=protein
	y=y[X.columns]
	y=y.loc[protein_list,]
	# Add header to cell
	X.columns=list(map(lambda x: header_list[i]+'-'+x, X.columns.tolist()))
	y.columns=list(map(lambda x: header_list[i]+'-'+x, y.columns.tolist()))
	if i==0:
		X_final=X
		y_final=y
	else:
		X_final=pd.concat([X_final,X], axis=1)
		y_final=pd.concat([y_final,y], axis=1)

# Normalize y
shared_proteins=y_final.apply(lambda x: not x.isna().any(),axis=1)
y=y_final.apply(lambda x: np.log((x+1.0)/stats.gmean(x[shared_proteins]+1.0)), axis=0)
del(y_final)
# Use direct UMI count is hard to transfer across experiment
# Let's try normalize the data with seurat like method
X=X_final.apply(lambda x: np.log((x*10000.0/sum(x))+1.0), axis=0)
del(X_final)
# random cell order
X=X.T
X=X.sample(frac=1,random_state=4905)
# separate test data from train data
if interactive:
	X_test=X.sample(n=918,random_state=4905)# Need change after test
else:
	X_test=X.sample(n=5187,random_state=4905)# Need change after test

y_test=y[X_test.index]

X=X.drop(X_test.index)
y=y.drop(columns=y_test.columns)
y=y[X.index]


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
n_batches=32


class Net(nn.Module):
    def __init__(self):
        super(Net, self).__init__()
        #self.fc1 = nn.Linear(14505, 1000)
        self.fc1 = nn.Linear(X.shape[1], 1000)
        self.fc2 = nn.Linear(1000, 256)
        self.fc3 = nn.ModuleDict({})
        for p in protein_list:
        	self.fc3[p]=nn.Linear(256, 64)
        self.fc4 = nn.ModuleDict({})
        for p in protein_list:
        	self.fc4[p]=nn.Linear(64, 1)
    def forward(self, x):
        x = F.relu(self.fc1(x))
        x = F.relu(self.fc2(x))
        outputs={}
        for p in protein_list:
        	outputs[p]=self.fc4[p](F.relu(self.fc3[p](x)))
        return outputs

net = Net()
if repi==4:
	net.load_state_dict(torch.load('../model_optimal/training_04262020model_rep'+str(repi)+'_ep29'))
else:
	net.load_state_dict(torch.load('../model_optimal/training_04262020model_rep'+str(repi)+'_ep54'))


criterion = nn.MSELoss()
optimizer = optim.Adam(net.parameters(), lr=0.001,amsgrad=True, weight_decay=0.001)
# optimizer = optim.Adagrad(net.parameters(), lr=lr_vec[repi], lr_decay=0.001)

max_epochs=200

train_loss=pd.DataFrame(np.zeros(shape=(len(protein_list),max_epochs)),index=protein_list)
test_loss=pd.DataFrame(np.zeros(shape=(len(protein_list),max_epochs)),index=protein_list)

# Init early stop
patience=30
best_score=None
Dy=len(protein_list)
estop_counter=pd.Series(np.zeros(Dy),index=protein_list)
early_stop=pd.Series([False]*Dy,index=protein_list)

for epoch in range(max_epochs):
	if all(early_stop):
		break
	running_loss=pd.Series(np.zeros(Dy),index=protein_list)
	for i in range(int(y.shape[0]/n_batches)):
		# Local batches and labels
		local_X, local_y = X[i*n_batches:min((i+1)*n_batches,X.shape[0]-1),], y[i*n_batches:min((i+1)*n_batches,y.shape[0]-1),]
		# zero the parameter gradients
		optimizer.zero_grad()
		# forward + backward + optimize
		outputs_dict = net(local_X)
		loss=None
		loss_count=0.0
		for p in protein_list:
			notNaN=(local_y[:,protein_list.index(p):(protein_list.index(p)+1)]==local_y[:,protein_list.index(p):(protein_list.index(p)+1)])
			loss_p=criterion(outputs_dict[p][notNaN],local_y[:,protein_list.index(p):(protein_list.index(p)+1)][notNaN])
			if not torch.isnan(loss_p):
				loss_count+=1.0
				running_loss[p]+=loss_p.item()
				if loss is None:
					loss=loss_p
				else:
					loss=loss+loss_p
		loss.backward()
		optimizer.step()
		if(i==(int(y.shape[0]/n_batches)-1)):
			train_loss.iloc[:,epoch]=(running_loss / 150)
		if i % 150 == 149:    # print every mini-batches
			print('[%d, %5d] loss: %.3f' % (epoch + 1, i + 1, sum(running_loss / 150)))
			running_loss=pd.Series(np.zeros(Dy),index=protein_list)
			sys.stdout.flush()
	test_outputs = net(X_test)
	test_outputs = [test_outputs[p] for p in protein_list]
	test_outputs=torch.transpose(torch.stack(test_outputs),0,1).view(X_test.shape[0],-1)
	test_loss_i=pd.Series([criterion(test_outputs[:,pi][y_test[:,pi]==y_test[:,pi]], y_test[:,pi][y_test[:,pi]==y_test[:,pi]]).item() for pi in range(Dy)],index=protein_list)
	test_loss.iloc[:,epoch]=test_loss_i
	if epoch % 10 == 9:
		f,ax=plt.subplots(figsize=(6,6))
		ax.scatter(y_test.detach().numpy(),test_outputs.detach().numpy())
		ax.plot([-2,5],[-2,5],ls='--',c='.3')
		#ax.text(3,-2,'correlation: '+str(np.corrcoef(test_outputs.detach().numpy().flatten(),y_test.detach().numpy().flatten())[1,0]))
		df = pd.DataFrame({"y_pred":test_outputs.detach().numpy().flatten(),'y_truth':y_test.detach().numpy().flatten()})
		ax.text(3,-2,'correlation: '+str(round(df.corr().values[1,0],4)))
		fig = ax.get_figure()
		fig.savefig(loc+'figure_rep'+str(repi)+'_ep'+str(epoch)+'.pdf')
		sys.stdout.flush()
		plt.close(fig)
	if epoch % 5 == 4:
		torch.save(net.state_dict(), loc+'model_rep'+str(repi)+'_ep'+str(epoch))
	# Implement early stopping
	if best_score is None:
		best_score=test_loss_i
	else:
		for p in protein_list:
			if test_loss_i[p]>(best_score[p]-0.001) and (not early_stop[p]):
				estop_counter[p]+=1
				if estop_counter[p]>=patience:
					early_stop[p]=True
			else:
				best_score[p]=test_loss_i[p]
				estop_counter[p]=0
		print(estop_counter)

print('Finished Training')


torch.save(net.state_dict(), loc+'model_rep'+str(repi)+'_ep'+str(epoch))

train_loss.index=['train_'+p for p in protein_list]
test_loss.index=['test_'+p for p in protein_list]
log=pd.concat([train_loss,test_loss])

log.to_csv(loc+'log_rep'+str(repi)+'.csv')

