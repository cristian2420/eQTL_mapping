import sys
import peer
import scipy as SP
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

#########
ifile = sys.argv[1] #files with gene counts
confounders = int(sys.argv[2]) #number of peer factors to find
odir = sys.argv[3] #output directory
mstvar = sys.argv[4] #get genes with most variance
n_genes = sys.argv[5] #number of genes with most variance
if n_genes != 'NA':
    n_genes = int(n_genes)
#######
#ifile='/mnt/bioadhoc-temp/Groups/vd-vijay/Cristian/DICE_LungCancer/eQTL_pipeline/data/expression/input/tumor/CD8/0/gene_expression.txt'
#
expr = pd.read_table(ifile)
#####if no peer factors, it would be set to 1 just to run the analysis. THe pipeline will not take it into consideration
if confounders == 0:
    confounders = 1

####
donors =  list(expr.columns)
genes =  list(expr.index.values)
###
if mstvar == 'True':
    sort_var_genes = expr.var(axis=1).sort_values(ascending=False)
    genes_use = sort_var_genes.axes[0].tolist()[0:n_genes]
    expr = expr[expr.index.isin(genes_use)]
    expr = expr.to_numpy()
else:
    expr = expr.to_numpy()

expr = expr.transpose()

model = peer.PEER()
model.setPhenoMean(expr)
model.getPhenoMean().shape
model.setNk(confounders)
model.getNk()
#model.setNmax_iterations()
model.update()

#####
factors = model.getX()
factors.shape
weights = model.getW()
weights.shape
precision = model.getAlpha()
precision.shape
residuals = model.getResiduals()
residuals.shape


#### plot
df1 = pd.DataFrame(precision, columns = ['Alpha'])
df1['Factor'] = range(1,len(precision)+1)
df1['variance_perc'] = 1/df1['Alpha']
fig, (ax1, ax2) = plt.subplots(1, 2,figsize=(8,10))
fig.suptitle('PEER Factors:' + str(confounders))
ax1.scatter(df1['Factor'],precision)
ax1.plot(df1['Factor'], precision)
ax2.scatter(df1['Factor'], 1/precision)
ax2.plot(df1['Factor'], 1/precision)

ax1.set_xlabel("Factors")
ax1.set_ylabel("Alpha")
ax2.set_xlabel("Factors")
ax2.set_ylabel("Inverse of variance")
# sns.set()
# sns.lineplot( x="Factor", y="Alpha", markers=True, data = df1, ax=axs[0])
# sns.lineplot( x="Factor", y="variance_perc", markers=True, data = df1, ax=axs[1])
plt.savefig(odir + '/PF_stats.pdf', format = 'pdf')


####save factors
factors = pd.DataFrame(factors)
factors.columns = ["Factor" + str(el) for el in range(1, confounders+1)]
factors.index = donors

factors.to_csv(odir + '/factors.txt', sep = ',')
df1.to_csv(odir + '/PF_stats.txt', sep = '\t', index=False)
