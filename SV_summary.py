#!/usr/bin/env python
# coding: utf-8

# In[111]:


import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
import seaborn as sns
import math
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


# In[2]:


jointdf = pd.read_csv("C:/PhD/Rotations/Rotation_1/data/SV2/bedpe/combine/jointdf.bedpe",sep="\t",header=None)
jointdf.columns = ["Chrom1","start1","end1","Chrom2","start2","end2","type","score","samples","Description"]
print(jointdf.type.unique())
jointdf.head()
#order data
jointdf["Chrom1"] = pd.Categorical(jointdf["Chrom1"],categories = ["chr" + n for n in map(str,range(1,23))+["X","Y"]])
ethrank = ['Finnish','SOUTHERN HAN CHINESE','PUERTO RICAN','UTAH/MORMON','YORUBA/Nigeria','HAN CHINESE/China',
           'JAPANESE/Japan','USA/MEXICAN','USA/AFRICAN','ITALY/TOSCANI','Caucasian']
jointdf["Description"] = pd.Categorical(jointdf["Description"],categories = ethrank)
jointdf = jointdf.sort_values(by=["Description","Chrom1","start1"],ascending=True)


# ## General Visualization

# ### Number of SV by population

# In[49]:


jointdf = pd.read_csv("C:/PhD/Rotations/Rotation_1/data/SV2/bedpe/combine/jointdf.bedpe",sep="\t",header=None)
jointdf.columns = ["Chrom1","start1","end1","Chrom2","start2","end2","type","score","samples","Description"]
a = jointdf.groupby(["Description","samples"]).Chrom1.count().reset_index()
b = jointdf.groupby(["Description"]).samples.nunique().reset_index()
a = pd.merge(a,b,on="Description",how="outer")
a.columns = ["Description","samples","n_event","n_samples_population"]
a = a.sort_values(["n_samples_population","Description"],ascending=False)
ethranka = jointdf.groupby(["Description"]).samples.nunique().reset_index().sort_values(["samples"],ascending=False)
a["Description"] = pd.Categorical(a["Description"],categories=ethranka.loc[:,"Description"])


# In[4]:


plt.figure(figsize=(16,10))
g = sns.boxplot(x="Description",y="n_event",data = a, hue = "n_samples_population",dodge=False,palette="Blues")
g.axes.set_title("Number of structural variants by population",fontsize=30)
g.set_ylabel("Number of SV in each sample",fontsize=25)
g.set_xlabel("Population",fontsize=25)
g.tick_params(labelsize=20)
plt.xticks(rotation=45)
plt.setp(g.get_legend().get_texts(), fontsize=20) # for legend text
plt.setp(g.get_legend().get_title(), fontsize=20)
plt.tight_layout()
plt.show()
# g.figure.savefig("C:/PhD/Rotations/Rotation_1/writing/plots/pdf/ncontig_population.pdf")


# ### Components of SV by population

# In[88]:


jointdf = pd.read_csv("C:/PhD/Rotations/Rotation_1/data/SV2/bedpe/combine/jointdf.bedpe",sep="\t",header=None)
jointdf.columns = ["Chrom1","start1","end1","Chrom2","start2","end2","type","score","samples","Description"]
a = jointdf.groupby(["Description","type"]).start1.count().reset_index()
b = jointdf.groupby(["Description"]).start1.count().reset_index()
a = pd.merge(a,b,on="Description",how="outer")
a.columns = ["Description","type","n_SV","total_SV"]
a.loc[:,"SV_perc"] = a.loc[:,"n_SV"] / a.loc[:,"total_SV"]
a = a.sort_values(["Description","type"])
a.loc[:,"SV_sum_perc"] = a.groupby(["Description"]).SV_perc.cumsum()
a = a.sort_values(["Description","SV_sum_perc"],ascending = False)


# In[97]:


plt.figure(figsize=(16,10))
g = sns.barplot(x="Description",y="SV_sum_perc",data = a, hue = "type",dodge=False,palette="GnBu") #GnBu_d is reverse
g.axes.set_title("Component of structural variants by population",fontsize=30)
g.set_ylabel("Percentage of SV by type",fontsize=25)
g.set_xlabel("Population",fontsize=25)
g.tick_params(labelsize=20)
g.legend(frameon=False,title="Type",loc='center right', bbox_to_anchor=(1.4, 0.5), ncol=1)
plt.xticks(rotation=45)
plt.setp(g.get_legend().get_texts(), fontsize=20) # for legend text
plt.setp(g.get_legend().get_title(), fontsize=20)
plt.tight_layout()
plt.show()


# ### Number of SV - PCA

# In[125]:


sumdf = jointdf.groupby(["Description","samples","type"]).start1.count().reset_index()
sumdf2 = sumdf.pivot("samples","type","start1").fillna(0)

#PCA
X = sumdf2.iloc[:,:].values
scaler = StandardScaler()
scaler.fit(X)
X=scaler.transform(X)
pca = PCA(n_components = 3)
pc = pca.fit_transform(X)
pcdf = pd.concat([pd.Series(list(sumdf2.index)), pd.DataFrame(pc)], axis=1)
pcdf.columns = ["samples","PC1","PC2","PC3"]

print(pca.explained_variance_ratio_)

#ethnicity information
ethdf = pd.read_excel("C:/PhD/Rotations/Rotation_1/70_samples.xlsx")
ethdf = ethdf.iloc[:,[0,1]]
ethdf.columns = ["samples","ethnicity"]
ethdf["samples"] = ethdf.loc[:,"samples"].str.split("_",n=1,expand=True).loc[:,1]
pcdf = pd.merge(pcdf,ethdf,on="samples",how="outer")
pcdf2 = pd.merge(pcdf,sumdf2,on="samples",how="outer")

#Visualization
plt.figure(figsize=(16,10))
g = sns.scatterplot("PC1","PC2",hue="ethnicity",data=pcdf2,palette="Paired",s=100).set_title("PCA based on number of SV, \n colored by population",size=20)
plt.show()

