#!/usr/bin/env python
# coding: utf-8

# In[1]:


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


# ## STATISTICS

# In[3]:


print(jointdf[(jointdf["Chrom2"]=="chr22") & (jointdf["type"]=="inverted-reciprocal")].start1.count())
print(jointdf[(jointdf["Chrom2"]=="chr22") & (jointdf["type"]=="reciprocal")].start1.count())
print(np.mean(jointdf[(jointdf["Chrom2"]=="chr22") & (jointdf["type"]=="inverted-reciprocal")].start2))
# jointdf[(jointdf["Chrom2"]=="chr4") & (jointdf["type"]=="reciprocal")].reset_index()
# jointdf[(jointdf["Chrom1"]=="chr4") & (jointdf["type"]=="inverted-reciprocal")]


# ### Mean and sd of #SV in each population

# In[34]:


a = jointdf.groupby(["Description","samples"]).Chrom1.count().reset_index()
a1 = a.groupby(["Description"]).Chrom1.mean().reset_index().sort_values(by="Chrom1")
a1.columns = ["Description","mean"]
a2 = a.groupby(["Description"]).Chrom1.std().reset_index().sort_values(by="Chrom1")
a2.columns = ["Description","sd"]
a3 = pd.merge(a1,a2,on="Description")
a3


# ### Number of structural variants

# In[26]:


a = jointdf.groupby(["type","samples"]).Chrom1.count().reset_index()
a1 = a.groupby(["type"]).Chrom1.mean().reset_index().sort_values(by="Chrom1")
a1.columns = ["type","mean"]
a2 = a.groupby(["type"]).Chrom1.std().reset_index().sort_values(by="Chrom1")
a2.columns = ["type","sd"]
a3 = pd.merge(a1,a2,on="type")
a3


# ## DATA TABLE

# In[29]:


datatable = jointdf.groupby(["Description","samples","type"]).Chrom1.count().reset_index()
datatable = datatable.sort_values(["Description","samples","type","Chrom1"])
datatable.columns = ["Description","samples","type","n_SV"]
datatable.to_csv("C:/PhD/Rotations/Rotation_1/writing_valor2/table/SV_summary.csv")


# In[30]:


a = datatable.groupby(["Description","type"]).n_SV.mean().reset_index()


# In[32]:


a.pivot(index="Description",columns="type").fillna(0).to_csv("C:/PhD/Rotations/Rotation_1/writing_valor2/table/SV_summary2.csv")


# In[11]:


datatable2 = datatable.groupby(["Description","type"]).n_SV.mean()


# ## GENERAL VISUALIZATION

# ### Number of SV by population

# In[171]:


jointdf = pd.read_csv("C:/PhD/Rotations/Rotation_1/data/SV2/bedpe/combine/jointdf.bedpe",sep="\t",header=None)
jointdf.columns = ["Chrom1","start1","end1","Chrom2","start2","end2","type","score","samples","Description"]
a = jointdf.groupby(["Description","samples"]).Chrom1.count().reset_index()
b = jointdf.groupby(["Description"]).samples.nunique().reset_index()
a = pd.merge(a,b,on="Description",how="outer")
a.columns = ["Description","samples","n_event","n_samples_population"]
a = a.sort_values(["n_samples_population","Description"],ascending=False)
ethranka = jointdf.groupby(["Description"]).samples.nunique().reset_index().sort_values(["samples"],ascending=False)
a["Description"] = pd.Categorical(a["Description"],categories=ethranka.loc[:,"Description"])


# In[173]:


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
g.figure.savefig("C:/PhD/Rotations/Rotation_1/writing_valor2/plots/nSV_population.pdf")


# ### Components of SV by population

# In[174]:


jointdf = pd.read_csv("C:/PhD/Rotations/Rotation_1/data/SV2/bedpe/combine/jointdf.bedpe",sep="\t",header=None)
jointdf.columns = ["Chrom1","start1","end1","Chrom2","start2","end2","type","score","samples","Description"]
a = jointdf.groupby(["Description","type"]).start1.count().reset_index()
b = jointdf.groupby(["Description"]).start1.count().reset_index()
a = pd.merge(a,b,on="Description",how="outer")
a.columns = ["Description","type","n_SV","total_SV"]
a["type"] = pd.Categorical(a["type"],categories=["inversion","deletion","translocation","inverted-translocation",
                                                 "reciprocal","inverted-reciprocal","duplication","inverted-duplication"],ordered=True)
a.loc[:,"SV_perc"] = a.loc[:,"n_SV"] / a.loc[:,"total_SV"]
a = a.sort_values(["Description","type"],ascending=False)
a.loc[:,"SV_sum_perc"] = a.groupby(["Description"]).SV_perc.cumsum()


# In[175]:


plt.figure(figsize=(16,10))
g = sns.barplot(x="Description",y="SV_sum_perc",data = a, hue = "type",dodge=False,palette="tab20") #GnBu_d is reverse
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
g.figure.savefig("C:/PhD/Rotations/Rotation_1/writing_valor2/plots/components_SV_population.pdf")


# ### Number of SV - PCA

# In[7]:


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
# plt.figure(figsize=(16,10))
plt.figure(figsize=(15,6))
g = sns.scatterplot("PC1","PC2",hue="ethnicity",data=pcdf2,palette="Paired",s=100)
g.axes.set_title("PCA based on number of SV, \n colored by population",size=20)
g.legend(frameon=False,title="Population")
plt.show()


# ## INVERSION AND DELETION
# 
# ### Subset preparation

# In[177]:


sub = jointdf[jointdf["type"].isin(["inversion","deletion"])]
sub.loc[:,"start"] = jointdf[['start1','end1']].min(axis=1)
sub.loc[:,"end"] = jointdf[['start2','end2']].min(axis=1)
sub.loc[:,"length"] = sub.loc[:,"end"] - sub.loc[:,"start"]
print(all(sub["Chrom1"]==sub["Chrom2"]))
sub.to_csv("C:/PhD/Rotations/Rotation_1/data/SV2/bedpe/combine/inv_del.bedpe",sep="\t",index=False,header=False)


# ### Process with bedtools:
# 
# ```pairToBed -a inv_del.bedpe -b /mnt/c/PhD/Rotations/Rotation_1/scripts/hg_ucsc.bed > inv_del_gene_annotation.bed```
# 
# ### Read back gene annotations for Inversion and Deletion

# In[192]:


hgbed = pd.read_csv("C:/PhD/Rotations/Rotation_1/data/SV2/bedpe/combine/inv_del_gene_annotation.bed",
                    sep='\t', lineterminator='\n',header=None)
# hgbed = hgbed.loc[:,[11,12,13,14,15,6,0,1,2,3,4,5,16]]
hgbed = hgbed.loc[:,[8,9,13,14,15,6,0,1,2,3,4,5,16]]
hgbed.columns = ["samples","ethnicity","Chrom1","position","end","type",
                "chra","starta","enda","chrb","startb","endb","transcription"]
hgbed["length"] = hgbed["end"] - hgbed["position"]
hgbed["Chrom1"] = pd.Categorical(hgbed["Chrom1"],categories = ["chr" + n for n in map(str,range(1,23))+["X","Y"]])
hgbed["chra"] = pd.Categorical(hgbed["chra"],categories = ["chr" + n for n in map(str,range(1,23))+["X","Y"]])
hgbed["chrb"] = pd.Categorical(hgbed["chrb"],categories = ["chr" + n for n in map(str,range(1,23))+["X","Y"]])
hgbed = hgbed.sort_values(by=["samples","Chrom1","position","chra","chrb"])
hgbed["ethnicity"] = hgbed["ethnicity"].replace("\r","",regex=True)
hgbed.to_csv("C:/PhD/Rotations/Rotation_1/writing_valor2/table/inv_del_gene.csv",index=False)


# ## HOTSPOTS

# In[5]:


jointdf.type.unique()


# ### Hotspot for inversion and deletion
# 
# Inversion and deletion hotspots are defined as in the same Mb location, more than 30 out of 68 samples have the event. 
# 
# We identified 7 inversion hotspots distributed on Chrom1, 5, 7, 14, but no deletions that have more than 30 samples happen within the same Mb unit. 

# In[22]:


sub = jointdf[jointdf["type"].isin(["inversion","deletion"])]
sub.loc[:,"start_location"] =  (sub[['start1', 'end1']].mean(axis=1)/1000000).round()
sub.loc[:,"end_location"] =  (sub[['start2', 'end2']].mean(axis=1)/1000000).round()


# In[61]:


start_summary = sub.groupby(["type","Chrom1","start_location"]).samples.nunique().reset_index()
end_summary = sub.groupby(["type","Chrom2","end_location"]).samples.nunique().reset_index()

deletion_summary_start = start_summary[start_summary["type"]=="deletion"]
deletion_summary_end = end_summary[end_summary["type"]=="deletion"]

inversion_summary_start = start_summary[start_summary["type"]=="inversion"]
inversion_summary_end = end_summary[end_summary["type"]=="inversion"]

# deletion_summary_start.sort_values(["samples"],ascending = False).head()

inversion_hotspot = inversion_summary_start[inversion_summary_start["samples"]>=30]
deletion_hotspot = deletion_summary_start[deletion_summary_start["samples"]>=30]


# In[142]:


inversion_top.head()


# In[95]:


#7 inversion hotspots
print(inversion_hotspot.sort_values(["samples"],ascending = False))

#population distribution of first hotspot (only 1 sample from Utah/Mormon does not have the inversion)
inversion_top = sub[(sub["type"]=="inversion") & (sub["start_location"]==104.0) & (sub["Chrom1"]=="chr1")]
a = inversion_top.groupby(["Chrom1","start_location","Description"]).samples.nunique().reset_index()
nsample_df = jointdf.groupby(["Description"]).samples.nunique().reset_index()
nsample_df.columns = ["Description","nsamples"]
a = pd.merge(a,nsample_df,on="Description",how="outer")
a

# sns.scatterplot(data=inversion_top,x=(sub[['start1', 'end1']].mean(axis=1)),y=(sub[['start2', 'end2']].mean(axis=1)),hue="Description",palette="Paired")


# ### Hotspot for inter-chromosome events
# 
# Consistent with previous definition of hotspots, events happen at same Mb unit for >=30 out of 68 samples are defined as hotspots.
# 
# We identified 3 and 5 hotspot at original locations for reciprocal and inverted-reciprocal respectively, and 4 hotspots at target locations for reciprocal and 6 for inverted-reciprocal. Among these hotspots, location 44 Mb on Chr22 and location 156 Mb on Chrom 7 are the most active spots for both original and target events, while location 7 Mb on Chr4 is a target hotspot only. 
# 
# Specifically, location 44Mb on Chr22 is identified as hotspot for both reciprocal and inverted-reciprocal, for both original and target locations. 
# 
# For other events, there're no hotspots found.

# #### General analysis

# In[105]:


sub = jointdf[~jointdf["type"].isin(["inversion","deletion"])]
sub.loc[:,"from_location"] =  (sub[['start1', 'end1']].mean(axis=1)/1000000).round()
sub.loc[:,"to_location"] =  (sub[['start2', 'end2']].mean(axis=1)/1000000).round()


# In[121]:


# hotspot for original location
from_summary = sub.groupby(["type","Chrom1","from_location"]).samples.nunique().reset_index()
# hotspot for new location
to_summary = sub.groupby(["type","Chrom2","to_location"]).samples.nunique().reset_index()

from_hotspot = from_summary[from_summary["samples"]>=30].sort_values(["type","samples"],ascending=False)
to_hotspot = to_summary[to_summary["samples"]>=30].sort_values(["type","samples"],ascending = False)


# #### Check top hotspot for (inv-)reciprocal: chr22 location 44mb

# In[139]:


reciprocal_from = sub[(sub["Chrom1"]=="chr22") & (sub["from_location"]==44.0)]
reciprocal_to = sub[(sub["Chrom2"]=="chr22") & (sub["to_location"]==44.0)]
reciprocal_from.to_csv("C:/PhD/Rotations/Rotation_1/plots/valor2/all/tables/reciprocal_inv_from_hotspot.csv",index=False)
reciprocal_to.to_csv("C:/PhD/Rotations/Rotation_1/plots/valor2/all/tables/reciprocal_inv_to_hotspot.csv",index=False)


# In[ ]:

