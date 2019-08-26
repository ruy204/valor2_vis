#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
import seaborn as sns


# ## Merge all bedpe files

# In[2]:


folder_location = "C:/PhD/Rotations/Rotation_1/data/SV2/bedpe/file_download/"


# In[3]:


#1. read in all bedpe files
bedpelist = list()
for a in os.listdir(folder_location):
    if a.endswith(".bedpe"):
        afile = pd.read_csv(folder_location + a,
                           sep="\t",lineterminator='\n',names=range(11))
        afile["samples"] = a.split("-")[0]
        afile.columns = ["Chrom1","start1","end1","Chrom2","start2","end2","type","score","number1","number2","number3","samples"]
        bedpelist.append(afile)

alldf = bedpelist[0]
for element in bedpelist[1:]:
    alldf = alldf.append(element)
print(alldf.shape)
print(alldf.samples.unique())

#2. add population information
popdf = pd.read_excel("C:/PhD/Rotations/Rotation_1/70_samples.xlsx")
popdf.columns = ["sample_names"]+list(popdf.columns.values[1:])
popdf["samples"] = popdf["sample_names"].str.split("_",1).str[1]
subpopdf = popdf[["samples","Description"]]
jointdf = pd.merge(alldf, subpopdf, on='samples')
print(jointdf.shape)

jointdf["Chrom1"] = pd.Categorical(jointdf["Chrom1"],categories = ["chr" + n for n in map(str,range(1,23))+["X","Y"]])
jointdf["Chrom2"] = pd.Categorical(jointdf["Chrom2"],categories = ["chr" + n for n in map(str,range(1,23))+["X","Y"]])
jointdf = jointdf.loc[:,["Chrom1","start1","end1","Chrom2","start2","end2","type","score","samples","Description"]]

jointdf.to_csv("C:/PhD/Rotations/Rotation_1/data/SV2/bedpe/combine/jointdf.bedpe",sep="\t",index=False,header=False)


# In[4]:


jointdf.head()


# In[ ]:





# In[ ]:




