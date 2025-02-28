#!/usr/bin/env python
# coding: utf-8

# In[72]:


#!/usr/bin/env python

import os
import pandas as pd

# Set the working directory
os.chdir('/scratch/richards/chen-yang.su/01.pQTL/01.aric/03.cis_pQTL/')
os.system('mkdir -p /scratch/richards/chen-yang.su/01.pQTL/01.aric/03.cis_pQTL/all_cis_pQTL/')
os.system('mkdir -p /scratch/richards/chen-yang.su/01.pQTL/01.aric/03.cis_pQTL/all_strict_cis/')

# Read the protein data
all_proteins = pd.read_csv('/scratch/richards/chen-yang.su/01.pQTL/01.aric/cis_trans_pQTL_identification/all_cis_trans_res.tsv', sep='\t')


# In[9]:


# Read the dictionary file for SeqId and EntrezGeneSymbol
sl = pd.read_csv('/scratch/richards/chen-yang.su/data/01.SL_proteins/SL_proteins.txt', sep='\t')

# Separate SeqId into seqid1 and seqid2, then unite them with '_'
sl[['seqid1', 'seqid2']] = sl['SeqId'].str.split('-', expand=True)
sl['seqid'] = sl['seqid1'] + '_' + sl['seqid2']

# Select seqid and EntrezGeneSymbol, and filter out rows with missing EntrezGeneSymbol
sl = sl[['seqid', 'EntrezGeneSymbol']].dropna()


# In[15]:


# Create a list of retained seqids: 4595 seqids
retained_list = list(set(sl['seqid']) & set(all_proteins['seqid']))


# In[46]:


# print(retained_list)


# In[22]:


# example
#ith_seqid = '10000_28'
##ith_protein_extracted = all_proteins[all_proteins['seqid'] == ith_seqid]
#ith_protein_extracted = all_proteins[(all_proteins['seqid'] == ith_seqid) & (all_proteins['cis_trans'] == 'cis')]
#print(ith_protein_extracted)


# In[51]:


def extract_and_save_protein_data(ith_seqid):
    ith_seqid = str(ith_seqid)  # Ensure ith_seqid is treated as a string
    # Filter rows where "cis_trans" is "cis"
    cis_rows = all_proteins[(all_proteins['seqid'] == ith_seqid) & (all_proteins['cis_trans'] == 'cis')]
    
    # Check if there are any rows that meet the condition
    if not cis_rows.empty:
        ith_filename = cis_rows['protein'].unique()[0] + '.tsv'
        cis_rows.to_csv('/scratch/richards/chen-yang.su/01.pQTL/01.aric/03.cis_pQTL/all_cis_pQTL/cis.' + ith_filename, sep='\t', header=True, index=False)
    else:
        print(f"No 'cis' rows found for seqid: {ith_seqid}. Skipping...")


# In[52]:


# example (NEED QUOTES to ensure it's a string, not variable!)
# extract_and_save_protein_data(10000_28)
# extract_and_save_protein_data('10037_98') # with cis


# In[50]:


# ith_seqid = '10037_98'
# cis_rows = all_proteins[(all_proteins['seqid'] == ith_seqid) & (all_proteins['cis_trans'] == 'cis')]
# cis_rows


# In[53]:


# Iterate over retained seqids and extract/save protein data
for ith_seqid in retained_list:
    extract_and_save_protein_data(str(ith_seqid)) # need str() to ensure it's a string, not variable!


# In[67]:


#####################
# From here, I'll create a strict cis tsv
#####################
# Filter rows with 'cis_trans' equal to 'cis'
cis_variants = all_proteins[all_proteins['cis_trans'] == 'cis']


# This groups by seqid ----

# # Group by 'rs_id' and count unique 'seqid' values
# variant_counts = cis_variants.groupby('rs_id')['seqid'].nunique()
# 
# # Filter variants with two or more seqids
# variants_with_two_or_more_seqids = variant_counts[variant_counts >= 2]



# Group by gene symbol name instead ----

# Group by 'rs_id' and count unique 'gene' values
variant_counts = cis_variants.groupby('rs_id')['gene'].nunique()

# Filter variants with two or more genes
variants_with_two_or_more_genes = variant_counts[variant_counts >= 2]


# Save the DataFrame to a TSV file
variants_with_two_or_more_genes.to_csv('/scratch/richards/chen-yang.su/01.pQTL/01.aric/03.cis_pQTL/variants_removed_by_strict_filter.tsv', sep='\t', header=True, index=True)

# Get the list of variants with two or more genes
selected_variants = variants_with_two_or_more_genes.index.tolist()


# In[69]:


# remove "selected variants" from all_proteins
# Create a boolean condition to select rows not in 'selected_variants'
condition = ~all_proteins['rs_id'].isin(selected_variants)

# Use the condition to filter the DataFrame
all_proteins_strict = all_proteins[condition]


# In[73]:


# define function to save each protein's strict cis result
def extract_and_save_strict_protein_data(ith_seqid):
    ith_seqid = str(ith_seqid)  # Ensure ith_seqid is treated as a string
    # Filter rows where "cis_trans" is "cis"
    cis_rows = all_proteins_strict[(all_proteins_strict['seqid'] == ith_seqid) & (all_proteins_strict['cis_trans'] == 'cis')]
    
    # Check if there are any rows that meet the condition
    if not cis_rows.empty:
        ith_filename = cis_rows['protein'].unique()[0] + '.tsv'
        cis_rows.to_csv('/scratch/richards/chen-yang.su/01.pQTL/01.aric/03.cis_pQTL/all_strict_cis/cis.' + ith_filename, sep='\t', header=True, index=False)
    else:
        print(f"No 'cis' rows found for seqid: {ith_seqid}. Skipping...")


# In[74]:


# Iterate over retained seqids and extract/save strict cis protein data
for ith_seqid in retained_list:
    extract_and_save_strict_protein_data(str(ith_seqid)) # need str() to ensure it's a string, not variable!


# In[ ]:




