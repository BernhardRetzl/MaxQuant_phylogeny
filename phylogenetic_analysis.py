import re
from ete3 import NCBITaxa
ncbi = NCBITaxa()


tax_list = list()
with open('proteinGroups.txt') as file:
    for line in file:
        line = line.strip().split('\t')
        uni_prot_info = line[7]
        taxonomy_infos = re.findall(r'OX=\d+', uni_prot_info)
        for tax_info in taxonomy_infos:
            tax_list.append(tax_info)

tax_list = list(set(tax_list))

tax_dict = dict()
for element in tax_list:
    tax_dict[element] = list()


with open('proteinGroups.txt') as file:
    for line in file:
        line = line.strip().split('\t')
        uni_prot_info = line[7]
        taxonomy_infos = re.findall(r'OX=\d+', uni_prot_info)
        peptide_sequence = line[27]
        peptide_sequence = ''.join(peptide_sequence.split(';'))
        peptide_dummy = '-'*len(peptide_sequence)
        for item in tax_dict:
            if item in taxonomy_infos:
                tax_dict[item].append(peptide_sequence)
            else:
                tax_dict[item].append(peptide_dummy)

for item in tax_dict:
    tax_dict[item] = ''.join(tax_dict[item])


with open('alignment.txt', 'wt') as out_file:
    for item in tax_dict:
        out_file.write('>'+item+'\n')
        out_file.write(tax_dict[item]+'\n')



length_dict = dict()
for item in tax_dict:
    sequence_length = ''.join(tax_dict[item].split('-'))
    length_dict[item] = len(sequence_length)

length_dict_items = list(length_dict.items())
length_dict_items.sort(key=lambda x: x[1], reverse=True)


for item in length_dict_items:
    ox = int(item[0].split('=')[-1])
    species = ncbi.get_taxid_translator([ox, ])
    print(species, item[1])



for item in length_dict_items:
    ox = int(item[0].split('=')[-1])
    species = ncbi.get_taxid_translator([ox, ])
    print(species, item[1])


seq_1 = tax_dict['OX=9606']
seq_2 = tax_dict['OX=435590']

sim = 0
not_sim = 0
for pos in range(len(seq_1)):
    if seq_2[pos]!= '-':
        if seq_1[pos] == seq_2[pos]:
            sim += 1
        else:
            not_sim += 1



import pandas as pd
from scipy.spatial.distance import pdist, squareform
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import linkage, dendrogram


data = {'OX': list(tax_dict.keys()), 'OX_2': list(tax_dict.keys()), 'Sequence': list(tax_dict.values())}
df = pd.DataFrame(data)
df.set_index('OX', inplace=True)

def calculate_similarity(seq1, seq2):
    sim = 0
    for pos in range(len(seq1)):
        if seq2[pos] != '-':
            if seq1[pos] == seq2[pos]:
                sim += 1
    return sim


def similarity_func(u, v):
    a = u
    seq_1 = u[1]
    len_seq_1 = length_dict[u[0]]
    seq_2 = v[1]
    len_seq_2 = length_dict[v[0]]
    return calculate_similarity(seq_1, seq_2)/min(len_seq_1, len_seq_2)

dists = pdist(df, similarity_func)
similarity_matrix = pd.DataFrame(squareform(dists), columns=df.index, index=df.index)
similarity_matrix = similarity_matrix.astype(float)

import numpy as np

# np.fill_diagonal(similarity_matrix.values, np.nan)

fig, ax = plt.subplots(figsize=(40,40))
sns.heatmap(similarity_matrix, annot=True, cmap='gnuplot2_r', linewidths=.5, ax=ax, annot_kws={"size": 6})





plt.show()

import seaborn as sns
sns.clustermap(similarity_matrix, cmap='gnuplot2_r', figsize=(30,30), annot=True, fmt='.1g')
plt.show()